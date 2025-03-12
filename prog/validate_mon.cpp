// Copyright 2018 Global Phasing Ltd.
//
// Part of gemmi-validate that does extra validation for cif files from
// the Refmac monomer dictionary ("gemmi validate --monomer").

#include "validate_mon.h"
#include <stdio.h>
#include "gemmi/chemcomp.hpp"     // for ChemComp
#include "gemmi/mmcif.hpp"        // for make_residue_from_chemcomp_block
#include "gemmi/topo.hpp"         // for Topo
#include "gemmi/calculate.hpp"    // for find_best_plane

namespace cif = gemmi::cif;
using gemmi::Restraints;
using gemmi::Topo;
using gemmi::ChemComp;

namespace {

// some rules for the number of bonds (currently only for H and P)
void check_valency(const ChemComp& cc) {
  for (const ChemComp::Atom& atom : cc.atoms) {
    if (cc.atoms.size() == 1)
      continue;
    float valency = 0.0f;
    for (const Restraints::Bond& bond : cc.rt.bonds)
      if (bond.id1 == atom.id || bond.id2 == atom.id)
        valency += order_of_bond_type(bond.type);
    bool ok = valency >= 0.5f;
    valency -= atom.charge;
    if (atom.is_hydrogen()) {
      ok = std::round(valency) == 1.0;
    } else if (atom.el == gemmi::El::P) {
      ok = (valency == 3.0f || valency == 5.0f || valency == 5.5f);
    }
    if (!ok)
      printf("%s [valency] %s (%s) has bond order %g\n", cc.name.c_str(),
             atom.id.c_str(), element_name(atom.el), valency);
  }
}

void check_bond_angle_consistency(const ChemComp& cc) {
  const std::string tag = cc.name + " [restr]";
  for (const Restraints::Angle& angle : cc.rt.angles) {
    if (!cc.rt.are_bonded(angle.id1, angle.id2) ||
        !cc.rt.are_bonded(angle.id2, angle.id3))
      printf("%s angle %s with non-bonded atoms\n", tag.c_str(),
             angle.str().c_str());
    if (angle.value < 20)
      printf("%s angle %s with low value: %g\n", tag.c_str(),
             angle.str().c_str(), angle.value);
  }
  for (const Restraints::Torsion& tor : cc.rt.torsions) {
    if (!cc.rt.are_bonded(tor.id1, tor.id2) ||
        !cc.rt.are_bonded(tor.id2, tor.id3) ||
        !cc.rt.are_bonded(tor.id3, tor.id4))
      printf("%s torsion %s with non-bonded atoms\n", tag.c_str(),
             tor.str().c_str());
  }
}

template <typename T>
bool check_esd(const std::string& name, const T* restr) {
  if (restr->esd <= 0.) {
    printf("%s [esd] %s %s has non-positive esd: %g\n", name.c_str(),
           restr->what(), restr->str().c_str(), restr->esd);
    return false;
  }
  return true;
}

void print_outliers(const Topo& topo, const std::string& name, double z_score) {
  for (const Topo::Bond& t : topo.bonds) {
    if (!check_esd(name, t.restr))
      continue;
    double value = t.calculate();
    if (std::abs(value - t.restr->value) > z_score * t.restr->esd)
      printf("%s [atom.xyz] bond %s should be %g (esd %g) but is %.2f\n", name.c_str(),
             t.restr->str().c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Angle& t : topo.angles) {
    if (!check_esd(name, t.restr))
      continue;
    double value = gemmi::deg(t.calculate());
    if (gemmi::angle_abs_diff(value, t.restr->value) > z_score * t.restr->esd)
      printf("%s [atom.xyz] angle %s should be %g (esd %g) but is %.2f\n", name.c_str(),
             t.restr->str().c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Torsion& t : topo.torsions) {
    if (!check_esd(name, t.restr))
      continue;
    double value = gemmi::deg(t.calculate());
    double full = 360. / std::max(1, t.restr->period);
    if (gemmi::angle_abs_diff(value, t.restr->value, full) > z_score * t.restr->esd)
      printf("%s [atom.xyz] torsion %s should be %g (period %d, esd %g) but is %.2f\n",
             name.c_str(),
             t.restr->str().c_str(), t.restr->value, t.restr->period, t.restr->esd, value);
  }
  for (const Topo::Chirality& t : topo.chirs) {
    double value = t.calculate();
    if (t.restr->is_wrong(value))
      printf("%s [atom.xyz] chir %s should be %s but is %.2f\n", name.c_str(),
             t.restr->str().c_str(), gemmi::chirality_to_string(t.restr->sign),
             value);
  }
  for (const Topo::Plane& t : topo.planes) {
    if (!check_esd(name, t.restr))
      continue;
    auto coeff = find_best_plane(t.atoms);
    for (const gemmi::Atom* atom : t.atoms) {
      double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
      if (dist > z_score * t.restr->esd)
        printf("%s [atom.xyz] plane %s has atom %s in a distance %.2f\n", name.c_str(),
               t.restr->str().c_str(), atom->name.c_str(), dist);
    }
  }
}

struct HeavyAtom {
  const ChemComp::Atom* atom;
  std::vector<std::string> hydrogens;
  HeavyAtom(const ChemComp::Atom* atom_) : atom(atom_) {}
  bool operator<(const HeavyAtom& o) const { return atom->id < o.atom->id; }
};

struct SortedAtoms {
  std::vector<HeavyAtom> heavys;
  std::map<std::string, const Restraints::Bond*> bond_map;
};

SortedAtoms sorted_heavy_atoms(const ChemComp& cc) {
  SortedAtoms sa;
  std::map<std::string, int> hydrogen_names;
  for (const ChemComp::Atom& a : cc.atoms) {
    if (a.is_hydrogen())
      hydrogen_names.emplace(a.id, 0);
    else
      sa.heavys.emplace_back(&a);
  }
  std::sort(sa.heavys.begin(), sa.heavys.end());
  auto process_h = [&](const std::string& h_name, const std::string& parent_name) -> bool {
    auto h = hydrogen_names.find(h_name);
    if (h == hydrogen_names.end())
      return false;
    if (h->second != 0)
      gemmi::fail("2+ bonds for hydrogen ", h_name, " in ", cc.name);
    h->second = 1;
    auto parent = std::find_if(sa.heavys.begin(), sa.heavys.end(),
                               [&](const HeavyAtom& a) { return a.atom->id == parent_name; });
    if (parent == sa.heavys.end())
      gemmi::fail("missing parent atom for hydrogen ", h_name, " in ", cc.name);
    parent->hydrogens.push_back(h_name);
    return true;
  };

  for (const Restraints::Bond& bond : cc.rt.bonds) {
    if (process_h(bond.id2.atom, bond.id1.atom) ||
        process_h(bond.id1.atom, bond.id2.atom)) {
      if (bond.type != gemmi::BondType::Single)
        printf("%s [ccd] bond %s (hydrogen atom) is not SINGle\n",
               cc.name.c_str(), bond.str().c_str());
    } else {
      auto result = sa.bond_map.emplace(bond.lexicographic_str(), &bond);
      if (!result.second)
        printf("%s [ccd:bond]   duplicated bond %s\n", cc.name.c_str(),
               result.first->first.c_str());
    }
  }
  for (const auto& h : hydrogen_names)
    if (h.second == 0)
      gemmi::fail("hydrogen atom without any bond: ", h.first, " in ", cc.name);
  for (HeavyAtom& heavy : sa.heavys)
    std::sort(heavy.hydrogens.begin(), heavy.hydrogens.end());
  return sa;
}

// both arguments are sorted vectors
bool is_subset(const std::vector<std::string>& subset,
               const std::vector<std::string>& superset) {
  size_t i = 0, j = 0;
  for (;;) {
    if (i == subset.size())
      return true;
    if (j == superset.size())
      return false;
    if (subset[i] == superset[j])
      ++i;
    ++j;
  }
  gemmi::unreachable();
}

void check_consistency_with_ccd(const ChemComp& lib, const cif::Block& ccd_block,
                                bool verbose) {
  const char* name = lib.name.c_str();
  const ChemComp ccd = gemmi::make_chemcomp_from_block(ccd_block);

  // compare heavy atoms
  bool same_atoms = true;
  SortedAtoms lib_sa = sorted_heavy_atoms(lib);
  SortedAtoms ccd_sa = sorted_heavy_atoms(ccd);
  if (lib_sa.heavys.size() == ccd_sa.heavys.size()) {
    for (size_t i = 0; i != lib_sa.heavys.size(); ++i) {
      if (lib_sa.heavys[i].atom->id == ccd_sa.heavys[i].atom->id) {
        const HeavyAtom& lib_heavy = lib_sa.heavys[i];
        const HeavyAtom& ccd_heavy = ccd_sa.heavys[i];
        const char* atom_id = lib_heavy.atom->id.c_str();
        if (lib_heavy.atom->el != ccd_heavy.atom->el)
          printf("%s [ccd] different element for %s\n", name, atom_id);
        if (lib_heavy.hydrogens.size() != ccd_heavy.hydrogens.size())
          printf("%s [ccd] different protonation of %s, #H=%zu (%zu in CCD)\n",
                 name, atom_id, lib_heavy.hydrogens.size(), ccd_heavy.hydrogens.size());
        const auto* shorter = &lib_heavy.hydrogens;
        const auto* longer = &ccd_heavy.hydrogens;
        bool ccd_less = (longer->size() < shorter->size());
        if (ccd_less)
          std::swap(shorter, longer);
        if (is_subset(*shorter, *longer)) {
          if (ccd_less)
            for (const std::string& h : lib_heavy.hydrogens)
              for (const HeavyAtom& heavy : ccd_sa.heavys) {
                if (&heavy != &ccd_heavy && gemmi::in_vector(h, heavy.hydrogens))
                  printf("%s [ccd] wrong parent for hydrogen %s: %s (%s in CCD)\n",
                         name, h.c_str(), atom_id, heavy.atom->id.c_str());
              }
        } else {
          printf("%s [ccd] different names of hydrogens on %s: %s  vs  %s\n",
                 name, atom_id,
                 gemmi::join_str(lib_heavy.hydrogens, ' ').c_str(),
                 gemmi::join_str(ccd_heavy.hydrogens, ' ').c_str());
        }
      } else {
        printf("%s [ccd] different names of heavy atoms\n", name);
        same_atoms = false;
        break;
      }
    }
  } else {
    printf("%s [ccd] different number of heavy atoms: %zu (%zu in CCD)\n", name,
           lib_sa.heavys.size(), ccd_sa.heavys.size());
    same_atoms = false;
  }
  if (!same_atoms && verbose) {
    auto getter = [](const HeavyAtom& a) { return a.atom->id; };
    printf("%s [ccd]   %s\n", name, gemmi::join_str(lib_sa.heavys, ' ', getter).c_str());
    printf("%s [ccd]   %s\n", name, gemmi::join_str(ccd_sa.heavys, ' ', getter).c_str());
    cif::Table audit = const_cast<cif::Block&>(ccd_block)
                         .find("_pdbx_chem_comp_audit.", {"action_type", "date"});
    int audit_len = audit.length();
    if (audit_len > 0) {
      cif::Table::Row last_row = audit[audit_len-1];
      printf("%s [ccd]   Last modification: %s  %s\n", name,
             last_row[1].c_str(), last_row.str(0).c_str());
    } else {
      printf("%s [ccd]   missing _pdbx_chem_comp_audit in CCD\n", name);
    }
  }

  // check bonds between heavy atoms
  for (const auto& lib_bond : lib_sa.bond_map) {
    const std::string& bond_str = lib_bond.first;
    auto ccd_iter = ccd_sa.bond_map.find(bond_str);
    if (ccd_iter == ccd_sa.bond_map.end()) {
      printf("%s [ccd:bond]   extra bond %s\n", name, bond_str.c_str());
      continue;
    }
    if (lib_bond.second->type != ccd_iter->second->type && verbose)
      printf("%s [ccd:bond]   %s bond type is: %s (%s in CCD)\n", name,
             bond_str.c_str(),
             gemmi::bond_type_to_string(lib_bond.second->type),
             gemmi::bond_type_to_string(ccd_iter->second->type));
    ccd_sa.bond_map.erase(ccd_iter);
  }
  for (const auto& ccd_iter : ccd_sa.bond_map)
    printf("%s [ccd:bond]   missing bond %s\n", name, ccd_iter.first.c_str());
}

}  // anonymous namespace

void check_monomer(const cif::Block& block, double z_score) {
  gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
  check_valency(cc);
  check_bond_angle_consistency(cc);
  if (z_score != +INFINITY) {
    // check consistency of _chem_comp_atom.x/y/z with restraints
    gemmi::Residue res = gemmi::make_residue_from_chemcomp_block(block, gemmi::ChemCompModel::Xyz);
    Topo topo;
    topo.apply_restraints(cc.rt, res, nullptr, gemmi::Asu::Same, '\0', '\0', false);
    print_outliers(topo, cc.name, z_score);
  }
}

void compare_monomer_with_ccd(const cif::Block& lib_block, const cif::Block& ccd_block,
                              bool verbose) {
  check_consistency_with_ccd(gemmi::make_chemcomp_from_block(lib_block), ccd_block, verbose);
}
