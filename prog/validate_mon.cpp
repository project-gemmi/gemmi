// Copyright 2018 Global Phasing Ltd.
//
// Part of gemmi-validate that does extra validation for cif files from
// the Refmac monomer dictionary ("gemmi validate --monomer").

#include "validate_mon.h"
#include <stdio.h>
#include <exception>
#include "gemmi/chemcomp.hpp"     // for ChemComp
#include "gemmi/chemcomp_xyz.hpp" // for make_residue_from_chemcomp_block
#include "gemmi/topo.hpp"         // for Topo
#include "gemmi/calculate.hpp"    // for find_best_plane

namespace cif = gemmi::cif;
using gemmi::Restraints;
using gemmi::Topo;

namespace {

// some rules for the number of bonds (currently only for H and P)
void check_valency(const gemmi::ChemComp& cc) {
  for (const gemmi::ChemComp::Atom& atom : cc.atoms) {
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

void check_bond_angle_consistency(const gemmi::ChemComp& cc) {
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

template<typename T, typename Compare>
std::vector<const T*> sorted_pointers(const std::vector<T>& vec, Compare comp, bool h) {
  std::vector<const T*> pp;
  for (const T& x : vec)
    if (h == x.is_hydrogen())
      pp.push_back(&x);
  std::sort(pp.begin(), pp.end(), comp);
  return pp;
}

void check_consistency_with_ccd(const gemmi::ChemComp& lib, const cif::Block& ccd_block,
                                bool verbose) {
  using gemmi::ChemComp;
  const char* name = lib.name.c_str();
  const ChemComp ccd = gemmi::make_chemcomp_from_block(ccd_block);
  auto cmp = [](const ChemComp::Atom* a, const ChemComp::Atom* b) { return a->id < b->id; };
  for (int h = 0; h <= 1; ++h) {  // bool : { false, true }
    bool same_atoms = true;
    const char* kind = !h ? "heavy atom" : "hydrogen";
    std::vector<const ChemComp::Atom*> lib_atoms = sorted_pointers(lib.atoms, cmp, !!h);
    std::vector<const ChemComp::Atom*> ccd_atoms = sorted_pointers(ccd.atoms, cmp, !!h);
    if (lib_atoms.size() == ccd_atoms.size()) {
      for (size_t i = 0; i != lib_atoms.size(); ++i) {
        if (lib_atoms[i]->id == ccd_atoms[i]->id) {
          if (lib_atoms[i]->el != ccd_atoms[i]->el)
            printf("%s [ccd] different element for %s\n", name, lib_atoms[i]->id.c_str());
        } else {
          printf("%s [ccd] different %ss names\n", name, kind);
          same_atoms = false;
          break;
        }
      }
    } else {
      printf("%s [ccd] different number of %ss: %zu (%zu in CCD)\n", name, kind,
             lib_atoms.size(), ccd_atoms.size());
      same_atoms = false;
    }
    if (!same_atoms && verbose) {
      auto getter = [](const ChemComp::Atom* a) { return a->id; };
      printf("%s [ccd]   %s\n", name, gemmi::join_str(lib_atoms, ' ', getter).c_str());
      printf("%s [ccd]   %s\n", name, gemmi::join_str(ccd_atoms, ' ', getter).c_str());
      cif::Table tab = const_cast<cif::Block&>(ccd_block)
                           .find("_pdbx_chem_comp_audit.", {"action_type", "date"});
      int audit_len = tab.length();
      if (audit_len > 0) {
        cif::Table::Row last_row = tab[audit_len-1];
        printf("%s [ccd]   Last modification: %s  %s\n", name,
               last_row[1].c_str(), last_row.str(0).c_str());
      } else {
        printf("%s [ccd]   missing _pdbx_chem_comp_audit in CCD\n", name);
      }
    }
  }

  // check bonds
  std::map<std::string, const Restraints::Bond*> ccd_bonds;
  for (const Restraints::Bond& bond : ccd.rt.bonds)
    ccd_bonds.emplace(bond.lexicographic_str(), &bond);
  std::string bond_str;
  for (const Restraints::Bond& bond : lib.rt.bonds) {
    std::string prev_bond_str = bond_str;
    bond_str = bond.lexicographic_str();
    if (bond_str == prev_bond_str) {
      printf("%s [ccd:bond]   duplicated bond %s\n", name, bond_str.c_str());
      continue;
    }
    auto ccd_iter = ccd_bonds.find(bond_str);
    if (ccd_iter == ccd_bonds.end()) {
      printf("%s [ccd:bond]   extra bond %s\n", name, bond_str.c_str());
    } else {
      if (bond.type != ccd_iter->second->type)
        printf("%s [ccd:bond]   %s bond type is: %s (%s in CCD)\n", name,
               bond_str.c_str(),
               gemmi::bond_type_to_string(bond.type),
               gemmi::bond_type_to_string(ccd_iter->second->type));
      ccd_bonds.erase(ccd_iter);
    }
  }
  for (const auto& it : ccd_bonds)
    printf("%s [ccd:bond]   missing bond %s\n", name, it.first.c_str());
}

}  // anonymous namespace

void check_monomer_doc(const cif::Document& doc, bool normal_checks, double z_score,
                       const std::map<std::string, cif::Block>& ccd_map, bool verbose) {
  for (const cif::Block& block : doc.blocks) {
    if (block.name == "comp_list")
      continue;
    gemmi::ChemComp cc;
    try {
      cc = gemmi::make_chemcomp_from_block(block);
    } catch (const std::exception& e) {
      fprintf(stderr, "Failed to interpret %s from %s:\n %s\n",
              block.name.c_str(), doc.source.c_str(), e.what());
    }
    if (normal_checks) {
      check_valency(cc);
      check_bond_angle_consistency(cc);
    }
    if (z_score != +INFINITY) {
      // check consistency of _chem_comp_atom.x/y/z with restraints
      gemmi::Residue res = gemmi::make_residue_from_chemcomp_block(block,
                                                  gemmi::ChemCompModel::Xyz);
      Topo topo;
      topo.apply_restraints(cc.rt, res, nullptr, gemmi::Asu::Same, '\0', '\0', false);
      print_outliers(topo, cc.name, z_score);
    }
    if (!ccd_map.empty()) {
      auto it = ccd_map.find(cc.name);
      if (it != ccd_map.end())
        check_consistency_with_ccd(cc, it->second, verbose);
      else
        printf("%s [ccd] monomer not found in provided CCD file(s)\n", cc.name.c_str());
    }
  }
}
