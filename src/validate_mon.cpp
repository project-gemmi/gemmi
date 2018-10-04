// Copyright 2018 Global Phasing Ltd.
//
// Handles of "gemmi validate --monomer".
// Part of gemmi-validate that does extra validation for cif files from
// Refmac monomer dictionary.

#include <stdio.h>
#include <stdexcept>
#include "gemmi/cifdoc.hpp"
#include "gemmi/chemcomp.hpp"
#include "gemmi/calculate.hpp" // for calculate_angle

namespace cif = gemmi::cif;
using gemmi::Restraints;
using gemmi::Topo;

// some rules for the number of bonds (currently only for H and P)
static void check_valency(const gemmi::ChemComp& cc) {
  const std::string tag = cc.name + " [valency]";
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
      printf("%s %s (%s) has bond order %g\n", tag.c_str(),
             atom.id.c_str(), element_name(atom.el), valency);
  }
}

static std::string repr(const Restraints::Bond& bond) {
  return bond.id1.atom + "-" + bond.id2.atom;
}
static std::string repr(const Restraints::Angle& angle) {
  return angle.id1.atom + "-" + angle.id2.atom + "-" + angle.id3.atom;
}
static std::string repr(const Restraints::Torsion& tor) {
  return tor.id1.atom + "-" + tor.id2.atom + "-" +
         tor.id3.atom + "-" + tor.id4.atom;
}
static std::string repr(const Restraints::Chirality& chir) {
  return chir.id_ctr.atom + "," + chir.id1.atom + "," +
         chir.id2.atom + "," + chir.id3.atom;
}

static std::string repr(const Restraints::Plane& plane) {
  return gemmi::join_str(plane.ids, ',',
                         [](const Restraints::AtomId& a) { return a.atom; });
}

static void check_bond_angle_consistency(const gemmi::ChemComp& cc) {
  const std::string tag = cc.name + " [restr]";
  for (const Restraints::Angle& angle : cc.rt.angles) {
    if (!cc.rt.are_bonded(angle.id1, angle.id2) ||
        !cc.rt.are_bonded(angle.id2, angle.id3))
      printf("%s angle %s with non-bonded atoms\n", tag.c_str(),
             repr(angle).c_str());
    if (angle.value < 20)
      printf("%s angle %s with low value: %g\n", tag.c_str(),
             repr(angle).c_str(), angle.value);
  }
  for (const Restraints::Torsion& tor : cc.rt.torsions) {
    if (!cc.rt.are_bonded(tor.id1, tor.id2) ||
        !cc.rt.are_bonded(tor.id2, tor.id3) ||
        !cc.rt.are_bonded(tor.id3, tor.id4))
      printf("%s torsion %s with non-bonded atoms\n", tag.c_str(),
             repr(tor).c_str());
  }
}

gemmi::Residue read_chem_comp_atom_as_residue(const cif::Block& block) {
  gemmi::Residue res;
  res.name = block.name;
  cif::Table table = const_cast<cif::Block&>(block).find(
      "_chem_comp_atom.", {"atom_id", "type_symbol", "x", "y", "z"});
  res.atoms.resize(table.length());
  int n = 0;
  for (auto row : table) {
    gemmi::Atom& atom = res.atoms[n++];
    atom.name = row.str(0);
    atom.element = gemmi::Element(row.str(1));
    atom.pos = gemmi::Position(cif::as_number(row[2]),
                               cif::as_number(row[3]),
                               cif::as_number(row[4]));
  }
  return res;
}

static double angle_abs_diff(double a, double b) {
  double d = std::abs(a - b);
  return d > 180 ? std::abs(d - 360.) : d;
}

void print_outliers(const Topo& topo, const char* tag) {
  const double esd_mult = 2.0;
  for (const Topo::Bond& t : topo.bonds) {
    double value = t.calculate();
    if (std::abs(value - t.restr->value) > esd_mult * t.restr->esd)
      printf("%s bond %s should be %g (esd %g) but is %.2f\n", tag,
             repr(*t.restr).c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Angle& t : topo.angles) {
    double value = gemmi::deg(t.calculate());
    if (angle_abs_diff(value, t.restr->value) > esd_mult * t.restr->esd)
      printf("%s angle %s should be %g (esd %g) but is %.2f\n", tag,
             repr(*t.restr).c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Torsion& t : topo.torsions) {
    double value = gemmi::deg(t.calculate());
    if (angle_abs_diff(value, t.restr->value) > esd_mult * t.restr->esd)
      printf("%s torsion %s should be %g (esd %g) but is %.2f\n", tag,
             repr(*t.restr).c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Chirality& t : topo.chirs) {
    double value = t.calculate();
    if ((t.restr->chir == Restraints::Chirality::Type::Positive && value < 0) ||
        (t.restr->chir == Restraints::Chirality::Type::Negative && value > 0))
      printf("%s chir %s should be %s but is %.2f\n", tag,
             repr(*t.restr).c_str(), gemmi::chirality_to_string(t.restr->chir),
             value);
  }
  for (const Topo::Plane& t : topo.planes) {
    auto coeff = find_best_plane(t.atoms);
    for (const gemmi::Atom* atom : t.atoms) {
      double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
      if (dist > esd_mult * t.restr->esd)
        printf("%s plane %s has atom %s in a distance %.2f\n", tag,
               repr(*t.restr).c_str(), atom->name.c_str(), dist);
    }
  }
}

static void check_xyz_consistency(const gemmi::ChemComp& cc,
                                  gemmi::Residue& res) {
  Topo topo;
  topo.chains.resize(1);
  Topo::ResInfo ri(&res);
  ri.chemcomp = cc;
  topo.apply_internal_restraints_to_residue(ri);
  topo.chains[0].residues.push_back(ri);
  print_outliers(topo, (cc.name + " [atom.xyz]").c_str());
}

void check_monomer_doc(const cif::Document& doc) {
  for (const cif::Block& block : doc.blocks)
    if (block.name != "comp_list") {
      try {
        gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
        check_valency(cc);
        check_bond_angle_consistency(cc);
        // check consistency of _chem_comp_atom.x/y/z with restraints
        gemmi::Residue res = read_chem_comp_atom_as_residue(block);
        check_xyz_consistency(cc, res);
      } catch (const std::exception& e) {
        fprintf(stderr, "Failed to interpret %s from %s:\n %s\n",
                block.name.c_str(), doc.source.c_str(), e.what());
      }
    }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
