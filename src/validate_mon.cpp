// Copyright 2018 Global Phasing Ltd.
//
// Handles of "gemmi validate --monomer".
// Part of gemmi-validate that does extra validation for cif files from
// Refmac monomer dictionary.

#include <iostream>
#include <stdexcept>
#include "gemmi/cifdoc.hpp"
#include "gemmi/chemcomp.hpp"
#include "gemmi/calculate.hpp" // for calculate_angle

namespace cif = gemmi::cif;
using gemmi::Restraints;
using gemmi::Topo;

// some rules for the number of bonds (currently only for H and P)
static void check_valency(const gemmi::ChemComp& cc) {
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
      std::cout << cc.name << ": " << atom.id << " (" << element_name(atom.el)
                << ") has bond order " << valency << std::endl;
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

static void check_bond_angle_consistency(const gemmi::ChemComp& cc) {
  for (const Restraints::Angle& angle : cc.rt.angles) {
    if (!cc.rt.are_bonded(angle.id1, angle.id2) ||
        !cc.rt.are_bonded(angle.id2, angle.id3))
      std::cout << cc.name << ": angle " << repr(angle) << " not bonded"
                << std::endl;
    if (angle.value < 20)
      std::cout << cc.name << ": angle " << repr(angle)
                << " with low value: " << angle.value << std::endl;
  }
  for (const Restraints::Torsion& tor : cc.rt.torsions) {
    if (!cc.rt.are_bonded(tor.id1, tor.id2) ||
        !cc.rt.are_bonded(tor.id2, tor.id3) ||
        !cc.rt.are_bonded(tor.id3, tor.id4))
      std::cout << cc.name << ": torsion " << repr(tor) << " not bonded"
                << std::endl;
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

// check _chem_comp_atom.x/y/z
static void check_xyz_consistency(const gemmi::ChemComp& cc,
                                  const cif::Block& block) {
  gemmi::Residue res = read_chem_comp_atom_as_residue(block);
  Topo topo;
  topo.chains.resize(1);
  Topo::ResInfo ri(&res);
  ri.chemcomp = cc;
  topo.apply_internal_restraints_to_residue(ri);
  topo.chains[0].residues.push_back(ri);
  for (const Topo::Bond& t : topo.bonds) {
    double value = t.calculate();
    if (std::abs(value - t.restr->value) > t.restr->esd)
      std::cout << " [xyz bond] " << repr(*t.restr) << " should be "
                << t.restr->value << " (esd " << t.restr->esd << ") but is "
                << value << std::endl;
  }
  for (const Topo::Angle& t : topo.angles) {
    double value = gemmi::deg(t.calculate());
    if (angle_abs_diff(value, t.restr->value) > t.restr->esd)
      std::cout << " [xyz angle] " << repr(*t.restr) << " should be "
                << t.restr->value << " (esd " << t.restr->esd << ") but is "
                << value << std::endl;
  }
  for (const Topo::Torsion& t : topo.torsions) {
    double value = gemmi::deg(t.calculate());
    if (angle_abs_diff(value, t.restr->value) > t.restr->esd)
      std::cout << " [xyz torsion] " << repr(*t.restr) << " should be "
                << t.restr->value << " (esd " << t.restr->esd << ") but is "
                << value << std::endl;
  }
  // TODO chiral, plane
}

void check_monomer_doc(const cif::Document& doc) {
  for (const cif::Block& block : doc.blocks)
    if (block.name != "comp_list") {
      try {
        gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
        check_valency(cc);
        check_bond_angle_consistency(cc);
        check_xyz_consistency(cc, block);
      } catch (const std::exception& e) {
        std::cerr << "Failed to interpret " << block.name << " from "
                  << doc.source << ":\n " << e.what() << std::endl;
      }
    }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
