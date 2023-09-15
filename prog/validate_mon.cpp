// Copyright 2018 Global Phasing Ltd.
//
// Part of gemmi-validate that does extra validation for cif files from
// the Refmac monomer dictionary ("gemmi validate --monomer").

#include <stdio.h>
#include <stdexcept>
#include "gemmi/cifdoc.hpp"
#include "gemmi/chemcomp.hpp"     // for ChemComp
#include "gemmi/chemcomp_xyz.hpp" // for make_residue_from_chemcomp_block
#include "gemmi/topo.hpp"         // for Topo
#include "gemmi/calculate.hpp"    // for find_best_plane

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

static void check_bond_angle_consistency(const gemmi::ChemComp& cc) {
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

void print_outliers(const Topo& topo, const char* tag) {
  const double esd_mult = 2.0;
  for (const Topo::Bond& t : topo.bonds) {
    double value = t.calculate();
    if (std::abs(value - t.restr->value) > esd_mult * t.restr->esd)
      printf("%s bond %s should be %g (esd %g) but is %.2f\n", tag,
             t.restr->str().c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Angle& t : topo.angles) {
    double value = gemmi::deg(t.calculate());
    if (gemmi::angle_abs_diff(value, t.restr->value) > esd_mult * t.restr->esd)
      printf("%s angle %s should be %g (esd %g) but is %.2f\n", tag,
             t.restr->str().c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Torsion& t : topo.torsions) {
    double value = gemmi::deg(t.calculate());
    double full = 360. / std::max(1, t.restr->period);
    if (gemmi::angle_abs_diff(value, t.restr->value, full) > esd_mult * t.restr->esd)
      printf("%s torsion %s should be %g (esd %g) but is %.2f\n", tag,
             t.restr->str().c_str(), t.restr->value, t.restr->esd, value);
  }
  for (const Topo::Chirality& t : topo.chirs) {
    double value = t.calculate();
    if (t.restr->is_wrong(value))
      printf("%s chir %s should be %s but is %.2f\n", tag,
             t.restr->str().c_str(), gemmi::chirality_to_string(t.restr->sign),
             value);
  }
  for (const Topo::Plane& t : topo.planes) {
    auto coeff = find_best_plane(t.atoms);
    for (const gemmi::Atom* atom : t.atoms) {
      double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
      if (dist > esd_mult * t.restr->esd)
        printf("%s plane %s has atom %s in a distance %.2f\n", tag,
               t.restr->str().c_str(), atom->name.c_str(), dist);
    }
  }
}

void check_monomer_doc(const cif::Document& doc) {
  for (const cif::Block& block : doc.blocks)
    if (block.name != "comp_list") {
      try {
        gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
        check_valency(cc);
        check_bond_angle_consistency(cc);
        // check consistency of _chem_comp_atom.x/y/z with restraints
        gemmi::Residue res = gemmi::make_residue_from_chemcomp_block(block,
                                                    gemmi::ChemCompModel::Xyz);
        Topo topo;
        topo.apply_restraints(cc.rt, res, nullptr, gemmi::Asu::Same, '\0', '\0', false);
        print_outliers(topo, (cc.name + " [atom.xyz]").c_str());
      } catch (const std::exception& e) {
        fprintf(stderr, "Failed to interpret %s from %s:\n %s\n",
                block.name.c_str(), doc.source.c_str(), e.what());
      }
    }
}
