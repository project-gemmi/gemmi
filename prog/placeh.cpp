// Copyright 2018 Global Phasing Ltd.
//
// for testing riding_h.hpp

#include <gemmi/riding_h.hpp>
#include <stdio.h>
#include <gemmi/cif.hpp>   // for read_file
#include <gemmi/fail.hpp>  // for fail

using namespace gemmi;

void print_restraint_summary(const std::string& id, const ChemComp& cc) {
  printf("%-5s %-5s ", cc.name.c_str(), id.c_str());
  fflush(stdout);
  std::string heavy_atom;

  for (const Restraints::Bond& bond : cc.rt.bonds) {
    const Restraints::AtomId* other_end = nullptr;
    if (bond.id1 == id)
      other_end = &bond.id2;
    else if (bond.id2 == id)
      other_end = &bond.id1;
    if (other_end) {
      if (!heavy_atom.empty())
        fail("H atom with 2+ bonds");
      heavy_atom = other_end->atom;
    }
  }
  if (heavy_atom.empty())
    fail("non-bonded H");
  printf("%-4s ", heavy_atom.c_str());

  int angle_count = 0;
  for (const Restraints::Angle& angle : cc.rt.angles) {
    const Restraints::AtomId* other_end = nullptr;
    if (angle.id1 == id)
      other_end = &angle.id3;
    else if (angle.id3 == id)
      other_end = &angle.id1;
    if (other_end) {
      if (angle.id2 != heavy_atom)
        fail("_chem_comp_angle.atom_id_2 is not H's heavy atom.");
      //printf("%.1f deg to %s, ", angle.value, other_end->atom.c_str());
      if (!cc.get_atom(other_end->atom).is_hydrogen()) {
        ++angle_count;
      }
    }
  }

  int tor_count = 0;
  for (const Restraints::Torsion& tor : cc.rt.torsions) {
    //if (tor.period > 1) continue;
    const Restraints::AtomId* other[3] = {nullptr, nullptr, nullptr};
    if (tor.id1 == id) {
      other[0] = &tor.id2;
      other[1] = &tor.id3;
      other[2] = &tor.id4;
    } else if (tor.id4 == id) {
      other[0] = &tor.id3;
      other[1] = &tor.id2;
      other[2] = &tor.id1;
    }
    if (other[0]) {
      if (other[0]->atom != heavy_atom)
        fail("_chem_comp_tor atom next to H is not H's heavy atom.");
      if (!cc.get_atom(other[2]->atom).is_hydrogen())
        ++tor_count;
    }
  }

  int chir_count = 0;
  for (const Restraints::Chirality& chir : cc.rt.chirs) {
    if (chir.id1 == id || chir.id2 == id || chir.id3 == id) {
      if (chir.id_ctr != heavy_atom)
        fail("_chem_comp_chir atom next to H is not H's heavy atom.");
      ++chir_count;
    }
  }

  int plane_count = 0;
  for (const Restraints::Plane& plane : cc.rt.planes) {
    const std::vector<Restraints::AtomId>& ids = plane.ids;
    if (std::find(ids.begin(), ids.end(), id) != ids.end()) {
      if (std::find(ids.begin(), ids.end(), heavy_atom) == ids.end())
        fail("H in _chem_comp_plane without its heavy atom.");
      ++plane_count;
    }
  }

  printf("%d angles, %d torsions, %d chiralities, %d planes\n",
         angle_count, tor_count, chir_count, plane_count);
}

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    cif::Document doc = cif::read_file(argv[i]);
    for (const cif::Block& block : doc.blocks)
      if (block.name != "comp_list") {
        ChemComp cc = make_chemcomp_from_block(block);
        for (const ChemComp::Atom& atom : cc.atoms)
          if (atom.el == El::H)
            try {
              print_restraint_summary(atom.id, cc);
            } catch (std::runtime_error& e) {
              fprintf(stderr, "%s %s:%s\n",
                      block.name.c_str(), atom.id.c_str(), e.what());
            }
      }
  }
}
