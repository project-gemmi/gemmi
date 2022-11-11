// Copyright 2022 Global Phasing Ltd.
//
// Create cif::Block with monomer library _chem_comp* categories
// from struct ChemComp.

#ifndef GEMMI_TO_CHEMCOMP_HPP_
#define GEMMI_TO_CHEMCOMP_HPP_

#include "chemcomp.hpp"  // for ChemComp
#include "sprintf.hpp"   // for to_str
#include "calculate.hpp" // for calculate_angle

namespace gemmi {

inline void add_chemcomp_to_block(const ChemComp& cc, cif::Block& block) {
  {
    cif::Table tab = block.find_or_add("_chem_comp_atom.",
        {"comp_id", "atom_id", "type_symbol", "type_energy", "charge"});
    for (const ChemComp::Atom& a : cc.atoms)
      tab.append_row({cc.name, a.id, a.el.name(), cif::quote(a.chem_type),
                      std::to_string(iround(a.charge))});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_bond.",
        {"comp_id", "atom_id_1", "atom_id_2", "type", "aromatic",
         "value_dist", "value_dist_esd", "value_dist_nucleus", "value_dist_nucleus_esd"});
    for (const Restraints::Bond& a : cc.rt.bonds)
      tab.append_row({cc.name, a.id1.atom, a.id2.atom, bond_type_to_string(a.type),
                      std::string(1, a.aromatic ? 'y' : 'n'),
                      to_str(a.value), to_str(a.esd),
                      to_str(a.value_nucleus), to_str(a.esd_nucleus)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_angle.",
        {"comp_id", "atom_id_1", "atom_id_2", "atom_id_3",
         "value_angle", "value_angle_esd"});
    for (const Restraints::Angle& a : cc.rt.angles)
      tab.append_row({cc.name, a.id1.atom, a.id2.atom, a.id3.atom,
                      to_str(a.value), to_str(a.esd)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_tor.",
        {"comp_id", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
         "value_angle", "value_angle_esd", "period"});
    for (const Restraints::Torsion& a : cc.rt.torsions)
      tab.append_row({cc.name, a.label, a.id1.atom, a.id2.atom, a.id3.atom, a.id4.atom,
                      to_str(a.value), to_str(a.esd), std::to_string(a.period)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_chir.",
        {"comp_id", "id", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3",
         "volume_sign"});
    for (const Restraints::Chirality& a : cc.rt.chirs) {
      std::string label = "chir_" + std::to_string(tab.length() + 1);
      tab.append_row({cc.name, label, a.id_ctr.atom, a.id1.atom, a.id2.atom, a.id3.atom,
                      chirality_to_string(a.sign)});
    }
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_plane_atom.",
        {"comp_id", "plane_id", "atom_id", "dist_esd"});
    for (const Restraints::Plane& p : cc.rt.planes)
      for (const Restraints::AtomId& atom_id : p.ids)
        tab.append_row({cc.name, p.label, atom_id.atom, to_str(p.esd)});
  }
}

inline ChemComp make_chemcomp_with_restraints(const Residue& res) {
  ChemComp cc;
  cc.name = res.name;
  cc.type_or_group = "?";  // cf. is_ad_hoc()
  cc.group = ChemComp::Group::Null;
  // add atoms
  cc.atoms.reserve(res.atoms.size());
  for (const Atom& a : res.atoms) {
    Element el = a.element == El::X ? Element(El::N) : a.element;
    const std::string& chem_type = el.uname();
    cc.atoms.push_back(ChemComp::Atom{a.name, el, float(a.charge), chem_type});
  }
  // prepare pairs of atoms
  struct Pair {
    size_t n1, n2;
    double dist;
  };
  std::vector<Pair> pairs;
  // first heavy atoms only
  for (size_t i = 0; i != res.atoms.size(); ++i) {
    const Atom& at1 = res.atoms[i];
    if (at1.is_hydrogen())
      continue;
    float r1 = at1.element.covalent_r();
    for (size_t j = i+1; j != res.atoms.size(); ++j) {
      const Atom& at2 = res.atoms[j];
      if (at2.is_hydrogen())
        continue;
      double d2 = at1.pos.dist_sq(at2.pos);
      float r2 = at2.element.covalent_r();
      double dmax = std::max(2.0, 1.3 * std::max(r1, r2));
      if (d2 < sq(dmax))
        pairs.push_back(Pair{i, j, std::sqrt(d2)});
    }
  }
  // now each hydrogen with the nearest heavy atom
  for (size_t i = 0; i != res.atoms.size(); ++i) {
    const Atom& at1 = res.atoms[i];
    if (at1.is_hydrogen()) {
      size_t nearest = (size_t)-1;
      double min_d2 = sq(2.5);
      for (size_t j = 0; j != res.atoms.size(); ++j) {
        const Atom& at2 = res.atoms[j];
        if (!at2.is_hydrogen()) {
          double d2 = at1.pos.dist_sq(at2.pos);
          if (d2 < min_d2) {
            min_d2 = d2;
            nearest = j;
          }
        }
      }
      if (nearest != (size_t)-1) {
        pairs.push_back(Pair{nearest, i, std::sqrt(min_d2)});
      }
    }
  }

  // add bonds
  for (const Pair& p : pairs) {
    Restraints::Bond bond;
    bond.id1 = Restraints::AtomId{1, res.atoms[p.n1].name};
    bond.id2 = Restraints::AtomId{1, res.atoms[p.n2].name};
    bond.type = BondType::Unspec;
    bond.aromatic = false;
    double rounded_dist = 0.001 * std::round(1000 * p.dist);
    bond.value = bond.value_nucleus = rounded_dist;
    bond.esd = bond.esd_nucleus = 0.02;
    cc.rt.bonds.push_back(bond);
  }
  // add angles
  struct Triple {
    size_t n1, n2, n3;
  };
  std::vector<Triple> triples;
  for (size_t i = 0; i != pairs.size(); ++i)
    for (size_t j = i+1; j != pairs.size(); ++j) {
      if (pairs[i].n1 == pairs[j].n1)
        triples.push_back(Triple{pairs[i].n2, pairs[i].n1, pairs[j].n2});
      else if (pairs[i].n1 == pairs[j].n2)
        triples.push_back(Triple{pairs[i].n2, pairs[i].n1, pairs[j].n1});
      else if (pairs[i].n2 == pairs[j].n1)
        triples.push_back(Triple{pairs[i].n1, pairs[i].n2, pairs[j].n2});
      else if (pairs[i].n2 == pairs[j].n2)
        triples.push_back(Triple{pairs[i].n1, pairs[i].n2, pairs[j].n1});
    }
  for (const Triple& triple : triples) {
    Restraints::Angle angle;
    angle.id1 = Restraints::AtomId{1, res.atoms[triple.n1].name};
    angle.id2 = Restraints::AtomId{1, res.atoms[triple.n2].name};
    angle.id3 = Restraints::AtomId{1, res.atoms[triple.n3].name};
    double angle_rad = calculate_angle(res.atoms[triple.n1].pos,
                                       res.atoms[triple.n2].pos,
                                       res.atoms[triple.n3].pos);
    angle.value = 0.01 * std::round(100 * deg(angle_rad));
    angle.esd = 3.0;
    cc.rt.angles.push_back(angle);
  }
  return cc;
}

} // namespace gemmi
#endif
