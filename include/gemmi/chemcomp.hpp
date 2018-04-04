// Copyright 2018 Global Phasing Ltd.
//
// Chemical component represents a monomer from the Refmac monomer library,
// or from PDB CCD.

#ifndef GEMMI_CHEMCOMP_HPP_
#define GEMMI_CHEMCOMP_HPP_

#include <cassert>
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"  // for Element
#include "numb.hpp"  // for as_number

namespace gemmi {

struct ChemComp {
  enum BondType { Single, Double, Triple, Aromatic, Deloc, Metal };
  enum ChiralityType { Positive, Negative, Both };
  struct Atom {
    std::string id;
    Element el;
    std::string chem_type;
  };
  struct Bond {
    std::string id1, id2;
    BondType type;
    bool aromatic;
    double value;
    double esd;
  };
  struct Angle {
    std::string id1, id2, id3;
    double value;
    double esd;
  };
  struct Torsion {
    std::string label;
    std::string id1, id2, id3, id4;
    double value;
    double esd;
    int period;
  };
  struct Chirality {
    std::string id_ctr, id1, id2, id3;
    ChiralityType chir;
  };
  struct Plane {
    std::vector<std::string> ids;
    double esd = 0.02;
  };

  std::string name;
  std::string group;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::map<std::string, Plane> planes;

  const Atom& get_atom(const std::string& atom_id) const {
    for (const Atom& a : atoms)
      if (a.id == atom_id)
        return a;
    throw std::out_of_range("Chemical componenent " + name + " has no atom "
                            + atom_id);
  }

  Bond& get_bond(const std::string& a1, const std::string& a2) {
    for (Bond& b : bonds)
      if ((b.id1 == a1 && b.id2 == a2) || (b.id1 == a2 && b.id2 == a1))
        return b;
    throw std::out_of_range("Chemical componenent " + name + " has no bond "
                            + a1 + "-" + a2);
  }
  const Bond& get_bond(const std::string& a1, const std::string& a2) const {
    return const_cast<ChemComp*>(this)->get_bond(a1, a2);
  }

  Angle& get_angle(const std::string& atom1, const std::string& atom2,
                   const std::string& atom3) {
    for (Angle& a : angles)
      if (a.id2 == atom2 && ((a.id1 == atom1 && a.id3 == atom3) ||
                             (a.id1 == atom3 && a.id3 == atom1)))
        return a;
    throw std::out_of_range("Chemical componenent " + name + " has no angle "
                            + atom1 + "-" + atom2 + "-" + atom3);
  }
  const Angle& get_angle(const std::string& a1, const std::string& a2,
                         const std::string& a3) const {
    return const_cast<ChemComp*>(this)->get_angle(a1, a2, a3);
  }

  template<typename T>
  void reorder_atoms(std::vector<T>& alist) const {
    for (const T& a : alist)
      get_atom(a.name); // check that all atoms exists in _chem_comp_atom
    std::vector<T> ordered;
    ordered.reserve(alist.size());
    for (const Atom& cca : atoms)
      for (const T& a : alist)
        if (a.name == cca.id)
          ordered.push_back(a);
    assert(alist.size() == ordered.size());
    alist.swap(ordered);
  }

  double chiral_abs_volume(const Chirality& ch) const {
    double mult = get_bond(ch.id_ctr, ch.id1).value *
                  get_bond(ch.id_ctr, ch.id2).value *
                  get_bond(ch.id_ctr, ch.id3).value;
    double x = 1;
    double y = 2;
    for (double a : {get_angle(ch.id1, ch.id_ctr, ch.id2).value,
                     get_angle(ch.id2, ch.id_ctr, ch.id3).value,
                     get_angle(ch.id3, ch.id_ctr, ch.id1).value}) {
      constexpr double deg2rad = 3.1415926535897932384626433832795029 / 180.0;
      double cosine = a == 90. ? 0. : std::cos(deg2rad * a);
      x -= cosine * cosine;
      y *= cosine;
    }
    return mult * std::sqrt(x + y);
  }
};

inline double calculate_chiral_volume(const Position& actr, const Position& a1,
                                      const Position& a2, const Position& a3) {
  return (a1 - actr).dot((a2 - actr).cross(a3 - actr));
}

inline ChemComp::BondType bond_type_from_string(const std::string& s) {
  if (istarts_with(s, "sing"))
    return ChemComp::Single;
  if (istarts_with(s, "doub"))
    return ChemComp::Double;
  if (istarts_with(s, "trip"))
    return ChemComp::Triple;
  if (istarts_with(s, "arom"))
    return ChemComp::Aromatic;
  if (istarts_with(s, "delo"))
    return ChemComp::Deloc;
  if (istarts_with(s, "metal"))
    return ChemComp::Metal;
  if (s == "1.5")
    return ChemComp::Deloc;
  throw std::out_of_range("Unexpected bond type: " + s);
}

// it doesn't handle crossN types from the monomer library
inline ChemComp::ChiralityType chirality_from_string(const std::string& s) {
  switch (s[0] | 0x20) {
    case 'p': return ChemComp::Positive;
    case 'n': return ChemComp::Negative;
    case 'b': return ChemComp::Both;
    default: throw std::out_of_range("Unexpected chirality: " + s);
  }
}

inline
ChemComp make_chemcomp_from_cif(const std::string& name, cif::Document doc) {
  ChemComp cc;
  cc.name = name;
  cif::Block* block = doc.find_block("comp_" + name);
  if (!block)
    block = doc.find_block(name);
  if (!block)
    throw std::runtime_error("data_comp_" + name + " not in the cif file");
  cif::Block* aux_block = doc.find_block("comp_list");
  cif::Column group_col = block->find_values("_chem_comp.group");
  if (!group_col && aux_block)
    group_col = aux_block->find_values("_chem_comp.group");
  if (group_col)
    cc.group = group_col.str(0);
  for (auto row : block->find("_chem_comp_atom.",
                              {"atom_id", "type_symbol", "type_energy"}))
    cc.atoms.emplace_back(ChemComp::Atom{row.str(0),
                                         Element(row.str(1)),
                                         row.str(2)});
  for (auto row : block->find("_chem_comp_bond.",
                              {"atom_id_1", "atom_id_2",
                               "type", "?aromatic",
                               "value_dist", "value_dist_esd"}))
    cc.bonds.emplace_back(ChemComp::Bond{row.str(0),
                                         row.str(1),
                                         bond_type_from_string(row[2]),
                                         row.has2(3) && row[3] == "y",
                                         cif::as_number(row[4]),
                                         cif::as_number(row[5])});
  for (auto row : block->find("_chem_comp_angle.",
                              {"atom_id_1", "atom_id_2", "atom_id_3",
                               "value_angle", "value_angle_esd"}))
    cc.angles.emplace_back(ChemComp::Angle{row.str(0), row.str(1), row.str(2),
                                           cif::as_number(row[3]),
                                           cif::as_number(row[4])});
  for (auto row : block->find("_chem_comp_tor.",
                              {"id", "atom_id_1", "atom_id_2",
                                     "atom_id_3", "atom_id_4",
                               "value_angle", "value_angle_esd", "period"}))
    cc.torsions.emplace_back(ChemComp::Torsion{
        row.str(0), row.str(1), row.str(2), row.str(3), row.str(4),
        cif::as_number(row[5]), cif::as_number(row[6]), cif::as_int(row[7])});
  for (auto row : block->find("_chem_comp_chir.",
                              {"atom_id_centre", "atom_id_1",
                               "atom_id_2", "atom_id_3", "volume_sign"}))
    cc.chirs.emplace_back(ChemComp::Chirality{
        row.str(0), row.str(1), row.str(2), row.str(3),
        chirality_from_string(row[4])});
  for (auto row : block->find("_chem_comp_plane_atom.",
                              {"plane_id", "atom_id" /*, "dist_esd"*/}))
    // at the moment dist_esd is ignored by Refmac and it is assumed 0.02
    cc.planes[row.str(0)].ids.emplace_back(row.str(1));
  return cc;
}


} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
