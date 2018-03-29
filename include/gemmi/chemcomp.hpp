// Copyright 2018 Global Phasing Ltd.
//
// Chemical component represents a monomer from the Refmac monomer library,
// or from PDB CCD.

#ifndef GEMMI_CHEMCOMP_HPP_
#define GEMMI_CHEMCOMP_HPP_

#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"  // for Element
#include "numb.hpp"  // for as_number

namespace gemmi {

struct ChemComp {
  enum BondType { Single, Double, Triple, Aromatic, Deloc, Metal };
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
    std::string id1, id2, id3, id4;
    double value;
    double esd;
    int period;
  };
  struct Chirality {
    std::string id_ctr, id1, id2, id3;
    bool positive;
  };
  struct Plane {
    std::vector<std::string> ids;
    double esd;
  };

  std::string name;
  std::string group;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  const Atom& get_atom(const std::string& atom_id) const {
    for (const Atom& a : atoms)
      if (a.id == atom_id)
        return a;
    throw std::out_of_range("Chemical componenent " + name + " has no atom "
                            + atom_id);
  }
};

ChemComp::BondType bond_type_from_string(const std::string& s) {
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
  return cc;
}


} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
