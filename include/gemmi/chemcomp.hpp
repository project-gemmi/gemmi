// Copyright 2018 Global Phasing Ltd.
//
// Chemical component represents a monomer from the Refmac monomer library,
// or from PDB CCD.

#ifndef GEMMI_CHEMCOMP_HPP_
#define GEMMI_CHEMCOMP_HPP_

#include <cassert>
#include <map>
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"  // for Element
#include "numb.hpp"  // for as_number
#include "util.hpp"  // for istarts_with
#include "model.hpp" // for Residue, Atom

namespace gemmi {

struct Restraints {
  struct AtomId {
    int comp;
    std::string atom;

    bool operator==(const AtomId& o) const {
      return comp == o.comp && atom == o.atom;
    }
    bool operator!=(const AtomId& o) const { return !operator==(o); }

    const Atom* get_from(const Residue& res, const Residue* res2) const {
      if (comp == 1 || res2 == nullptr)
        return res.find_atom(atom);
      else if (comp == 2)
        return res2->find_atom(atom);
      throw std::out_of_range("Unexpected component ID");
    }
  };

  struct Bond {
    enum Type { Unspec, Single, Double, Triple, Aromatic, Deloc, Metal };
    AtomId id1, id2;
    Type type;
    bool aromatic;
    double value;
    double esd;
  };

  struct Angle {
    AtomId id1, id2, id3;
    double value;
    double esd;
  };

  struct Torsion {
    std::string label;
    AtomId id1, id2, id3, id4;
    double value;
    double esd;
    int period;
  };

  struct Chirality {
    enum Type { Positive, Negative, Both };
    AtomId id_ctr, id1, id2, id3;
    Type chir;
  };

  struct Plane {
    std::vector<AtomId> ids;
    double esd = 0.0;
  };

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::map<std::string, Plane> planes;

  Bond& get_bond(const AtomId& a1, const AtomId& a2) {
    for (Bond& b : bonds)
      if ((b.id1 == a1 && b.id2 == a2) || (b.id1 == a2 && b.id2 == a1))
        return b;
    throw std::out_of_range("Bond restraint not found: " +
                            a1.atom + "-" + a2.atom);
  }
  const Bond& get_bond(const AtomId& a1, const AtomId& a2) const {
    return const_cast<Restraints*>(this)->get_bond(a1, a2);
  }

  Angle& get_angle(const AtomId& a1, const AtomId& a2, const AtomId& a3) {
    for (Angle& ang : angles)
      if (ang.id2 == a2 && ((ang.id1 == a1 && ang.id3 == a3) ||
                            (ang.id1 == a3 && ang.id3 == a1)))
        return ang;
    throw std::out_of_range("Angle restraint not found: " +
                            a1.atom + "-" + a2.atom + "-" + a3.atom);
  }
  const Angle& get_angle(const AtomId& a1, const AtomId& a2,
                         const AtomId& a3) const {
    return const_cast<Restraints*>(this)->get_angle(a1, a2, a3);
  }

  double chiral_abs_volume(const Restraints::Chirality& ch) const {
    double mult = get_bond(ch.id_ctr, ch.id1).value *
                  get_bond(ch.id_ctr, ch.id2).value *
                  get_bond(ch.id_ctr, ch.id3).value;
    double x = 1;
    double y = 2;
    for (double a : {get_angle(ch.id1, ch.id_ctr, ch.id2).value,
                     get_angle(ch.id2, ch.id_ctr, ch.id3).value,
                     get_angle(ch.id3, ch.id_ctr, ch.id1).value}) {
      constexpr double deg2rad = pi() / 180.0;
      double cosine = a == 90. ? 0. : std::cos(deg2rad * a);
      x -= cosine * cosine;
      y *= cosine;
    }
    return mult * std::sqrt(x + y);
  }
};

struct ChemComp {
  struct Atom {
    std::string id;
    Element el;
    std::string chem_type;
  };

  std::string name;
  std::string group;
  std::vector<Atom> atoms;
  Restraints rt;

  const Atom& get_atom(const std::string& atom_id) const {
    for (const Atom& a : atoms)
      if (a.id == atom_id)
        return a;
    throw std::out_of_range("Chemical componenent " + name + " has no atom "
                            + atom_id);
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
};

struct ChemLink {
  std::string id;
  std::string name;
  std::string comp[2];
  std::string mod[2];
  std::string group[2];
  Restraints rt;
};

struct ChemMod {
  std::string id;
  std::string name;
  std::string comp_id;
  std::string group_id;
};

inline Restraints::Bond::Type bond_type_from_string(const std::string& s) {
  if (istarts_with(s, "sing"))
    return Restraints::Bond::Single;
  if (istarts_with(s, "doub"))
    return Restraints::Bond::Double;
  if (istarts_with(s, "trip"))
    return Restraints::Bond::Triple;
  if (istarts_with(s, "arom"))
    return Restraints::Bond::Aromatic;
  if (istarts_with(s, "metal"))
    return Restraints::Bond::Metal;
  if (istarts_with(s, "delo") || s == "1.5")
    return Restraints::Bond::Deloc;
  if (cif::is_null(s))
    return Restraints::Bond::Unspec;
  throw std::out_of_range("Unexpected bond type: " + s);
}

inline std::string bond_type_to_string(Restraints::Bond::Type btype) {
  switch (btype) {
    case Restraints::Bond::Unspec: return ".";
    case Restraints::Bond::Single: return "single";
    case Restraints::Bond::Double: return "double";
    case Restraints::Bond::Triple: return "triple";
    case Restraints::Bond::Aromatic: return "aromatic";
    case Restraints::Bond::Deloc: return "deloc";
    case Restraints::Bond::Metal: return "metal";
  }
}

// it doesn't handle crossN types from the monomer library
inline Restraints::Chirality::Type chirality_from_string(const std::string& s) {
  switch (s[0] | 0x20) {
    case 'p': return Restraints::Chirality::Positive;
    case 'n': return Restraints::Chirality::Negative;
    case 'b': return Restraints::Chirality::Both;
    default: throw std::out_of_range("Unexpected chirality: " + s);
  }
}

inline ChemComp make_chemcomp_from_cif(const std::string& name,
                                       const cif::Document& doc) {
  ChemComp cc;
  cc.name = name;
  const cif::Block* block_ = doc.find_block("comp_" + name);
  if (!block_)
    block_ = doc.find_block(name);
  if (!block_)
    throw std::runtime_error("data_comp_" + name + " not in the cif file");
  cif::Block* block = const_cast<cif::Block*>(block_);
  cif::Block* aux_block = const_cast<cif::Block*>(doc.find_block("comp_list"));
  cif::Column group_col = block->find_values("_chem_comp.group");
  if (!group_col && aux_block)
    group_col = aux_block->find_values("_chem_comp.group");
  if (group_col)
    cc.group = group_col.str(0);
  for (auto row : block->find("_chem_comp_atom.",
                              {"atom_id", "type_symbol", "type_energy"}))
    cc.atoms.push_back({row.str(0), Element(row.str(1)), row.str(2)});
  for (auto row : block->find("_chem_comp_bond.",
                              {"atom_id_1", "atom_id_2",
                               "type",
                               "?aromatic",
                               "value_dist", "value_dist_esd"}))
    cc.rt.bonds.push_back({{1, row.str(0)}, {1, row.str(1)},
                           bond_type_from_string(row[2]),
                           row.has2(3) && row[3][0] == 'y',
                           cif::as_number(row[4]), cif::as_number(row[5])});
  for (auto row : block->find("_chem_comp_angle.",
                              {"atom_id_1", "atom_id_2", "atom_id_3",
                               "value_angle", "value_angle_esd"}))
    cc.rt.angles.push_back({{1, row.str(0)}, {1, row.str(1)}, {1, row.str(2)},
                            cif::as_number(row[3]), cif::as_number(row[4])});
  for (auto row : block->find("_chem_comp_tor.",
                              {"id",
                               "atom_id_1", "atom_id_2",
                               "atom_id_3", "atom_id_4",
                               "value_angle", "value_angle_esd",
                               "period"}))
    cc.rt.torsions.push_back({row.str(0),
                              {1, row.str(1)}, {1, row.str(2)},
                              {1, row.str(3)}, {1, row.str(4)},
                              cif::as_number(row[5]), cif::as_number(row[6]),
                              cif::as_int(row[7])});
  for (auto row : block->find("_chem_comp_chir.",
                              {"atom_id_centre",
                               "atom_id_1", "atom_id_2", "atom_id_3",
                               "volume_sign"}))
    cc.rt.chirs.push_back({{1, row.str(0)},
                           {1, row.str(1)}, {1, row.str(2)}, {1, row.str(3)},
                           chirality_from_string(row[4])});
  for (auto row : block->find("_chem_comp_plane_atom.",
                              {"plane_id", "atom_id" , "dist_esd"})) {
    Restraints::Plane& plane = cc.rt.planes[row.str(0)];
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[2]);
    plane.ids.push_back({1, row.str(1)});
  }
  return cc;
}

inline Restraints read_link_restraints(const cif::Block& block_) {
  auto read_aid = [](cif::Table::Row& row, int n) {
    return Restraints::AtomId{cif::as_int(row[n]), row.str(n+1)};
  };
  Restraints rt;
  cif::Block& block = const_cast<cif::Block&>(block_);
  for (auto row : block.find("_chem_link_bond.",
                             {"atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "type",
                              "value_dist", "value_dist_esd"}))
    rt.bonds.push_back({read_aid(row, 0), read_aid(row, 2),
                        bond_type_from_string(row[4]), false,
                        cif::as_number(row[5]), cif::as_number(row[6])});
  for (auto row : block.find("_chem_link_angle.",
                             {"atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "value_angle", "value_angle_esd"}))
    rt.angles.push_back({read_aid(row, 0), read_aid(row, 2), read_aid(row, 4),
                         cif::as_number(row[6]), cif::as_number(row[7])});
  for (auto row : block.find("_chem_link_tor.",
                             {"id",
                              "atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "atom_4_comp_id", "atom_id_4",
                              "value_angle", "value_angle_esd",
                              "period"}))
    rt.torsions.push_back({row.str(0),
                           read_aid(row, 1), read_aid(row, 3),
                           read_aid(row, 5), read_aid(row, 7),
                           cif::as_number(row[9]), cif::as_number(row[10]),
                           cif::as_int(row[11])});
  for (auto row : block.find("_chem_link_chir.",
                             {"atom_centre_comp_id", "atom_id_centre",
                              "atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "volume_sign"}))
    rt.chirs.push_back({read_aid(row, 0), read_aid(row, 2),
                        read_aid(row, 4), read_aid(row, 6),
                        chirality_from_string(row[8])});
  for (auto row : block.find("_chem_link_plane.",
                             {"plane_id", "atom_comp_id", "atom_id",
                              "dist_esd"})) {
    Restraints::Plane& plane = rt.planes[row.str(0)];
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[3]);
    plane.ids.push_back(read_aid(row, 1));
  }
  return rt;
}

inline std::map<std::string,ChemLink> read_chemlinks(cif::Document& doc) {
  std::map<std::string, gemmi::ChemLink> links;
  const cif::Block* list_block = doc.find_block("link_list");
  if (!list_block)
    throw std::runtime_error("data_link_list not in the cif file");
  for (auto row : const_cast<cif::Block*>(list_block)->find("_chem_link.",
                                 {"id", "name",
                                  "comp_id_1", "mod_id_1", "group_comp_1",
                                  "comp_id_2", "mod_id_2", "group_comp_2"})) {
    ChemLink link;
    link.id = row.str(0);
    link.name = row.str(1);
    for (int i : {0, 1}) {
      link.comp[i] = row.str(2 + i * 3);
      link.mod[i] = row.str(3 + i * 3);
      link.group[i] = row.str(4 + i * 3);
    }
    const cif::Block* block = doc.find_block("link_" + link.id);
    if (!block)
      throw std::runtime_error("inconsisted data_link_list");
    link.rt = read_link_restraints(*block);
    links.emplace(link.id, link);
  }
  return links;
}

inline std::map<std::string,ChemMod> read_chemmods(cif::Document& doc) {
  std::map<std::string, gemmi::ChemMod> mods;
  const cif::Block* list_block = doc.find_block("mod_list");
  if (!list_block)
    throw std::runtime_error("data_mod_list not in the cif file");
  for (auto row : const_cast<cif::Block*>(list_block)->find("_chem_mod.",
                                 {"id", "name", "comp_id", "group_id"})) {
    ChemMod mod;
    mod.id = row.str(0);
    mod.name = row.str(1);
    mod.comp_id = row.str(2);
    mod.group_id = row.str(3);
    const cif::Block* block = doc.find_block("mod_" + mod.id);
    if (!block)
      throw std::runtime_error("inconsisted data_mod_list");
    //TODO read actual modifications
    mods.emplace(mod.id, mod);
  }
  return mods;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
