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

    bool operator==(const std::string& name) const { return atom == name; }
    bool operator!=(const std::string& name) const { return atom != name; }

    const Atom* get_from(const Residue& res1, const Residue* res2,
                         char altloc) const {
      const Residue* residue;
      if (comp == 1 || res2 == nullptr)
        residue = &res1;
      else if (comp == 2)
        residue = res2;
      else
        throw std::out_of_range("Unexpected component ID");
      const Atom* ret = residue->find_atom(atom, altloc);
      return ret && ret->flag != 'M' ? ret : nullptr;
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
    std::string label;
    std::vector<AtomId> ids;
    double esd;
  };

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  bool empty() const {
    return bonds.empty() && angles.empty() && torsions.empty() &&
           chirs.empty() && planes.empty();
  }

  template<typename T>
  std::vector<Bond>::iterator find_bond(const T& a1, const T& a2) {
    return std::find_if(bonds.begin(), bonds.end(), [&](const Bond& b) {
        return (b.id1 == a1 && b.id2 == a2) || (b.id1 == a2 && b.id2 == a1);
    });
  }
  const Bond& get_bond(const AtomId& a1, const AtomId& a2) const {
    auto it = const_cast<Restraints*>(this)->find_bond(a1, a2);
    if (it == bonds.end())
      fail("Bond restraint not found: " + a1.atom + "-" + a2.atom);
    return *it;
  }

  template<typename T>
  std::vector<Angle>::iterator find_angle(const T& a, const T& b, const T& c) {
    return std::find_if(angles.begin(), angles.end(), [&](const Angle& ang) {
        return ang.id2 == b && ((ang.id1 == a && ang.id3 == c) ||
                                (ang.id1 == c && ang.id3 == a));
    });
  }
  const Angle& get_angle(const AtomId& a, const AtomId& b, const AtomId& c)
                                                                        const {
    auto it = const_cast<Restraints*>(this)->find_angle(a, b, c);
    if (it == angles.end())
      fail("Angle restraint not found: " + a.atom + "-" + b.atom + "-"+c.atom);
    return *it;
  }

  template<typename T>
  std::vector<Torsion>::iterator find_torsion(const T& a, const T& b,
                                              const T& c, const T& d) {
    return std::find_if(torsions.begin(), torsions.end(),
                        [&](const Torsion& t) {
        return (t.id1 == a && t.id2 == b && t.id3 == c && t.id4 == d) ||
               (t.id1 == d && t.id2 == c && t.id3 == b && t.id4 == a);
    });
  }

  template<typename T>
  std::vector<Chirality>::iterator find_chir(const T& ctr, const T& a,
                                             const T& b, const T& c) {
    return std::find_if(chirs.begin(), chirs.end(), [&](const Chirality& t) {
        return t.id_ctr == ctr && ((t.id1 == a && t.id2 == b && t.id3 == c) ||
                                   (t.id1 == b && t.id2 == c && t.id3 == a) ||
                                   (t.id1 == c && t.id2 == a && t.id3 == b));
    });
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

  std::vector<Plane>::iterator get_plane(const std::string& label) {
    return std::find_if(planes.begin(), planes.end(),
                        [&label](const Plane& p) { return p.label == label; });
  }

  Plane& get_or_add_plane(const std::string& label) {
    std::vector<Plane>::iterator it = get_plane(label);
    if (it != planes.end())
      return *it;
    planes.push_back(Plane{label, {}, 0.0});
    return planes.back();
  }
};

struct ChemComp {
  struct Atom {
    std::string id;
    Element el;
    std::string chem_type;

    gemmi::Atom to_full_atom() const {
      gemmi::Atom atom;
      atom.name = id;
      atom.flag = 'M'; // refmac convention for "modelled" missing atoms
      atom.occ = 0.0f;
      atom.b_iso = 0.0f;
      atom.element = el;
      return atom;
    }
  };

  std::string name;
  std::string group;
  std::vector<Atom> atoms;
  Restraints rt;

  std::vector<Atom>::iterator find_atom(const std::string& atom_id) {
    return std::find_if(atoms.begin(), atoms.end(),
                        [&](const Atom& a) { return a.id == atom_id; });
  }
  std::vector<Atom>::const_iterator find_atom(const std::string& atom_id) const{
    return const_cast<ChemComp*>(this)->find_atom(atom_id);
  }

  const Atom& get_atom(const std::string& atom_id) const {
    auto it = find_atom(atom_id);
    if (it == atoms.end())
      fail("Chemical componenent " + name + " has no atom " + atom_id);
    return *it;
  }
};

struct ChemLink {
  std::string id;
  std::string name;
  std::string comp[2];
  std::string mod[2];
  std::string group[2];
  Restraints rt;

  bool matches(const ChemLink& other) const {
    if (rt.bonds.empty() || other.rt.bonds.empty())
      return false;
    const Restraints::Bond& bond1 = rt.bonds[0];
    const Restraints::Bond& bond2 = other.rt.bonds[0];
    // TODO: handle all combinations?
    if (comp[0] == other.comp[0] && comp[1] == other.comp[1]) {
      return bond1.id1 == bond2.id1 && bond1.id2 == bond2.id2;
    } else if (comp[0] == other.comp[1] && comp[1] == other.comp[0]) {
      return bond1.id1.atom == bond2.id2.atom &&
             bond1.id2.atom == bond2.id1.atom;
    }
    return false;
  }
};

struct ChemMod {
  struct AtomMod {
    int func;
    std::string old_id;
    std::string new_id;
    Element el;
    std::string chem_type;
  };

  std::string id;
  std::string name;
  std::string comp_id;
  std::string group_id;
  std::vector<AtomMod> atom_mods;
  Restraints rt;

  void apply_to(ChemComp& cc) const;
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
    if (row[4][0] != 'c') // ignore crossN
      cc.rt.chirs.push_back({{1, row.str(0)},
                             {1, row.str(1)}, {1, row.str(2)}, {1, row.str(3)},
                             chirality_from_string(row[4])});
  for (auto row : block->find("_chem_comp_plane_atom.",
                              {"plane_id", "atom_id" , "dist_esd"})) {
    Restraints::Plane& plane = cc.rt.get_or_add_plane(row.str(0));
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
    if (row[4][0] != 'c') // ignore crossN
      rt.chirs.push_back({read_aid(row, 0), read_aid(row, 2),
                          read_aid(row, 4), read_aid(row, 6),
                          chirality_from_string(row[8])});
  for (auto row : block.find("_chem_link_plane.",
                             {"plane_id", "atom_comp_id", "atom_id",
                              "dist_esd"})) {
    Restraints::Plane& plane = rt.get_or_add_plane(row.str(0));
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
    fail("data_link_list not in the cif file");
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

// Helper function. str is one of "add", "delete", "change".
inline int chem_mod_type(const std::string& str) {
  char c = str[0] | 0x20;
  if (c != 'a' && c != 'd' && c != 'c')
    fail("Unexpected value of _chem_mod_*.function: " + str);
  return c;
}

inline Restraints read_restraint_modifications(const cif::Block& block_) {
  Restraints rt;
  cif::Block& block = const_cast<cif::Block&>(block_);
  for (auto row : block.find("_chem_mod_bond.",
                             {"function", "atom_id_1", "atom_id_2",
                              "new_type",
                              "new_value_dist", "new_value_dist_esd"}))
    rt.bonds.push_back({{chem_mod_type(row[0]), row.str(1)}, {1, row.str(2)},
                        bond_type_from_string(row[3]), false,
                        cif::as_number(row[4]), cif::as_number(row[5])});
  for (auto row : block.find("_chem_mod_angle.",
                             {"function", "atom_id_1",
                              "atom_id_2", "atom_id_3",
                              "new_value_angle", "new_value_angle_esd"}))
    rt.angles.push_back({{chem_mod_type(row[0]), row.str(1)},
                         {1, row.str(2)}, {1, row.str(3)},
                         cif::as_number(row[4]), cif::as_number(row[5])});
  for (auto row : block.find("_chem_mod_tor.",
                              {"function", "id", "atom_id_1",
                               "atom_id_2", "atom_id_3", "atom_id_4",
                               "new_value_angle", "new_value_angle_esd",
                               "new_period"}))
    rt.torsions.push_back({row.str(1), {chem_mod_type(row[0]), row.str(2)},
                           {1, row.str(3)}, {1, row.str(4)}, {1, row.str(5)},
                           cif::as_number(row[6]), cif::as_number(row[7]),
                           cif::as_int(row[8])});
  for (auto row : block.find("_chem_mod_chir.",
                             {"function", "atom_id_centre", "atom_id_1",
                              "atom_id_2", "atom_id_3",
                              "new_volume_sign"}))
    rt.chirs.push_back({{1, row.str(1)}, {chem_mod_type(row[0]), row.str(2)},
                        {1, row.str(3)}, {1, row.str(4)},
                        chirality_from_string(row[5])});
  for (auto row : block.find("_chem_mod_plane_atom.",
                             {"function", "plane_id", "atom_id" ,
                              "new_dist_esd"})) {
    Restraints::Plane& plane = rt.get_or_add_plane(row.str(1));
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[3]);
    plane.ids.push_back({chem_mod_type(row[0]), row.str(2)});
  }
  return rt;
}

inline std::map<std::string, ChemMod> read_chemmods(cif::Document& doc) {
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
    for (auto ra : const_cast<cif::Block*>(block)->find("_chem_mod_atom.",
                                {"function", "atom_id", "new_atom_id",
                                 "new_type_symbol", "new_type_energy"}))
      mod.atom_mods.push_back({chem_mod_type(ra[0]), ra.str(1), ra.str(2),
                               Element(ra.str(3)), ra.str(4)});
    mod.rt = read_restraint_modifications(*block);
    mods.emplace(mod.id, mod);
  }
  return mods;
}

namespace impl {
template <typename T>
T& add_or_set(std::vector<T>& items, typename std::vector<T>::iterator it,
              const T& x) {
  if (it == items.end()) {
    items.push_back(x);
    return items.back();
  }
  *it = x;
  return *it;
}
} // namespace impl

inline void ChemMod::apply_to(ChemComp& chemcomp) const {
  // _chem_mod_atom
  for (const AtomMod& mod : atom_mods) {
    auto it = chemcomp.find_atom(mod.old_id);
    switch (mod.func) {
      case 'a':
        if (chemcomp.find_atom(mod.new_id) == chemcomp.atoms.end())
          chemcomp.atoms.push_back({mod.new_id, mod.el, mod.chem_type});
        break;
      case 'd':
        if (it != chemcomp.atoms.end()) {
          chemcomp.atoms.erase(it);
          // delete restraints containing mod.old_id
          const std::string& old = mod.old_id;
          vector_remove_if(chemcomp.rt.bonds, [&](const Restraints::Bond& b) {
              return b.id1 == old || b.id2 == old;
          });
          vector_remove_if(chemcomp.rt.angles, [&](const Restraints::Angle& a) {
              return a.id1 == old || a.id2 == old || a.id3 == old;
          });
          vector_remove_if(chemcomp.rt.torsions,
              [&](const Restraints::Torsion& t) {
                return t.id1 == old || t.id2 == old || t.id3 == old ||
                       t.id4 == old;
          });
          vector_remove_if(chemcomp.rt.chirs,
              [&](const Restraints::Chirality& c) {
                return c.id_ctr == old || c.id1 == old || c.id2 == old ||
                       c.id3 == old;
          });
          for (Restraints::Plane& plane : chemcomp.rt.planes)
            vector_remove_if(plane.ids, [&](const Restraints::AtomId& a) {
                return a.atom == old;
            });
        }
        break;
      case 'c':
        if (it != chemcomp.atoms.end()) {
          if (!mod.new_id.empty())
            it->id = mod.new_id;
          if (mod.el != El::X)
            it->el = mod.el;
          if (!mod.chem_type.empty())
            it->chem_type = mod.chem_type;
        }
        break;
    }
  }

  // _chem_mod_bond
  for (const Restraints::Bond& mod : rt.bonds) {
    auto it = chemcomp.rt.find_bond(mod.id1.atom, mod.id2.atom);
    switch (mod.id1.comp) {
      case 'a':
        impl::add_or_set(chemcomp.rt.bonds, it, mod).id1.comp = 1;
        break;
      case 'd':
        if (it != chemcomp.rt.bonds.end())
          chemcomp.rt.bonds.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.bonds.end()) {
          if (mod.type != Restraints::Bond::Unspec)
            it->type = mod.type;
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
        }
        break;
    }
  }

  // _chem_mod_angle
  for (const Restraints::Angle& mod : rt.angles) {
    auto it = chemcomp.rt.find_angle(mod.id1.atom, mod.id2.atom, mod.id3.atom);
    switch (mod.id1.comp) {
      case 'a':
        impl::add_or_set(chemcomp.rt.angles, it, mod).id1.comp = 1;
        break;
      case 'd':
        if (it != chemcomp.rt.angles.end())
          chemcomp.rt.angles.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.angles.end()) {
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
        }
        break;
    }
  }

  // _chem_mod_tor
  for (const Restraints::Torsion& mod : rt.torsions) {
    auto it = chemcomp.rt.find_torsion(mod.id1.atom, mod.id2.atom,
                                       mod.id3.atom, mod.id4.atom);
    switch (mod.id1.comp) {
      case 'a':
        impl::add_or_set(chemcomp.rt.torsions, it, mod).id1.comp = 1;
        break;
      case 'd':
        if (it != chemcomp.rt.torsions.end())
          chemcomp.rt.torsions.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.torsions.end()) {
          if (!mod.label.empty())
            it->label = mod.label;
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
          if (mod.period != -1)
            it->period = mod.period;
        }
        break;
    }
  }

  // _chem_mod_chir
  for (const Restraints::Chirality& mod : rt.chirs) {
    auto it = chemcomp.rt.find_chir(mod.id_ctr.atom, mod.id1.atom,
                                    mod.id2.atom, mod.id3.atom);
    switch (mod.id1.comp) {
      case 'a':
        impl::add_or_set(chemcomp.rt.chirs, it, mod).id1.comp = 1;
        break;
      case 'd':
        if (it != chemcomp.rt.chirs.end())
          chemcomp.rt.chirs.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.chirs.end())
          it->chir = mod.chir;
        break;
    }
  }

  // _chem_mod_plane_atom
  for (const Restraints::Plane& mod : rt.planes)
    for (const Restraints::AtomId& atom_id : mod.ids)
      if (atom_id.comp == 'a') {
        Restraints::Plane& plane = chemcomp.rt.get_or_add_plane(mod.label);
        if (plane.esd == 0.0 && !std::isnan(mod.esd))
          plane.esd = mod.esd;
        auto it = std::find(plane.ids.begin(), plane.ids.end(), atom_id.atom);
        if (it == plane.ids.end())
          plane.ids.push_back({1, atom_id.atom});
      } else if (atom_id.comp == 'd') {
        auto it = chemcomp.rt.get_plane(mod.label);
        if (it != chemcomp.rt.planes.end()) {
          auto item = std::find(it->ids.begin(), it->ids.end(), atom_id.atom);
          if (item != it->ids.end())
            it->ids.erase(item);
        }
      }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
