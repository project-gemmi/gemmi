// Copyright 2018 Global Phasing Ltd.
//
// Monomer library - (Refmac) restraints dictionary,
// which is made of monomers (chemical components), links and modifications.

#ifndef GEMMI_MONLIB_HPP_
#define GEMMI_MONLIB_HPP_

#include <map>
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"       // for Element
#include "numb.hpp"       // for as_number
#include "fail.hpp"       // for fail, unreachable
#include "model.hpp"      // for Residue, Atom
#include "chemcomp.hpp"   // for ChemComp

namespace gemmi {

typedef cif::Document (*read_cif_func)(const std::string&);

inline void add_distinct_altlocs(const Residue& res, std::string& altlocs) {
  for (const Atom& atom : res.atoms)
    if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
      altlocs += atom.altloc;
}

inline bool atom_match_with_alias(const std::string& atom_id, const std::string& atom,
                                  const ChemComp::Aliasing* aliasing) {
  if (aliasing)
    if (const std::string* real_id = aliasing->name_from_alias(atom_id))
      return *real_id == atom;
  return atom_id == atom;
}

struct ChemLink {
  struct Side {
    using Group = ChemComp::Group;
    std::string comp;
    std::string mod;
    Group group = Group::Null;
    bool matches_group(Group res) const {
      if (group == Group::Null)
        return false;
      return res == group || (group == Group::Peptide && ChemComp::is_peptide_group(res))
                          || (group == Group::DnaRna && ChemComp::is_nucleotide_group(res));
    }
    int specificity() const {
      if (!comp.empty())
        return 3;
      return group == Group::PPeptide || group == Group::MPeptide ? 1 : 0;
    }
  };
  std::string id;
  std::string name;
  Side side1;
  Side side2;
  Restraints rt;
  cif::Block block;  // temporary, until we have ChemLink->Block function

  /// If multiple ChemLinks match a bond, the one with highest scores should be used.
  int calculate_score(const Residue& res1, const Residue* res2,
                      char alt, char alt2,
                      const ChemComp::Aliasing* aliasing1,
                      const ChemComp::Aliasing* aliasing2) const;
};

struct GEMMI_DLL ChemMod {
  struct AtomMod {
    int func;
    std::string old_id;
    std::string new_id;
    Element el;
    float charge;
    std::string chem_type;
  };

  std::string id;
  std::string name;
  std::string comp_id;
  std::string group_id;
  std::vector<AtomMod> atom_mods;
  Restraints rt;
  cif::Block block;  // temporary, until we have ChemMod->Block function

  void apply_to(ChemComp& chemcomp, ChemComp::Group alias_group) const;
};

struct EnerLib {
  enum class RadiusType {Vdw, Vdwh, Ion};
  struct Atom {
    Element element;
    char hb_type;
    double vdw_radius;
    double vdwh_radius;
    double ion_radius;
    int valency;
    int sp;
  };
  struct Bond {
    std::string atom_type_2;
    BondType type;
    double length;
    double value_esd;
  };

  EnerLib() {}
  void read(const cif::Document& doc) {
    cif::Block& block = const_cast<cif::Block&>(doc.blocks[0]);
    for (const auto& row : block.find("_lib_atom.",
                    {"type", "hb_type", "vdw_radius", "vdwh_radius",
                     "ion_radius", "element", "valency", "sp"}))
      atoms.emplace(row[0], Atom{Element(row[5]), row[1][0], cif::as_number(row[2]),
                                 cif::as_number(row[3]), cif::as_number(row[4]),
                                 cif::as_int(row[6], -1), cif::as_int(row[7], -1)});
    for (const auto& row : block.find("_lib_bond.",
                    {"atom_type_1", "atom_type_2", "type", "length", "value_esd"}))
      bonds.emplace(row.str(0), Bond{row.str(1), bond_type_from_string(row[2]),
                                 cif::as_number(row[3]), cif::as_number(row[4])});
  }
  std::map<std::string, Atom> atoms; // type->Atom
  std::multimap<std::string, Bond> bonds; // atom_type_1->Bond
};

struct GEMMI_DLL MonLib {
  std::string monomer_dir;
  std::string lib_version;
  EnerLib ener_lib;
  std::map<std::string, ChemComp> monomers;
  std::map<std::string, ChemLink> links;
  std::map<std::string, ChemMod> modifications;
  std::map<std::string, ChemComp::Group> cc_groups;

  const ChemLink* get_link(const std::string& link_id) const {
    auto link = links.find(link_id);
    return link != links.end() ? &link->second : nullptr;
  }
  const ChemMod* get_mod(const std::string& name) const {
    auto modif = modifications.find(name);
    return modif != modifications.end() ? &modif->second : nullptr;
  }

  // Returns the most specific link and a flag that is true
  // if the order is comp2-comp1 in the link definition.
  std::tuple<const ChemLink*, bool, const ChemComp::Aliasing*, const ChemComp::Aliasing*>
  match_link(const Residue& res1, const std::string& atom1, char alt1,
             const Residue& res2, const std::string& atom2, char alt2,
             double min_bond_sq=0) const {
    const ChemLink* best_link = nullptr;
    bool inverted = false;
    const ChemComp::Aliasing* aliasing1 = nullptr;
    const ChemComp::Aliasing* aliasing2 = nullptr;
    const ChemComp::Aliasing* aliasing1_final = nullptr;
    const ChemComp::Aliasing* aliasing2_final = nullptr;
    int best_score = -1000;
    for (auto& ml : links) {
      const ChemLink& link = ml.second;
      if (link.rt.bonds.empty() || starts_with(link.name, "auto-"))
        continue;
      // for now we don't have link definitions with >1 bonds
      const Restraints::Bond& bond = link.rt.bonds[0];
      if (sq(bond.value) < min_bond_sq)
        continue;
      if (link_side_matches_residue(link.side1, res1.name, &aliasing1) &&
          link_side_matches_residue(link.side2, res2.name, &aliasing2) &&
          atom_match_with_alias(bond.id1.atom, atom1, aliasing1) &&
          atom_match_with_alias(bond.id2.atom, atom2, aliasing2)) {
        int score = link.calculate_score(res1, &res2, alt1, alt2, aliasing1, aliasing2);
        if (score > best_score) {
          best_link = &link;
          best_score = score;
          aliasing1_final = aliasing1;
          aliasing2_final = aliasing2;
          inverted = false;
        }
      }
      if (link_side_matches_residue(link.side1, res2.name, &aliasing2) &&
          link_side_matches_residue(link.side2, res1.name, &aliasing1) &&
          atom_match_with_alias(bond.id1.atom, atom2, aliasing2) &&
          atom_match_with_alias(bond.id2.atom, atom1, aliasing1)) {
        int score = link.calculate_score(res2, &res1, alt2, alt1, aliasing2, aliasing1);
        if (score > best_score) {
          best_link = &link;
          best_score = score;
          aliasing1_final = aliasing1;
          aliasing2_final = aliasing2;
          inverted = true;
        }
      }
    }
    return std::make_tuple(best_link, inverted, aliasing1_final, aliasing2_final);
  }

  void add_monomer_if_present(const cif::Block& block) {
    if (block.has_tag("_chem_comp_atom.atom_id")) {
      ChemComp cc = make_chemcomp_from_block(block);
      if (cc.group == ChemComp::Group::Null) {
        auto it = cc_groups.find(cc.name);
        if (it != cc_groups.end())
          cc.group = it->second;
      }
      std::string name = cc.name;
      monomers.emplace(name, std::move(cc));
    }
  }

  bool link_side_matches_residue(const ChemLink::Side& side,
                                 const std::string& res_name,
                                 ChemComp::Aliasing const** aliasing) const {
    assert(aliasing);
    *aliasing = nullptr;
    if (!side.comp.empty())
      return side.comp == res_name;
    auto it = monomers.find(res_name);
    if (it != monomers.end()) {
      if (side.matches_group(it->second.group))
        return true;
      for (const ChemComp::Aliasing& a : it->second.aliases)
        if (side.matches_group(a.group)) {
          *aliasing = &a;
          return true;
        }
    }
    return false;
  }

  /// Returns path to the monomer cif file (the file may not exist).
  std::string path(const std::string& code) {
      return monomer_dir + relative_monomer_path(code);
  }

  static std::string relative_monomer_path(const std::string& code) {
    std::string path;
    if (!code.empty()) {
      path += lower(code[0]);
      path += '/';  // works also on Windows
      path += code;
      // On Windows several names are reserved (CON, PRN, AUX, ...), see
      // https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file
      // The workaround in CCP4 monomer library is to use CON_CON.cif, etc.
      if (code.size() == 3)
        switch (ialpha3_id(code.c_str())) {
          case ialpha3_id("AUX"):
          case ialpha3_id("COM"):
          case ialpha3_id("CON"):
          case ialpha3_id("LPT"):
          case ialpha3_id("PRN"):
            path += '_';
            path += code;
        }
      path += ".cif";
    }
    return path;
  }

  void read_monomer_doc(const cif::Document& doc);

  void read_monomer_cif(const std::string& path_, read_cif_func read_cif) {
    const cif::Document& doc = (*read_cif)(path_);
    if (!doc.blocks.empty() && doc.blocks[0].name == "lib")
      if (const std::string* ver = doc.blocks[0].find_value("_lib.version"))
        lib_version = *ver;
    read_monomer_doc(doc);
  }

  void set_monomer_dir(const std::string& monomer_dir_) {
    monomer_dir = monomer_dir_;
    if (monomer_dir.back() != '/' && monomer_dir.back() != '\\')
      monomer_dir += '/';
  }

  /// Read mon_lib_list.cif, ener_lib.cif and required monomers.
  /// Returns true if all requested monomers were added.
  bool read_monomer_lib(const std::string& monomer_dir_,
                        const std::vector<std::string>& resnames,
                        read_cif_func read_cif,
                        std::string* error=nullptr) {
    if (monomer_dir_.empty())
      fail("read_monomer_lib: monomer_dir not specified.");
    set_monomer_dir(monomer_dir_);

    read_monomer_cif(monomer_dir + "list/mon_lib_list.cif", read_cif);
    ener_lib.read((*read_cif)(monomer_dir + "ener_lib.cif"));

    bool ok = true;
    for (const std::string& name : resnames) {
      if (monomers.find(name) != monomers.end())
        continue;
      try {
        const cif::Document& doc = (*read_cif)(path(name));
        read_monomer_doc(doc);
      } catch (std::system_error& err) {
        if (error) {
          if (err.code().value() == ENOENT)
            cat_to(*error, "Monomer not in the library: ", name, ".\n");
          else
            cat_to(*error, "Failed to read ", name, ": ", err.what(), ".\n");
        }
        ok = false;
      } catch (std::runtime_error& err) {
        if (error)
          cat_to(*error, "Failed to read ", name, ": ", err.what(), ".\n");
        ok = false;
      }
    }
    return ok;
  }

  double find_ideal_distance(const const_CRA& cra1, const const_CRA& cra2) const;
  std::string update_old_atom_names(Structure& st) const;
};

// to be deprecated
inline MonLib read_monomer_lib(const std::string& monomer_dir,
                               const std::vector<std::string>& resnames,
                               read_cif_func read_cif,
                               const std::string& libin="",
                               bool ignore_missing=false) {
  MonLib monlib;
  if (!libin.empty())
    monlib.read_monomer_cif(libin, read_cif);
  std::string error;
  bool ok = monlib.read_monomer_lib(monomer_dir, resnames, read_cif, &error);
  if (!ignore_missing && !ok)
    fail(error + "Please create definitions for missing monomers.");
  return monlib;
}

} // namespace gemmi
#endif
