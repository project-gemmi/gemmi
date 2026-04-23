// Copyright 2018 Global Phasing Ltd.
//
// Monomer library - (Refmac) restraints dictionary,
// which consists of monomers (chemical components), links, and modifications.

#ifndef GEMMI_MONLIB_HPP_
#define GEMMI_MONLIB_HPP_

#include <map>
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"       // for Element
#include "ener_lib.hpp"   // for EnerLib
#include "fail.hpp"       // for fail, unreachable
#include "model.hpp"      // for Residue, Atom
#include "chemcomp.hpp"   // for ChemComp
#include "logger.hpp"     // for Logger

namespace gemmi {

/// @brief Check if an atom ID matches a canonical atom name, resolving aliases.
/// @param atom_id Atom identifier to check (may be aliased)
/// @param atom Canonical atom name
/// @param aliasing Optional aliasing rules to resolve atom_id
/// @return true if atom_id matches atom (directly or via aliasing)
inline bool atom_match_with_alias(const std::string& atom_id, const std::string& atom,
                                  const ChemComp::Aliasing* aliasing) {
  if (aliasing)
    if (const std::string* real_id = aliasing->name_from_alias(atom_id))
      return *real_id == atom;
  return atom_id == atom;
}

/// @brief Chemical link definition (bond, angle, dihedral between residues).
struct GEMMI_DLL ChemLink {
  /// @brief Specification of one side of a chemical link.
  struct Side {
    using Group = ChemComp::Group;
    std::string comp;          ///< Specific chemical component name, or empty for group-based matching
    std::string mod;           ///< Chemical modification identifier
    Group group = Group::Null; ///< Group type for general matching (peptide, nucleotide, etc.)

    /// @brief Check if this side matches a given chemical group.
    /// @param res Group to test against
    /// @return true if the side specification matches the group
    bool matches_group(Group res) const {
      if (group == Group::Null)
        return false;
      return res == group || (group == Group::Peptide && ChemComp::is_peptide_group(res))
                          || (group == Group::DnaRna && ChemComp::is_nucleotide_group(res));
    }

    /// @brief Calculate specificity score for matching priority.
    /// @return Higher scores indicate more specific matches (specific component > group-based)
    int specificity() const {
      if (!comp.empty())
        return 3;
      return group == Group::PPeptide || group == Group::MPeptide ? 1 : 0;
    }
  };

  std::string id;    ///< Link identifier
  std::string name;  ///< Link name
  Side side1;        ///< First residue specification
  Side side2;        ///< Second residue specification
  Restraints rt;     ///< Restraints (bonds, angles, dihedrals, etc.)
  cif::Block block;  ///< Temporary CIF block storage

  /// @brief Calculate matching score for this link between two residues.
  /// If multiple ChemLinks match a bond, the one with highest score should be used.
  /// @param res1 First residue
  /// @param res2 Second residue (nullptr if not available)
  /// @param alt First residue alternate location indicator
  /// @param alt2 Second residue alternate location indicator
  /// @param aliasing1 Aliasing rules for first residue
  /// @param aliasing2 Aliasing rules for second residue
  /// @return Numeric score indicating match quality (higher is better)
  int calculate_score(const Residue& res1, const Residue* res2,
                      char alt, char alt2,
                      const ChemComp::Aliasing* aliasing1,
                      const ChemComp::Aliasing* aliasing2) const;
};

/// @brief Chemical modification (alteration to a chemical component).
struct GEMMI_DLL ChemMod {
  /// @brief Modification to a single atom.
  struct AtomMod {
    int func;                ///< Modification function code
    std::string old_id;      ///< Original atom identifier
    std::string new_id;      ///< New atom identifier
    Element el;              ///< New element
    float charge;            ///< New formal charge
    std::string chem_type;   ///< New chemical type
  };

  std::string id;                   ///< Modification identifier
  std::string name;                 ///< Modification name
  std::string comp_id;              ///< Target chemical component
  std::string group_id;             ///< Group identifier
  std::vector<AtomMod> atom_mods;   ///< Atom modifications to apply
  Restraints rt;                    ///< Modified restraints
  cif::Block block;                 ///< Temporary CIF block storage

  /// @brief Apply this modification to a chemical component.
  /// @param chemcomp Chemical component to modify (in-place)
  /// @param alias_group Optional group alias to apply
  void apply_to(ChemComp& chemcomp, ChemComp::Group alias_group) const;
};

/// @brief Monomer library with chemical components, links, and modifications.
/// Stores the (Refmac) restraints dictionary including monomers, chemical links,
/// and modifications, along with atomic energy parameters.
struct GEMMI_DLL MonLib {
  std::string monomer_dir;                          ///< Directory containing monomer CIF files
  std::map<std::string, ChemComp> monomers;         ///< Chemical components indexed by name
  std::map<std::string, ChemLink> links;            ///< Chemical links indexed by ID
  std::map<std::string, ChemMod> modifications;     ///< Chemical modifications indexed by name
  std::map<std::string, ChemComp::Group> cc_groups; ///< Component group assignments
  EnerLib ener_lib;                                 ///< Energy library with atomic properties

  /// @brief Find a chemical link by identifier.
  /// @param link_id Link identifier
  /// @return Pointer to ChemLink, or nullptr if not found
  const ChemLink* get_link(const std::string& link_id) const {
    auto link = links.find(link_id);
    return link != links.end() ? &link->second : nullptr;
  }

  /// @brief Find a chemical modification by name.
  /// @param name Modification name
  /// @return Pointer to ChemMod, or nullptr if not found
  const ChemMod* get_mod(const std::string& name) const {
    auto modif = modifications.find(name);
    return modif != modifications.end() ? &modif->second : nullptr;
  }

  /// @brief Find the most specific chemical link between two residues and atoms.
  /// Returns the most specific link and a flag indicating if the residue order
  /// is inverted (comp2-comp1) in the link definition.
  /// @param res1 First residue
  /// @param atom1 Atom name in first residue
  /// @param alt1 Alternate location indicator for first atom
  /// @param res2 Second residue
  /// @param atom2 Atom name in second residue
  /// @param alt2 Alternate location indicator for second atom
  /// @param min_bond_sq Minimum squared bond length to accept
  /// @return Tuple of (link, inverted_flag, aliasing1, aliasing2);
  ///         link is nullptr if no match found
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
    for (const auto& ml : links) {
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
        // NOLINTNEXTLINE(readability-suspicious-call-argument)
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

  /// @brief Add a chemical component from a CIF block if it contains atom definitions.
  /// @param block CIF block containing chemical component data
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

  /// @brief Check if a link side specification matches a residue.
  /// @param side Link side specification to test
  /// @param res_name Residue name
  /// @param aliasing Output parameter: aliasing rules if matched via alias, nullptr otherwise
  /// @return true if side matches res_name (exactly or via group/alias)
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

  /// @brief Returns path to the monomer CIF file (the file may not exist).
  /// @param code Chemical component code
  /// @return Full file path constructed from monomer_dir and code
  std::string path(const std::string& code) const {
    return monomer_dir + relative_monomer_path(code);
  }

  /// @brief Get relative file path for a monomer within a standard directory structure.
  /// @param code Chemical component code
  /// @return Relative file path (e.g., "m/monomers/m_code.cif")
  static std::string relative_monomer_path(const std::string& code);

  /// @brief Read monomer library data from a CIF document.
  /// @param doc CIF document containing chemical components, links, and/or modifications
  void read_monomer_doc(const cif::Document& doc);

  /// @brief Read monomer library data from a CIF file.
  /// @param path_ File path to read
  void read_monomer_cif(const std::string& path_);

  /// @brief Set the directory for monomer CIF files.
  /// @param monomer_dir_ Directory path (trailing slash is optional and auto-added)
  void set_monomer_dir(const std::string& monomer_dir_) {
    monomer_dir = monomer_dir_;
    if (!monomer_dir.empty() && monomer_dir.back() != '/' && monomer_dir.back() != '\\')
      monomer_dir += '/';
  }

  /// @brief Read mon_lib_list.cif, ener_lib.cif and required monomers.
  /// @param monomer_dir_ Directory containing monomer library files
  /// @param resnames List of chemical component names to load
  /// @param logger Logger for diagnostic messages
  /// @return true if all requested monomers were added
  bool read_monomer_lib(const std::string& monomer_dir_,
                        const std::vector<std::string>& resnames,
                        const Logger& logger);

  /// @brief Find ideal bond distance from library for two atoms.
  /// @param cra1 First atom (chain, residue, atom reference)
  /// @param cra2 Second atom (chain, residue, atom reference)
  /// @return Ideal bond distance, or 0 if not found
  double find_ideal_distance(const const_CRA& cra1, const const_CRA& cra2) const;

  /// @brief Update old atom names in structure using alias information.
  /// @param st Structure to update (modified in-place)
  /// @param logger Logger for diagnostic messages
  void update_old_atom_names(Structure& st, const Logger& logger) const;
};

/// @brief Free function wrapper to read monomer library.
/// @deprecated Use MonLib::read_monomer_lib() method instead.
/// @param monomer_dir Directory containing monomer library files
/// @param resnames List of chemical component names to load
/// @param libin Optional path to additional library CIF file
/// @param ignore_missing If true, silently ignore missing components; if false, throw exception
/// @return Populated MonLib instance
inline MonLib read_monomer_lib(const std::string& monomer_dir,
                               const std::vector<std::string>& resnames,
                               const std::string& libin="",
                               bool ignore_missing=false) {
  MonLib monlib;
  if (!libin.empty())
    monlib.read_monomer_cif(libin);
  std::string error;
  Logger logger;
  if (!ignore_missing)
    logger.callback = [&error](const std::string& s) { cat_to(error, s, '\n'); };
  bool ok = monlib.read_monomer_lib(monomer_dir, resnames, logger);
  if (!ignore_missing && !ok) {
    error += "Please create definitions for missing monomers.";
    fail(error);
  }
  return monlib;
}

} // namespace gemmi
#endif
