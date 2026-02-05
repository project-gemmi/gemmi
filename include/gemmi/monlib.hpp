//! @file
//! @brief Monomer library - (Refmac) restraints dictionary.
//!
//! The monomer library consists of monomers (chemical components), links, and modifications.
//! Provides classes for managing chemical component definitions, covalent links between
//! residues, chemical modifications, and energy library parameters used in refinement.

// Copyright 2018 Global Phasing Ltd.

#ifndef GEMMI_MONLIB_HPP_
#define GEMMI_MONLIB_HPP_

#include <map>
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"       // for Element
#include "fail.hpp"       // for fail, unreachable
#include "model.hpp"      // for Residue, Atom
#include "chemcomp.hpp"   // for ChemComp
#include "logger.hpp"     // for Logger

namespace gemmi {

//! @brief Check if atom name matches, considering aliases.
//! @param atom_id Atom identifier from link/restraint definition
//! @param atom Atom name from structure
//! @param aliasing Optional aliasing table
//! @return True if names match (possibly through aliasing)
//!
//! Helper function for matching atom names in link definitions when aliases are present.
inline bool atom_match_with_alias(const std::string& atom_id, const std::string& atom,
                                  const ChemComp::Aliasing* aliasing) {
  if (aliasing)
    if (const std::string* real_id = aliasing->name_from_alias(atom_id))
      return *real_id == atom;
  return atom_id == atom;
}

//! @brief Chemical link definition between residues.
//!
//! Defines covalent bonds and associated restraints between two residues or
//! chemical components. Links can specify particular residue types or match
//! groups (peptides, nucleotides, etc.). Used for peptide bonds, disulfide
//! bridges, glycosidic bonds, and other inter-residue connections.
struct GEMMI_DLL ChemLink {
  //! @brief One side of a chemical link.
  //!
  //! Specifies which residue type(s) can participate on this side of the link.
  //! Can match specific component names, modifications, or chemical groups.
  struct Side {
    using Group = ChemComp::Group;
    std::string comp;  //!< Specific component name (e.g., "ALA", "GLY")
    std::string mod;   //!< Modification name
    Group group = Group::Null;  //!< Chemical group (peptide, nucleotide, etc.)

    //! @brief Check if this side matches a residue group.
    //! @param res Residue group to check
    //! @return True if the residue group matches this side's criteria
    bool matches_group(Group res) const {
      if (group == Group::Null)
        return false;
      return res == group || (group == Group::Peptide && ChemComp::is_peptide_group(res))
                          || (group == Group::DnaRna && ChemComp::is_nucleotide_group(res));
    }
    //! @brief Calculate specificity level of this side.
    //! @return Specificity score (3 for specific comp, 1 for L/D-peptide, 0 for generic group)
    //!
    //! Higher specificity indicates more restrictive matching criteria.
    int specificity() const {
      if (!comp.empty())
        return 3;
      return group == Group::PPeptide || group == Group::MPeptide ? 1 : 0;
    }
  };

  std::string id;    //!< Link identifier (e.g., "TRANS", "CIS", "SS")
  std::string name;  //!< Descriptive name
  Side side1;        //!< First side of the link
  Side side2;        //!< Second side of the link
  Restraints rt;     //!< Geometric restraints for the link
  cif::Block block;  //!< Temporary: until we have ChemLink->Block function

  //! @brief Calculate matching score for this link with given residues.
  //! @param res1 First residue
  //! @param res2 Second residue (may be nullptr)
  //! @param alt Alternate location for first residue
  //! @param alt2 Alternate location for second residue
  //! @param aliasing1 Aliasing for first residue
  //! @param aliasing2 Aliasing for second residue
  //! @return Score (higher is better match)
  //!
  //! If multiple ChemLinks match a bond, the one with highest score should be used.
  int calculate_score(const Residue& res1, const Residue* res2,
                      char alt, char alt2,
                      const ChemComp::Aliasing* aliasing1,
                      const ChemComp::Aliasing* aliasing2) const;
};

//! @brief Chemical modification definition.
//!
//! Describes how to modify a chemical component by adding, deleting, or changing
//! atoms. Used for post-translational modifications, methylation, phosphorylation,
//! etc. Contains atom modifications and updated restraints.
struct GEMMI_DLL ChemMod {
  //! @brief Single atom modification within a ChemMod.
  //!
  //! Specifies an operation to perform on an atom (add, delete, change).
  struct AtomMod {
    int func;                 //!< Function code (1=add, 2=delete, 3=change, etc.)
    std::string old_id;       //!< Original atom name
    std::string new_id;       //!< New atom name (for rename)
    Element el;               //!< Element type
    float charge;             //!< Atomic charge
    std::string chem_type;    //!< Chemical type
  };

  std::string id;                    //!< Modification identifier
  std::string name;                  //!< Descriptive name
  std::string comp_id;               //!< Component to modify
  std::string group_id;              //!< Chemical group
  std::vector<AtomMod> atom_mods;    //!< Atom modifications to apply
  Restraints rt;                     //!< Updated restraints after modification
  cif::Block block;                  //!< Temporary: until we have ChemMod->Block function

  //! @brief Apply this modification to a chemical component.
  //! @param chemcomp Chemical component to modify (modified in place)
  //! @param alias_group Alias group for the component
  //!
  //! Performs all atom modifications and updates restraints in the given ChemComp.
  void apply_to(ChemComp& chemcomp, ChemComp::Group alias_group) const;
};

//! @brief Energy library with atomic parameters for refinement.
//!
//! Contains atom types with van der Waals radii, hydrogen bonding types,
//! and ideal bond lengths. Read from ener_lib.cif in the monomer library.
struct EnerLib {
  //! @brief Type of atomic radius to use.
  enum class RadiusType {Vdw, Vdwh, Ion};

  //! @brief Atomic parameters for energy calculations.
  struct Atom {
    Element element;       //!< Chemical element
    char hb_type;          //!< Hydrogen bonding type
    double vdw_radius;     //!< Van der Waals radius
    double vdwh_radius;    //!< VDW radius including hydrogen
    double ion_radius;     //!< Ionic radius
    int valency;           //!< Valency
    int sp;                //!< Hybridization (sp, sp2, sp3)
  };

  //! @brief Ideal bond parameters.
  struct Bond {
    std::string atom_type_2;  //!< Second atom type
    BondType type;            //!< Bond type (single, double, etc.)
    double length;            //!< Ideal bond length in Angstroms
    double value_esd;         //!< Estimated standard deviation
  };

  EnerLib() {}

  //! @brief Read energy library from CIF document.
  //! @param doc CIF document containing energy library data
  void read(const cif::Document& doc);

  std::map<std::string, Atom> atoms;       //!< Atom types: type -> Atom
  std::multimap<std::string, Bond> bonds;  //!< Bonds: atom_type_1 -> Bond
};

//! @brief Monomer library (Refmac restraints dictionary).
//!
//! Central repository for chemical component definitions, links, modifications,
//! and energy parameters. Provides methods to load monomer definitions from files,
//! match links between residues, and find ideal geometric parameters.
struct GEMMI_DLL MonLib {
  std::string monomer_dir;                         //!< Path to monomer library directory
  std::map<std::string, ChemComp> monomers;        //!< Chemical components by name
  std::map<std::string, ChemLink> links;           //!< Link definitions by ID
  std::map<std::string, ChemMod> modifications;    //!< Modifications by ID
  std::map<std::string, ChemComp::Group> cc_groups;//!< Component name -> group mapping
  EnerLib ener_lib;                                //!< Energy library parameters

  //! @brief Get link definition by ID.
  //! @param link_id Link identifier (e.g., "TRANS", "SS")
  //! @return Pointer to ChemLink, or nullptr if not found
  const ChemLink* get_link(const std::string& link_id) const {
    auto link = links.find(link_id);
    return link != links.end() ? &link->second : nullptr;
  }

  //! @brief Get modification definition by name.
  //! @param name Modification identifier
  //! @return Pointer to ChemMod, or nullptr if not found
  const ChemMod* get_mod(const std::string& name) const {
    auto modif = modifications.find(name);
    return modif != modifications.end() ? &modif->second : nullptr;
  }

  //! @brief Find the best matching link for a bond between two residues.
  //! @param res1 First residue
  //! @param atom1 Atom name in first residue
  //! @param alt1 Alternate location for first residue
  //! @param res2 Second residue
  //! @param atom2 Atom name in second residue
  //! @param alt2 Alternate location for second residue
  //! @param min_bond_sq Minimum bond length squared (for filtering)
  //! @return Tuple of (best link, inverted flag, aliasing1, aliasing2)
  //!
  //! Returns the most specific link and a flag that is true
  //! if the order is comp2-comp1 in the link definition (inverted).
  //! Also returns aliasing tables if aliases were used for matching.
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

  //! @brief Get full path to monomer CIF file.
  //! @param code Monomer/residue code
  //! @return Full path (file may not exist)
  //!
  //! Returns path to the monomer cif file (the file may not exist).
  //! Combines monomer_dir with the standard relative path for the code.
  std::string path(const std::string& code) const {
    return monomer_dir + relative_monomer_path(code);
  }

  //! @brief Get relative path for a monomer file.
  //! @param code Monomer/residue code
  //! @return Relative path within monomer library (e.g., "a/ALA.cif")
  static std::string relative_monomer_path(const std::string& code);

  //! @brief Read monomer definitions from CIF document.
  //! @param doc CIF document containing monomer definitions
  void read_monomer_doc(const cif::Document& doc);

  //! @brief Read single monomer CIF file.
  //! @param path_ Path to monomer CIF file
  void read_monomer_cif(const std::string& path_);

  //! @brief Set monomer library directory path.
  //! @param monomer_dir_ Path to monomer library root directory
  //!
  //! Ensures path ends with directory separator (/ or \).
  void set_monomer_dir(const std::string& monomer_dir_) {
    monomer_dir = monomer_dir_;
    if (!monomer_dir.empty() && monomer_dir.back() != '/' && monomer_dir.back() != '\\')
      monomer_dir += '/';
  }

  //! @brief Read complete monomer library for specified residues.
  //! @param monomer_dir_ Path to monomer library directory
  //! @param resnames List of residue names to load
  //! @param logger Logger for messages
  //! @return True if all requested monomers were added
  //!
  //! Read mon_lib_list.cif, ener_lib.cif and required monomers.
  //! Returns true if all requested monomers were added.
  bool read_monomer_lib(const std::string& monomer_dir_,
                        const std::vector<std::string>& resnames,
                        const Logger& logger);

  //! @brief Find ideal distance between two atoms.
  //! @param cra1 First atom (Chain/Residue/Atom)
  //! @param cra2 Second atom (Chain/Residue/Atom)
  //! @return Ideal distance in Angstroms
  double find_ideal_distance(const const_CRA& cra1, const const_CRA& cra2) const;

  //! @brief Update old atom names to current standard.
  //! @param st Structure to update (modified in place)
  //! @param logger Logger for messages
  void update_old_atom_names(Structure& st, const Logger& logger) const;
};

//! @deprecated Use MonLib::read_monomer_lib() method instead
//! @brief Read monomer library (deprecated function).
//! @param monomer_dir Path to monomer library directory
//! @param resnames List of residue names to load
//! @param libin Optional additional CIF file to read
//! @param ignore_missing If true, don't fail on missing monomers
//! @return Populated MonLib
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
