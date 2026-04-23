// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include <map>           // for multimap
#include <memory>        // for unique_ptr
#include <unordered_map> // for unordered_map
#include "chemcomp.hpp"  // for ChemComp
#include "monlib.hpp"    // for MonLib
#include "model.hpp"     // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral
#include "logger.hpp"    // for Logger

namespace gemmi {

/// @brief Specification for how to modify hydrogen atoms during topology preparation.
enum class HydrogenChange {
  NoChange,         ///< Leave hydrogen atoms as they are in the input model.
  Shift,            ///< Move hydrogen atoms to standard positions based on geometry.
  Remove,           ///< Remove all hydrogen atoms from the model.
  ReAdd,            ///< Remove all hydrogen atoms and then re-add them at standard positions.
  ReAddButWater,    ///< Remove and re-add hydrogen atoms, except in water molecules.
  ReAddKnown        ///< Re-add only hydrogen atoms that are known in the monomer library.
};

/// @brief Topology of restraints from a monomer library applied to a crystallographic model.
///
/// Non-copyable due to internal atom pointers set up during apply_restraints()
/// that reference ResInfo chemical component restraint data.
struct GEMMI_DLL Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  /// @brief A bond restraint between two atoms.
  ///
  /// Holds a reference to the restraint specification and the two atoms,
  /// and provides methods to calculate the bond distance and z-score.
  struct Bond {
    /// @brief Pointer to the restraint specification from the monomer library.
    const Restraints::Bond* restr;
    /// @brief The two atoms involved in the bond.
    std::array<Atom*, 2> atoms;
    /// @brief Asymmetric unit relationship between the atoms.
    Asu asu;
    /// @brief Calculate the bond distance in Angstroms.
    /// @return Distance between atoms, or NAN if atoms are in different asymmetric units.
    double calculate() const {
      return asu != Asu::Different ? atoms[0]->pos.dist(atoms[1]->pos) : NAN;
    }
    /// @brief Calculate the z-score for a given distance.
    /// @param d The distance value.
    /// @return Z-score: (distance - ideal_value) / esd
    double calculate_z_(double d) const { return std::abs(d - restr->value) / restr->esd; }
    /// @brief Calculate the z-score for the current bond distance.
    /// @return Z-score of the observed distance relative to the restraint.
    double calculate_z() const { return calculate_z_(calculate()); }
  };
  /// @brief An angle restraint between three atoms.
  ///
  /// Holds a reference to the restraint specification and the three atoms,
  /// and provides methods to calculate the angle and z-score.
  struct Angle {
    /// @brief Pointer to the restraint specification from the monomer library.
    const Restraints::Angle* restr;
    /// @brief The three atoms involved in the angle (atom[0]-atom[1]-atom[2]).
    std::array<Atom*, 3> atoms;
    /// @brief Calculate the angle value.
    /// @return Angle in radians.
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    /// @brief Calculate the z-score for the current angle.
    /// @return Z-score of the observed angle relative to the restraint.
    double calculate_z() const { return angle_z(calculate(), *restr); }
  };
  /// @brief A torsion (dihedral) angle restraint between four atoms.
  ///
  /// Holds a reference to the restraint specification and the four atoms,
  /// and provides methods to calculate the dihedral angle and z-score.
  struct Torsion {
    /// @brief Pointer to the restraint specification from the monomer library.
    const Restraints::Torsion* restr;
    /// @brief The four atoms involved in the torsion (atoms[0]-atoms[1]-atoms[2]-atoms[3]).
    std::array<Atom*, 4> atoms;
    /// @brief Calculate the dihedral angle value.
    /// @return Dihedral angle in radians.
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    /// @brief Calculate the z-score for the current dihedral angle.
    ///
    /// The z-score accounts for the periodicity of the torsion restraint.
    /// @return Z-score of the observed dihedral relative to the restraint.
    double calculate_z() const {
      return angle_z(calculate(), *restr, 360. / std::max(1, restr->period));
    }
  };
  /// @brief A chirality restraint on the stereochemistry around a chiral center.
  ///
  /// Holds a reference to the restraint specification and the four atoms
  /// (center and three substituents), and provides methods to calculate the
  /// chiral volume and z-score.
  struct Chirality {
    /// @brief Pointer to the restraint specification from the monomer library.
    const Restraints::Chirality* restr;
    /// @brief The four atoms: atoms[0] is the chiral center, atoms[1-3] are substituents.
    std::array<Atom*, 4> atoms;
    /// @brief Calculate the chiral volume.
    /// @return The signed chiral volume (a scalar triple product).
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    /// @brief Calculate the z-score for the chiral volume.
    /// @param ideal_abs_vol Ideal absolute value of the chiral volume.
    /// @param esd Standard deviation of the restraint.
    /// @return Z-score: absolute deviation from ideal value divided by esd.
    double calculate_z(double ideal_abs_vol, double esd) const {
      double calc = calculate();
      if (restr->sign == ChiralityType::Negative ||
          (restr->sign == ChiralityType::Both && calc < 0))
        ideal_abs_vol *= -1;
      return std::abs(calc - ideal_abs_vol) / esd;
    }
    /// @brief Check whether the chirality is correct.
    /// @return True if the chirality sign matches the restraint specification.
    bool check() const { return !restr->is_wrong(calculate()); }
  };
  /// @brief A planar restraint on a group of atoms.
  ///
  /// Holds a reference to the restraint specification and the atoms that
  /// should lie in a plane.
  struct Plane {
    /// @brief Pointer to the restraint specification from the monomer library.
    const Restraints::Plane* restr;
    /// @brief The atoms that should lie in the plane.
    std::vector<Atom*> atoms;
    /// @brief Check whether an atom is part of this plane restraint.
    /// @param atom The atom to check.
    /// @return True if the atom is in the atoms vector.
    bool has(const Atom* atom) const {
      return in_vector(const_cast<Atom*>(atom), atoms);
    }
  };

  /// @brief Type of restraint rule.
  enum class RKind {
    Bond,      ///< Bond distance restraint.
    Angle,     ///< Angle restraint.
    Torsion,   ///< Torsion (dihedral) angle restraint.
    Chirality, ///< Chirality restraint.
    Plane      ///< Planarity restraint.
  };

  /// @brief A reference to a restraint rule.
  ///
  /// Identifies which type of restraint and its index in the corresponding
  /// vector (bonds, angles, torsions, chirs, or planes) in Topo.
  struct Rule {
    /// @brief The kind of restraint.
    RKind rkind;
    /// @brief Index in the respective vector (bonds, angles, torsions, chirs, or planes).
    size_t index;
  };

  /// @brief A link between two residues with associated restraints.
  ///
  /// Describes a covalent link (such as a peptide bond or a disulfide bridge)
  /// between two residues, including the restraint rules and bonding information.
  struct Link {
    /// @brief Link name from the monomer library (e.g., "PLNK", "disulf").
    std::string link_id;
    /// @brief Pointer to the first residue.
    Residue* res1 = nullptr;
    /// @brief Pointer to the second residue.
    Residue* res2 = nullptr;
    /// @brief Restraint rules applied by this link.
    std::vector<Rule> link_rules;
    /// @brief Alternate location indicator for res1.
    char alt1 = '\0';
    /// @brief Alternate location indicator for res2.
    char alt2 = '\0';
    /// @brief Asymmetric unit relationship between res1 and res2.
    ///
    /// Used only in Links in ChainInfo::extras.
    Asu asu = Asu::Any;
    /// @brief Helper field for CISPEP record generation in output.
    bool is_cis = false;

    /// @brief Cached atom name ID for res1 (used in find_polymer_link).
    int atom1_name_id = 0;
    /// @brief Cached atom name ID for res2 (used in find_polymer_link).
    int atom2_name_id = 0;

    /// @brief Pointer to aliasing information for res1.
    ///
    /// Points to a vector element in ChemComp::aliases.
    /// The pointer remains valid even if a ChemComp is moved.
    const ChemComp::Aliasing* aliasing1 = nullptr;
    /// @brief Pointer to aliasing information for res2.
    ///
    /// Points to a vector element in ChemComp::aliases.
    /// The pointer remains valid even if a ChemComp is moved.
    const ChemComp::Aliasing* aliasing2 = nullptr;

    /// @brief Calculate the pointer difference between residues.
    ///
    /// Only valid for polymer links where res1 and res2 are in the same Chain.
    /// @return Signed distance in residues (res1 - res2).
    std::ptrdiff_t res_distance() const { return res1 - res2; }
  };

  /// @brief A chemical modification applied to a residue.
  ///
  /// Describes a ChemMod from the monomer library and the specific atom
  /// group (aliasing) to which it applies.
  struct Mod {
    /// @brief ID of the ChemMod from the monomer library dictionary.
    std::string id;
    /// @brief Atom group alias to which the modification applies.
    ChemComp::Group alias;
    /// @brief Alternate location indicator ('\0' for all conformers).
    char altloc;

    /// @brief Check equality between two modifications.
    /// @param o The other modification to compare.
    /// @return True if id, alias, and altloc are identical.
    bool operator==(const Mod& o) const {
      return id == o.id && alias == o.alias && altloc == o.altloc;
    }
  };

  /// @brief Final chemical component with modifications applied.
  ///
  /// Represents a ChemComp with all modifications already applied.
  struct FinalChemComp {
    /// @brief Alternate location indicator for which these restraints apply.
    char altloc;
    /// @brief Pointer to the ChemComp with modifications applied.
    const ChemComp* cc;
  };

  /// @brief Information about a residue in the topology.
  ///
  /// Contains the residue, its chemical composition (with modifications),
  /// link information, and hydrogen bonding data.
  struct ResInfo {
    /// @brief Pointer to the residue in the model.
    Residue* res;
    /// @brief Links to previous residue(s).
    ///
    /// In case of microheterogeneity, there may be 2 or more previous residues.
    std::vector<Link> prev;
    /// @brief Chemical modifications applied to this residue.
    std::vector<Mod> mods;
    /// @brief Pointer to the original ChemComp from MonLib::monomers.
    const ChemComp* orig_chemcomp = nullptr;
    /// @brief ChemComps with modifications applied, per conformer.
    std::vector<FinalChemComp> chemcomps;
    /// @brief Restraint rules applied to this residue.
    std::vector<Rule> monomer_rules;
    /// @brief Two hydrogen-bonded donors with lowest energy (for DSSP).
    std::array<Topo::ResInfo*, 2> donors = {nullptr, nullptr};
    /// @brief Two hydrogen-bonded acceptors with lowest energy (for DSSP).
    std::array<Topo::ResInfo*, 2> acceptors = {nullptr, nullptr};
    /// @brief Energies of the two donor hydrogen bonds.
    std::array<double, 2> donor_energies = {0.0, 0.0};
    /// @brief Energies of the two acceptor hydrogen bonds.
    std::array<double, 2> acceptor_energies = {0.0, 0.0};

    /// @brief Constructor.
    /// @param r The residue to associate with this ResInfo.
    ResInfo(Residue* r) : res(r) {}
    /// @brief Add a modification to this residue.
    /// @param m The ID of the modification (from MonLib).
    /// @param aliasing Pointer to the aliasing information, or nullptr.
    /// @param altloc Alternate location indicator ('\0' for all conformers).
    void add_mod(const std::string& m, const ChemComp::Aliasing* aliasing, char altloc) {
      if (!m.empty()) {
        auto alias_group = aliasing ? aliasing->group : ChemComp::Group::Null;
        Mod mod{m, alias_group, altloc};
        if (!in_vector(mod, mods))
          mods.push_back(mod);
      }
    }

    /// @brief Get the final ChemComp for a specific conformer.
    /// @param altloc Alternate location indicator.
    /// @return Reference to the ChemComp for the specified conformer, or the first one if not found.
    const ChemComp& get_final_chemcomp(char altloc) const {
      if (chemcomps.size() == 1)
        return *chemcomps[0].cc;
      assert(!chemcomps.empty());
      for (const FinalChemComp& it : chemcomps)
        if (it.altloc == altloc)
          return *it.cc;
      return *chemcomps[0].cc;
    }
  };

  /// @brief Information about a sub-chain (continuous polymer segment).
  struct ChainInfo {
    /// @brief Reference to the full Chain.
    const Chain& chain_ref;
    /// @brief Name of the sub-chain.
    std::string subchain_name;
    /// @brief Entity ID from the PDB ENTITY_POLY record.
    std::string entity_id;
    /// @brief Whether this sub-chain is a polymer.
    bool polymer;
    /// @brief Type of polymer (protein, DNA, RNA, etc.).
    PolymerType polymer_type;
    /// @brief Residue information for each residue in the sub-chain.
    std::vector<ResInfo> res_infos;

    /// @brief Constructor.
    /// @param subchain The residue span for this sub-chain.
    /// @param chain The full Chain.
    /// @param ent Pointer to the Entity, or nullptr.
    ChainInfo(ResidueSpan& subchain, const Chain& chain, const Entity* ent);
    /// @brief Iterator type for ResInfo.
    using iterator = std::vector<ResInfo>::iterator;
    /// @brief Find the end of a residue group.
    ///
    /// Residues belong to the same group if they have the same group_key().
    /// @param b Iterator to the start of the group.
    /// @return Iterator to the first residue not in the same group.
    iterator group_end(iterator b) const {
      auto e = b + 1;
      while (e != res_infos.end() && e->res->group_key() == b->res->group_key())
        ++e;
      return e;
    }
  };

  /// @brief Check whether an atom is part of a structure.
  /// @tparam T A structure with an atoms member (e.g., Bond, Angle, etc.).
  /// @param a The atom to search for.
  /// @param t The structure to search in.
  /// @return Index of the atom in t.atoms, or -1 if not found.
  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  /// @brief Logger for warnings and informational messages.
  Logger logger{};
  /// @brief Internal flag for apply_restraints().
  bool only_bonds = false;
  /// @brief Information about each sub-chain in the model.
  std::vector<ChainInfo> chain_infos;
  /// @brief Extra links not bound to specific chains.
  std::vector<Link> extras;

  /// @brief Bond restraints applied to the model.
  std::vector<Bond> bonds;
  /// @brief Angle restraints applied to the model.
  std::vector<Angle> angles;
  /// @brief Torsion restraints applied to the model.
  std::vector<Torsion> torsions;
  /// @brief Chirality restraints applied to the model.
  std::vector<Chirality> chirs;
  /// @brief Planarity restraints applied to the model.
  std::vector<Plane> planes;

  /// @brief Index of bonds by atom.
  ///
  /// Maps each atom to the bonds it is part of.
  std::multimap<const Atom*, Bond*> bond_index;
  /// @brief Index of angles by center atom.
  ///
  /// Maps each atom to the angles where it is the center atom (atoms[1]).
  std::multimap<const Atom*, Angle*> angle_index;
  /// @brief Index of torsions by middle atoms.
  ///
  /// Maps atoms[1] and atoms[2] to the torsions containing them.
  std::multimap<const Atom*, Torsion*> torsion_index;
  /// @brief Index of planes by atom.
  ///
  /// Maps each atom to the planes it is part of.
  std::multimap<const Atom*, Plane*> plane_index;

  /// @brief Find the ResInfo for a residue.
  /// @param res The residue to search for.
  /// @return Pointer to the ResInfo, or nullptr if not found.
  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chain_infos)
      for (ResInfo& ri : ci.res_infos)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  /// @brief Get the first bond restraint in a link.
  /// @param link The link to search.
  /// @return Pointer to the first Bond in link.link_rules, or nullptr if none.
  Bond* first_bond_in_link(const Link& link) {
    for (const Rule& rule : link.link_rules)
      if (rule.rkind == RKind::Bond)
        return &bonds[rule.index];
    return nullptr;
  }

  /// @brief Find a bond restraint between two atoms.
  /// @param a First atom.
  /// @param b Second atom.
  /// @return Pointer to the bond restraint, or nullptr if no bond is restrained.
  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    auto range = bond_index.equal_range(a);
    for (auto i = range.first; i != range.second; ++i) {
      const Bond* bond = i->second;
      if ((bond->atoms[0] == b && bond->atoms[1] == a) ||
          (bond->atoms[1] == b && bond->atoms[0] == a))
        return bond->restr;
    }
    return nullptr;
  }

  /// @brief Find an angle restraint between three atoms.
  /// @param a First atom.
  /// @param b Center atom.
  /// @param c Third atom.
  /// @return Pointer to the angle restraint, or nullptr if no angle is restrained.
  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    auto range = angle_index.equal_range(b);
    for (auto i = range.first; i != range.second; ++i) {
      const Angle* ang = i->second;
      if ((ang->atoms[0] == a && ang->atoms[2] == c) ||
          (ang->atoms[0] == c && ang->atoms[2] == a))
        return ang->restr;
    }
    return nullptr;
  }

  /// @brief Get the chirality restraint for a chiral center.
  /// @param ctr The chiral center atom.
  /// @return Pointer to the Chirality restraint, or nullptr if none.
  const Chirality* get_chirality(const Atom* ctr) const {
    for (const Chirality& chir : chirs)
      if (chir.atoms[0] == ctr)
        return &chir;
    return nullptr;
  }

  /// @brief Get the ideal absolute chiral volume for a chirality restraint.
  /// @param ch The chirality restraint.
  /// @return Ideal absolute value of the chiral volume.
  double ideal_chiral_abs_volume(const Chirality &ch) const;

  /// @brief Apply restraints from a Restraints object to a residue.
  /// @param rt The restraints specification.
  /// @param res The residue to apply restraints to.
  /// @param res2 Second residue (for inter-residue restraints), or nullptr.
  /// @param asu Asymmetric unit relationship between atoms.
  /// @param altloc1 Alternate location indicator for atoms in res.
  /// @param altloc2 Alternate location indicator for atoms in res2.
  /// @param require_alt If true, only apply restraints matching the altloc.
  /// @return Vector of restraint rules applied.
  std::vector<Rule> apply_restraints(const Restraints& rt,
                                     Residue& res, Residue* res2, Asu asu,
                                     char altloc1, char altloc2, bool require_alt);
  /// @brief Apply restraints from a Link.
  /// @param link The Link to apply.
  /// @param monlib The monomer library (used to get ChemComp information).
  void apply_restraints_from_link(Link& link, const MonLib& monlib);

  /// @brief Initialize the topology from a Structure and MonLib.
  ///
  /// This method populates the internal topology state from the model and monomer library.
  ///
  /// @param st The Structure (non-const to assign link_id to connections).
  /// @param model0 The Model (non-const to store pointers to residues).
  /// @param monlib The MonLib (non-const; may be modified by addition of extra links).
  /// @param ignore_unknown_links If true, skip links not in the library.
  ///
  /// @note After this step, do not add or remove residues from the model,
  /// as Topo holds internal pointers to them.
  /// @note The monlib may be modified by the addition of extra links from the model.
  void initialize_refmac_topology(Structure& st, Model& model0,
                                  MonLib& monlib, bool ignore_unknown_links=false);

  /// @brief Apply all restraints from the monomer library to the model.
  ///
  /// This populates the bonds, angles, torsions, chirs, and planes vectors.
  ///
  /// @param monlib The monomer library (used only for link information).
  ///
  /// @note This step stores pointers to gemmi::Atom's from model0,
  /// so after this step do not add or remove atoms from the model.
  void apply_all_restraints(const MonLib& monlib);

  /// @brief Prepare the atom-to-restraint indices.
  ///
  /// Populates bond_index, angle_index, torsion_index, and plane_index
  /// for efficient lookups.
  void create_indices();

  /// @brief Find a polymer link between two atoms.
  ///
  /// Searches for a matching Link in ResInfo::prev lists.
  ///
  /// @param a1 First atom address.
  /// @param a2 Second atom address.
  /// @return Pointer to the Link, or nullptr if no matching link is found.
  Link* find_polymer_link(const AtomAddress& a1, const AtomAddress& a2);

  /// @brief Generate CISPEP records in the Structure based on cis peptide bonds.
  /// @param st The Structure (non-const to add CISPEP records).
  void set_cispeps_in_structure(Structure& st);

private:
  // storage for link restraints modified by aliases
  std::vector<std::unique_ptr<Restraints>> rt_storage;
  // cache for ChemComps after applying modifications
  std::unordered_map<std::string, std::unique_ptr<ChemComp>> cc_cache;
  // storage for ad-hoc ChemComps (placeholders for those missing in MonLib)
  std::vector<std::unique_ptr<ChemComp>> cc_storage;

  void setup_connection(Connection& conn, Model& model0, MonLib& monlib,
                        bool ignore_unknown_links);
};

/// @brief Prepare the topology for a Structure.
///
/// Creates and initializes a Topo object for the given Structure and model,
/// with optional hydrogen atom adjustments.
///
/// @param st The Structure.
/// @param monlib The monomer library.
/// @param model_index Index of the model to use (typically 0).
/// @param h_change Specification for how to modify hydrogen atoms.
/// @param reorder If true, reorder atoms in the model.
/// @param logger Logger for warnings and informational messages.
/// @param ignore_unknown_links If true, skip links not in the library.
/// @param use_cispeps If true, identify and mark cis peptide bonds.
/// @return A unique_ptr to the initialized Topo.
GEMMI_DLL std::unique_ptr<Topo>
prepare_topology(Structure& st, MonLib& monlib, size_t model_index,
                 HydrogenChange h_change, bool reorder,
                 const Logger& logger={}, bool ignore_unknown_links=false,
                 bool use_cispeps=false);

/// @brief Create a ChemComp with restraints for a residue.
///
/// Generates a ChemComp with default restraints (bonds, angles, torsions, chiral volumes)
/// inferred from the residue's atom geometry.
///
/// @param res The residue to create a ChemComp for.
/// @return A unique_ptr to the generated ChemComp.
GEMMI_DLL std::unique_ptr<ChemComp> make_chemcomp_with_restraints(const Residue& res);

/// @brief Find atoms in the model that are missing from the restraints.
///
/// Identifies atoms that are present in the model but not found in the
/// corresponding chemical component definitions.
///
/// @param topo The topology with applied restraints.
/// @param including_hydrogen If true, include hydrogen atoms in the search.
/// @return Vector of addresses of missing atoms.
GEMMI_DLL std::vector<AtomAddress> find_missing_atoms(const Topo& topo,
                                                      bool including_hydrogen=false);

} // namespace gemmi
#endif
