//! @file
//! @brief Generating biological assemblies by applying operations from struct Assembly to a Model.
//!
//! Provides functions for generating biological assemblies from crystallographic
//! asymmetric units by applying symmetry operations. Also includes utilities for
//! chain naming and renaming when creating assembly copies.

// Copyright 2020 Global Phasing Ltd.

#ifndef GEMMI_ASSEMBLY_HPP_
#define GEMMI_ASSEMBLY_HPP_

#include "model.hpp"  // for Model
#include "util.hpp"   // for in_vector
#include "logger.hpp" // for Logger

namespace gemmi {

//! @brief Strategy for naming copied chains in assemblies.
enum class HowToNameCopiedChain {
  Short,      //!< Use short 1-2 letter names (A, B, ..., Z, a, ..., z, AA, AB, ...)
  AddNumber,  //!< Add numeric postfix to original name (A becomes A1, A2, ...)
  Dup         //!< Keep duplicate names (not recommended for most uses)
};

//! @brief Utility for generating unique chain names when creating assemblies.
//!
//! Tracks used chain names and generates new unique names according to the
//! specified naming strategy. Prevents name collisions when copying chains.
struct ChainNameGenerator {
  HowToNameCopiedChain how;  //!< Naming strategy to use
  std::vector<std::string> used_names;  //!< Set of already-used chain names

  //! @brief Construct with naming strategy only.
  //! @param how_ Naming strategy
  ChainNameGenerator(HowToNameCopiedChain how_) : how(how_) {}

  //! @brief Construct with model and naming strategy.
  //! @param model Model whose chain names to register as used
  //! @param how_ Naming strategy
  //!
  //! Registers all existing chain names from the model unless using Dup strategy.
  ChainNameGenerator(const Model& model, HowToNameCopiedChain how_) : how(how_) {
    if (how != HowToNameCopiedChain::Dup)
      for (const Chain& chain : model.chains)
        used_names.push_back(chain.name);
  }

  //! @brief Try to register a name as used.
  //! @param name Chain name to register
  //! @return True if name was available and registered, false if already used
  bool try_add(const std::string& name) {
    if (in_vector(name, used_names))
      return false;
    used_names.push_back(name);
    return true;
  }

  //! @brief Generate a short unique chain name.
  //! @param preferred Preferred name (tried first)
  //! @return Unique 1- or 2-letter chain name
  //! @throws std::runtime_error if all 1- and 2-letter names are exhausted
  //!
  //! Tries the preferred name first. If unavailable, generates names in order:
  //! A-Z, a-z, 0-9, then two-letter combinations AA, AB, ..., 99.
  std::string make_short_name(const std::string& preferred) {
    static const char symbols[] = {
      'A','B','C','D','E','F','G','H','I','J','K','L','M',
      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
      'a','b','c','d','e','f','g','h','i','j','k','l','m',
      'n','o','p','q','r','s','t','u','v','w','x','y','z',
      '0','1','2','3','4','5','6','7','8','9'
    };
    if (try_add(preferred))
      return preferred;
    std::string name(1, 'A');
    for (char symbol : symbols) {
      name[0] = symbol;
      if (try_add(name))
        return name;
    }
    name += 'A';
    for (char symbol1 : symbols) {
      name[0] = symbol1;
      for (char symbol2 : symbols) {
        name[1] = symbol2;
        if (try_add(name))
          return name;
      }
    }
    fail("run out of 1- and 2-letter chain names");
  }

  //! @brief Generate a name with numeric postfix.
  //! @param base Base chain name
  //! @param n Starting number to try
  //! @return Unique name with numeric postfix (e.g., "A1", "A2")
  //!
  //! Tries base+n, then increments n until a unique name is found.
  std::string make_name_with_numeric_postfix(const std::string& base, int n) {
    std::string name = base;
    name += std::to_string(n);
    while (!try_add(name)) {
      name.resize(base.size());
      name += std::to_string(++n);
    }
    return name;
  }

  //! @brief Generate a new name according to the naming strategy.
  //! @param old Original chain name
  //! @param n Copy number (used for AddNumber strategy)
  //! @return Unique chain name
  //!
  //! Dispatches to the appropriate naming method based on the strategy.
  std::string make_new_name(const std::string& old, int n) {
    switch (how) {
      case HowToNameCopiedChain::Short: return make_short_name(old);
      case HowToNameCopiedChain::AddNumber: return make_name_with_numeric_postfix(old, n);
      case HowToNameCopiedChain::Dup: return old;
    }
    unreachable();
  }
};

//! @brief Ensure chain has a unique name within the model.
//! @param model Model containing the chain
//! @param chain Chain whose name to make unique
//!
//! Renames the chain if necessary to avoid name collision with other chains
//! in the model. Uses short naming strategy (A, B, ..., AA, AB, ...).
inline void ensure_unique_chain_name(const Model& model, Chain& chain) {
  ChainNameGenerator namegen(HowToNameCopiedChain::Short);
  for (const Chain& ch : model.chains)
    if (&ch != &chain)
      namegen.try_add(ch.name);
  chain.name = namegen.make_short_name(chain.name);
}

//! @brief Generate biological assembly by applying assembly operations.
//! @param assembly Assembly definition with operators
//! @param model Source model (asymmetric unit)
//! @param how Chain naming strategy for copied chains
//! @param logging Logger for messages
//! @return New model containing the assembled structure
//!
//! Creates a new model by applying the transformation operators defined in the
//! assembly to selected chains. Chains are renamed according to the specified
//! naming strategy to avoid collisions.
GEMMI_DLL Model make_assembly(const Assembly& assembly, const Model& model,
                              HowToNameCopiedChain how, const Logger& logging);

//! @brief Create pseudo-assembly representing the unit cell.
//! @param cell Unit cell with symmetry operations
//! @return Assembly that expands to the full unit cell (P1)
//!
//! Generates an Assembly object that applies all symmetry operations to expand
//! the asymmetric unit to the full crystallographic unit cell.
inline Assembly pseudo_assembly_for_unit_cell(const UnitCell& cell) {
  Assembly assembly("unit_cell");
  std::vector<Assembly::Operator> operators(cell.images.size() + 1);
  // operators[0] stays as identity
  for (size_t i = 1; i != operators.size(); ++i) {
    const FTransform& op = cell.images[i-1];
    operators[i].transform = cell.orth.combine(op.combine(cell.frac));
  }
  assembly.generators.push_back({{"(all)"}, {}, operators});
  return assembly;
}

//! @brief Transform structure to specified biological assembly.
//! @param st Structure to transform (modified in place)
//! @param assembly_name Name of assembly to generate, or "unit_cell" for P1 expansion
//! @param how Chain naming strategy for copied chains
//! @param logging Logger for messages
//! @param keep_spacegroup If true, preserves space group and unit cell metadata
//! @param merge_dist Maximum distance for merging overlapping atoms (Angstroms)
//!
//! If called with assembly_name="unit_cell" changes structure to unit cell (P1).
//! The keep_spacegroup parameter preserves space group and unit cell - is it needed?
GEMMI_DLL void transform_to_assembly(Structure& st, const std::string& assembly_name,
                                     HowToNameCopiedChain how, const Logger& logging,
                                     bool keep_spacegroup=false, double merge_dist=0.2);

//! @brief Expand model by applying NCS (non-crystallographic symmetry) operations.
//! @param model Source model to expand
//! @param ncs Vector of NCS operators to apply
//! @param how Chain naming strategy for NCS copies
//! @return New model with NCS copies added
//!
//! Creates copies of chains by applying non-crystallographic symmetry operators.
//! Useful for generating complete NCS assemblies from the asymmetric unit.
GEMMI_DLL Model expand_ncs_model(const Model& model, const std::vector<NcsOp>& ncs,
                                 HowToNameCopiedChain how);

//! @brief Merge overlapping equivalent atoms from different chains.
//! @param model Model to process (modified in place)
//! @param cell Unit cell for distance calculations
//! @param max_dist Maximum distance for considering atoms equivalent (Angstroms)
//! @param compare_serial If true, also compare atom serial numbers
//!
//! Searches and merges overlapping equivalent atoms from different chains.
//! To be used after expand_ncs() and make_assembly() to clean up duplicate atoms
//! that may appear at special positions or assembly interfaces.
GEMMI_DLL void merge_atoms_in_expanded_model(Model& model, const UnitCell& cell,
                                             double max_dist=0.2, bool compare_serial=true);

//! @brief Shorten all chain names in structure to 1-2 letters.
//! @param st Structure to modify
//!
//! Renames all chains using short naming convention (A, B, ..., AA, AB, ...).
//! Useful for compatibility with older PDB format limitations.
GEMMI_DLL void shorten_chain_names(Structure& st);

//! @brief Expand structure by applying NCS operations.
//! @param st Structure to expand (modified in place)
//! @param how Chain naming strategy for NCS copies
//! @param merge_dist Maximum distance for merging overlapping atoms (Angstroms)
//!
//! Applies NCS operations defined in the structure to all models, creating
//! NCS copies of chains. Optionally merges overlapping atoms.
GEMMI_DLL void expand_ncs(Structure& st, HowToNameCopiedChain how, double merge_dist=0.2);

//! @brief Split chains by segment ID into separate chains.
//! @param model Model to process (modified in place)
//! @param how Chain naming strategy (Dup adds segment name to chain name)
//!
//! Creates separate chains for each segment ID within a chain. When using
//! HowToNameCopiedChain::Dup, the segment name is added to the chain name.
//! Useful for handling PDB files with segment IDs.
GEMMI_DLL void split_chains_by_segments(Model& model, HowToNameCopiedChain how);

} // namespace gemmi
#endif
