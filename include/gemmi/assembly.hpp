// Copyright 2020 Global Phasing Ltd.
//
// Generating biological assemblies by applying operations
// from struct Assembly to a Model.
// Includes chain (re)naming utilities.

#ifndef GEMMI_ASSEMBLY_HPP_
#define GEMMI_ASSEMBLY_HPP_

#include "model.hpp"  // for Model
#include "util.hpp"   // for in_vector
#include "logger.hpp" // for Logger

namespace gemmi {

/// Naming strategy for chains copied during assembly or NCS expansion.
///
/// @brief Specifies how to name copied chains to ensure uniqueness.
enum class HowToNameCopiedChain {
  /// Use 1–2 character chain names cycling through A–Z, a–z, 0–9
  Short,
  /// Append a numeric suffix to the original chain name
  AddNumber,
  /// Keep the original chain name (may result in duplicates)
  Dup
};

/// @brief Utility to generate unique chain names for assembly and NCS expansion.
///
/// Manages a set of used chain names and provides methods to generate new,
/// unique names according to the specified naming strategy.
struct ChainNameGenerator {
  /// The naming strategy to apply
  HowToNameCopiedChain how;
  /// List of already-assigned chain names to avoid duplicates
  std::vector<std::string> used_names;

  /// Initialize with a naming strategy.
  ///
  /// @param how_ The chain naming strategy to use
  ChainNameGenerator(HowToNameCopiedChain how_) : how(how_) {}

  /// Initialize with chain names pre-populated from a Model.
  ///
  /// @param model The model to extract existing chain names from
  /// @param how_ The chain naming strategy to use
  ChainNameGenerator(const Model& model, HowToNameCopiedChain how_) : how(how_) {
    if (how != HowToNameCopiedChain::Dup)
      for (const Chain& chain : model.chains)
        used_names.push_back(chain.name);
  }

  /// Register a chain name if not already used.
  ///
  /// @param name The name to register
  /// @return true if the name was successfully added, false if already in used_names
  bool try_add(const std::string& name) {
    if (in_vector(name, used_names))
      return false;
    used_names.push_back(name);
    return true;
  }

  /// Generate a unique 1–2 character name starting from preferred.
  ///
  /// Cycles through A–Z, a–z, 0–9 to find an available name.
  ///
  /// @param preferred The preferred 1-character name to try first
  /// @return A unique short name
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

  /// Generate a name by appending a numeric postfix to a base name.
  ///
  /// Increments the numeric suffix until a unique name is found.
  ///
  /// @param base The base chain name
  /// @param n The starting number to append
  /// @return A unique name of the form "base" + number
  std::string make_name_with_numeric_postfix(const std::string& base, int n) {
    std::string name = base;
    name += std::to_string(n);
    while (!try_add(name)) {
      name.resize(base.size());
      name += std::to_string(++n);
    }
    return name;
  }

  /// Generate a new name according to the configured naming strategy.
  ///
  /// Dispatches to make_short_name(), make_name_with_numeric_postfix(),
  /// or returns the original name based on the how field.
  ///
  /// @param old The original chain name
  /// @param n The numeric suffix (used only for AddNumber strategy)
  /// @return A new unique chain name
  std::string make_new_name(const std::string& old, int n) {
    switch (how) {
      case HowToNameCopiedChain::Short: return make_short_name(old);
      case HowToNameCopiedChain::AddNumber: return make_name_with_numeric_postfix(old, n);
      case HowToNameCopiedChain::Dup: return old;
    }
    unreachable();
  }
};

/// Rename a chain to ensure its name is unique within the model.
///
/// Uses the Short naming strategy to generate a 1–2 character name
/// that does not conflict with other chains in the model.
///
/// @param model The model containing the chain
/// @param chain The chain to rename
inline void ensure_unique_chain_name(const Model& model, Chain& chain) {
  ChainNameGenerator namegen(HowToNameCopiedChain::Short);
  for (const Chain& ch : model.chains)
    if (&ch != &chain)
      namegen.try_add(ch.name);
  chain.name = namegen.make_short_name(chain.name);
}

/// Apply all generators in an assembly to a model to generate a biological unit.
///
/// Creates a new Model containing copies of the input model transformed
/// according to each generator in the assembly. Chain names are made unique
/// according to the specified strategy.
///
/// @param assembly The assembly record containing transformation operators
/// @param model The input model to transform
/// @param how The strategy for naming copied chains
/// @param logging Logger for warning messages
/// @return A new Model containing the assembled biological unit
GEMMI_DLL Model make_assembly(const Assembly& assembly, const Model& model,
                              HowToNameCopiedChain how, const Logger& logging);

/// Create a pseudo-assembly representing all unit-cell images.
///
/// Generates an Assembly whose operators produce all crystallographic
/// images of the structure within the unit cell. Useful for packing displays
/// and structure visualization.
///
/// @param cell The unit cell with precomputed image transformations
/// @return An Assembly object whose generators apply all unit-cell symmetry
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

/// Modify a Structure in-place to contain only a named assembly.
///
/// Replaces all models in the structure with the expanded assembly,
/// optionally preserving the space group. If assembly_name="unit_cell",
/// converts the structure to the unit cell (P1).
///
/// @param st The structure to modify in-place
/// @param assembly_name The name of the assembly to apply
/// @param how The strategy for naming copied chains
/// @param logging Logger for warning messages
/// @param keep_spacegroup If true, preserve the original space group and unit cell
/// @param merge_dist Distance threshold for merging nearby atoms (default 0.2 Å)
GEMMI_DLL void transform_to_assembly(Structure& st, const std::string& assembly_name,
                                     HowToNameCopiedChain how, const Logger& logging,
                                     bool keep_spacegroup=false, double merge_dist=0.2);


/// Expand a model by applying NCS operations.
///
/// Generates copies of the input model by applying each NCS (non-crystallographic symmetry)
/// operation, returning a new Model with all copies combined.
///
/// @param model The input model to expand
/// @param ncs Vector of NCS operations to apply
/// @param how The strategy for naming copied chains
/// @return A new Model containing the original and all NCS-transformed copies
GEMMI_DLL Model expand_ncs_model(const Model& model, const std::vector<NcsOp>& ncs,
                                 HowToNameCopiedChain how);

/// Merge atoms within a distance threshold after NCS or assembly expansion.
///
/// Searches for overlapping atoms from different chains that are equivalent
/// under the unit cell symmetry and merges them. Typically used after
/// expand_ncs() and make_assembly() to eliminate duplicate atoms.
///
/// @param model The model to modify in-place
/// @param cell The unit cell (used for distance calculations)
/// @param max_dist Maximum distance for merging atoms (default 0.2 Å)
/// @param compare_serial If true, only merge atoms with matching serial numbers (default true)
GEMMI_DLL void merge_atoms_in_expanded_model(Model& model, const UnitCell& cell,
                                             double max_dist=0.2, bool compare_serial=true);


/// Rename all chains to use the shortest possible unique names.
///
/// Replaces chain names with 1–2 character names (A–Z, a–z, 0–9) while
/// maintaining uniqueness within the structure.
///
/// @param st The structure to modify in-place
GEMMI_DLL void shorten_chain_names(Structure& st);

/// Expand NCS operations in all models of a structure in-place.
///
/// Applies all NCS (non-crystallographic symmetry) operations to every model
/// in the structure, generating copies and optionally merging nearby atoms.
///
/// @param st The structure to expand in-place
/// @param how The strategy for naming copied chains
/// @param merge_dist Distance threshold for merging atoms (default 0.2 Å)
GEMMI_DLL void expand_ncs(Structure& st, HowToNameCopiedChain how, double merge_dist=0.2);

/// Split chains at residue segment boundaries into separate chains.
///
/// Divides each chain into multiple chains at segment breaks, with new chain
/// names generated according to the strategy. If how=HowToNameCopiedChain::Dup,
/// the segment name is appended to the chain name.
///
/// @param model The model to modify in-place
/// @param how The strategy for naming split chains
GEMMI_DLL void split_chains_by_segments(Model& model, HowToNameCopiedChain how);

/// Find all crystallographic symmetry images within a distance radius.
///
/// Searches for non-identity symmetry images of the structure around a given
/// position. Uses only the first model and is intended for viewer applications.
///
/// @param st The structure to query
/// @param pos The center position for the search
/// @param radius The search radius in Ångströms
/// @return Vector of NearestImage objects describing the nearby symmetry images
GEMMI_DLL std::vector<NearestImage> get_nearby_sym_ops(const Structure& st,
                                                       const Position& pos,
                                                       double radius);

/// Extract a symmetry image of the structure as a new Structure object.
///
/// Creates and returns a copy of the input structure transformed to the
/// coordinates of the specified symmetry image.
///
/// @param st The input structure
/// @param image The symmetry image specification
/// @return A new Structure object containing the transformed copy
GEMMI_DLL Structure get_sym_image(const Structure& st, const NearestImage& image);

} // namespace gemmi
#endif
