// Copyright Global Phasing Ltd.
//
// FlatStructure, FlatAtom

#ifndef GEMMI_FLAT_HPP_
#define GEMMI_FLAT_HPP_

#include <vector>
#include "model.hpp"

namespace gemmi {

/// @brief A flat representation of a single atom for efficient array storage.
/// Stores atomic properties in a contiguous, columnar-friendly format suitable
/// for NumPy arrays and bulk data processing.
struct FlatAtom {
  /// @brief Atom name (e.g., "CA", "CB").
  char atom_name[8] {};
  /// @brief Residue name (e.g., "ALA", "GLY").
  char residue_name[8] = {};
  /// @brief Chain identifier.
  char chain_id[8] = {};
  /// @brief Subchain identifier.
  char subchain[8] = {};
  /// @brief Entity identifier.
  char entity_id[8] = {};
  /// @brief Sequence identifier (residue number with insertion code).
  SeqId seq_id;
  /// @brief Cartesian position in Angstroms.
  Position pos;
  /// @brief Occupancy (0.0 to 1.0).
  float occ = 1.0f;
  /// @brief Isotropic B-factor (arbitrary default 20.0 Angstrom^2).
  float b_iso = 20.0f;
  /// @brief Alternate location indicator (0 if not set).
  char altloc = '\0';
  /// @brief Atom type flag ('A' = ATOM, 'H' = HETATM, 0 = unspecified).
  char het_flag = '\0';
  /// @brief Entity type classification.
  EntityType entity_type = EntityType::Unknown;
  /// @brief Chemical element.
  Element element = El::X;
  /// @brief Formal charge (-8 to +8).
  signed char charge = 0;
  /// @brief Anisotropic thermal parameters.
  SMat33<float> aniso = {0, 0, 0, 0, 0, 0};
  /// @brief Model number.
  int model_num;
  /// @brief Atom serial number.
  int serial = 0;
  /// @brief Selection flag.
  bool selected = false;

  /// @brief Generate a formatted atom identifier string.
  /// @return String containing chain, residue, and atom identifiers.
  std::string atom_str() const {
    ResidueId resid{seq_id, "", residue_name};
    return gemmi::atom_str(chain_id, resid, atom_name, altloc);
  }
};

/// @brief A flat representation of a structure as a contiguous table of atoms.
/// Converts hierarchical Structure into a single table of FlatAtom entries,
/// enabling efficient access patterns and integration with array-based tools.
struct GEMMI_DLL FlatStructure {
  /// @brief Empty structure holding only metadata (name, cell, spacegroup).
  Structure empty_st;
  /// @brief Table of flat atoms from the structure.
  std::vector<FlatAtom> table;
  /// @brief If true, string fields are interpreted as numeric when accessed as NumPy arrays.
  bool strings_as_numbers = true;

  /// @brief Convert a macromolecular structure to flat representation.
  /// @param st The structure to flatten.
  FlatStructure(const Structure& st);

  /// @brief Reconstruct a hierarchical structure from the flat representation.
  /// @return A Structure with chains, residues, and atoms rebuilt from the table.
  Structure generate_structure();
};

} // namespace gemmi

#endif
