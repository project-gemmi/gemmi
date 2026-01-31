//! @file
//! @brief Flat (table-like) structure representation for fast access.
//!
//! FlatStructure, FlatAtom

// Copyright Global Phasing Ltd.
//
// FlatStructure, FlatAtom

#ifndef GEMMI_FLAT_HPP_
#define GEMMI_FLAT_HPP_

#include <vector>
#include "model.hpp"

namespace gemmi {

//! @brief Single atom in flat representation with fixed-size string fields.
//!
//! All atom properties stored inline for efficient memory layout and NumPy compatibility.
struct FlatAtom {
  char atom_name[8] {};  //!< Atom name (padded with nulls)
  char residue_name[8] = {};  //!< Residue name
  char chain_id[8] = {};  //!< Chain identifier
  char subchain[8] = {};  //!< Subchain identifier
  char entity_id[8] = {};  //!< Entity identifier
  SeqId seq_id;  //!< Sequence ID with insertion code
  Position pos;  //!< Cartesian coordinates
  float occ = 1.0f;  //!< Occupancy
  float b_iso = 20.0f;  //!< Isotropic B-factor (arbitrary default value)
  char altloc = '\0';  //!< Alternate location (0 if not set)
  char het_flag = '\0';  //!< 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  EntityType entity_type = EntityType::Unknown;  //!< Entity classification
  Element element = El::X;  //!< Element type
  signed char charge = 0;  //!< Formal charge [-8, +8]
  SMat33<float> aniso = {0, 0, 0, 0, 0, 0};  //!< Anisotropic displacement parameters
  int model_num;  //!< Model number
  int serial = 0;  //!< Atom serial number
  bool selected = false;  //!< Selection flag

  //! @brief Get atom identifier string.
  //! @return String like "A/123/CA" for chain/seqid/atom
  std::string atom_str() const {
    ResidueId resid{seq_id, "", residue_name};
    return gemmi::atom_str(chain_id, resid, atom_name, altloc);
  }
};

//! @brief Flat table representation of entire structure.
//!
//! Converts hierarchical Structure to flat table of atoms for efficient processing
//! and NumPy array compatibility in Python bindings.
struct GEMMI_DLL FlatStructure {
  Structure empty_st;  //!< Template structure for metadata
  std::vector<FlatAtom> table;  //!< Flat table of all atoms
  bool strings_as_numbers = true;  //!< Encode strings as numbers for NumPy (Python bindings)

  //! @brief Construct from hierarchical structure.
  //! @param st Structure to flatten
  // Structure <-> FlatStructure
  FlatStructure(const Structure& st);
  //! @brief Reconstruct hierarchical structure from flat table.
  //! @return Hierarchical Structure
  Structure generate_structure();
};

} // namespace gemmi

#endif
