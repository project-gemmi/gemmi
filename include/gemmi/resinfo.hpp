//! @file
//! @brief List of common residues with basic classification data.

// Copyright 2018 Global Phasing Ltd.
//
// List of common residues with basic data.

#ifndef GEMMI_RESINFO_HPP_
#define GEMMI_RESINFO_HPP_

#include <cstdint>  // for uint8_t
#include <string>
#include <vector>
#include "fail.hpp"

namespace gemmi {

//! @brief Classification of residue types.
//!
//! Simple approximate classification:
//! - AA: Standard amino acid
//! - AAD: D-amino acid
//! - PAA: Proline-like amino acid
//! - MAA: Methylated amino acid
//! - RNA, DNA: Nucleic acids
//! - HOH: Water or heavy water
//! - PYR: Pyranose (from refmac dictionary)
//! - KET: Ketopyranose (from refmac dictionary)
//! - BUF: Crystallization buffer agent (from PISA agents.dat)
//! - ELS: Everything else (ligand)
enum class ResidueKind : unsigned char {
  // when changing this list update check_polymer_type()
  UNKNOWN=0, AA, AAD, PAA, MAA, RNA, DNA, BUF, HOH, PYR, KET, ELS
};

//! @brief Information about a residue type.
//!
//! Contains classification, linking info, one-letter code, and molecular weight.
struct ResidueInfo {
  char name[8];  //!< Residue name (e.g., "ALA", "GLY")
  ResidueKind kind;  //!< Residue classification
  std::uint8_t linking_type;  //!< 0=n/a, 1=peptide, 2=nucleic, 3=both
  char one_letter_code;  //!< One-letter code (uppercase if standard)
  std::uint8_t hydrogen_count;  //!< Approx H count for mass estimation
  float weight;  //!< Molecular weight

  //! @brief Check if residue info was found.
  //! @return True if kind is not UNKNOWN
  bool found() const { return kind != ResidueKind::UNKNOWN; }

  //! @brief Check if residue is water.
  //! @return True if kind is HOH
  bool is_water() const { return kind == ResidueKind::HOH; }

  //! @brief Check if residue is DNA.
  //! @return True if kind is DNA
  bool is_dna() const { return kind == ResidueKind::DNA; }

  //! @brief Check if residue is RNA.
  //! @return True if kind is RNA
  bool is_rna() const { return kind == ResidueKind::RNA; }

  //! @brief Check if residue is nucleic acid.
  //! @return True if DNA or RNA
  bool is_nucleic_acid() const { return is_dna() || is_rna(); }

  //! @brief Check if residue is amino acid.
  //! @return True if any AA type (AA, AAD, PAA, MAA)
  bool is_amino_acid() const {
    return kind == ResidueKind::AA || kind == ResidueKind::AAD ||
           kind == ResidueKind::PAA || kind == ResidueKind::MAA;
  }

  //! @brief Check if residue is buffer or water.
  //! @return True if HOH or BUF
  bool is_buffer_or_water() const {
    return kind == ResidueKind::HOH || kind == ResidueKind::BUF;
  }

  //! @brief Check if residue is standard (uppercase one-letter code).
  //!
  //! PDB format has non-standard residues (modified AA) marked as HETATM.
  //! @return True if one-letter code is uppercase
  bool is_standard() const { return (one_letter_code & 0x20) == 0; }

  //! @brief Get FASTA code.
  //! @return One-letter code if standard, otherwise 'X'
  char fasta_code() const { return is_standard() ? one_letter_code : 'X'; }

  //! @brief Check if residue is peptide linking.
  //! @return True if linking_type has bit 0 set
  bool is_peptide_linking() const { return (linking_type & 1); }

  //! @brief Check if residue is nucleic acid linking.
  //! @return True if linking_type has bit 1 set
  bool is_na_linking() const { return (linking_type & 2); }
};

//! @brief Get residue info by table index.
//! @param idx Index in internal residue table
//! @return Reference to ResidueInfo
GEMMI_DLL ResidueInfo& get_residue_info(size_t idx);

//! @brief Find index of residue in internal table.
//! @param name Residue name to search for
//! @return Table index, or unknown_tabulated_residue_idx() if not found
GEMMI_DLL size_t find_tabulated_residue_idx(const std::string& name);

//! @brief Get index for unknown residue type.
//! @return Index value representing unknown residue
constexpr size_t unknown_tabulated_residue_idx() { return 367; };

//! @brief Find residue info by name.
//! @param name Residue name (e.g., "ALA", "GLY")
//! @return Reference to ResidueInfo (unknown if not found)
GEMMI_DLL ResidueInfo& find_tabulated_residue(const std::string& name);

//! @brief Expand one-letter code to three-letter residue name.
//! @param c One-letter code (e.g., 'A', 'C', 'G')
//! @param kind Residue type (AA, RNA, or DNA)
//! @return Three-letter code (e.g., "ALA", "CYS") or nullptr if invalid
inline const char* expand_one_letter(char c, ResidueKind kind) {
  static const char* names =
    // amino-acids (all letters but J are used)
    "ALA\0ASX\0CYS\0ASP\0GLU\0PHE\0GLY\0HIS\0ILE\0\0   LYS\0LEU\0MET\0"  // A-M
    "ASN\0PYL\0PRO\0GLN\0ARG\0SER\0THR\0SEC\0VAL\0TRP\0UNK\0TYR\0GLX\0"  // N-Z
    // DNA
    "DA\0 \0\0  DC\0 \0\0  \0\0  \0\0  DG\0 \0\0  DI\0 \0\0  \0\0  \0\0  \0\0  "   // A-M
    "DN\0 \0\0  \0\0  \0\0  \0\0  \0\0  DT\0 DU\0 \0\0  \0\0  \0\0  \0\0  \0\0  "; // N-Z
  c &= ~0x20;
  const char* ret = nullptr;
  if (c >= 'A' && c <= 'Z') {
    ret = &names[4 * (c - 'A')];
    if (kind == ResidueKind::AA) {
      // ret is already set
    } else if (kind == ResidueKind::DNA) {
      ret += 4 * 26;
    } else if (kind == ResidueKind::RNA && c != 'T') {
      ret += 4 * 26 + 1;
    } else {
      ret = nullptr;
    }
  }
  return (ret && *ret) ? ret : nullptr;
}

//! @brief Expand one-letter sequence to three-letter residue names.
//! @param seq Sequence string (e.g., "ACGT" or "MKFL")
//! @param kind Residue type (AA, RNA, or DNA)
//! @return Vector of three-letter codes
GEMMI_DLL std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
                                                              ResidueKind kind);

//! @deprecated Use expand_one_letter(c, ResidueKind::AA) instead
inline const char* expand_protein_one_letter(char c) {
  return expand_one_letter(c, ResidueKind::AA);
}

//! @deprecated Use expand_one_letter_sequence() instead
GEMMI_DLL std::vector<std::string> expand_protein_one_letter_string(const std::string& s);


} // namespace gemmi
#endif
