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

struct ResidueInfo {
  // Simple approximate classification.
  // AA - aminoacid
  // AAD - D-aminoacid
  // PAA - proline-like aminoacid
  // MAA - methylated aminoacid
  // RNA, DNA - nucleic acids
  // HOH - water or heavy water (OH, H3O, D3O are not included here)
  // PYR - pyranose according to the refmac dictionary
  // KET - ketopyranose according to the refmac dictionary
  // BUF - agent from crystallization buffer according to PISA agents.dat
  // ELS - something else (ligand).
  enum Kind : unsigned char {
    // when changing this list update check_polymer_type()
    UNKNOWN=0, AA, AAD, PAA, MAA, RNA, DNA, BUF, HOH, PYR, KET, ELS
  };
  Kind kind;
  // linking type: 0=n/a, 1=peptide-linking, 2=nucl.-linking, 3=(2|1)
  std::uint8_t linking_type;
  // one-letter code or space (uppercase iff it is a standard residues)
  char one_letter_code;
  // rough count of hydrogens used to estimate mass with implicit hydrogens
  std::uint8_t hydrogen_count;
  // molecular weight
  float weight;

  bool found() const { return kind != UNKNOWN; }
  bool is_water() const { return kind == HOH; }
  bool is_dna() const { return kind == DNA; }
  bool is_rna() const { return kind == RNA; }
  bool is_nucleic_acid() const { return is_dna() || is_rna(); }
  bool is_amino_acid() const {
    return kind == AA || kind == AAD || kind == PAA || kind == MAA;
  }
  bool is_buffer_or_water() const { return kind == HOH || kind == BUF; }
  // PDB format has non-standard residues (modified AA) marked as HETATM.
  bool is_standard() const { return (one_letter_code & 0x20) == 0; }
  char fasta_code() const { return is_standard() ? one_letter_code : 'X'; }
  bool is_peptide_linking() const { return (linking_type & 1); }
  bool is_na_linking() const { return (linking_type & 2); }
};

GEMMI_DLL ResidueInfo find_tabulated_residue(const std::string& name);

inline const char* expand_protein_one_letter(char c) {
  static const char* data =
    "ALA\0ASX\0CYS\0ASP\0GLU\0PHE\0GLY\0HIS\0ILE\0\0   LYS\0LEU\0MET\0" // A-M
    "ASN\0PYL\0PRO\0GLN\0ARG\0SER\0THR\0SEC\0VAL\0TRP\0UNK\0TYR\0GLX";  // N-Z
  c &= ~0x20;
  if (c < 'A' || c > 'Z' || c == 'J')
    return nullptr;
  return &data[4 * (c - 'A')];
}

inline std::vector<std::string> expand_protein_one_letter_string(const std::string& s) {
  std::vector<std::string> r;
  r.reserve(s.size());
  for (char c : s) {
    const char* three_letters = expand_protein_one_letter(c);
    if (!three_letters)
      fail("unexpected letter in protein sequence: ", c);
    r.emplace_back(three_letters, 3);
  }
  return r;
}

} // namespace gemmi
#endif
