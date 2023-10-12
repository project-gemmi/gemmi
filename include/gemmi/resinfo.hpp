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
enum class ResidueKind : unsigned char {
  // when changing this list update check_polymer_type()
  UNKNOWN=0, AA, AAD, PAA, MAA, RNA, DNA, BUF, HOH, PYR, KET, ELS
};

struct ResidueInfo {
  ResidueKind kind;
  // linking type: 0=n/a, 1=peptide-linking, 2=nucl.-linking, 3=(2|1)
  std::uint8_t linking_type;
  // one-letter code or space (uppercase iff it is a standard residues)
  char one_letter_code;
  // rough count of hydrogens used to estimate mass with implicit hydrogens
  std::uint8_t hydrogen_count;
  // molecular weight
  float weight;

  bool found() const { return kind != ResidueKind::UNKNOWN; }
  bool is_water() const { return kind == ResidueKind::HOH; }
  bool is_dna() const { return kind == ResidueKind::DNA; }
  bool is_rna() const { return kind == ResidueKind::RNA; }
  bool is_nucleic_acid() const { return is_dna() || is_rna(); }
  bool is_amino_acid() const {
    return kind == ResidueKind::AA || kind == ResidueKind::AAD ||
           kind == ResidueKind::PAA || kind == ResidueKind::MAA;
  }
  bool is_buffer_or_water() const {
    return kind == ResidueKind::HOH || kind == ResidueKind::BUF;
  }
  // PDB format has non-standard residues (modified AA) marked as HETATM.
  bool is_standard() const { return (one_letter_code & 0x20) == 0; }
  char fasta_code() const { return is_standard() ? one_letter_code : 'X'; }
  bool is_peptide_linking() const { return (linking_type & 1); }
  bool is_na_linking() const { return (linking_type & 2); }
};

GEMMI_DLL ResidueInfo find_tabulated_residue(const std::string& name);

/// kind can be AA, RNA or DNA
inline const char* expand_one_letter(char c, ResidueKind kind) {
  static const char* names =
    // amino-acids (all letters but J are used)
    "ALA\0ASX\0CYS\0ASP\0GLU\0PHE\0GLY\0HIS\0ILE\0\0   LYS\0LEU\0MET\0"  // A-M
    "ASN\0PYL\0PRO\0GLN\0ARG\0SER\0THR\0SEC\0VAL\0TRP\0UNK\0TYR\0GLX\0"  // N-Z
    // DNA
    "DA\0 \0b  DC\0 \0d  \0e  \0f  DG\0 \0h  DI\0 \0j  \0k  \0l  \0m  "  // A-M
    "DN\0 \0o  \0p  \0q  \0r  \0s  DT\0 DU\0 \0v  \0w  \0x  \0y  \0z  "; // N-Z
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

/// kind can be AA, RNA or DNA
inline std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
                                                           ResidueKind kind) {
  std::vector<std::string> r;
  r.reserve(seq.size());
  auto kind_str = [&]() {
    switch (kind) {
      case ResidueKind::AA: return "peptide";
      case ResidueKind::RNA: return "RNA";
      case ResidueKind::DNA: return "DNA";
      default: return "unknown";
    }
  };
  for (size_t i = 0; i != seq.size(); ++i) {
    char c = seq[i];
    if (c == '(') { // special case, e.g. (MSE)
      size_t start = i + 1;
      i = seq.find(')', start);
      if (i == std::string::npos)
        fail("unmatched '(' in sequence");
      r.emplace_back(seq, start, i - start);
    } else {
      const char* str = expand_one_letter(c, kind);
      if (str == nullptr)
        fail("unexpected letter in ", kind_str(), " sequence: ", c);
      r.emplace_back(str);
    }
  }
  return r;
}

// deprecated
inline const char* expand_protein_one_letter(char c) {
  return expand_one_letter(c, ResidueKind::AA);
}
// deprecated
inline std::vector<std::string> expand_protein_one_letter_string(const std::string& s) {
  return expand_one_letter_sequence(s, ResidueKind::AA);
}


} // namespace gemmi
#endif
