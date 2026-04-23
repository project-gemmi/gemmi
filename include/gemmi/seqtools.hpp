// Copyright Global Phasing Ltd.
//
// Functions for working with sequences (other than alignment).

#ifndef GEMMI_SEQTOOLS_HPP_
#define GEMMI_SEQTOOLS_HPP_

#include "resinfo.hpp"   // for find_tabulated_residue
#include "metadata.hpp"  // for Entity::first_mon, PolymerType

namespace gemmi {

/// @brief Get the molecular weight of water (H2O).
/// @return Water weight in atomic mass units
constexpr double h2o_weight() { return 2 * 1.00794 + 15.9994; }

/// @brief Calculate the molecular weight of a polymer sequence.
/// @details Sums the weights of individual residues and subtracts water molecules
/// for the peptide/nucleic acid bonds formed.
/// @param seq Vector of residue names (3-letter codes)
/// @param unknown Weight to use for unknown residues (default 100.0)
/// @return Molecular weight in atomic mass units
inline double calculate_sequence_weight(const std::vector<std::string>& seq,
                                        double unknown=100.) {
  double weight = 0.;
  for (const std::string& item : seq) {
    size_t idx = find_tabulated_residue_idx(Entity::first_mon(item));
    if (idx == unknown_tabulated_residue_idx())
      weight += unknown;
    else
      weight += get_residue_info(idx).weight;
    //ResidueInfo res_info = find_tabulated_residue(Entity::first_mon(item));
    //weight += res_info.found() ? res_info.weight : unknown;
  }
  return weight - (seq.size() - 1) * h2o_weight();
}

/// @brief Convert a sequence to single-letter FASTA code.
/// @param seq Vector of residue names (3-letter codes)
/// @return String of single-letter codes (X for unknown residues)
inline std::string one_letter_code(const std::vector<std::string>& seq) {
  std::string r;
  for (const std::string& item : seq)
    r += find_tabulated_residue(Entity::first_mon(item)).fasta_code();
  return r;
}

/// @brief Convert sequence to PDBx format with non-standard residues in parentheses.
/// @details Returns the format used in _entity_poly.pdbx_seq_one_letter_code,
/// in which non-standard amino acids/nucleotides are represented by CCD codes
/// in parentheses, e.g. AA(MSE)H.
/// @param seq Vector of residue names (3-letter codes)
/// @param kind Type of residue (AA, DNA, RNA) to filter
/// @return String with single-letter codes for standard residues and (CCD) for non-standard
inline std::string pdbx_one_letter_code(const std::vector<std::string>& seq,
                                        ResidueKind kind) {
  std::string r;
  for (const std::string& item : seq) {
    std::string code = Entity::first_mon(item);
    const ResidueInfo ri = find_tabulated_residue(code);
    if (ri.is_standard() && ri.kind == kind)
      r += ri.one_letter_code;
    else
      cat_to(r, '(', code, ')');
  }
  return r;
}

/// @brief Convert polymer type to residue kind.
/// @details Used with expand_one_letter_sequence() to determine which
/// single-letter codes to expect for a polymer.
/// @param ptype Polymer type
/// @return Residue kind (AA, DNA, or RNA)
/// @throws Fails with error if polymer type is Unknown
inline ResidueKind sequence_kind(PolymerType ptype) {
  if (is_polypeptide(ptype))
    return ResidueKind::AA;
  if (ptype == PolymerType::Dna)
    return ResidueKind::DNA;
  if (ptype == PolymerType::Rna || ptype == PolymerType::DnaRnaHybrid)
    return ResidueKind::RNA;
  if (ptype == PolymerType::Unknown)
    fail("sequence_kind(): unknown polymer type");
  return ResidueKind::AA;
}

} // namespace gemmi
#endif
