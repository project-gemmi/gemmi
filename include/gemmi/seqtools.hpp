//! @file
//! @brief Sequence manipulation utilities (weight calculation, code conversion).
//!
//! Functions for working with sequences (other than alignment).

// Copyright Global Phasing Ltd.
//
// Functions for working with sequences (other than alignment).

#ifndef GEMMI_SEQTOOLS_HPP_
#define GEMMI_SEQTOOLS_HPP_

#include "resinfo.hpp"   // for find_tabulated_residue
#include "metadata.hpp"  // for Entity::first_mon, PolymerType

namespace gemmi {

//! @brief Get molecular weight of water.
//! @return Weight of H2O in daltons
constexpr double h2o_weight() { return 2 * 1.00794 + 15.9994; }

//! @brief Calculate molecular weight of a polymer sequence.
//! @param seq Vector of residue names
//! @param unknown Weight to use for unknown residues (default 100 Da)
//! @return Molecular weight in daltons, accounting for peptide/phosphodiester bonds
//!
//! Subtracts water molecules lost during condensation polymerization.
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

//! @brief Convert sequence to single-letter codes.
//! @param seq Vector of residue names (3-letter codes)
//! @return String of single-letter codes
//!
//! Uses standard FASTA single-letter codes for amino acids and nucleotides.
inline std::string one_letter_code(const std::vector<std::string>& seq) {
  std::string r;
  for (const std::string& item : seq)
    r += find_tabulated_residue(Entity::first_mon(item)).fasta_code();
  return r;
}

//! @brief Convert sequence to PDBx one-letter code format.
//! @param seq Vector of residue names
//! @param kind Expected residue type (AA, DNA, RNA)
//! @return String with standard codes and non-standard in parentheses
//!
//! Returns the format used in _entity_poly.pdbx_seq_one_letter_code,
//! in which non-standard amino acids/nucleotides are represented by CCD codes
//! in parenthesis, e.g. AA(MSE)H.
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

//! @brief Determine residue kind from polymer type.
//! @param ptype Polymer type (e.g., PeptideL, Dna, Rna)
//! @return Corresponding residue kind (AA, DNA, or RNA)
//! @throws Error if polymer type is Unknown
//!
//! Used with expand_one_letter_sequence()
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
