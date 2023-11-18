// Copyright Global Phasing Ltd.
//
// Functions for working with sequences (other than alignment).

#ifndef GEMMI_SEQTOOLS_HPP_
#define GEMMI_SEQTOOLS_HPP_

#include "resinfo.hpp"   // for find_tabulated_residue
#include "metadata.hpp"  // for Entity::first_mon, PolymerType

namespace gemmi {

constexpr double h2o_weight() { return 2 * 1.00794 + 15.9994; }

inline double calculate_sequence_weight(const std::vector<std::string>& seq,
                                        double unknown=0.) {
  double weight = 0.;
  for (const std::string& item : seq) {
    ResidueInfo res_info = find_tabulated_residue(Entity::first_mon(item));
    weight += res_info.found() ? res_info.weight : unknown;
  }
  return weight - (seq.size() - 1) * h2o_weight();
}

inline std::string one_letter_code(const std::vector<std::string>& seq) {
  std::string r;
  for (const std::string& item : seq)
    r += find_tabulated_residue(Entity::first_mon(item)).fasta_code();
  return r;
}

/// Returns the format used in _entity_poly.pdbx_seq_one_letter_code,
/// in which non-standard amino acids/nucleotides are represented by CCD codes
/// in parenthesis, e.g. AA(MSE)H.
inline std::string pdbx_one_letter_code(const std::vector<std::string>& seq,
                                        ResidueKind kind) {
  std::string r;
  for (const std::string& item : seq) {
    std::string code = Entity::first_mon(item);
    ResidueInfo ri = find_tabulated_residue(code);
    if (ri.is_standard() && ri.kind == kind)
      r += ri.one_letter_code;
    else
      cat_to(r, '(', code, ')');
  }
  return r;
}

/// used with expand_one_letter_sequence()
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
