// Copyright 2023 Global Phasing Ltd.

#include <gemmi/align.hpp>
#include <gemmi/polyheur.hpp> // for check_polymer_type
#include <gemmi/resinfo.hpp>  // for expand_one_letter_sequence, find_tabulated_residue
#include <gemmi/seqtools.hpp> // for sequence_kind

namespace gemmi {

namespace {

using Seq = std::vector<std::string>;  // three-letter-code sequence

// Check if a residue belongs to the given polymer type.
// Returns false for cations, ligands, and other non-polymer residues
// that may be present in the chain but are not part of the sequence.
bool is_compatible_residue(const std::string& name, PolymerType ptype) {
  ResidueInfo info = find_tabulated_residue(name);
  if (!info.found())
    return false;
  if (is_polypeptide(ptype))
    return info.is_amino_acid() || info.is_peptide_linking();
  if (is_polynucleotide(ptype))
    return info.is_nucleic_acid() || info.is_na_linking();
  return false;
}

struct ScoredMatch {
  const Seq* seq = nullptr; // only better scores are accepted
  double score = -5000;
};

// Find the best matching FASTA sequence for a polymer span, filtering out
// residues that are not compatible with the polymer type (e.g., cations
// like Ca or Zn in a protein chain).
ScoredMatch find_best_matching_sequence(const ConstResidueSpan& polymer,
                                        PolymerType ptype,
                                        const std::vector<Seq>& sequences) {
  const AlignmentScoring* scoring = AlignmentScoring::partial_model();
  // Build a filtered model sequence and gap-opening penalties,
  // excluding non-compatible residues (cations, ligands, etc.)
  std::vector<std::string> filtered_seq;
  std::vector<int> target_gapo;
  const Residue* prev = nullptr;
  for (const Residue& res : polymer.first_conformer()) {
    if (!is_compatible_residue(res.name, ptype))
      continue;
    if (filtered_seq.empty()) {
      target_gapo.push_back(0);  // free gap opening at start
    } else if (prev) {
      bool connected = are_connected3(*prev, res, ptype);
      target_gapo.push_back(connected ? scoring->bad_gapo : scoring->good_gapo);
    }
    filtered_seq.push_back(res.name);
    prev = &res;
  }
  if (filtered_seq.empty())
    return {};
  target_gapo.push_back(0);  // free gap after end of chain

  ScoredMatch best;
  for (const Seq& seq : sequences) {
    AlignmentResult result = align_string_sequences(seq, filtered_seq,
                                                     target_gapo, scoring);
    if (result.score > best.score) {
      best.score = result.score;
      best.seq = &seq;
    }
  }
  return best;
}

} // anonymous namespace

void assign_best_sequences(Structure& st, const std::vector<std::string>& fasta_sequences) {
  if (st.models.empty() || fasta_sequences.empty())
    return;
  // Resolve Unknown polymer types from the chain content before alignment.
  // This handles cases where setup_entities() wasn't called or detection
  // was incomplete.
  for (Entity& ent : st.entities)
    if (ent.entity_type == EntityType::Polymer &&
        ent.polymer_type == PolymerType::Unknown)
      for (const std::string& sub : ent.subchains) {
        ConstResidueSpan span = st.models[0].get_subchain(sub);
        if (!span.empty()) {
          ent.polymer_type = check_polymer_type(span);
          break;
        }
      }
  // Fasta sequence can be protein (A->ALA), RNA (A->A) or DNA (A->DA),
  // we treat the three separately.
  std::array<bool,16> present_polymer_types{};
  for (Entity& ent : st.entities)
    if (ent.entity_type == EntityType::Polymer)
      present_polymer_types[(int)ent.polymer_type] = true;
  PolymerType ptypes[3] = { PolymerType::PeptideL, PolymerType::Rna, PolymerType::Dna };
  for (int i = 0; i < 3; ++i) {
    PolymerType ptype = ptypes[i];
    if (!present_polymer_types[(int)ptype])
      continue;
    std::vector<Seq> sequences;
    sequences.reserve(fasta_sequences.size());
    for (const std::string& seq : fasta_sequences) {
      try {
        sequences.push_back(expand_one_letter_sequence(seq, sequence_kind(ptype)));
      } catch (std::runtime_error&) {}
    }
    for (Entity& ent : st.entities) {
      // Only match entities whose polymer type is compatible
      if (ent.entity_type != EntityType::Polymer || ent.polymer_type != ptype)
        continue;
      ScoredMatch best;
      // Try all subchains to handle copies with different truncations.
      // Pick the FASTA giving the best alignment across all copies.
      for (const std::string& subchain_name : ent.subchains) {
        ConstResidueSpan polymer = st.models[0].get_subchain(subchain_name);
        if (polymer.empty())
          continue;
        ScoredMatch match = find_best_matching_sequence(polymer, ptype, sequences);
        if (match.seq && match.score > best.score)
          best = match;
      }
      if (best.seq)
        ent.full_sequence = *best.seq;
    }
  }
}

} // namespace gemmi
