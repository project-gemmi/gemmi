// Copyright 2023 Global Phasing Ltd.

#include <gemmi/align.hpp>
#include <gemmi/resinfo.hpp>  // for expand_one_letter_sequence
#include <gemmi/seqtools.hpp> // for sequence_kind

namespace gemmi {

namespace {

using Seq = std::vector<std::string>;  // three-letter-code sequence

const Seq* find_best_matching_sequence(const ConstResidueSpan& polymer,
                                       PolymerType ptype,
                                       const std::vector<Seq>& sequences) {
  double best_score = -5000; // only better scores are accepted
  const Seq* best_seq = nullptr;
  for (const Seq& seq : sequences) {
    AlignmentResult result = align_sequence_to_polymer(seq, polymer, ptype);
    if (result.score > best_score) {
      best_score = result.score;
      best_seq = &seq;
    }
  }
  return best_seq;
}

} // anonymous namespace

void assign_best_sequences(Structure& st, const std::vector<std::string>& fasta_sequences) {
  if (st.models.empty() && fasta_sequences.empty())
    return;
  // Fasta sequence can be protein (A->ALA), RNA (A->A) or DNA (A->DA),
  // we treat the three separately.
  std::array<bool,16> present_polymer_types{};
  for (Entity& ent : st.entities)
    if (ent.entity_type == EntityType::Polymer)
      present_polymer_types[(int)ent.polymer_type] = true;
  PolymerType ptypes[3] = { PolymerType::PeptideL, PolymerType::Rna, PolymerType::Dna };
  for (int i = 0; i < 3; ++i) {
    PolymerType ptype = ptypes[i];
    if (present_polymer_types[(int)ptype]) {
      std::vector<Seq> sequences;
      sequences.reserve(fasta_sequences.size());
      for (const std::string& seq : fasta_sequences) {
        try {
          sequences.push_back(expand_one_letter_sequence(seq, sequence_kind(ptype)));
        } catch (std::runtime_error&) {}
      }
      for (Entity& ent : st.entities) {
        if (ent.entity_type == EntityType::Polymer)
          for (std::string& subchain_name : ent.subchains) {
            ConstResidueSpan polymer = st.models[0].get_subchain(subchain_name);
            const Seq* best_seq = find_best_matching_sequence(polymer, ptype, sequences);
            if (best_seq) {
              ent.full_sequence = *best_seq;
              break;
            }
          }
      }
    }
  }
}

} // namespace gemmi
