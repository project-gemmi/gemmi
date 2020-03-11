// Copyright 2019 Global Phasing Ltd.
//
// Heuristic method to assign _atom_site.label_seq_id.
// Involves sequence alignment.

#ifndef GEMMI_LABELSEQ_HPP_
#define GEMMI_LABELSEQ_HPP_

#include "model.hpp"
#include "seqalign.hpp"  // for align_sequences
#include "polyheur.hpp"  // for are_connected3

namespace gemmi {

inline std::vector<bool> prepare_free_gapo(const ConstResidueSpan& polymer,
                                           PolymerType polymer_type) {
  std::vector<bool> gaps;
  gaps.reserve(polymer.size());
  gaps.push_back(true); // free gap opening at the beginning of sequence
  auto first_conformer = polymer.first_conformer();
  auto res = first_conformer.begin();
  for (auto next_res = res; ++next_res != first_conformer.end(); res = next_res)
    gaps.push_back(!are_connected3(*res, *next_res, polymer_type));
  return gaps;
}

// pre: !!polymer
inline AlignmentResult align_sequence_to_polymer(
                                     const std::vector<std::string>& full_seq,
                                     const ConstResidueSpan& polymer,
                                     PolymerType polymer_type,
                                     const AlignmentScoring& scoring) {
  std::map<std::string, std::uint8_t> encoding;
  for (const Residue& res : polymer)
    encoding.emplace(res.name, encoding.size());
  for (const std::string& mon_list : full_seq)
    encoding.emplace(Entity::first_mon(mon_list), encoding.size());

  std::vector<std::uint8_t> encoded_full_seq(full_seq.size());
  for (size_t i = 0; i != full_seq.size(); ++i)
    encoded_full_seq[i] = encoding.at(Entity::first_mon(full_seq[i]));

  std::vector<std::uint8_t> encoded_model_seq;
  encoded_model_seq.reserve(polymer.size());
  for (const Residue& res : polymer.first_conformer())
    encoded_model_seq.push_back(encoding.at(res.name));

  return align_sequences(encoded_full_seq, encoded_model_seq,
                         prepare_free_gapo(polymer, polymer_type),
                         encoding.size(), scoring);
}

inline bool seqid_matches_seqres(const ConstResidueSpan& polymer,
                                 const Entity& ent) {
  for (const Residue& res : polymer.first_conformer()) {
    size_t seqid = (size_t) *res.seqid.num;
    if (res.seqid.has_icode() ||
        seqid >= ent.full_sequence.size() ||
        Entity::first_mon(ent.full_sequence[seqid]) != res.name)
      return false;
  }
  return true;
}

inline void assign_label_seq_id(ResidueSpan& polymer, const Entity& ent) {
  if (seqid_matches_seqres(polymer, ent)) {
    for (Residue& res : polymer)
      res.label_seq = res.seqid.num;
    return;
  }
  AlignmentScoring scoring;
  AlignmentResult result = align_sequence_to_polymer(ent.full_sequence, polymer,
                                                     ent.polymer_type, scoring);
  auto res_group = polymer.first_conformer().begin();
  int id = 1;
  for (AlignmentResult::Item item : result.cigar) {
    switch (item.op()) {
      case 'I':
        id += item.len();
        break;
      case 'D':  // leaving label_seq as it is
        for (uint32_t i = 0; i < item.len(); ++i)
          res_group++;
        break;
      case 'M':  // not checking for mismatches
        for (uint32_t i = 0; i < item.len(); ++i, ++id)
          for (Residue* res = &*res_group++; res != &*res_group; ++res)
            res->label_seq = id;
        break;
    }
  }
}

inline void clear_label_seq_id(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        res.label_seq = Residue::OptionalNum();
}

inline void assign_label_seq_id(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains) {
      ResidueSpan polymer = chain.get_polymer();
      if (const Entity* ent = st.get_entity_of(polymer))
        assign_label_seq_id(polymer, *ent);
    }
}

} // namespace gemmi
#endif
