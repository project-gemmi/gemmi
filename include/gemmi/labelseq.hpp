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

struct AlignmentScoring {
  int match = 1;
  int mismatch = -1;
  int gapo = 1;
  int gape = 1;
};

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
inline Alignment align_polymer(const ConstResidueSpan& polymer,
                               const Entity& ent,
                               const AlignmentScoring& scoring) {
  std::map<std::string, std::uint8_t> encoding;
  for (const Residue& res : polymer)
    encoding.emplace(res.name, 0);
  for (const std::string& mon_list : ent.full_sequence)
    encoding.emplace(Entity::first_mon(mon_list), 0);
  std::vector<std::string> all_monomers;
  size_t n_mon = encoding.size();
  all_monomers.reserve(n_mon);
  for (auto& item : encoding) {
    item.second = (std::uint8_t) all_monomers.size();
    all_monomers.push_back(item.first);
  }
  std::vector<std::uint8_t> model_seq;
  model_seq.reserve(polymer.size());
  auto first_conformer = polymer.first_conformer();
  for (const Residue& res : first_conformer)
    model_seq.push_back(encoding.at(res.name));
  std::vector<std::uint8_t> full_seq;
  full_seq.reserve(ent.full_sequence.size());
  for (const std::string& mon_list : ent.full_sequence)
    full_seq.push_back(encoding.at(Entity::first_mon(mon_list)));
  std::vector<std::int8_t> score_matrix(n_mon * n_mon, scoring.mismatch);
  for (size_t i = 0; i != n_mon; ++i)
    score_matrix[i * n_mon + i] = scoring.match;
  return align_sequences(full_seq.size(), full_seq.data(),
                         model_seq.size(), model_seq.data(),
                         prepare_free_gapo(polymer, ent.polymer_type),
                         n_mon, score_matrix.data(),
                         scoring.gapo, scoring.gape);
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
  Alignment result = align_polymer(polymer, ent, scoring);
  auto res_group = polymer.first_conformer().begin();
  int id = 1;
  for (Alignment::Item item : result.cigar) {
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
