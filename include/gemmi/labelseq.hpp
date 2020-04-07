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

inline AlignmentScoring prepare_blosum62_scoring() {
  AlignmentScoring s;
  s.match = 1;
  s.mismatch = -4;
  s.gapo = -10;  // BLAST uses BLOSUM-62 with gap cost (10,1)
  s.gape = -1;
  s.score_matrix = {
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
   -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
   -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
   -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
   -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
   -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
   -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
   -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,
  };
  s.matrix_encoding = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  };
  return s;
}

inline std::vector<bool> prepare_free_gapo(const ConstResidueSpan& polymer,
                                           PolymerType polymer_type) {
  std::vector<bool> gaps;
  gaps.reserve(polymer.size());
  gaps.push_back(true); // free gap opening at the beginning of sequence
  if (!is_polypeptide(polymer_type) && !is_polynucleotide(polymer_type))
    return gaps;
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
  for (const std::string& res_name : scoring.matrix_encoding)
    encoding.emplace(res_name, (std::uint8_t)encoding.size());
  for (const Residue& res : polymer)
    encoding.emplace(res.name, (std::uint8_t)encoding.size());
  for (const std::string& mon_list : full_seq)
    encoding.emplace(Entity::first_mon(mon_list), (std::uint8_t)encoding.size());
  if (encoding.size() > 255)
    return AlignmentResult();

  std::vector<std::uint8_t> encoded_full_seq(full_seq.size());
  for (size_t i = 0; i != full_seq.size(); ++i)
    encoded_full_seq[i] = encoding.at(Entity::first_mon(full_seq[i]));

  std::vector<std::uint8_t> encoded_model_seq;
  encoded_model_seq.reserve(polymer.size());
  for (const Residue& res : polymer.first_conformer())
    encoded_model_seq.push_back(encoding.at(res.name));

  return align_sequences(encoded_full_seq, encoded_model_seq,
                         prepare_free_gapo(polymer, polymer_type),
                         (std::uint8_t)encoding.size(), scoring);
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

inline void assign_label_seq_id(ResidueSpan& polymer, const Entity* ent) {
  // sequence not known
  if (!ent || ent->full_sequence.empty()) {
    int n = 1;
    SeqId prev;
    for (Residue& res : polymer) {
      res.label_seq = n;
      if (prev != res.seqid)
        ++n;
      prev = res.seqid;
    }
    return;
  }

  // exact match - common case that doesn't require alignment
  if (seqid_matches_seqres(polymer, *ent)) {
    for (Residue& res : polymer)
      res.label_seq = res.seqid.num;
    return;
  }

  // sequence alignment
  AlignmentScoring scoring;
  AlignmentResult result =
    align_sequence_to_polymer(ent->full_sequence, polymer, ent->polymer_type,
                              scoring);
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

inline void assign_label_seq_id(Structure& st, bool force) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      if (ResidueSpan polymer = chain.get_polymer())
        if (force || !polymer.front().label_seq || !polymer.back().label_seq) {
          const Entity* ent = st.get_entity_of(polymer);
          assign_label_seq_id(polymer, ent);
        }
}

inline void setup_for_mmcif(Structure& st) {
  setup_entities(st);
  assign_label_seq_id(st, false);
}

} // namespace gemmi
#endif
