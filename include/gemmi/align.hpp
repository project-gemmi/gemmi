//! @file
//! @brief Sequence alignment, label_seq_id assignment, and structure superposition.
//!
//! Provides functions for aligning sequences to coordinate models, assigning
//! label_seq_id based on SEQRES, and calculating structure superpositions using
//! QCP algorithm with optional iterative trimming.

// Copyright 2020 Global Phasing Ltd.

#ifndef GEMMI_ALIGN_HPP_
#define GEMMI_ALIGN_HPP_

#include "model.hpp"
#include "seqalign.hpp"  // for align_sequences
#include "qcp.hpp"       // for superpose_positions
#include "polyheur.hpp"  // for are_connected3

namespace gemmi {

// Sequence alignment and label_seq_id assignment

//! @brief Prepare gap opening penalties for sequence alignment.
//! @param polymer Polymer chain from coordinate model
//! @param polymer_type Type of polymer (peptide, DNA, RNA)
//! @param scoring Alignment scoring parameters (nullptr for default)
//! @return Vector of gap opening penalties for each position
//!
//! Helper function for sequence alignment. Returns position-specific gap opening
//! penalties based on connectivity between residues. Gaps are penalized less
//! where residues are not connected (chain breaks).
inline std::vector<int> prepare_target_gapo(const ConstResidueSpan& polymer,
                                            PolymerType polymer_type,
                                            const AlignmentScoring* scoring=nullptr) {
  if (!scoring)
    scoring = AlignmentScoring::partial_model();
  std::vector<int> gaps;
  gaps.reserve(polymer.size());
  gaps.push_back(0); // free gap opening at the beginning of sequence
  if (is_polypeptide(polymer_type) || is_polynucleotide(polymer_type)) {
    auto first_conformer = polymer.first_conformer();
    auto res = first_conformer.begin();
    for (auto next_res = res; ++next_res != first_conformer.end(); res = next_res) {
      bool connected = are_connected3(*res, *next_res, polymer_type);
      gaps.push_back(connected ? scoring->bad_gapo : scoring->good_gapo);
    }
    gaps.push_back(0); // free gap after the end of chain
  }
  return gaps;
}

//! @brief Align a full sequence (SEQRES) to a polymer chain from coordinates.
//! @param full_seq Full sequence from SEQRES or entity
//! @param polymer Polymer chain from coordinate model
//! @param polymer_type Type of polymer (peptide, DNA, RNA)
//! @param scoring Alignment scoring parameters (nullptr for default)
//! @return Alignment result with CIGAR string
//!
//! Performs sequence alignment between the full sequence (from SEQRES or entity)
//! and the residues present in the coordinate model. Uses position-specific gap
//! penalties to account for chain breaks in the model.
inline AlignmentResult align_sequence_to_polymer(
                                     const std::vector<std::string>& full_seq,
                                     const ConstResidueSpan& polymer,
                                     PolymerType polymer_type,
                                     const AlignmentScoring* scoring=nullptr) {
  if (!polymer)
    return AlignmentResult();
  std::map<std::string, std::uint8_t> encoding;
  if (!scoring)
    scoring = AlignmentScoring::partial_model();
  for (const std::string& res_name : scoring->matrix_encoding)
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
                         prepare_target_gapo(polymer, polymer_type, scoring),
                         (std::uint8_t)encoding.size(), *scoring);
}

//! @brief Check for exact match between model sequence and full sequence (SEQRES).
//! @param polymer Polymer chain from coordinate model
//! @param ent Entity with full sequence
//! @return True if sequences match exactly and seqid numbering is sequential
//!
//! Fast path for the common case where model sequence exactly matches SEQRES
//! and residue numbering is sequential (1, 2, 3, ...) with no insertion codes.
inline bool seqid_matches_seqres(const ConstResidueSpan& polymer,
                                 const Entity& ent) {
  if (ent.full_sequence.size() != polymer.size())
    return false;
  int idx = 0;
  for (const Residue& res : polymer) {
    if (ent.full_sequence[idx] != res.name ||
        ++idx != *res.seqid.num || res.seqid.has_icode())
      return false;
  }
  return true;
}

//! @brief Clear all sequences from structure entities.
//! @param st Structure to modify
//!
//! Removes full_sequence, database references, and SIFTS UniProt accessions
//! from all entities in the structure.
inline void clear_sequences(Structure& st) {
  for (Entity& ent : st.entities) {
    ent.full_sequence.clear();
    ent.dbrefs.clear();
    ent.sifts_unp_acc.clear();
  }
}

//! @brief Assign best matching sequences from FASTA to structure entities.
//! @param st Structure to populate with sequences
//! @param fasta_sequences Vector of FASTA sequence strings
//!
//! Aligns provided FASTA sequences to the coordinate models and assigns the
//! best matches to each entity's full_sequence.
GEMMI_DLL
void assign_best_sequences(Structure& st, const std::vector<std::string>& fasta_sequences);

//! @brief Assign label_seq_id to polymer using sequence alignment.
//! @param polymer Polymer chain to assign label_seq_id
//! @param ent Entity with full sequence (may be nullptr)
//! @param force If true, assign label_seq even if full sequence unknown (assumes no gaps)
//!
//! Uses sequence alignment (model to SEQRES) to assign label_seq.
//! The force parameter allows assignment even without full sequence, assuming no gaps.
inline void assign_label_seq_to_polymer(ResidueSpan& polymer,
                                        const Entity* ent, bool force) {
  AlignmentResult result;

  // sequence not known
  if (!ent || ent->full_sequence.empty()) {
    if (!force)
      return;
    PolymerType ptype = get_or_check_polymer_type(ent, polymer);
    const Residue* prev = nullptr;
    for (const Residue& res : polymer.first_conformer()) {
      if (prev && !are_connected3(*prev, res, ptype))
        result.push_cigar(1, 1);  // assume a single insertion
      result.push_cigar(0, 1);
      prev = &res;
    }

  // exact match - common case that doesn't require alignment
  } else if (seqid_matches_seqres(polymer, *ent)) {
    result.push_cigar(0, (int)ent->full_sequence.size());

  // sequence alignment
  } else {
    PolymerType ptype = get_or_check_polymer_type(ent, polymer);
    result = align_sequence_to_polymer(ent->full_sequence, polymer, ptype);
  }

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

//! @brief Clear all label_seq_id values from structure.
//! @param st Structure to modify
//!
//! Resets label_seq to unset for all residues in all models.
inline void clear_label_seq_id(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        res.label_seq = Residue::OptionalNum();
}

//! @brief Assign label_seq_id to all polymers in structure.
//! @param st Structure to process
//! @param force If true, assign even without full sequence (assumes no gaps)
//!
//! Iterates over all polymer chains in all models and assigns label_seq_id
//! based on sequence alignment to SEQRES. Only processes chains where
//! label_seq is not already assigned at both ends.
inline void assign_label_seq_id(Structure& st, bool force) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      if (ResidueSpan polymer = chain.get_polymer())
        if (!polymer.front().label_seq || !polymer.back().label_seq) {
          const Entity* ent = st.get_entity_of(polymer);
          assign_label_seq_to_polymer(polymer, ent, force);
        }
}


// superposition

//! @brief Atom selection for structure superposition.
enum class SupSelect {
  CaP,       //!< Only CA (amino acids) or P (nucleotides) atoms
  MainChain, //!< Only main chain atoms (N,CA,C,O for peptides; P,O5',C5',C4',C3',O3' for nucleotides)
  All        //!< All atoms (matched by name and element)
};

//! @brief Prepare position vectors for superposition calculation.
//! @param pos1 Output: positions from fixed structure
//! @param pos2 Output: positions from movable structure
//! @param fixed Fixed polymer chain
//! @param movable Movable polymer chain
//! @param ptype Polymer type
//! @param sel Atom selection strategy
//! @param altloc Alternate location to use ('\0' for default)
//! @param ca_offsets Optional output: indices of CA/P positions in pos1/pos2
//!
//! Aligns sequences and extracts matching atom positions from both structures.
//! Only atoms present in both structures at aligned positions are included.
inline void prepare_positions_for_superposition(std::vector<Position>& pos1,
                                                std::vector<Position>& pos2,
                                                ConstResidueSpan fixed,
                                                ConstResidueSpan movable,
                                                PolymerType ptype,
                                                SupSelect sel,
                                                char altloc='\0',
                                                std::vector<int>* ca_offsets=nullptr) {
  AlignmentResult result = align_sequence_to_polymer(fixed.extract_sequence(),
                                                     movable, ptype,
                                                     AlignmentScoring::blosum62());
  auto it1 = fixed.first_conformer().begin();
  auto it2 = movable.first_conformer().begin();
  std::vector<AtomNameElement> used_atoms;
  bool is_na = is_polynucleotide(ptype);
  const AtomNameElement* ca_p = nullptr;
  if (sel == SupSelect::CaP) {
    used_atoms.push_back({is_na ? "P" : "CA", is_na ? El::P : El::C});
  } else if (sel == SupSelect::MainChain) {
    used_atoms = get_mainchain_atoms(ptype);
    ca_p = &used_atoms[is_na ? 0 : 1];
  }
  for (AlignmentResult::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      int ca_offset = -1;
      if (op == 'M' && it1->name == it2->name) {
        if (!used_atoms.empty()) {
          for (const AtomNameElement& ane : used_atoms) {
            const Atom* a1 = it1->find_atom(ane.atom_name, altloc, ane.el);
            const Atom* a2 = it2->find_atom(ane.atom_name, altloc, ane.el);
            if (a1 && a2) {
              if (&ane == ca_p)
                ca_offset = (int)pos1.size();
              pos1.push_back(a1->pos);
              pos2.push_back(a2->pos);
            }
          }
        } else {
          for (const Atom& a1 : it1->atoms)
            if (a1.altloc_matches(altloc))
              if (const Atom* a2 = it2->find_atom(a1.name, altloc, a1.element)) {
                pos1.push_back(a1.pos);
                pos2.push_back(a2->pos);
              }
        }
      }
      if (op == 'M' || op == 'I') {
        ++it1;
        if (ca_offsets)
          ca_offsets->push_back(ca_offset);
      }
      if (op == 'M' || op == 'D')
        ++it2;
    }
  }
}

//! @brief Calculate current RMSD between two structures without superposition.
//! @param fixed Fixed polymer chain
//! @param movable Movable polymer chain
//! @param ptype Polymer type
//! @param sel Atom selection strategy
//! @param altloc Alternate location to use ('\0' for default)
//! @return Result with RMSD and atom count
//!
//! Calculates RMSD between aligned structures in their current positions,
//! without performing any superposition. Useful for checking RMSD after
//! applying a known transformation.
inline SupResult calculate_current_rmsd(ConstResidueSpan fixed,
                                        ConstResidueSpan movable,
                                        PolymerType ptype,
                                        SupSelect sel,
                                        char altloc='\0') {
  std::vector<Position> pos1, pos2;
  prepare_positions_for_superposition(pos1, pos2, fixed, movable, ptype, sel, altloc);
  SupResult r;
  r.count = pos1.size();
  double sd = 0;
  for (size_t i = 0; i != pos1.size(); ++i)
    sd += pos1[i].dist_sq(pos2[i]);
  r.rmsd = std::sqrt(sd / r.count);
  return r;
}

//! @brief Calculate optimal superposition between two polymer chains.
//! @param fixed Fixed polymer chain (reference)
//! @param movable Movable polymer chain (to be transformed)
//! @param ptype Polymer type
//! @param sel Atom selection strategy
//! @param trim_cycles Number of iterative trimming cycles (0 for no trimming)
//! @param trim_cutoff RMSD multiplier for trimming threshold
//! @param altloc Alternate location to use ('\0' for default)
//! @return Superposition result with transformation and RMSD
//!
//! Calculates optimal superposition using QCP algorithm. With trim_cycles > 0,
//! iteratively removes outlier atom pairs (distance > trim_cutoff × RMSD) and
//! recalculates superposition. Useful for improving alignment quality.
inline SupResult calculate_superposition(ConstResidueSpan fixed,
                                         ConstResidueSpan movable,
                                         PolymerType ptype,
                                         SupSelect sel,
                                         int trim_cycles=0,
                                         double trim_cutoff=2.0,
                                         char altloc='\0') {
  std::vector<Position> pos1, pos2;
  prepare_positions_for_superposition(pos1, pos2, fixed, movable, ptype, sel, altloc);
  const double* weights = nullptr;
  size_t len = pos1.size();
  SupResult sr = superpose_positions(pos1.data(), pos2.data(), len, weights);

  for (int n = 0; n < trim_cycles; ++n) {
    double max_dist_sq = sq(trim_cutoff * sr.rmsd);
    size_t p = 0;
    for (size_t i = 0; i != len; ++i) {
      Vec3 m2 = sr.transform.apply(pos2[i]);
      if (m2.dist_sq(pos1[i]) <= max_dist_sq) {
        if (i != p) {
          pos1[p] = pos1[i];
          pos2[p] = pos2[i];
        }
        ++p;
      }
    }
    if (p == len)
      break;
    len = p;
    if (len < 3)
      fail("in calculate_superposition(): only ", std::to_string(len),
           " atoms after trimming");
    sr = superpose_positions(pos1.data(), pos2.data(), len, weights);
  }

  return sr;
}

//! @brief Calculate local superpositions in moving window around each residue.
//! @param fixed Fixed polymer chain (reference)
//! @param movable Movable polymer chain (to be transformed)
//! @param ptype Polymer type
//! @param radius Window radius around each CA/P atom (Angstroms)
//! @return Vector of superposition results, one per residue in fixed
//!
//! Returns superpositions for all residues in fixed.first_conformer(),
//! performed by superposing backbone in radius=10.0 from residue's Ca.
//! Useful for analyzing local structural flexibility or identifying rigid regions.
inline std::vector<SupResult> calculate_superpositions_in_moving_window(
                                      ConstResidueSpan fixed,
                                      ConstResidueSpan movable,
                                      PolymerType ptype,
                                      double radius=10.0) {
  const double radius_sq = radius * radius;
  std::vector<Position> pos1, pos2;
  char altloc = '\0';
  SupSelect sel = SupSelect::MainChain;
  std::vector<int> ca_offsets;
  prepare_positions_for_superposition(pos1, pos2, fixed, movable, ptype,
                                      sel, altloc, &ca_offsets);
  const double* weights = nullptr;
  std::vector<SupResult> result;
  for (int offset : ca_offsets) {
    if (offset == -1) {
      result.push_back(SupResult{NAN, 0, {}, {}, {}});
      continue;
    }
    const Position& ca_pos = pos1[offset];
    int a = offset;
    while (a > 0 && ca_pos.dist_sq(pos1[a-1]) < radius_sq)
      --a;
    int b = offset;
    while (b+1 < (int)pos1.size() && ca_pos.dist_sq(pos1[b+1]) < radius_sq)
      ++b;
    result.push_back(superpose_positions(&pos1[a], &pos2[a], b-a+1, weights));
  }
  return result;
}

} // namespace gemmi
#endif
