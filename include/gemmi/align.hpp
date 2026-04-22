// Copyright 2020 Global Phasing Ltd.
//
// Sequence alignment, label_seq_id assignment, structure superposition.

#ifndef GEMMI_ALIGN_HPP_
#define GEMMI_ALIGN_HPP_

#include "model.hpp"
#include "seqalign.hpp"  // for align_sequences
#include "qcp.hpp"       // for superpose_positions
#include "polyheur.hpp"  // for are_connected3

namespace gemmi {

/// @name Sequence Alignment and label_seq_id Assignment
/// @{

/// @brief Compute per-position gap-open penalties for sequence alignment.
/// Accounts for gaps and discontinuities in the observed polymer sequence.
/// @param polymer Polymer span to analyze
/// @param polymer_type Type of polymer (peptide or nucleic acid)
/// @param scoring Optional AlignmentScoring parameters; uses default partial_model() if nullptr
/// @return Vector of gap-opening penalties for each position
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

/// @brief Align a full sequence (SEQRES/Entity) to observed residues in a polymer.
/// @param full_seq Full sequence from SEQRES or Entity
/// @param polymer Polymer chain span (observed residues)
/// @param polymer_type Type of polymer (peptide or nucleic acid)
/// @param scoring Optional AlignmentScoring parameters; uses default partial_model() if nullptr
/// @return AlignmentResult with CIGAR string describing the alignment
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

/// @brief Check if model sequence exactly matches SEQRES numbering (fast path).
/// @param polymer Polymer chain span
/// @param ent Entity with full_sequence from SEQRES
/// @return True if label_seq values already match Entity::full_sequence exactly
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

/// @brief Clear full_sequence and database references from all entities.
/// @param st Structure to clear
inline void clear_sequences(Structure& st) {
  for (Entity& ent : st.entities) {
    ent.full_sequence.clear();
    ent.dbrefs.clear();
    ent.sifts_unp_acc.clear();
  }
}

/// @brief Match FASTA sequences to entities and set Entity::full_sequence.
/// @param st Structure to update
/// @param fasta_sequences Vector of FASTA sequence strings
GEMMI_DLL
void assign_best_sequences(Structure& st, const std::vector<std::string>& fasta_sequences);

/// @brief Assign label_seq_id to residues using alignment to Entity full_sequence.
/// @param polymer Polymer chain span to assign label_seq to
/// @param ent Entity with full_sequence from SEQRES (nullptr allowed if force=true)
/// @param force If true, assign label_seq sequentially even without Entity information
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

/// @brief Reset all Residue::label_seq to unset state.
/// @param st Structure to clear
inline void clear_label_seq_id(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        res.label_seq = Residue::OptionalNum();
}

/// @brief Assign label_seq_id to all polymer chains in a structure.
/// @param st Structure to process
/// @param force If true, assign label_seq sequentially even without Entity information
inline void assign_label_seq_id(Structure& st, bool force) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      if (ResidueSpan polymer = chain.get_polymer())
        if (!polymer.front().label_seq || !polymer.back().label_seq) {
          const Entity* ent = st.get_entity_of(polymer);
          assign_label_seq_to_polymer(polymer, ent, force);
        }
}


/// @}
/// @name Structure Superposition
/// @{

/// @brief Atom selection mode for superposition calculations.
enum class SupSelect {
  /// Use only Cα atoms (for peptides) or P atoms (for nucleic acids)
  CaP,
  /// Use backbone atoms (N, CA, C, O for peptides; P, O5', C5', C4', C3', O3' for nucleic acids)
  MainChain,
  /// Use all atoms
  All
};

/// @brief Extract matching atom positions from two polymer spans for superposition.
/// @param pos1 Output vector of positions from fixed chain
/// @param pos2 Output vector of positions from movable chain
/// @param fixed Fixed polymer chain (reference structure)
/// @param movable Movable polymer chain (structure to be superposed)
/// @param ptype Polymer type (peptide or nucleic acid)
/// @param sel Atom selection mode (CaP, MainChain, or All)
/// @param altloc Alternate location code to select; '\0' means all conformers
/// @param ca_offsets Optional output vector storing indices of Cα/P atoms in pos1/pos2
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

/// @brief Calculate RMSD between two polymer spans without transformation.
/// @param fixed Fixed polymer chain (reference structure)
/// @param movable Movable polymer chain (structure to compare)
/// @param ptype Polymer type (peptide or nucleic acid)
/// @param sel Atom selection mode (CaP, MainChain, or All)
/// @param altloc Alternate location code to select; '\0' means all conformers
/// @return SupResult with RMSD and atom count (rotation/translation matrices are empty)
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

/// @brief Calculate optimal superposition using QCP (quaternion characteristic polynomial).
/// @param fixed Fixed polymer chain (reference structure)
/// @param movable Movable polymer chain (structure to be superposed)
/// @param ptype Polymer type (peptide or nucleic acid)
/// @param sel Atom selection mode (CaP, MainChain, or All)
/// @param trim_cycles Number of iterative trimming cycles (0 = no trimming)
/// @param trim_cutoff Outlier distance threshold in Å for trimming (default 2.0)
/// @param altloc Alternate location code to select; '\0' means all conformers
/// @return SupResult with rotation matrix, translation vector, and RMSD
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

/// @brief Compute local superpositions using a sliding window around each Cα/P.
/// Calculates superpositions by sliding a sphere of given radius around each Cα/P position.
/// @param fixed Fixed polymer chain (reference structure)
/// @param movable Movable polymer chain (structure to be superposed)
/// @param ptype Polymer type (peptide or nucleic acid)
/// @param radius Window radius in Ångströms around each Cα/P (default 10.0)
/// @return Vector of SupResult, one per residue in fixed.first_conformer()
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

/// @}

} // namespace gemmi
#endif
