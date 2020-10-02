// Copyright 2017-2018 Global Phasing Ltd.
//
// Heuristic methods for working with chains and polymers.
// Includes also a few well-defined functions, such as removal of hydrogens.

#ifndef GEMMI_POLYHEUR_HPP_
#define GEMMI_POLYHEUR_HPP_

#include <vector>
#include <set>
#include "model.hpp"
#include "resinfo.hpp"   // for find_tabulated_residue
#include "seqalign.hpp"  // for align_sequences
#include "util.hpp"      // for vector_remove_if

namespace gemmi {

// A simplistic classification. It may change in the future.
// It returns PolymerType which corresponds to _entity_poly.type,
// but here we use only PeptideL, Rna, Dna, DnaRnaHybrid and Unknown.
inline PolymerType check_polymer_type(const ConstResidueSpan& polymer) {
  if (polymer.size() < 2)
    return PolymerType::Unknown;
  size_t counts[ResidueInfo::ELS+1] = {0};
  size_t aa = 0;
  size_t na = 0;
  for (const Residue& r : polymer)
    if (r.entity_type == EntityType::Unknown ||
        r.entity_type == EntityType::Polymer) {
      ResidueInfo info = find_tabulated_residue(r.name);
      if (info.found())
        counts[info.kind]++;
      else if (r.get_ca())
        ++aa;
      else if (r.get_p())
        ++na;
    }
  aa += counts[ResidueInfo::AA] + counts[ResidueInfo::AAD] +
        counts[ResidueInfo::PAA] + counts[ResidueInfo::MAA];
  na += counts[ResidueInfo::RNA] + counts[ResidueInfo::DNA];
  if (aa == polymer.size() || (aa > 10 && 2 * aa > polymer.size()))
    return counts[ResidueInfo::AA] >= counts[ResidueInfo::AAD]
           ? PolymerType::PeptideL : PolymerType::PeptideD;
  if (na == polymer.size() || (na > 10 && 2 * na > polymer.size())) {
    if (counts[ResidueInfo::DNA] == 0)
      return PolymerType::Rna;
    else if (counts[ResidueInfo::RNA] == 0)
      return PolymerType::Dna;
    else
      return PolymerType::DnaRnaHybrid;
  }
  return PolymerType::Unknown;
}

inline double calculate_sequence_weight(const std::vector<std::string>& seq,
                                        double unknown=0.) {
  double weight = 0.;
  for (const std::string& item : seq) {
    ResidueInfo res_info = find_tabulated_residue(Entity::first_mon(item));
    weight += res_info.found() ? res_info.weight : unknown;
  }
  return weight;
}

inline std::string one_letter_code(const std::vector<std::string>& seq) {
  std::string r;
  for (const std::string& item : seq)
    r += find_tabulated_residue(Entity::first_mon(item)).fasta_code();
  return r;
}

inline std::string one_letter_code(const ConstResidueSpan& polymer) {
  std::string r;
  for (const Residue& res : polymer.first_conformer())
    r += find_tabulated_residue(res.name).fasta_code();
  return r;
}

inline bool is_polymer_residue(const Residue& res, PolymerType ptype) {
  ResidueInfo info = find_tabulated_residue(res.name);
  // If a standard residue is HETATM we assume that it is in the buffer.
  if (info.found() && info.is_standard() && res.het_flag == 'H')
    return false;
  switch (ptype) {
    case PolymerType::PeptideL:
    case PolymerType::PeptideD:
      // here we don't mind mixing D- and L- peptides
      return info.found() ? info.is_amino_acid() : !!res.get_ca();
    case PolymerType::Dna:
      return info.found() ? info.is_dna() : !!res.get_p();
    case PolymerType::Rna:
      return info.found() ? info.is_rna() : !!res.get_p();
    case PolymerType::DnaRnaHybrid:
      return info.found() ? info.is_nucleic_acid() : !!res.get_p();
    default:
      return false;
  }
}

inline bool are_connected(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype)) {
    const Atom* a1 = r1.get_c();
    const Atom* a2 = r2.get_n();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
  }
  if (is_polynucleotide(ptype)) {
    const Atom* a1 = r1.get_o3prim();
    const Atom* a2 = r2.get_p();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
  }
  return false;
}

// are_connected2() is less exact, but requires only CA (or P) atoms.
inline bool are_connected2(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype)) {
    const Atom* a1 = r1.get_ca();
    const Atom* a2 = r2.get_ca();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(5.0);
  }
  if (is_polynucleotide(ptype)) {
    const Atom* a1 = r1.get_p();
    const Atom* a2 = r2.get_p();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(7.5);
  }
  return false;
}

// are_connected3() = are_connected() + fallback to are_connected2()
inline bool are_connected3(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype)) {
    if (const Atom* a1 = r1.get_c())
      if (const Atom* a2 = r2.get_n())
        return a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
    if (const Atom* a1 = r1.get_ca())
      if (const Atom* a2 = r2.get_ca())
        return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(5.0);
  } else if (is_polynucleotide(ptype)) {
    if (const Atom* a1 = r1.get_o3prim())
      if (const Atom* a2 = r2.get_p())
        return a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
    if (const Atom* a1 = r1.get_p())
      if (const Atom* a2 = r2.get_p())
        return a1->pos.dist_sq(a2->pos) < sq(7.5);
  }
  return false;
}

inline std::string make_one_letter_sequence(const ConstResidueSpan& polymer) {
  std::string seq;
  const Residue* prev = nullptr;
  PolymerType ptype = check_polymer_type(polymer);
  for (const Residue& residue : polymer.first_conformer()) {
    ResidueInfo info = find_tabulated_residue(residue.name);
    if (prev && !are_connected2(*prev, residue, ptype))
      seq += '-';
    seq += (info.one_letter_code != ' ' ? info.one_letter_code : 'X');
    prev = &residue;
  }
  return seq;
}

inline bool has_subchains_assigned(const Chain& chain) {
  return std::all_of(chain.residues.begin(), chain.residues.end(),
                     [](const Residue& r) { return !r.subchain.empty(); });
}

inline void add_entity_types(Chain& chain, bool overwrite) {
  PolymerType ptype = check_polymer_type(chain.whole());
  auto it = chain.residues.begin();
  for (; it != chain.residues.end(); ++it)
    if (overwrite || it->entity_type == EntityType::Unknown) {
      if (!is_polymer_residue(*it, ptype))
        break;
      it->entity_type = EntityType::Polymer;
    } else if (it->entity_type != EntityType::Polymer) {
      break;
    }
  for (; it != chain.residues.end(); ++it)
    if (overwrite || it->entity_type == EntityType::Unknown)
      it->entity_type = it->is_water() ? EntityType::Water
                                       : EntityType::NonPolymer;
}

inline void add_entity_types(Structure& st, bool overwrite) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      add_entity_types(chain, overwrite);
}

// The subchain field in the residue is where we store_atom_site.label_asym_id
// from mmCIF files. As of 2018 wwPDB software splits author's chains
// (auth_asym_id) into label_asym_id units:
// * linear polymer,
// * non-polymers (each residue has different separate label_asym_id),
// * and waters.
// Refmac/makecif is doing similar thing but using different naming and
// somewhat different rules (it was written in 1990's before PDBx/mmCIF).
//
// Here we use naming and rules different from both wwPDB and makecif.
inline void assign_subchain_names(Chain& chain) {
  for (Residue& res : chain.residues) {
    res.subchain = chain.name;
    switch (res.entity_type) {
      case EntityType::Polymer:    res.subchain += "poly";          break;
      case EntityType::NonPolymer: res.subchain += res.seqid.str(); break;
      case EntityType::Water:      res.subchain += "wat";           break;
      case EntityType::Branched:  // FIXME
      case EntityType::Unknown: break; // should not happen
    }
  }
}

inline void assign_subchains(Structure& st, bool force) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      if (force || !has_subchains_assigned(chain)) {
        add_entity_types(chain, false);
        assign_subchain_names(chain);
      }
}

inline void ensure_entities(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (ResidueSpan& sub : chain.subchains()) {
        Entity* ent = st.get_entity_of(sub);
        if (!ent) {
          EntityType etype = sub[0].entity_type;
          std::string name;
          if (etype == EntityType::Polymer)
            name = chain.name;
          else if (etype == EntityType::NonPolymer)
            name = sub[0].name + "!";
          else if (etype == EntityType::Water)
            name = "water";
          if (!name.empty()) {
            ent = &impl::find_or_add(st.entities, name);
            ent->entity_type = etype;
            ent->subchains.push_back(sub.subchain_id());
          }
        }
        // ensure we have polymer_type set where needed
        if (ent && ent->entity_type == EntityType::Polymer &&
            ent->polymer_type == PolymerType::Unknown)
          ent->polymer_type = check_polymer_type(sub);
      }
}

inline bool operator==(const Entity::DbRef& a, const Entity::DbRef& b) {
  return a.db_name == b.db_name &&
         a.id_code == b.id_code &&
         a.isoform == b.isoform &&
         a.seq_begin == b.seq_begin && a.seq_end == b.seq_end &&
         a.db_begin == b.db_begin && a.db_end == b.db_end;
}

inline void deduplicate_entities(Structure& st) {
  for (auto i = st.entities.begin(); i != st.entities.end(); ++i)
    if (!i->full_sequence.empty())
      for (auto j = i + 1; j != st.entities.end(); ++j)
        if (j->polymer_type == i->polymer_type &&
            j->full_sequence == i->full_sequence &&
            j->dbrefs == i->dbrefs) {
          vector_move_extend(i->subchains, std::move(j->subchains));
          st.entities.erase(j--);
        }
}

inline void setup_entities(Structure& st) {
  assign_subchains(st, false);
  ensure_entities(st);
  deduplicate_entities(st);
}

// Sequence alignment and label_seq_id assignment

// helper function for sequence alignment
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

// check for exact match between model sequence and full sequence (SEQRES)
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

// Uses sequence alignment (model to SEQRES) to assign label_seq.
// force: assign label_seq even if full sequence is not known (assumes no gaps)
inline void assign_label_seq_to_polymer(ResidueSpan& polymer,
                                        const Entity* ent, bool force) {
  AlignmentResult result;

  // sequence not known
  if (!ent || ent->full_sequence.empty()) {
    if (!force)
      return;
    PolymerType ptype = (ent && ent->polymer_type != PolymerType::Unknown
                         ? ent->polymer_type
                         : check_polymer_type(polymer));
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
    AlignmentScoring scoring;
    result = align_sequence_to_polymer(ent->full_sequence, polymer,
                                       ent->polymer_type, scoring);
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
        if (!polymer.front().label_seq || !polymer.back().label_seq) {
          const Entity* ent = st.get_entity_of(polymer);
          assign_label_seq_to_polymer(polymer, ent, force);
        }
}
// end of sequence alignment and label_seq_id business


// Remove alternative conformations.
template<class T> void remove_alternative_conformations(T& obj) {
  for (auto& child : obj.children())
    remove_alternative_conformations(child);
}
template<> inline void remove_alternative_conformations(Chain& chain) {
  std::set<SeqId> seqids;
  for (size_t i = 0; i < chain.residues.size(); ) {
    if (seqids.insert(chain.residues[i].seqid).second)
      ++i;
    else
      chain.residues.erase(chain.residues.begin() + i);
  }
  for (Residue& residue : chain.residues) {
    std::set<std::string> names;
    for (size_t i = 0; i < residue.atoms.size(); ) {
      Atom& atom = residue.atoms[i];
      atom.altloc = '\0';
      if (names.insert(atom.name).second)
        ++i;
      else
        residue.atoms.erase(residue.atoms.begin() + i);
    }
  }
}

// Remove hydrogens.
template<class T> void remove_hydrogens(T& obj) {
  for (auto& child : obj.children())
    remove_hydrogens(child);
}
template<> inline void remove_hydrogens(Residue& res) {
  vector_remove_if(res.atoms, [](const Atom& a) {
    return a.element == El::H || a.element == El::D;
  });
}

// Remove waters. It may leave empty chains.
template<class T> void remove_waters(T& obj) {
  for (auto& child : obj.children())
    remove_waters(child);
}
template<> inline void remove_waters(Chain& ch) {
  vector_remove_if(ch.residues,
                   [](const Residue& res) { return res.is_water(); });
}

// Remove ligands and waters. It may leave empty chains.
inline void remove_ligands_and_waters(Chain& ch) {
  PolymerType ptype = check_polymer_type(ch.whole());
  vector_remove_if(ch.residues, [&](const Residue& res) {
      if (res.entity_type == EntityType::Unknown) {
        // TODO: check connectivity
        return !is_polymer_residue(res, ptype);
      }
      return res.entity_type != EntityType::Polymer;
  });
}

inline void remove_ligands_and_waters(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      remove_ligands_and_waters(chain);
}

// Remove empty chains.
inline void remove_empty_chains(Model& m) {
  vector_remove_if(m.chains,
                   [](const Chain& chain) { return chain.residues.empty(); });
}
inline void remove_empty_chains(Structure& st) {
  for (Model& model : st.models)
    remove_empty_chains(model);
}

// Trim to alanine. Returns true if trimmed, false if it's (likely) not AA.
inline bool trim_to_alanine(Residue& res) {
  static const std::pair<std::string, El> ala_atoms[6] = {
    {"N", El::N}, {"CA", El::C}, {"C", El::C}, {"O", El::O}, {"CB", El::C},
    {"OXT", El::O}
  };
  if (res.get_ca() == nullptr)
    return false;
  vector_remove_if(res.atoms, [](const Atom& a) {
      for (const auto& name_el : ala_atoms)
        if (a.name == name_el.first && a.element == name_el.second)
          return false;
      return true;
  });
  return true;
}

inline void trim_to_alanine(Chain& chain) {
  for (Residue& res : chain.residues)
    trim_to_alanine(res);
}

} // namespace gemmi
#endif
