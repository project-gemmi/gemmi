// Copyright 2017-2018 Global Phasing Ltd.
//
// Heuristic methods for working with chains and polymers.
// Includes also a few well-defined functions, such as removal of waters.

#ifndef GEMMI_POLYHEUR_HPP_
#define GEMMI_POLYHEUR_HPP_

#include <vector>
#include "model.hpp"
#include "resinfo.hpp"   // for find_tabulated_residue
#include "util.hpp"      // for vector_remove_if

namespace gemmi {

// A simplistic classification. It may change in the future.
// It returns PolymerType which corresponds to _entity_poly.type,
// but here we use only PeptideL, Rna, Dna, DnaRnaHybrid and Unknown.
GEMMI_DLL PolymerType check_polymer_type(const ConstResidueSpan& span,
                                         bool ignore_entity_type=false);

inline PolymerType get_or_check_polymer_type(const Entity* ent,
                                             const ConstResidueSpan& polymer) {
  if (ent && ent->polymer_type != PolymerType::Unknown)
    return ent->polymer_type;
  return check_polymer_type(polymer);
}

struct AtomNameElement { std::string atom_name; El el; };

inline std::vector<AtomNameElement> get_mainchain_atoms(PolymerType ptype) {
  if (is_polynucleotide(ptype))
    return {{"P", El::P}, {"O5'", El::O}, {"C5'", El::C},
            {"C4'", El::C}, {"O4'", El::O}, {"C3'", El::C}, {"O3'", El::O},
            {"C2'", El::C}, {"O2'", El::O}, {"C1'", El::C}};
  return {{"N", El::N}, {"CA", El::C}, {"C", El::C}, {"O", El::O}};
}

/// distance-based check for peptide bond
inline bool in_peptide_bond_distance(const Atom* a1, const Atom* a2) {
  return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
}
inline bool have_peptide_bond(const Residue& r1, const Residue& r2) {
  return in_peptide_bond_distance(r1.get_c(), r2.get_n());
}

/// distance-based check for phosphodiester bond between nucleotide
inline bool in_nucleotide_bond_distance(const Atom* a1, const Atom* a2) {
  return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
}
inline bool have_nucleotide_bond(const Residue& r1, const Residue& r2) {
  return in_nucleotide_bond_distance(r1.get_o3prim(), r2.get_p());
}

/// check C-N or O3'-P distance
inline bool are_connected(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype))
    return have_peptide_bond(r1, r2);
  if (is_polynucleotide(ptype))
    return have_nucleotide_bond(r1, r2);
  return false;
}

/// are_connected2() is less exact, but requires only CA (or P) atoms.
inline bool are_connected2(const Residue& r1, const Residue& r2, PolymerType ptype) {
  auto this_or_first = [](const Atom* a, const Residue& r, El el) -> const Atom* {
    if (a || r.atoms.empty())
      return a;
    if (const Atom* b = r.find_by_element(el))
      return b;
    return &r.atoms.front();
  };
  if (is_polypeptide(ptype)) {
    const Atom* a1 = this_or_first(r1.get_ca(), r1, El::C);
    const Atom* a2 = this_or_first(r2.get_ca(), r2, El::C);
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(5.0);
  }
  if (is_polynucleotide(ptype)) {
    const Atom* a1 = this_or_first(r1.get_p(), r1, El::P);
    const Atom* a2 = this_or_first(r2.get_p(), r2, El::P);
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(7.5);
  }
  return false;
}

/// are_connected3() = are_connected() + fallback to are_connected2()
inline bool are_connected3(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype)) {
    if (const Atom* a1 = r1.get_c())
      if (const Atom* a2 = r2.get_n())
        return a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
    if (const Atom* a1 = r1.get_ca())
      if (const Atom* a2 = r2.get_ca())
        return a1->pos.dist_sq(a2->pos) < sq(5.0);
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

/// Assigns entity_type=Polymer|NonPolymer|Water for each Residue (only
/// for residues with entity_type==Unknown, unless overwrite=true).
/// Determining where the polymer ends and ligands start is sometimes
/// arbitrary -- there can be a non-standard residue at the end that can
/// be regarded as as either the last residue or a linked ligand.
GEMMI_DLL void add_entity_types(Chain& chain, bool overwrite);
GEMMI_DLL void add_entity_types(Structure& st, bool overwrite);

/// Assigns entity_type=Unknown for all residues.
GEMMI_DLL void remove_entity_types(Structure& st);

/// Assigns Residue::entity_id based on Residue::subchain and Entity::subchains.
GEMMI_DLL void add_entity_ids(Structure& st, bool overwrite);

/// The subchain field in the residue is where we store_atom_site.label_asym_id
/// from mmCIF files. As of 2018 wwPDB software splits author's chains
/// (auth_asym_id) into label_asym_id units:
/// * linear polymer,
/// * non-polymers (each residue has different separate label_asym_id),
/// * and waters.
/// Refmac/makecif is doing similar thing but using different naming and
/// somewhat different rules (it was written in 1990's before PDBx/mmCIF).
///
/// Here we use naming and rules different from both wwPDB and makecif.
/// Note: call add_entity_types() first.
GEMMI_DLL void assign_subchain_names(Chain& chain, int& nonpolymer_counter);

GEMMI_DLL void assign_subchains(Structure& st, bool force, bool fail_if_unknown=true);

GEMMI_DLL void ensure_entities(Structure& st);

GEMMI_DLL void deduplicate_entities(Structure& st);

inline void setup_entities(Structure& st) {
  add_entity_types(st, /*overwrite=*/false);
  assign_subchains(st, /*force=*/false);
  ensure_entities(st);
  deduplicate_entities(st);
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
template<class T> void remove_ligands_and_waters(T& obj) {
  for (auto& child : obj.children())
    remove_ligands_and_waters(child);
}
template<> inline void remove_ligands_and_waters(Chain& ch) {
  vector_remove_if(ch.residues, [&](const Residue& res) {
      if (res.entity_type == EntityType::Unknown)
        fail("remove_ligands_and_waters(): missing entity_type in chain ", ch.name);
      return res.entity_type != EntityType::Polymer;
  });
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

// Functions for switching between long (>3 chars) residue names (CCD codes)
// and shortened ones that are compatible with the PDB format.
GEMMI_DLL
void change_ccd_code(Structure& st, const std::string& old, const std::string& new_);

GEMMI_DLL void shorten_ccd_codes(Structure& st);

inline void restore_full_ccd_codes(Structure& st) {
  for (const OldToNew& item : st.shortened_ccd_codes)
    change_ccd_code(st, item.new_, item.old);
  st.shortened_ccd_codes.clear();
}

} // namespace gemmi
#endif
