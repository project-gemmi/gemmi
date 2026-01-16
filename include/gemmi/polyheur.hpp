//! @file
//! @brief Heuristic methods for working with chains and polymers.
//!
//! Provides heuristic functions for polymer classification, connectivity checking,
//! and chain manipulation. Also includes well-defined functions such as removal of waters.

// Copyright 2017-2018 Global Phasing Ltd.

#ifndef GEMMI_POLYHEUR_HPP_
#define GEMMI_POLYHEUR_HPP_

#include <vector>
#include "model.hpp"
#include "util.hpp"      // for vector_remove_if

namespace gemmi {

//! @brief Heuristically determine polymer type from residue sequence.
//! @param span Residue span to analyze
//! @param ignore_entity_type If true, don't use entity type information
//! @return PolymerType (PeptideL, Rna, Dna, DnaRnaHybrid, or Unknown)
//!
//! A simplistic classification. It may change in the future.
//! It returns PolymerType which corresponds to _entity_poly.type,
//! but here we use only PeptideL, Rna, Dna, DnaRnaHybrid and Unknown.
GEMMI_DLL PolymerType check_polymer_type(const ConstResidueSpan& span,
                                         bool ignore_entity_type=false);

//! @brief Get polymer type from entity or check heuristically.
//! @param ent Entity (may be nullptr)
//! @param polymer Residue span
//! @return PolymerType from entity, or determined heuristically
inline PolymerType get_or_check_polymer_type(const Entity* ent,
                                             const ConstResidueSpan& polymer) {
  if (ent && ent->polymer_type != PolymerType::Unknown)
    return ent->polymer_type;
  return check_polymer_type(polymer);
}

//! @brief Atom name and element pair for mainchain atoms.
struct AtomNameElement { std::string atom_name; El el; };

//! @brief Get list of mainchain atoms for a polymer type.
//! @param ptype Polymer type
//! @return Vector of mainchain atom specifications
//!
//! Returns backbone atoms: N,CA,C,O for peptides; sugar-phosphate for nucleotides.
inline std::vector<AtomNameElement> get_mainchain_atoms(PolymerType ptype) {
  if (is_polynucleotide(ptype))
    return {{"P", El::P}, {"O5'", El::O}, {"C5'", El::C},
            {"C4'", El::C}, {"O4'", El::O}, {"C3'", El::C}, {"O3'", El::O},
            {"C2'", El::C}, {"O2'", El::O}, {"C1'", El::C}};
  return {{"N", El::N}, {"CA", El::C}, {"C", El::C}, {"O", El::O}};
}

//! @brief Distance-based check for peptide bond between two atoms.
//! @param a1 First atom (typically C)
//! @param a2 Second atom (typically N)
//! @return True if atoms are within peptide bond distance (<2.0 Å)
inline bool in_peptide_bond_distance(const Atom* a1, const Atom* a2) {
  return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
}

//! @brief Check if two residues have a peptide bond (C to N).
//! @param r1 First residue
//! @param r2 Second residue
//! @return True if peptide bond exists
inline bool have_peptide_bond(const Residue& r1, const Residue& r2) {
  return in_peptide_bond_distance(r1.get_c(), r2.get_n());
}

//! @brief Distance-based check for phosphodiester bond between nucleotides.
//! @param a1 First atom (typically O3')
//! @param a2 Second atom (typically P)
//! @return True if atoms are within nucleotide bond distance (<2.4 Å)
inline bool in_nucleotide_bond_distance(const Atom* a1, const Atom* a2) {
  return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
}

//! @brief Check if two residues have a phosphodiester bond (O3' to P).
//! @param r1 First residue
//! @param r2 Second residue
//! @return True if nucleotide bond exists
inline bool have_nucleotide_bond(const Residue& r1, const Residue& r2) {
  return in_nucleotide_bond_distance(r1.get_o3prim(), r2.get_p());
}

//! @brief Check if two residues are connected by a polymer bond.
//! @param r1 First residue
//! @param r2 Second residue
//! @param ptype Polymer type
//! @return True if residues are connected
//!
//! Check C-N distance for peptides or O3'-P distance for nucleotides.
inline bool are_connected(const Residue& r1, const Residue& r2, PolymerType ptype) {
  if (is_polypeptide(ptype))
    return have_peptide_bond(r1, r2);
  if (is_polynucleotide(ptype))
    return have_nucleotide_bond(r1, r2);
  return false;
}

//! @brief Check residue connectivity using only CA or P atoms.
//! @param r1 First residue
//! @param r2 Second residue
//! @param ptype Polymer type
//! @return True if residues appear connected
//!
//! are_connected2() is less exact than are_connected(), but requires only CA (or P) atoms.
//! Uses CA-CA distance (<5 Å) for peptides, P-P distance (<7.5 Å) for nucleotides.
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

//! @brief Check connectivity with fallback: are_connected() + are_connected2().
//! @param r1 First residue
//! @param r2 Second residue
//! @param ptype Polymer type
//! @return True if residues are connected
//!
//! are_connected3() = are_connected() + fallback to are_connected2().
//! First tries exact bond distance check, falls back to CA/P distance if needed.
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

GEMMI_DLL std::string make_one_letter_sequence(const ConstResidueSpan& polymer);

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

/// Determine ATOM/HETATM record type, based on Residue::entity_type
GEMMI_DLL char recommended_het_flag(const Residue& res);
/// R = recommended_het_flag(), other valid values are A, H and '\0'
template<class T> void assign_het_flags(T& obj, char flag='R') {
  for (auto& child : obj.children())
    assign_het_flags(child, flag);
}
template<> inline void assign_het_flags(Residue& res, char flag) {
  flag &= ~0x20; // uppercase letters, ' ' -> \0
  if (flag != 'R' && flag != '\0' && flag != 'A' && flag != 'H')
    fail("assign_het_flags(): the only allowed values are A, H, ' ' and R");
  res.het_flag = flag == 'R' ? recommended_het_flag(res) : flag;
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
GEMMI_DLL bool trim_to_alanine(Residue& res);

inline void trim_to_alanine(Chain& chain) {
  for (Residue& res : chain.residues)
    trim_to_alanine(res);
}

// Functions for switching between long (>3 chars) residue names (CCD codes)
// and shortened ones that are compatible with the PDB format.
GEMMI_DLL void shorten_ccd_codes(Structure& st);
GEMMI_DLL void restore_full_ccd_codes(Structure& st);

/// Modifies Entity::full_sequence. Uses only the first chain for each Entity.
GEMMI_DLL void add_microhetero_to_sequences(Structure& st, bool overwrite=false);

GEMMI_DLL void add_tls_group_ids(Structure& st);

} // namespace gemmi
#endif
