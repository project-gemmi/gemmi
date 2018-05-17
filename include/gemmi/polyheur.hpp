// Copyright 2017 Global Phasing Ltd.
//
// Heuristics for working with polymers.

#ifndef GEMMI_POLYHEUR_HPP_
#define GEMMI_POLYHEUR_HPP_

#include <vector>
#include "model.hpp"
#include "resinfo.hpp"

namespace gemmi {

// A simplistic classification. It may change in the future.
// It returns PolymerType which corresponds to _entity_poly.type,
// but here we use only PeptideL, Rna, Dna, DnaRnaHybrid and Unknown.
inline PolymerType check_polymer_type(const std::vector<Residue>& rr) {
  size_t counts[9] = {0};
  size_t aa = 0;
  size_t na = 0;
  for (const Residue& r : rr) {
    ResidueInfo info = find_tabulated_residue(r.name);
    if  (info.found()) {
      counts[info.kind]++;
    } else {
      if (r.get_ca())
        ++aa;
      else if (r.get_p())
        ++na;
    }
  }
  aa += counts[ResidueInfo::AA] + counts[ResidueInfo::AAD];
  na += counts[ResidueInfo::RNA] + counts[ResidueInfo::DNA];
  if (aa == rr.size() || (aa > 10 && 2 * aa > rr.size()))
    return counts[ResidueInfo::AA] >= counts[ResidueInfo::AAD]
           ? PolymerType::PeptideL : PolymerType::PeptideD;
  if (na == rr.size() || (na > 10 && 2 * na > rr.size())) {
    if (counts[ResidueInfo::DNA] == 0)
      return PolymerType::Rna;
    else if (counts[ResidueInfo::RNA] == 0)
      return PolymerType::Dna;
    else
      return PolymerType::DnaRnaHybrid;
  }
  return PolymerType::Unknown;
}

// TODO: use it in remove_ligands_and_waters()
inline bool is_polymer_residue(const Residue& res, PolymerType ptype) {
  ResidueInfo info = find_tabulated_residue(res.name);
  switch (ptype) {
    case PolymerType::PeptideL:
      return info.found() ? info.kind == ResidueInfo::AA : !!res.get_ca();
    case PolymerType::PeptideD:
      // D-peptide can contain AA in addition AAD (for example GLY)
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

inline bool are_connected(const Residue& r1, const Residue& r2,
                          PolymerType ptype) {
  if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD) {
    // similar to has_peptide_bond_to()
    const Atom* a1 = r1.get_c();
    const Atom* a2 = r2.get_n();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
  }
  if (ptype == PolymerType::Dna || ptype == PolymerType::Rna ||
      ptype == PolymerType::DnaRnaHybrid) {
    const Atom* a1 = r1.get_o3prim();
    const Atom* a2 = r2.get_p();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
  }
  return false;
}

// not a good check, but requires only CA (or P) atoms
inline bool are_connected2(const Residue& r1, const Residue& r2,
                           PolymerType ptype) {
  if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD) {
    const Atom* a1 = r1.get_ca();
    const Atom* a2 = r2.get_ca();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(5.0);
  }
  if (ptype == PolymerType::Dna || ptype == PolymerType::Rna ||
      ptype == PolymerType::DnaRnaHybrid) {
    const Atom* a1 = r1.get_p();
    const Atom* a2 = r2.get_p();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(7.5);
  }
  return false;
}

inline std::string make_one_letter_sequence(const std::vector<Residue>& rr,
                                            PolymerType ptype) {
  std::string seq;
  const Residue* prev = nullptr;
  for (const Residue& residue : rr)
    if (is_polymer_residue(residue, ptype)) {
      ResidueInfo info = find_tabulated_residue(residue.name);
      if (prev && !are_connected2(*prev, residue, ptype))
        seq += '-';
      seq += (info.one_letter_code != ' ' ? info.one_letter_code : 'X');
      prev = &residue;
    }
  return seq;
}

// returns a string such as AAL:GSHMTTPSHLSDRYEL
inline std::string extract_sequence_info(const Chain& chain) {
  PolymerType ptype = check_polymer_type(chain.residues);
  std::string info = polymer_type_abbr(ptype);
  info += ':';
  info += make_one_letter_sequence(chain.residues, ptype);
  return info;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
