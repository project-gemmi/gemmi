// Copyright 2017 Global Phasing Ltd.
//
// Heuristics for working with polymers.

#ifndef GEMMI_POLYHEUR_HPP_
#define GEMMI_POLYHEUR_HPP_

//#include <vector>     // for vector

#include "model.hpp"
//#include "util.hpp"
#include "resinfo.hpp"

namespace gemmi {

// A simplistic classification. It may change in the future.
// It returns PolymerType which corresponds to _entity_poly.type,
// but here we use only PeptideL, Rna, Dna, DnaRnaHybrid and Unknown.
inline PolymerType check_polymer_type(const std::vector<Residue>& rr) {
  size_t aa = 0;
  size_t dna = 0;
  size_t rna = 0;
  size_t na = 0;
  for (const Residue& r : rr)
    switch (find_tabulated_residue(r.name).kind) {
      case ResidueInfo::AA: ++aa; break;
      case ResidueInfo::DNA: ++dna; break;
      case ResidueInfo::RNA: ++rna; break;
      case ResidueInfo::UNKNOWN:
        if (r.get_ca())
          ++aa;
        else if (r.get_p())
          ++na;
        break;
      default: break;
    }
  if (aa == rr.size() || (aa > 10 && 2 * aa > rr.size()))
    return PolymerType::PeptideL;
  na += dna + rna;
  if (na == rr.size() || (na > 10 && 2 * na > rr.size())) {
    if (dna == 0)
      return PolymerType::Rna;
    else if (rna == 0)
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
    case PolymerType::PeptideD:
      // for now we don't differentiate L and D
      return info.found() ? info.is_amino() : (res.get_ca() != nullptr);
    case PolymerType::Dna:
      return info.found() ? info.is_dna() : (res.get_p() != nullptr);
    case PolymerType::Rna:
      return info.found() ? info.is_rna() : (res.get_p() != nullptr);
    case PolymerType::DnaRnaHybrid:
      return info.found() ? info.is_nucleic() : (res.get_p() != nullptr);
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
