// Copyright 2018 Global Phasing Ltd.
//
// EntityType, PolymerType <-> PDBx/mmCIF names
// PolymerType <-> 3-letter string

#ifndef GEMMI_ENTSTR_HPP_
#define GEMMI_ENTSTR_HPP_

#include "model.hpp"  // for EntityType, PolymerType

namespace gemmi {

inline std::string entity_type_to_string(EntityType entity_type) {
  switch (entity_type) {
    case EntityType::Polymer: return "polymer";
    case EntityType::NonPolymer: return "non-polymer";
    case EntityType::Water: return "water";
    default /*EntityType::Unknown*/: return "?";
  }
}

inline EntityType entity_type_from_string(const std::string& t) {
  if (t == "polymer")     return EntityType::Polymer;
  if (t == "non-polymer") return EntityType::NonPolymer;
  if (t == "water")       return EntityType::Water;
  return EntityType::Unknown;
}


inline std::string polymer_type_to_string(PolymerType polymer_type) {
  switch (polymer_type) {
    case PolymerType::PeptideL: return "polypeptide(L)";
    case PolymerType::PeptideD: return "polypeptide(D)";
    case PolymerType::Dna: return "polydeoxyribonucleotide";
    case PolymerType::Rna: return "polyribonucleotide";
    case PolymerType::DnaRnaHybrid:
      return "polydeoxyribonucleotide/polyribonucleotide hybrid";
    case PolymerType::SaccharideD: return "polysaccharide(D)";
    case PolymerType::SaccharideL: return "polysaccharide(L)";
    case PolymerType::Other: return "other";
    case PolymerType::Pna: return "peptide nucleic acid";
    case PolymerType::CyclicPseudoPeptide: return "cyclic-pseudo-peptide";
    default /*PolymerType::Unknown*/: return "?";
  }
}

inline PolymerType polymer_type_from_string(const std::string& t) {
  if (t == "polypeptide(L)")          return PolymerType::PeptideL;
  if (t == "polydeoxyribonucleotide") return PolymerType::Dna;
  if (t == "polyribonucleotide")      return PolymerType::Rna;
  if (t == "polydeoxyribonucleotide/polyribonucleotide hybrid")
                                      return PolymerType::DnaRnaHybrid;
  if (t == "polypeptide(D)")          return PolymerType::PeptideD;
  if (t == "polysaccharide(D)")       return PolymerType::SaccharideD;
  if (t == "other")                   return PolymerType::Other;
  if (t == "peptide nucleic acid")    return PolymerType::Pna;
  if (t == "cyclic-pseudo-peptide")   return PolymerType::CyclicPseudoPeptide;
  if (t == "polysaccharide(L)")       return PolymerType::SaccharideL;
  return PolymerType::Unknown;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
