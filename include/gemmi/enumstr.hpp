//! @file
//! @brief Conversions between enums and mmCIF strings.
//!
//! Converts between enums (EntityType, PolymerType, Connection::Type,
//! SoftwareItem::Classification) and mmCIF strings.

// Copyright 2018 Global Phasing Ltd.
//
// Converts between enums (EntityType, PolymerType, Connection::Type,
// SoftwareItem::Classification) and mmCIF strings.

#ifndef GEMMI_ENUMSTR_HPP_
#define GEMMI_ENUMSTR_HPP_

#include "metadata.hpp"  // for EntityType, PolymerType, SoftwareItem
#include "util.hpp"      // for iequal

namespace gemmi {

//! @brief Convert EntityType to mmCIF string.
//! @param entity_type Entity type enum
//! @return mmCIF string representation
inline const char* entity_type_to_string(EntityType entity_type) {
  switch (entity_type) {
    case EntityType::Polymer: return "polymer";
    case EntityType::Branched: return "branched";
    case EntityType::NonPolymer: return "non-polymer";
    case EntityType::Water: return "water";
    default /*EntityType::Unknown*/: return "?";
  }
}

//! @brief Convert mmCIF string to EntityType.
//! @param t mmCIF string
//! @return Entity type enum
inline EntityType entity_type_from_string(const std::string& t) {
  if (t == "polymer")     return EntityType::Polymer;
  if (t == "branched")    return EntityType::Branched;
  if (t == "non-polymer") return EntityType::NonPolymer;
  if (t == "water")       return EntityType::Water;
  return EntityType::Unknown;
}


//! @brief Convert PolymerType to mmCIF string.
//! @param polymer_type Polymer type enum
//! @return mmCIF string representation
inline const char* polymer_type_to_string(PolymerType polymer_type) {
  switch (polymer_type) {
    case PolymerType::PeptideL: return "polypeptide(L)";
    case PolymerType::PeptideD: return "polypeptide(D)";
    case PolymerType::Dna: return "polydeoxyribonucleotide";
    case PolymerType::Rna: return "polyribonucleotide";
    case PolymerType::DnaRnaHybrid:
      return "'polydeoxyribonucleotide/polyribonucleotide hybrid'";
    case PolymerType::SaccharideD: return "polysaccharide(D)";
    case PolymerType::SaccharideL: return "polysaccharide(L)";
    case PolymerType::Other: return "other";
    case PolymerType::Pna: return "'peptide nucleic acid'";
    case PolymerType::CyclicPseudoPeptide: return "cyclic-pseudo-peptide";
    default /*PolymerType::Unknown*/: return "?";
  }
}

//! @brief Convert mmCIF string to PolymerType.
//! @param t mmCIF string
//! @return Polymer type enum
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


//! @brief Convert Connection::Type to mmCIF string.
//! @param t Connection type enum
//! @return mmCIF string representation
inline const char* connection_type_to_string(Connection::Type t) {
  static constexpr const char* type_ids[] = {
    "covale", "disulf", "hydrog", "metalc", "."
  };
  return type_ids[t];
}

//! @brief Convert mmCIF string to Connection::Type.
//! @param t mmCIF string
//! @return Connection type enum
inline Connection::Type connection_type_from_string(const std::string& t) {
  for (int i = 0; i != Connection::Unknown; ++i)
    if (connection_type_to_string(Connection::Type(i)) == t)
      return Connection::Type(i);
  return Connection::Unknown;
}

//! @brief Convert SoftwareItem::Classification to string.
//! @param c Software classification enum
//! @return String representation
inline
std::string software_classification_to_string(SoftwareItem::Classification c) {
  switch (c) {
    case SoftwareItem::DataCollection: return "data collection";
    case SoftwareItem::DataExtraction: return "data extraction";
    case SoftwareItem::DataProcessing: return "data processing";
    case SoftwareItem::DataReduction:  return "data reduction";
    case SoftwareItem::DataScaling:    return "data scaling";
    case SoftwareItem::ModelBuilding:  return "model building";
    case SoftwareItem::Phasing:        return "phasing";
    case SoftwareItem::Refinement:     return "refinement";
    case SoftwareItem::Unspecified:    return "";
  }
  unreachable();
}

//! @brief Convert string to SoftwareItem::Classification.
//! @param str String representation
//! @return Software classification enum
inline SoftwareItem::Classification
software_classification_from_string(const std::string& str) {
  if (iequal(str, "data collection")) return SoftwareItem::DataCollection;
  if (iequal(str, "data extraction")) return SoftwareItem::DataExtraction;
  if (iequal(str, "data processing")) return SoftwareItem::DataProcessing;
  if (iequal(str, "data reduction"))  return SoftwareItem::DataReduction;
  if (iequal(str, "data scaling"))    return SoftwareItem::DataScaling;
  if (iequal(str, "model building"))  return SoftwareItem::ModelBuilding;
  if (iequal(str, "phasing"))         return SoftwareItem::Phasing;
  if (iequal(str, "refinement"))      return SoftwareItem::Refinement;
  return SoftwareItem::Unspecified;
}

} // namespace gemmi
#endif
