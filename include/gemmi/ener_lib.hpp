// Copyright 2026 Global Phasing Ltd.
//
// Data from CCP4 ener_lib.cif used by monlib and related helpers.

#ifndef GEMMI_ENER_LIB_HPP_
#define GEMMI_ENER_LIB_HPP_

#include <string>
#include <vector>

#include "gemmi/chemcomp.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/fail.hpp"  // for GEMMI_DLL

namespace gemmi {

/// @brief Energy library from CCP4 ener_lib.cif.
/// Stores atomic properties and ideal bond parameters used for structure validation.
struct GEMMI_DLL EnerLib {
  /// @brief Atom radius types.
  enum class RadiusType {
    Vdw,   ///< Van der Waals radius
    Vdwh,  ///< Van der Waals radius with hydrogen
    Ion    ///< Ionic radius
  };

  /// @brief Atomic properties indexed by atom type.
  struct Atom {
    Element element;    ///< Chemical element
    char hb_type;       ///< Hydrogen bond type
    double vdw_radius;  ///< Van der Waals radius
    double vdwh_radius; ///< Van der Waals radius (hydrogen atoms)
    double ion_radius;  ///< Ionic radius
    int valency;        ///< Valence
    int sp;             ///< sp hybridization state
  };

  /// @brief Ideal bond parameters.
  struct Bond {
    std::string atom_type_1;  ///< First atom type
    std::string atom_type_2;  ///< Second atom type
    BondType type;            ///< Bond type
    double length;            ///< Ideal bond length
    double value_esd;         ///< Standard deviation

    /// @brief Comparison operator for sorting by atom types.
    /// Sorts first by atom_type_1, then by atom_type_2.
    bool operator<(const Bond& o) const {
      if (atom_type_1 != o.atom_type_1)
        return atom_type_1 < o.atom_type_1;
      return atom_type_2 < o.atom_type_2;
    }
  };

  EnerLib() {}

  /// @brief Read energy library data from a CIF document.
  /// @param doc CIF document containing ener_lib data tables
  void read(const cif::Document& doc);

  std::map<std::string, Atom> atoms; ///< Atom properties indexed by type
  std::vector<Bond> bonds;           ///< Ideal bond parameters
};

/// @brief Compare bond with atom type string (for binary search).
/// @param lhs Bond to compare
/// @param rhs Atom type string
/// @return true if bond's first atom type is less than the string
inline bool operator<(const EnerLib::Bond& lhs, const std::string& rhs) {
  return lhs.atom_type_1 < rhs;
}

/// @brief Compare atom type string with bond (for binary search).
/// @param lhs Atom type string
/// @param rhs Bond to compare
/// @return true if string is less than bond's first atom type
inline bool operator<(const std::string& lhs, const EnerLib::Bond& rhs) {
  return lhs < rhs.atom_type_1;
}

}  // namespace gemmi

#endif  // GEMMI_ENER_LIB_HPP_
