// Copyright 2026 Global Phasing Ltd.
//
// Data from CCP4 ener_lib.cif used by monlib and related helpers.

#ifndef GEMMI_ENER_LIB_HPP_
#define GEMMI_ENER_LIB_HPP_

//:#include <map>
#include <string>
#include <vector>

#include "gemmi/chemcomp.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/fail.hpp"  // for GEMMI_DLL

namespace gemmi {

struct GEMMI_DLL EnerLib {
  enum class RadiusType {Vdw, Vdwh, Ion};
  struct Atom {
    Element element;
    char hb_type;
    double vdw_radius;
    double vdwh_radius;
    double ion_radius;
    int valency;
    int sp;
  };
  struct Bond {
    std::string atom_type_1;
    std::string atom_type_2;
    BondType type;
    double length;
    double value_esd;

    bool operator<(const Bond& o) const {
      if (atom_type_1 != o.atom_type_1)
        return atom_type_1 < o.atom_type_1;
      return atom_type_2 < o.atom_type_2;
    }
  };

  EnerLib() {}
  void read(const cif::Document& doc);
  std::map<std::string, Atom> atoms; // type->Atom
  std::vector<Bond> bonds;
};

inline bool operator<(const EnerLib::Bond& lhs, const std::string& rhs) {
  return lhs.atom_type_1 < rhs;
}

inline bool operator<(const std::string& lhs, const EnerLib::Bond& rhs) {
  return lhs < rhs.atom_type_1;
}

}  // namespace gemmi

#endif  // GEMMI_ENER_LIB_HPP_
