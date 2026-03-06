// Copyright 2026 Global Phasing Ltd.
//
// Data from CCP4 ener_lib.cif used by monlib and related helpers.

#ifndef GEMMI_ENER_LIB_HPP_
#define GEMMI_ENER_LIB_HPP_

#include <map>
#include <string>

#include "gemmi/chemcomp.hpp"
#include "gemmi/cifdoc.hpp"

namespace gemmi {

struct EnerLib {
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
    std::string atom_type_2;
    BondType type;
    double length;
    double value_esd;
  };

  EnerLib() {}
  void read(const cif::Document& doc);
  std::map<std::string, Atom> atoms; // type->Atom
  std::multimap<std::string, Bond> bonds; // atom_type_1->Bond
};

}  // namespace gemmi

#endif  // GEMMI_ENER_LIB_HPP_
