// Copyright 2026 Global Phasing Ltd.
//
// Data from CCP4 ener_lib.cif used by monlib and related helpers.

#include "gemmi/ener_lib.hpp"

#include <algorithm>

namespace gemmi {

void EnerLib::read(const cif::Document& doc) {
  cif::Block& block = const_cast<cif::Block&>(doc.blocks[0]);
  for (const auto& row : block.find("_lib_atom.",
                  {"type", "hb_type", "vdw_radius", "vdwh_radius",
                   "ion_radius", "element", "valency", "sp"}))
    atoms.emplace(row[0], Atom{Element(row[5]), row[1][0], cif::as_number(row[2]),
                               cif::as_number(row[3]), cif::as_number(row[4]),
                               cif::as_int(row[6], -1), cif::as_int(row[7], -1)});
  for (const auto& row : block.find("_lib_bond.",
                  {"atom_type_1", "atom_type_2", "type", "length", "value_esd"}))
    bonds.push_back(Bond{row.str(0), row.str(1), bond_type_from_string(row[2]),
                         cif::as_number(row[3]), cif::as_number(row[4])});
  std::sort(bonds.begin(), bonds.end());
}

}  // namespace gemmi
