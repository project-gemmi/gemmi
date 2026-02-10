//! @file
//! @brief Place riding hydrogens from restraints.
//!
//! Place hydrogens according to bond lengths and angles from monomer library.

// Copyright 2018-2022 Global Phasing Ltd.
//
// Place hydrogens according to bond lengths and angles from monomer library.

#ifndef GEMMI_RIDING_H_HPP_
#define GEMMI_RIDING_H_HPP_

#include <cmath>         // for isnan
#include "topo.hpp"      // for Topo

namespace gemmi {

//! @brief Place all hydrogens using riding positions.
//! @param topo Topology with restraints
GEMMI_DLL void place_hydrogens_on_all_atoms(Topo& topo);

//! @brief Adjust hydrogen-heavy atom distances to restraint values.
//! @param topo Topology with bonds
//! @param of Which distance to use from restraints
//! @param default_scale Scale factor if restraint distance is NaN
inline void adjust_hydrogen_distances(Topo& topo, Restraints::DistanceOf of,
                                      double default_scale=1.) {
  for (const Topo::Bond& t : topo.bonds) {
    assert(t.atoms[0] != nullptr && t.atoms[1] != nullptr);
    if (t.atoms[0]->is_hydrogen() || t.atoms[1]->is_hydrogen()) {
      Position u = t.atoms[1]->pos - t.atoms[0]->pos;
      double scale = t.restr->distance(of) / u.length();
      if (std::isnan(scale))
        scale = default_scale;
      if (t.atoms[1]->is_hydrogen())
        t.atoms[1]->pos = t.atoms[0]->pos + u * scale;
      else
        t.atoms[0]->pos = t.atoms[1]->pos - u * scale;
    }
  }
}

} // namespace gemmi
#endif
