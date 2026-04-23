/// @file reciproc.hpp
/// @brief Utility functions for working with reflections and reciprocal space.
///
/// Provides functions to enumerate unique reflections within a resolution range,
/// respecting space group symmetry and systematically absent reflections.

// Copyright 2020 Global Phasing Ltd.
//
// Reciprocal space helper functions.

#ifndef GEMMI_RECIPROC_HPP_
#define GEMMI_RECIPROC_HPP_

#include "symmetry.hpp"  // for SpaceGroup
#include "unitcell.hpp"  // for UnitCell

namespace gemmi {

/// @brief Iterate over all reflections within a resolution range, applying a function to each.
/// @tparam Func Callable taking a Miller index array
/// @param func Function to apply to each valid reflection
/// @param cell Unit cell parameters
/// @param spacegroup Space group for symmetry and systematic absences
/// @param dmin Minimum d-spacing (Angstroms, lower resolution limit)
/// @param dmax Maximum d-spacing (Angstroms, higher resolution limit); 0 = unlimited; INFINITY = exact high-res limit
/// @param unique If true, iterate only over asymmetric unit (one per symmetry equivalent); if false, all in range
///
/// Generates Miller indices (h, k, l) in the resolution range [dmax, dmin] (or [dmin, inf) if dmax=0).
/// Checks space group symmetries and excludes systematically absent reflections.
/// @note dmin should include a small margin (e.g., 1e-6 Angstrom) to avoid numerical boundary issues.
template<typename Func>
void for_all_reflections(Func func,
                         const UnitCell& cell, const SpaceGroup* spacegroup,
                         double dmin, double dmax=0., bool unique=true) {
  Miller lim = cell.get_hkl_limits(dmin);
  double inv_dmin2 = 1. / sq(dmin);
  double inv_dmax2 = 0.;
  if (dmax > 0)
    inv_dmax2 = dmax == INFINITY ? -1 : 1. / sq(dmax);
  ReciprocalAsu asu(spacegroup);
  GroupOps gops = spacegroup->operations();
  Miller hkl;
  for (hkl[0] = -lim[0]; hkl[0] <= lim[0]; ++hkl[0])
    for (hkl[1] = -lim[1]; hkl[1] <= lim[1]; ++hkl[1])
      for (hkl[2] = -lim[2]; hkl[2] <= lim[2]; ++hkl[2])
        if (!unique || asu.is_in(hkl)) {
          double inv_d2 = cell.calculate_1_d2(hkl);
          if (inv_d2 <= inv_dmin2 && inv_d2 > inv_dmax2 &&
              !gops.is_systematically_absent(hkl))
            func(hkl);
        }
}

/// @brief Count the number of unique (or all) reflections in a resolution range.
/// @param cell Unit cell parameters
/// @param spacegroup Space group
/// @param dmin Minimum d-spacing (Angstroms)
/// @param dmax Maximum d-spacing (Angstroms); 0 = no upper limit
/// @param unique If true, count only asymmetric unit; if false, count all
/// @return Total number of reflections satisfying the criteria
///
/// @note dmin should include a small margin for numerical accuracy.
inline int count_reflections(const UnitCell& cell, const SpaceGroup* spacegroup,
                             double dmin, double dmax=0., bool unique=true) {
  int counter = 0;
  for_all_reflections([&counter](const Miller&) { ++counter; },
                      cell, spacegroup, dmin, dmax, unique);
  return counter;
}

/// @brief Generate a vector of all reflections in a resolution range.
/// @param cell Unit cell parameters
/// @param spacegroup Space group
/// @param dmin Minimum d-spacing (Angstroms)
/// @param dmax Maximum d-spacing (Angstroms); 0 = no upper limit
/// @param unique If true, return only asymmetric unit; if false, all in range
/// @return Vector of Miller indices [h, k, l] sorted in enumeration order
///
/// Useful for generating a complete or unique list of expected reflections for comparison with data.
inline std::vector<Miller>
make_miller_vector(const UnitCell& cell, const SpaceGroup* spacegroup,
                   double dmin, double dmax=0., bool unique=true) {
  std::vector<Miller> hkls;
  for_all_reflections([&hkls](const Miller& hkl) { hkls.push_back(hkl); },
                      cell, spacegroup, dmin, dmax, unique);
  return hkls;
}


} // namespace gemmi
#endif
