// Copyright 2026 Global Phasing Ltd.
//
// Generate idealized ChemComp coordinates from restraint dictionaries.

#ifndef GEMMI_CHEMCOMP_XYZ_HPP_
#define GEMMI_CHEMCOMP_XYZ_HPP_

#include "chemcomp.hpp"

namespace gemmi {

/// @brief Generate idealized 3D coordinates for a chemical component.
/// Generates a deterministic idealized conformer by applying bond lengths,
/// angles, and torsion restraints in sequence. Modifies cc.atoms[*].xyz in-place.
/// @param cc ChemComp to generate coordinates for; atoms must be present.
/// @return Number of atoms assigned finite coordinates.
/// @note Atoms without restraints may remain uninitialized (NAN coordinates).
GEMMI_DLL int generate_chemcomp_xyz_from_restraints(ChemComp& cc);

/// @brief Refine chemical component coordinates against restraints.
/// Refines atom coordinates against bond and angle restraints using
/// Levenberg-Marquardt least squares optimization. Modifies cc.atoms[*].xyz in-place.
/// @param cc ChemComp with initial coordinates to refine.
/// @return Final weighted sum of squared residuals (WSSR) of the fit.
/// @note Requires atoms to have initial finite coordinates (e.g., from generate_chemcomp_xyz_from_restraints).
GEMMI_DLL double refine_chemcomp_xyz(ChemComp& cc);

} // namespace gemmi

#endif
