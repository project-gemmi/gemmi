// Copyright 2026 Global Phasing Ltd.
//
// Generate idealized ChemComp coordinates from restraint dictionaries.

#ifndef GEMMI_CHEMCOMP_XYZ_HPP_
#define GEMMI_CHEMCOMP_XYZ_HPP_

#include "chemcomp.hpp"

namespace gemmi {

/// Generate a deterministic idealized conformer from bond/angle/torsion
/// restraints. Returns the number of atoms assigned finite coordinates.
GEMMI_DLL int generate_chemcomp_xyz_from_restraints(ChemComp& cc);

/// Refine monomer coordinates against bond and angle restraints
/// using Levenberg-Marquardt least squares. Returns final WSSR.
GEMMI_DLL double refine_chemcomp_xyz(ChemComp& cc);

} // namespace gemmi

#endif
