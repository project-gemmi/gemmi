// Copyright 2026 Global Phasing Ltd.
//
// CCP4 energy-type assignment for ChemComp atoms.

#ifndef GEMMI_CCP4ENER_HPP_
#define GEMMI_CCP4ENER_HPP_

#include "chemcomp.hpp"

namespace gemmi {

// Assign CCP4 atom energy types (_chem_comp_atom.type_energy) in-place.
// This uses local graph analysis only and does not require AceDRG tables.
GEMMI_DLL void assign_chemcomp_ccp4_types(ChemComp& cc);

}  // namespace gemmi

#endif
