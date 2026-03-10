// Copyright 2026 Global Phasing Ltd.
//
// Chemical-adjustment helpers used by ace_cc pipeline.

#ifndef GEMMI_CC_ADJ_HPP_
#define GEMMI_CC_ADJ_HPP_

#include "gemmi/chemcomp.hpp"

namespace gemmi {

// Chemical normalization pass used by prepare_chemcomp() and general users.
// It normalizes functional groups (nitro, carboxylate, phosphate, etc.)
// to a consistent protonation and resonance state.
GEMMI_DLL void normalize_chemcomp(ChemComp& cc);

// Deprecated alias for normalize_chemcomp.
inline void apply_chemical_adjustments(ChemComp& cc) { normalize_chemcomp(cc); }

GEMMI_DLL bool add_n_terminal_h3(ChemComp& cc);
GEMMI_DLL void sync_n_terminal_h3_angles(ChemComp& cc);

} // namespace gemmi

#endif
