// Copyright 2026 Global Phasing Ltd.
//
// Chemical-adjustment helpers used by ace_cc pipeline.

#ifndef GEMMI_CC_ADJ_HPP_
#define GEMMI_CC_ADJ_HPP_

#include "gemmi/chemcomp.hpp"

namespace gemmi {

// Chemical normalization pass used by prepare_chemcomp().
void apply_chemical_adjustments(ChemComp& cc);
bool add_n_terminal_h3(ChemComp& cc);
void sync_n_terminal_h3_angles(ChemComp& cc);

} // namespace gemmi

#endif
