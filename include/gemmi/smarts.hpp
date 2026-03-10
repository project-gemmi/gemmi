// Copyright 2026 Global Phasing Ltd.
//
// Lightweight SMARTS-subset matcher for ChemComp.

#ifndef GEMMI_SMARTS_HPP_
#define GEMMI_SMARTS_HPP_

#include <string>
#include <vector>
#include "gemmi/chemcomp.hpp"

namespace gemmi {

// A match is a vector of atom indices in the ChemComp, corresponding
// to the atoms in the SMARTS pattern.
using SmartsMatch = std::vector<int>;

// Find all non-overlapping or overlapping matches of a SMARTS pattern.
// This implementation supports a small subset of SMARTS:
// - Atomic symbols: [C], [N], [O], etc. (or C, N, O without brackets)
// - Wildcards: *
// - Aromaticity: [c], [n], etc.
// - Constraints: H<n> (hydrogen count), X<n> (connectivity)
// - Bonds: - (single), = (double), ~ (any)
// - Branching: ( )
GEMMI_DLL std::vector<SmartsMatch> match_smarts(const ChemComp& cc, const std::string& pattern);

} // namespace gemmi

#endif
