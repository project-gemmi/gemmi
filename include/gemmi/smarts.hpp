// Copyright 2026 Global Phasing Ltd.
//
// Lightweight SMARTS-subset matcher for ChemComp.

#ifndef GEMMI_SMARTS_HPP_
#define GEMMI_SMARTS_HPP_

#include <string>
#include <vector>
#include "chemcomp.hpp"

namespace gemmi {

/// @brief A SMARTS pattern match represented as atom indices.
/// Each match is a vector of indices into the ChemComp's atom list,
/// corresponding in order to atoms in the SMARTS pattern.
using SmartsMatch = std::vector<int>;

/// @brief Find all matches of a SMARTS pattern in a chemical component.
/// @details
/// This implementation supports a small subset of SMARTS notation:
/// - Atomic symbols: [C], [N], [O], etc. (or C, N, O without brackets)
/// - Wildcards: * (matches any atom)
/// - Aromaticity: [c], [n], etc. (aromatic atoms)
/// - Constraints: H<n> (hydrogen count), X<n> (connectivity/degree)
/// - Bonds: - (single), = (double), ~ (any)
/// - Branching: ( ) for subgraph grouping
/// @param cc The chemical component to search.
/// @param pattern The SMARTS pattern string.
/// @return Vector of all matches found; may include overlapping matches.
GEMMI_DLL std::vector<SmartsMatch> match_smarts(const ChemComp& cc, const std::string& pattern);

} // namespace gemmi

#endif
