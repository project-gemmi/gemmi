// Copyright 2026 Global Phasing Ltd.
// Internal carborane processing for AceDRG-style pipeline.

#ifndef GEMMI_ACE_CARBORANE_HPP_
#define GEMMI_ACE_CARBORANE_HPP_

#include "chemcomp.hpp"
#include "ace_graph.hpp"
#include <set>
#include <string>

namespace gemmi {

bool has_carborane_seed(const ChemComp& cc, const AceBondAdjacency& adj);
bool is_carborane_mode_component(const ChemComp& cc, const AceBondAdjacency& adj);
void apply_carborane_mode(ChemComp& cc, bool no_angles);
void apply_mixed_carborane_mode(ChemComp& cc, bool no_angles,
                                const std::string& tables_dir);

} // namespace gemmi

#endif
