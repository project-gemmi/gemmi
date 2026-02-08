// Copyright 2025 Global Phasing Ltd.
//
// prepare_chemcomp() — restraint generation pipeline.

#ifndef GEMMI_ACE_CC_HPP_
#define GEMMI_ACE_CC_HPP_

#include <map>
#include <string>
#include "chemcomp.hpp"  // for ChemComp

namespace gemmi {

struct AcedrgTables;

/// Run the full restraint-generation pipeline on a ChemComp:
/// chemical-group adjustments, protonation, fill_restraints,
/// torsion/chirality/plane generation, and CCP4 type assignment.
/// \param atom_stereo  maps atom names to pdbx_stereo_config strings
///                     (needed for chirality generation).
GEMMI_DLL void prepare_chemcomp(ChemComp& cc, const AcedrgTables& tables,
                                const std::map<std::string, std::string>& atom_stereo = {});

} // namespace gemmi
#endif
