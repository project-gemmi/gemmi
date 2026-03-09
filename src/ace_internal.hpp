// Copyright 2026 Global Phasing Ltd.
// Shared internal helpers for AceDRG-style processing.

#ifndef GEMMI_SRC_ACE_INTERNAL_HPP_
#define GEMMI_SRC_ACE_INTERNAL_HPP_

#include "gemmi/chemcomp.hpp"
#include "gemmi/ace_cc.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>

namespace gemmi {

inline bool env_flag_enabled(const char* name) {
  if (const char* env = std::getenv(name))
    return env[0] != '\0' && env[0] != '0';
  return false;
}

inline bool option_flag(const PrepareOverride* setting, const char* env_name) {
  if (setting) {
    if (*setting == PrepareOverride::Enable)
      return true;
    if (*setting == PrepareOverride::Disable)
      return false;
  }
  return env_flag_enabled(env_name);
}

// Global/thread-local options state used by internal functions.
extern thread_local const PrepareChemcompOptions* active_prepare_options;

inline bool ace_strict_mode() {
  return option_flag(active_prepare_options ? &active_prepare_options->strict_mode : nullptr,
                     "GEMMI_ACE_STRICT");
}

inline bool ace_compat_mode() {
  return option_flag(active_prepare_options ? &active_prepare_options->compat_mode : nullptr,
                     "GEMMI_ACE_COMPAT");
}

inline bool ace_trace_mode() {
  return option_flag(active_prepare_options ? &active_prepare_options->trace_mode : nullptr,
                     "GEMMI_ACE_TRACE");
}

struct AceRuleStats {
  int atom_count = 0;
  int bond_count = 0;
  int angle_count = 0;
  int torsion_count = 0;
  int chir_count = 0;
  int plane_count = 0;
  double charge_sum = 0.0;
};

inline AceRuleStats collect_rule_stats(const ChemComp& cc) {
  AceRuleStats s;
  s.atom_count = static_cast<int>(cc.atoms.size());
  s.bond_count = static_cast<int>(cc.rt.bonds.size());
  s.angle_count = static_cast<int>(cc.rt.angles.size());
  s.torsion_count = static_cast<int>(cc.rt.torsions.size());
  s.chir_count = static_cast<int>(cc.rt.chirs.size());
  s.plane_count = static_cast<int>(cc.rt.planes.size());
  for (const auto& atom : cc.atoms)
    s.charge_sum += atom.charge;
  return s;
}

inline bool rule_stats_changed(const AceRuleStats& before, const AceRuleStats& after) {
  if (before.atom_count != after.atom_count) return true;
  if (before.bond_count != after.bond_count) return true;
  if (before.angle_count != after.angle_count) return true;
  if (before.torsion_count != after.torsion_count) return true;
  if (before.chir_count != after.chir_count) return true;
  if (before.plane_count != after.plane_count) return true;
  return std::fabs(before.charge_sum - after.charge_sum) > 1e-6;
}

inline void trace_phase_delta(const char* phase, const AceRuleStats& before, const AceRuleStats& after) {
  if (!ace_trace_mode())
    return;
  std::fprintf(stderr,
               "[ace-phase %s] atoms %+d bonds %+d angles %+d torsions %+d chirs %+d planes %+d charge %+0.3f\n",
               phase,
               after.atom_count - before.atom_count,
               after.bond_count - before.bond_count,
               after.angle_count - before.angle_count,
               after.torsion_count - before.torsion_count,
               after.chir_count - before.chir_count,
               after.plane_count - before.plane_count,
               after.charge_sum - before.charge_sum);
}

} // namespace gemmi

#endif
