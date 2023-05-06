// Copyright 2020 Global Phasing Ltd.
//
// Reindex merged or unmerged MTZ file.

#include <cstdio>
#include <cstring>  // for strpbrk
#include <algorithm>
#include <iostream>  // for cerr
#include <gemmi/mtz.hpp>
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#define GEMMI_PROG reindex
#include "options.h"

using std::fprintf;

namespace {

enum OptionIndex { Hkl=4, NoHistory, NoSort, Asu };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_MTZ OUTPUT_MTZ"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Hkl, 0, "", "hkl", Arg::Required,
    "  --hkl=OP  \tReindexing transform as triplet (e.g. k,h,-l)." },
  { NoHistory, 0, "", "no-history", Arg::None,
    "  --no-history  \tDo not add 'Reindexed with...' line to mtz HISTORY." },
  { NoSort, 0, "", "no-sort", Arg::None,
    "  --no-sort  \tDo not reorder reflections." },
  { Asu, 0, "", "asu", Arg::AsuChoice,
    "  --asu=ccp4|tnt  \tWrite merged data in CCP4 (default) or TNT ASU." },
  { NoOp, 0, "", "", Arg::None,
    "\nInput file can be gzipped." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);
  if (!p.options[Hkl] && !p.options[Asu]) {
    fprintf(stderr, "Specify transform with option --hkl\n");
    return 1;
  }
  try {
    std::string from_line = "From gemmi-reindex " GEMMI_VERSION;
    gemmi::Op op;
    if (p.options[Hkl]) {
      std::string hkl_arg = p.options[Hkl].arg;
      op = gemmi::parse_triplet(hkl_arg);
      if (std::strpbrk(hkl_arg.c_str(), "xyzabcXYZABC"))
        gemmi::fail("specify OP in terms of h, k and l");
      if (op.tran != gemmi::Op::Tran{{0, 0, 0}})
        gemmi::fail("reindexing operator should not have a translation");
      gemmi::cat_to(from_line, " with [", hkl_arg, ']');
    }

    gemmi::Mtz mtz;
    if (verbose) {
      fprintf(stderr, "Reading %s ...\n", input_path);
      mtz.warnings = &std::cerr;
    }
    mtz.read_input(gemmi::MaybeGzipped(input_path), true);

    if (p.options[Hkl])
      mtz.reindex(op, &std::cerr);

    if (mtz.is_merged()) {
      bool tnt_asu = false;
      if (p.options[Asu] && p.options[Asu].arg[0] == 't')
        tnt_asu = true;
      mtz.ensure_asu(tnt_asu);
    }

    if (!p.options[NoSort])
      mtz.sort();
    if (!p.options[NoHistory])
      mtz.history.emplace(mtz.history.begin(), from_line);
    if (verbose)
      fprintf(stderr, "Writing %s ...\n", output_path);
    mtz.write_to_file(output_path);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
