// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "input.h"
#include <cstdlib>  // for strtod
#include <cstdio>  // for fprintf
#define EXE_NAME "gemmi-map"
#include "options.h"

namespace mol = gemmi::mol;

enum OptionIndex { Unknown, Verbose, OutputMode };

static const option::Descriptor Usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] CCP4_MAP[...]\n" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

int main(int argc, char **argv) {
  OptParser parse;
  auto options = parse.simple_parse(argc, argv, Usage);
  bool verbose = options[Verbose];

  if (parse.nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    return 0;
  }

  try {
    for (int i = 0; i < parse.nonOptionsCount(); ++i) {
      const char* input = parse.nonOption(i);
      gemmi::Grid<> grid;
      if (verbose)
        std::fprintf(stderr, "Reading %s ...\n", input);
      grid.read_ccp4(input);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
