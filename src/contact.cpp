// Copyright 2018 Global Phasing Ltd.
//
// Searches for contacts -- neighbouring atoms.

#include <gemmi/cellmethod.hpp>
#include "input.h"
#define GEMMI_PROG contact
#include "options.h"
#include <stdio.h>
#include <algorithm>  // for min, max

using namespace gemmi;

enum OptionIndex { Verbose=3 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nSearches for contacts in a model (PDB or mmCIF)."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_contacts(const Structure& st, bool /*verbose*/) {
  CellMethod cm(st, 5);
  printf(" Cell grid: %d x %d x %d\n", cm.grid.nu, cm.grid.nv, cm.grid.nw);
  size_t min_count = SIZE_MAX, max_count = 0, total_count = 0;
  for (const auto& el : cm.grid.data) {
    min_count = std::min(min_count, el.size());
    max_count = std::max(max_count, el.size());
    total_count += el.size();
  }
  printf(" Items per cell: from %zu to %zu, average: %.2g\n",
         min_count, max_count, double(total_count) / cm.grid.data.size());
  // TODO
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  bool verbose = p.options[Verbose];
  if (p.nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    return 0;
  }
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.nonOption(i);
      if (is_pdb_code(input))
        input = expand_pdb_code_to_path_or_fail(input);
      if (i > 0)
        std::printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        std::printf("File: %s\n", input.c_str());
      Structure st = read_structure(input);
      print_contacts(st, verbose);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
