// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"
#include <stdio.h>

#define GEMMI_PROG sg
#include "options.h"

// enum OptionIndex { Verbose=3 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] SPACEGROUP[...]"
    "\nPrints information about the space group."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  //{ Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_symmetry_operations(const gemmi::GroupOps& ops) {
  printf("%zu x %zu symmetry operations:\n",
         ops.cen_ops.size(), ops.sym_ops.size());
  for (const gemmi::Op& op : ops)
    printf("    %s\n", op.triplet().c_str());
}

static void process_arg(const char* arg) {
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name(arg);
  if (sg == nullptr) {
    try {
      gemmi::GroupOps ops = gemmi::symops_from_hall(arg);
      sg = gemmi::find_spacegroup_by_ops(ops);
      if (sg == nullptr) {
        printf("Hall symbol: %s\n", arg);
        print_symmetry_operations(ops);
      }
    } catch (std::runtime_error&) {
    }
  }
  if (sg == nullptr) {
    fprintf(stderr, "Space group not found: %s\n", arg);
    return;
  }
  printf("Number: %d\n", sg->number);
  printf("CCP4 number: %d\n", sg->ccp4);
  printf("Hermann–Mauguin: %s\n", sg->hm);
  printf("Extended H-M: %s\n", sg->xhm().c_str());
  printf("Hall symbol: %s\n", sg->hall);
  printf("Point group: %s\n", sg->point_group_hm());
  gemmi::GroupOps ops = sg->operations();
  printf("Is centric: %s\n", ops.is_centric() ? "yes" : "no");
  print_symmetry_operations(ops);
  printf("\n");
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  for (int i = 0; i < p.nonOptionsCount(); ++i)
    process_arg(p.nonOption(i));
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
