// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"
#include <stdio.h>

void print_symmetry_operations(const gemmi::GroupOps& ops) {
  printf("%zu x %zu symmetry operations:\n",
         ops.cen_ops.size(), ops.sym_ops.size());
  for (const gemmi::Op& op : ops)
    printf("    %s\n", op.triplet().c_str());
}

void process_arg(const char* arg) {
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
  printf("extended H-M: %s\n", sg->xhm().c_str());
  printf("Hall symbol: %s\n", sg->hall);
  print_symmetry_operations(sg->operations());
  printf("\n");
}

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i)
    process_arg(argv[i]);
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
