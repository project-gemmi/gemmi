// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"
#include <stdio.h>

namespace sym = gemmi::sym;

void process_arg(const char* arg) {
  const sym::SpaceGroup* sg = sym::find_spacegroup_by_name(arg);
  if (sg == nullptr) {
    fprintf(stderr, "H-M name not found: %s\n", arg);
    return;
  }
  printf("Number: %d\n", sg->number);
  printf("CCP4 number: %d\n", sg->ccp4);
  printf("Hermann–Mauguin: %s\n", sg->hm);
  printf("extended H-M: %s\n", sg->xhm().c_str());
  printf("Hall symbol: %s\n", sg->hall);
  printf("Symmetry operations:\n");
  for (const sym::Op& op : sg->operations())
    printf("    %s\n", op.triplet().c_str());
  printf("\n");
}

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i)
    process_arg(argv[i]);
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
