// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"
#include "gemmi/grid.hpp"
#include "gemmi/asumask.hpp"  // for get_asu_mask
#include <cstdio>
#include <cstdlib>  // for atoi

#define GEMMI_PROG sg
#include "options.h"

using std::printf;

namespace {

enum OptionIndex { Asu=4 };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] SPACEGROUP[...]"
    "\nPrints information about the space group."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Asu, 0, "", "asu", Arg::Int,
    "  --asu=N  \tDraw ASU in NxNxN map grid and exit. Uses N(N+1) columns." },
  { 0, 0, 0, 0, 0, 0 }
};

void print_symmetry_operations(const gemmi::GroupOps& ops) {
  printf("%zu x %zu symmetry operations:\n",
         ops.cen_ops.size(), ops.sym_ops.size());
  for (gemmi::Op op : ops)
    printf("    %s\n", op.triplet().c_str());
}

void print_verbose_info(const char* hall) {
  using gemmi::Op;
  printf("The operations are generated from Hall symbol: %s\n", hall);
  gemmi::GroupOps ops = gemmi::generators_from_hall(hall);
  printf("%zu centering vector(s):\n", ops.cen_ops.size());
  for (const Op::Tran& cenop : ops.cen_ops)
    printf("    %s\n", Op{Op::identity().rot, cenop}.triplet().c_str());
  printf("%zu generator(s) of primitive symops (not counting identity):\n",
         ops.sym_ops.size() - 1);
  for (size_t i = 1; i < ops.sym_ops.size(); ++i)
    printf("    %s\n", ops.sym_ops[i].triplet().c_str());
  ops.add_missing_elements();
  printf("give %zu primitive symmetry operation(s):\n", ops.sym_ops.size());
  for (const Op& symop : ops.sym_ops)
    printf("    %s\n", symop.triplet().c_str());
}

void draw_asu(const gemmi::SpaceGroup* sg, int n) {
  gemmi::Grid<float> grid;
  grid.spacegroup = sg;
  grid.set_size(n, n, n);
  std::vector<std::int8_t> mask = gemmi::get_asu_mask(grid);
  int idx = 0;
  for (int w = 0; w != n; ++w) {
    for (int v = 0; v != n; ++v) {
      for (int u = 0; u != n; ++u, ++idx)
        std::putchar(mask[idx] == 0 ? '+' : '.');
      std::putchar(' ');
    }
    std::putchar('\n');
  }
}

const gemmi::SpaceGroup* find_spacegroup(const char* arg, bool verbose) {
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name(arg);
  if (sg == nullptr) {
    try {
      gemmi::GroupOps ops = gemmi::symops_from_hall(arg);
      sg = gemmi::find_spacegroup_by_ops(ops);
      if (sg == nullptr) {
        printf("Hall symbol: %s\n", arg);
        print_symmetry_operations(ops);
        if (verbose)
          print_verbose_info(arg);
      }
    } catch (std::runtime_error&) {
    }
  }
  return sg;
}

void print_info(const gemmi::SpaceGroup* sg, bool verbose) {
  printf("Number: %d\n", sg->number);
  bool is_reference = sg->is_reference_setting();
  printf("Is standard setting for this space group: %s\n",
         is_reference ? "yes" : "no");

  printf("Change-of-basis operator to standard setting: %s\n",
         sg->basisop_str());
  printf("CCP4 number: %d\n", sg->ccp4);
  printf("Hermann-Mauguin: %s\n", sg->hm);
  printf("Extended H-M: %s\n", sg->xhm().c_str());
  printf("Short name: %s\n", sg->short_name().c_str());
  printf("Hall symbol: %s\n", sg->hall);
  printf("Point group: %s\n", sg->point_group_hm());
  printf("Laue class: %s\n", sg->laue_str());
  printf("Crystal system: %s\n", sg->crystal_system_str());
  gemmi::GroupOps ops = sg->operations();
  printf("Is centrosymmetric: %s\n", ops.is_centrosymmetric() ? "yes" : "no");
  printf("Is enantiomorphic: %s\n", sg->is_enantiomorphic() ? "yes" : "no");
  std::array<int, 3> gf = ops.find_grid_factors();
  printf("Grid restrictions: NX=%dn NY=%dn NZ=%dn\n", gf[0], gf[1], gf[2]);
  for (bool tnt : {false, true})
    printf("Reciprocal space ASU (%s)%s: %s%s\n",
           tnt ? "TNT" : "CCP4",
           is_reference ? "" : " wrt. standard setting",
           tnt ? " " : "",
           gemmi::ReciprocalAsu(sg, tnt).condition_str());
  gemmi::AsuBrick brick = gemmi::find_asu_brick(sg);
  printf("Direct space ASU brick: %s\n", brick.str().c_str());
  print_symmetry_operations(ops);
  if (verbose)
    print_verbose_info(sg->hall);
  printf("\n");
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  bool verbose = p.options[Verbose];
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* arg = p.nonOption(i);
    const gemmi::SpaceGroup* sg = find_spacegroup(arg, verbose);
    if (sg == nullptr) {
      std::fprintf(stderr, "Space group not found: %s\n", arg);
      continue;
    }
    try {
      if (p.options[Asu])
        draw_asu(sg, std::atoi(p.options[Asu].arg));
      else
        print_info(sg, verbose);
    } catch (std::runtime_error& e) {
      std::fprintf(stderr, "ERROR: %s\n", e.what());
      return 1;
    }
  }
  return 0;
}
