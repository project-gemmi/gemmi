// Copyright 2018 Global Phasing Ltd.

// Microbenchmark of gemmi::UnitCell::find_nearest_image().

#include "gemmi/pdb.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/subcells.hpp"
#include <benchmark/benchmark.h>
#include <algorithm>

static const char* path;

static void read_pdb_file(benchmark::State& state) {
  while (state.KeepRunning()) {
    gemmi::Structure st = gemmi::read_pdb_file(path);
    benchmark::DoNotOptimize(st);
  }
}

static void find_atom_image(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(0).atoms.at(0).pos;
  while (state.KeepRunning()) {
    double sum_dist = 0;
    for (const Chain& ch1 : model.chains)
      for (const Residue& res1 : ch1.residues)
        for (const Atom& a1 : res1.atoms) {
          SymImage im = st.cell.find_nearest_image(ref, a1.pos, Asu::Any);
          sum_dist += im.dist_sq;
        }
    benchmark::DoNotOptimize(sum_dist);
  }
}

static void subcells_ctor(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  while (state.KeepRunning()) {
    SubCells sc(st.models.at(0), st.cell, 5.0);
    benchmark::DoNotOptimize(sc);
  }
}

static void subcells_find(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(2).atoms.at(0).pos;
  SubCells sc(st.models.at(0), st.cell, 5.0);
  while (state.KeepRunning()) {
    auto r = sc.find_atoms(ref, '\0', 4);
    benchmark::DoNotOptimize(r);
  }
}

static void subcells_for_each(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(2).atoms.at(0).pos;
  SubCells sc(st.models.at(0), st.cell, 5.0);
  while (state.KeepRunning()) {
    double sum = 0;
    sc.for_each(ref, '\0', 4, [&sum](SubCells::Mark&, float d) { sum += d; });
    benchmark::DoNotOptimize(sum);
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Call it with path to a pdb file as an argument.\n");
    return 1;
  }
  path = argv[argc-1];
  {
    gemmi::Structure st = gemmi::read_pdb_file(path);
    printf("PDB file: %s with %zu atom sites.\n",
           st.name.c_str(), count_atom_sites(st.models.at(0)));
  }
  benchmark::RegisterBenchmark("read_pdb_file", read_pdb_file);
  benchmark::RegisterBenchmark("find_atom_image", find_atom_image);
  benchmark::RegisterBenchmark("subcells_ctor", subcells_ctor);
  benchmark::RegisterBenchmark("subcells_find", subcells_find);
  benchmark::RegisterBenchmark("subcells_for_each", subcells_for_each);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}

/* Output from my desktop:

*/

// vim:sw=2:ts=2:et:path^=../include,../third_party
