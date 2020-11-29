// Copyright 2018 Global Phasing Ltd.

// Microbenchmark of gemmi::UnitCell::find_nearest_image().

#include "gemmi/model.hpp"
#include "gemmi/pdb.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/neighbor.hpp"
#include <benchmark/benchmark.h>

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

static void neighbor_search_ctor(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  while (state.KeepRunning()) {
    NeighborSearch ns(st.models.at(0), st.cell, 5.0);
    benchmark::DoNotOptimize(ns);
  }
}

static void neighbor_search_find(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(2).atoms.at(0).pos;
  NeighborSearch ns(st.models.at(0), st.cell, 5.0);
  while (state.KeepRunning()) {
    auto r = ns.find_atoms(ref, '\0', 4);
    benchmark::DoNotOptimize(r);
  }
}

static void neighbor_search_for_each(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(2).atoms.at(0).pos;
  NeighborSearch ns(st.models.at(0), st.cell, 5.0);
  while (state.KeepRunning()) {
    double sum = 0;
    ns.for_each(ref, '\0', 4,
                [&sum](NeighborSearch::Mark&, float d) { sum += d; });
    benchmark::DoNotOptimize(sum);
  }
}

static void calculate_box(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  while (state.KeepRunning()) {
    Box<Position> box = calculate_box(st);
    benchmark::DoNotOptimize(box);
  }
}

static void fractional_box(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  while (state.KeepRunning()) {
    Box<Fractional> box = calculate_fractional_box(st);
    benchmark::DoNotOptimize(box);
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
  benchmark::RegisterBenchmark("neighbor_search_ctor", neighbor_search_ctor);
  benchmark::RegisterBenchmark("neighbor_search_find", neighbor_search_find);
  benchmark::RegisterBenchmark("neighbor_search_for_each",
                               neighbor_search_for_each);
  benchmark::RegisterBenchmark("calculate_box", calculate_box);
  benchmark::RegisterBenchmark("fractional_box", fractional_box);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}

/* Output from my desktop:

*/
