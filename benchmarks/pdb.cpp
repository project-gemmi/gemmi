// Copyright 2018 Global Phasing Ltd.

// Microbenchmark of gemmi::UnitCell::find_nearest_image().

#include "gemmi/model.hpp"
#include "gemmi/pdb.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/neighbor.hpp"
#include "gemmi/select.hpp"  // count_atom_sites
#include "gemmi/remarks.hpp"  // read_metadata_from_remarks
#include <benchmark/benchmark.h>

static const char* path;

static void read_pdb_file(benchmark::State& state) {
  for (auto _ : state) {
    gemmi::Structure st = gemmi::read_pdb_file(path);
    benchmark::DoNotOptimize(st);
  }
}

static void read_pdb_remarks(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    st.meta = gemmi::Metadata();
    read_metadata_from_remarks(st);
    benchmark::DoNotOptimize(st.meta);
  }
}

static void find_atom_image(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(0).atoms.at(0).pos;
  for (auto _ : state) {
    double sum_dist = 0;
    for (const Chain& ch1 : model.chains)
      for (const Residue& res1 : ch1.residues)
        for (const Atom& a1 : res1.atoms) {
          NearestImage im = st.cell.find_nearest_image(ref, a1.pos, Asu::Any);
          sum_dist += im.dist_sq;
        }
    benchmark::DoNotOptimize(sum_dist);
  }
}

static void neighbor_search_ctor(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
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
  for (auto _ : state) {
    auto r = ns.find_atoms(ref, '\0', 0, 4);
    benchmark::DoNotOptimize(r);
  }
}

static void neighbor_search_for_each(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  const Model& model = st.models[0];
  Position ref = model.chains.at(0).residues.at(2).atoms.at(0).pos;
  NeighborSearch ns(st.models.at(0), st.cell, 5.0);
  for (auto _ : state) {
    double sum = 0;
    ns.for_each(ref, '\0', 4,
                [&sum](NeighborSearch::Mark&, double d) { sum += d; });
    benchmark::DoNotOptimize(sum);
  }
}

static void calculate_box(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    Box<Position> box = calculate_box(st);
    benchmark::DoNotOptimize(box);
  }
}

static void fractional_box(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    Box<Fractional> box = calculate_fractional_box(st);
    benchmark::DoNotOptimize(box);
  }
}

static bool has_hydrogen_with_levels(const gemmi::Structure& st) {
  for (const gemmi::Model& model : st.models)
    for (const gemmi::Chain& chain : model.chains)
      for (const gemmi::Residue& res : chain.residues)
        for (const gemmi::Atom& atom : res.atoms)
          if (atom.is_hydrogen())
            return true;
  return false;
}

static bool has_hydrogen_with_cra(const gemmi::Structure& st) {
  for (const gemmi::Model& model : st.models)
    for (gemmi::const_CRA cra : model.all())
      if (cra.atom->is_hydrogen())
        return true;
  return false;
}

static bool has_hydrogen_with_selection(const gemmi::Structure& st) {
  gemmi::Selection sel("[H,D]");
  return count_atom_sites(st, &sel) != 0;
}

static void has_hydrogen1(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    bool has_hydr = has_hydrogen_with_levels(st);
    benchmark::DoNotOptimize(has_hydr);
  }
}

static void has_hydrogen2(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    bool has_hydr = has_hydrogen_with_cra(st);
    benchmark::DoNotOptimize(has_hydr);
  }
}

static void has_hydrogen3(benchmark::State& state) {
  using namespace gemmi;
  Structure st = read_pdb_file(path);
  for (auto _ : state) {
    bool has_hydr = has_hydrogen_with_selection(st);
    benchmark::DoNotOptimize(has_hydr);
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
  benchmark::RegisterBenchmark("read_pdb_remarks", read_pdb_remarks);
  benchmark::RegisterBenchmark("find_atom_image", find_atom_image);
  benchmark::RegisterBenchmark("neighbor_search_ctor", neighbor_search_ctor);
  benchmark::RegisterBenchmark("neighbor_search_find", neighbor_search_find);
  benchmark::RegisterBenchmark("neighbor_search_for_each",
                               neighbor_search_for_each);
  benchmark::RegisterBenchmark("calculate_box", calculate_box);
  benchmark::RegisterBenchmark("fractional_box", fractional_box);
  benchmark::RegisterBenchmark("has_hydrogen1", has_hydrogen1);
  benchmark::RegisterBenchmark("has_hydrogen2", has_hydrogen2);
  benchmark::RegisterBenchmark("has_hydrogen3", has_hydrogen3);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}

/* Output from my desktop:

*/
