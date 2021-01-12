// Copyright 2017 Global Phasing Ltd.

// Microbenchmark of gemmi::Residue::has_standard_pdb_name().
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include resinfo.cpp $GB/src/libbenchmark.a -pthread

#include "gemmi/resinfo.hpp"
#include <benchmark/benchmark.h>
#include <stdlib.h>               // for rand

static const std::string residue_names[10] =
    { "GLY", "ASN", "ASN", "HOH", "LEU", "DC", "G", "GLU", "ALA", "SER" };

static void find_tabulated_residue_x10(benchmark::State& state) {
  std::string names[10];
  for (int i = 0; i != 10; ++i)
    names[i] = rand() % 1000 == 0 ? "GLN" : residue_names[i];
  while (state.KeepRunning()) {
    int nh = 0;
    for (int i = 0; i != 10; ++i)
      nh += gemmi::find_tabulated_residue(names[i]).hydrogen_count;
    benchmark::DoNotOptimize(nh);
  }
}

BENCHMARK(find_tabulated_residue_x10);
BENCHMARK_MAIN();

/* Output from my laptop:

------------------------------------------------------------------
Benchmark                           Time           CPU Iterations
------------------------------------------------------------------
find_tabulated_residue_x10         52 ns         52 ns   13458745

*/
