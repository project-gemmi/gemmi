// Copyright 2017 Global Phasing Ltd.

// Microbenchmark of functions from symmetry.hpp
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include bench_sym.cpp $GB/src/libbenchmark.a -pthread

#include <benchmark/benchmark.h>
#include <gemmi/symmetry.hpp>

namespace sym = gemmi::sym;

constexpr const char* TRIPLETS[10] = {
  "x,y,z", "-x,-y+1/2,z", "x+1/4,y+1/4,z", "x+1/4,y+1/4,-x+z-1/4",
  "1/2+X,1/2-Y,1/2+Z", "-Z,-X,Y", "Y+1/2,X+1/2,Z+1/2",
  "2/3+X,1/3+X-Y,5/6+Z", "-x+1/2,-y+1/2,z", " -x,-y, -z" };

constexpr const char* HALL_SYMBOLS[] = {
  "p 1", "C 2y", "I 2y", "C 2y (x,y,-x+z)",
  "P 2ac 2ab", "F 4d 2 3 -1cd", "-F 4ud 2vw 3 (x-1/8,y-1/8,z-1/8)"
};

static void bm_parse_triplet_all(benchmark::State& state) {
  while (state.KeepRunning())
    for (const char* triplet : TRIPLETS)
      benchmark::DoNotOptimize(sym::parse_triplet(triplet));
}

static void bm_make_triplet(benchmark::State& state) {
  int n = state.range(0);
  sym::Op op = sym::parse_triplet(TRIPLETS[n]);
  while (state.KeepRunning())
    benchmark::DoNotOptimize(op.triplet());
}

static void bm_make_triplet_all(benchmark::State& state) {
  sym::Op ops[10];
  for (int i = 0; i != 10; ++i)
    ops[i] = sym::parse_triplet(TRIPLETS[i]);
  while (state.KeepRunning())
    for (const sym::Op& op : ops)
      benchmark::DoNotOptimize(op.triplet());
}

static void bm_generators_from_hall(benchmark::State& state) {
  int n = state.range(0);
  const char* hall = HALL_SYMBOLS[n];
  while (state.KeepRunning())
    benchmark::DoNotOptimize(sym::generators_from_hall(hall));
}

static void bm_add_elements(benchmark::State& state) {
  int n = state.range(0);
  const char* hall = HALL_SYMBOLS[n];
  sym::GroupOps gops = sym::generators_from_hall(hall);
  while (state.KeepRunning()) {
    sym::GroupOps copy = gops;
    copy.add_missing_elements();
    benchmark::DoNotOptimize(copy);
  }
}

BENCHMARK(bm_parse_triplet_all);
//BENCHMARK(bm_make_triplet)->DenseRange(0, 9);
BENCHMARK(bm_make_triplet_all);
BENCHMARK(bm_generators_from_hall)->DenseRange(0, 6);
BENCHMARK(bm_add_elements)->DenseRange(0, 6);
BENCHMARK_MAIN();

// sym::parse_triplet(): 50-300 ns/triplet
// sym::make_triplet(): 40-120 ns/triplet
// sym::generators_from_hall(): 60-1000 ns
// GroupOps::add_missing_elements(): 140-90,000 ns

// vim:sw=2:ts=2:et:path^=../include,../third_party
