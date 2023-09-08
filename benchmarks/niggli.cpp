// Copyright 2022 Global Phasing Ltd.

// Microbenchmark of functions from cellred.hpp
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include niggli.cpp $GB/src/libbenchmark.a -pthread

#include <benchmark/benchmark.h>
#include <gemmi/cellred.hpp>

#if 0
static const gemmi::UnitCell cell(52.237, 55.222, 83.046, 90.000, 96.160, 90.000);
static const char centring = 'P';
#elif 0
static const gemmi::UnitCell cell(172.1, 172.1, 80.0, 90., 90., 120.);
static const char centring = 'R';
#else
static const gemmi::UnitCell cell(63.78, 63.86, 124.40, 90.0, 90.0, 90.0);
static const char centring = 'I';
#endif

void print_iterations() {
  printf("Niggli: %d iterations\n", gemmi::GruberVector(cell, centring).niggli_reduce());
  printf("Buerger: %d iterations\n", gemmi::GruberVector(cell, centring).buerger_reduce());
  printf("Selling: %d iterations\n", gemmi::SellingVector(cell, centring).reduce());
}

static void run_niggli(benchmark::State& state) {
  //print_iterations();
  for (auto _ : state) {
    gemmi::GruberVector gv(cell, centring, false);
    gv.niggli_reduce();
    benchmark::DoNotOptimize(gv.parameters());
  }
}

static void run_niggli_with_tracking(benchmark::State& state) {
  for (auto _ : state) {
    gemmi::GruberVector gv(cell, centring, true);
    gv.niggli_reduce();
    benchmark::DoNotOptimize(gv.parameters());
  }
}

static void run_buerger(benchmark::State& state) {
  for (auto _ : state) {
    gemmi::GruberVector gv(cell, centring);
    gv.buerger_reduce();
    benchmark::DoNotOptimize(gv.parameters());
  }
}

static void run_selling(benchmark::State& state) {
  for (auto _ : state) {
    gemmi::SellingVector sv(cell, centring);
    sv.reduce();
    benchmark::DoNotOptimize(sv.g6_parameters());
  }
}

BENCHMARK(run_niggli);
BENCHMARK(run_niggli_with_tracking);
BENCHMARK(run_buerger);
BENCHMARK(run_selling);
BENCHMARK_MAIN();
