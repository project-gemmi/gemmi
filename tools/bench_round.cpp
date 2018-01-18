// Copyright 2017 Global Phasing Ltd.

#include <cmath>
#include <benchmark/benchmark.h>

static double use_round(double x) { return std::round(x); }
static double use_floor(double x) { return std::floor(x + 0.5); }
static double use_nearbyint(double x) { return std::nearbyint(x); }
static double use_rint(double x) { return std::rint(x); }
static long use_lrint(double x) { return std::lrint(x); }

static long use2_round(double x) { return (long) std::round(x); }
static long use2_floor(double x) { return (long) std::floor(x + 0.5); }
static long use2_nearbyint(double x) { return (long) std::nearbyint(x); }
static long use2_rint(double x) { return (long) std::rint(x); }

double numbers[] = {
  1.331, 6.3, 4, -4.5, 99341.4106,
  0, -6.05823, 1e-9, 18.419, -12345.67
};

inline void run(benchmark::State& state, double(*func)(double)) {
  for (double& x : numbers) { // otherwise it would be all optimized out
    x = std::rand() % 10 == 0 ? -x : x;
  }
  while (state.KeepRunning())
    for (int x : numbers)
      benchmark::DoNotOptimize((*func)(x));
}

inline void run2(benchmark::State& state, long(*func)(double)) {
  for (double& x : numbers) { // otherwise it would be all optimized out
    x = std::rand() % 10 == 0 ? -x : x;
  }
  while (state.KeepRunning())
    for (int x : numbers)
      benchmark::DoNotOptimize((*func)(x));
}


BENCHMARK_CAPTURE(run, round, use_round);
BENCHMARK_CAPTURE(run, floor, use_floor);
BENCHMARK_CAPTURE(run, neabyint, use_nearbyint);
BENCHMARK_CAPTURE(run, rint, use_rint);
BENCHMARK_CAPTURE(run2, lrint, use_lrint);
BENCHMARK_CAPTURE(run2, round, use2_round);
BENCHMARK_CAPTURE(run2, floor, use2_floor);
BENCHMARK_CAPTURE(run2, neabyint, use2_nearbyint);
BENCHMARK_CAPTURE(run2, rint, use2_rint);
BENCHMARK_MAIN()

/*
-----------------------------------------------------
Benchmark              Time           CPU Iterations
-----------------------------------------------------
run/round              8 ns          8 ns   86728922
run/floor             25 ns         25 ns   27625341
run/neabyint           8 ns          8 ns   87853604
run/rint               8 ns          8 ns   88405267
run2/lrint            18 ns         18 ns   38545176
run2/round             5 ns          5 ns  139178181
run2/floor            18 ns         18 ns   38108709
run2/neabyint          5 ns          5 ns  135124876
run2/rint              5 ns          5 ns  143791413
*/

// vim:sw=2:ts=2:et:path^=../include,../third_party
