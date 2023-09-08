// Copyright 2017 Global Phasing Ltd.

#include <cmath>
#include <cstdlib> // for rand
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
  for (auto _ : state)
    for (double x : numbers)
      benchmark::DoNotOptimize((*func)(x));
}

inline void run2(benchmark::State& state, long(*func)(double)) {
  for (double& x : numbers) { // otherwise it would be all optimized out
    x = std::rand() % 10 == 0 ? -x : x;
  }
  for (auto _ : state)
    for (double x : numbers)
      benchmark::DoNotOptimize((*func)(x));
}


BENCHMARK_CAPTURE(run, round, use_round);
BENCHMARK_CAPTURE(run, floor, use_floor);
BENCHMARK_CAPTURE(run, nearbyint, use_nearbyint);
BENCHMARK_CAPTURE(run, rint, use_rint);
BENCHMARK_CAPTURE(run2, lrint, use_lrint);
BENCHMARK_CAPTURE(run2, round, use2_round);
BENCHMARK_CAPTURE(run2, floor, use2_floor);
BENCHMARK_CAPTURE(run2, nearbyint, use2_nearbyint);
BENCHMARK_CAPTURE(run2, rint, use2_rint);
BENCHMARK_MAIN();

/*
-----------------------------------------------------
Benchmark              Time           CPU Iterations
-----------------------------------------------------
run/round             26 ns         26 ns   23927270
run/floor             18 ns         18 ns   37953966
run/nearbyint         18 ns         18 ns   38603220
run/rint              13 ns         13 ns   53661537
run2/lrint            18 ns         18 ns   38754555
run2/round            27 ns         27 ns   26094020
run2/floor            18 ns         18 ns   38176320
run2/nearbyint        16 ns         16 ns   43588313
run2/rint             13 ns         13 ns   53491145
*/
