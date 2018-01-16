// Copyright 2017 Global Phasing Ltd.

#include <cmath>
#include <benchmark/benchmark.h>

static double use_round(double x) { return std::round(x); }
static double use_floor(double x) { return std::floor(x + 0.5); }
static double use_nearbyint(double x) { return std::nearbyint(x); }


inline void run(benchmark::State& state, double(*func)(double)) {
  double numbers[] = { 1.331, 6.3, 4, -4.5, 99341.4106,
                       0, -6.05823, 1e-9, 18.419, -12345.67 };
  for (double& x : numbers) { // otherwise it would be all optimized out
    x = std::rand() % 10 == 0 ? -x : x;
  }
  while (state.KeepRunning())
    for (int x : numbers)
      benchmark::DoNotOptimize((*func)(x));
}

BENCHMARK_CAPTURE(run, 1, use_round);
BENCHMARK_CAPTURE(run, 2, use_floor);
BENCHMARK_CAPTURE(run, 3, use_nearbyint);
BENCHMARK_MAIN()

/*
--------------------------------------------------
Benchmark           Time           CPU Iterations
--------------------------------------------------
run/1              32 ns         32 ns   21878391
run/2              17 ns         17 ns   41218368
run/3              14 ns         14 ns   49802527
run/4              21 ns         21 ns   33985958
run/5              11 ns         11 ns   64842324
*/

// vim:sw=2:ts=2:et:path^=../include,../third_party
