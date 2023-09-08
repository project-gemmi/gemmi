// Copyright 2017 Global Phasing Ltd.

// I couldn't decide how to implement positive modulo.
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include mod.cpp $GB/src/libbenchmark.a -pthread

#include <climits>
#include <cstdio>
#include <cstdlib>  // for rand
#include <benchmark/benchmark.h>

static const int n = 12;

static int positive_modulo1(int i) {
  return (i % n + n) % n;
}

static int positive_modulo2(int i) {
  return i >= 0 ? i % n : ((i + 1) % n + n - 1);
}

static int positive_modulo3(int i) {
  constexpr int BIG = INT_MAX / 2 / n;
  return (i + n*BIG) % n;
}

static int positive_modulo4(int i) {
  int j = i % n;
  return j < 0 ? j + n : j;
}

static int positive_modulo5(int i) {
  if (i >= n)
    return i % n;
  if (i < 0)
    return ((i+1) % n) + (n - 1);
  return i;
}

inline void run(benchmark::State& state, int(*func)(int)) {
  int divident[] = { 1, 0, 6, 4, -4, 6, 0, 18, 0, 0, -6, 0, 24, 3, 8, -2 };
  for (int& i : divident) { // otherwise it would be all optimized out
    i = std::rand() % 10 == 0 ? -i : i;
    //i += std::rand() % 1000 - 500;
  }
  for (int i = -30; i < 30; ++i)
    if ((*func)(i) != positive_modulo1(i))
      std::printf("ERROR: %d != %d modulo %d\n", (*func)(i), i, n);
  for (auto _ : state)
    for (int i : divident)
      benchmark::DoNotOptimize((*func)(i));
}

BENCHMARK_CAPTURE(run, 1, positive_modulo1);
BENCHMARK_CAPTURE(run, 2, positive_modulo2);
BENCHMARK_CAPTURE(run, 3, positive_modulo3);
BENCHMARK_CAPTURE(run, 4, positive_modulo4);
BENCHMARK_CAPTURE(run, 5, positive_modulo5);
BENCHMARK_MAIN();

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
