// Copyright 2017 Global Phasing Ltd.

// Microbenchmark of gemmi::SMat33::r_u_r().
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include smat33.cpp $GB/src/libbenchmark.a -pthread

#include "gemmi/math.hpp"
#include <benchmark/benchmark.h>
#include <stdlib.h>               // for rand

static void r_u_r1(benchmark::State& state) {
  gemmi::Vec3 rr[10] = {
    {-1., -2., -3.},
    { 1.,  1., -5.},
    {-1.,  6.,  3.},
    { 5., -2.,  9.},
    { 4.,  0., -3.},
    {-1.,  7., -8.},
    { 3., -2.,  3.},
    { 1., -4.,  0.},
    { 0.,  0.,  3.},
    { 1.,  2., -3.}
  };
  gemmi::SMat33<double> smat{
    0.000211, 0.000119, 0.0004575, -4.71e-05, 0.0001238, 3.91e-05
  };
  for (int i = 0; i != 10; ++i)
    if (rand() % 1000 == 0)
      rr[i].x = 10.;
  for (auto _ : state) {
    double sum = 0.;
    for (int i = 0; i != 10; ++i)
      sum += smat.r_u_r(rr[i]);
    benchmark::DoNotOptimize(sum);
  }
}

static void r_u_r2(benchmark::State& state) {
  std::array<int,3> rr[10] = {
    {-1, -2, -3},
    { 1,  1, -5},
    {-1,  6,  3},
    { 5, -2,  9},
    { 4,  0, -3},
    {-1,  7, -8},
    { 3, -2,  3},
    { 1, -4,  0},
    { 0,  0,  3},
    { 1,  2, -3}
  };
  for (int i = 0; i != 10; ++i)
    if (rand() % 1000 == 0)
      rr[i][0] = 10;
  gemmi::SMat33<double> smat{
    0.000211, 0.000119, 0.0004575, -4.71e-05, 0.0001238, 3.91e-05
  };
  for (auto _ : state) {
    double sum = 0.;
    for (int i = 0; i != 10; ++i)
      sum += smat.r_u_r(rr[i]);
    benchmark::DoNotOptimize(sum);
  }
}


BENCHMARK(r_u_r1);
BENCHMARK(r_u_r2);
BENCHMARK_MAIN();

/* Output from my laptop:

--------------------------------------------------
Benchmark           Time           CPU Iterations
--------------------------------------------------
r_u_r1             30 ns         30 ns   22750626
r_u_r2             39 ns         39 ns   18032077

*/
