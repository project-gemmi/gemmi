
#include <cstdio>
#include <string>
#include <vector>
#include <benchmark/benchmark.h>
#include <gemmi/atox.hpp>

static int std_stoi(const std::string& str) {
  return std::stoi(str);
}
static int to_int_true(const std::string& str) {
  return gemmi::string_to_int(str, true);
}

static int to_int_false(const std::string& str) {
  return gemmi::string_to_int(str, false);
}

void sequential(benchmark::State& state, int(*func)(const std::string&)) {
  std::vector<std::string> v(20000);
  for (size_t i = 0; i != v.size(); ++i)
    v[i] = std::to_string(i);
  for (size_t i = 0; i != v.size(); ++i)
    if ((*func)(v[i]) != (int)i)
      std::printf("ERROR: at %zu\n", i);
  for (auto _ : state)
    for (const std::string& s : v)
      benchmark::DoNotOptimize((*func)(s));
}


BENCHMARK_CAPTURE(sequential, std_stoi, std_stoi);
BENCHMARK_CAPTURE(sequential, to_int_false, to_int_false);
BENCHMARK_CAPTURE(sequential, to_int_true, to_int_true);
BENCHMARK_MAIN();
