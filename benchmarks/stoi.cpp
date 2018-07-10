
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>
#include <benchmark/benchmark.h>
#include <gemmi/pdb.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/stoi.hpp>

static int std_stoi(const std::string& str) {
  return std::stoi(str);
}
static int to_int_true(const std::string& str) {
  return gemmi::string_to_int(str, true);
}

static int to_int_false(const std::string& str) {
  return gemmi::string_to_int(str, false);
}

static int pdb_hpp(const std::string& str) {
  return gemmi::pdb_impl::read_int(str.c_str(), str.size());
}

static int cif_as_int(const std::string& str) {
  return gemmi::cif::as_int(str);
}

static int cif_as_int_noexcept(const std::string& str) {
  return gemmi::cif::as_int_noexcept(str, 0);
}

void sequential(benchmark::State& state, int(*func)(const std::string&)) {
  std::vector<std::string> v(20000);
  for (size_t i = 0; i != v.size(); ++i)
    v[i] = std::to_string(i);
  for (size_t i = 0; i != v.size(); ++i)
    if ((*func)(v[i]) != (int)i)
      std::printf("ERROR: at %zu\n", i);
  while (state.KeepRunning())
    for (const std::string& s : v)
      benchmark::DoNotOptimize((*func)(s));
}


BENCHMARK_CAPTURE(sequential, std_stoi, std_stoi);
BENCHMARK_CAPTURE(sequential, to_int_false, to_int_false);
BENCHMARK_CAPTURE(sequential, to_int_true, to_int_true);
BENCHMARK_CAPTURE(sequential, pdb_hpp, pdb_hpp);
BENCHMARK_CAPTURE(sequential, cif_as_int, cif_as_int);
BENCHMARK_CAPTURE(sequential, cif_as_int_noexcept, cif_as_int_noexcept);
BENCHMARK_MAIN();

// vim:sw=2:ts=2:et
