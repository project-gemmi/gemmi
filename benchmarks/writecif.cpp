// Copyright 2018 Global Phasing Ltd.

// Microbenchmark of CIF writing function.

#include "gemmi/to_cif.hpp"
#include "gemmi/read_cif.hpp"
#include <sstream>
#include <benchmark/benchmark.h>

static void write_cif(benchmark::State& state, const gemmi::cif::Document& doc,
                      gemmi::cif::Style options) {
  for (auto _ : state) {
    std::ostringstream os;
    gemmi::cif::write_cif_to_stream(os, doc, options);
    benchmark::DoNotOptimize(os);
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Call it with path to a cif file as an argument.\n");
    return 1;
  }
  const char* path = argv[argc-1];
  gemmi::cif::Document doc = gemmi::read_cif_gz(path);
  benchmark::RegisterBenchmark("write_cif", write_cif, doc, gemmi::cif::Style::Simple);
  benchmark::RegisterBenchmark("write_cif2", write_cif, doc, gemmi::cif::Style::Pdbx);
  benchmark::RegisterBenchmark("write_cif3", write_cif, doc, gemmi::cif::Style::Aligned);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}
