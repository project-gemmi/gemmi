// Copyright 2017 Global Phasing Ltd.

// Microbenchmark of gemmi::Residue::has_standard_pdb_name().
// Requires the google/benchmark library. It can be built manually:
// c++ -Wall -O2 -I../include -I$GB/include resinfo.cpp $GB/src/libbenchmark.a -pthread

#include "gemmi/resinfo.hpp"
#include "gemmi/seqtools.hpp"  // for calculate_sequence_weight
#include <benchmark/benchmark.h>
#include <stdlib.h>               // for rand

static const std::string residue_names[10] =
    { "GLY", "ASN", "ASN", "HOH", "LEU", "DC", "G", "GLU", "ALA", "SER" };

static void find_tabulated_residue_x10(benchmark::State& state) {
  std::string names[10];
  for (int i = 0; i != 10; ++i)
    names[i] = rand() % 1000 == 0 ? "GLN" : residue_names[i];
  for (auto _ : state) {
    int nh = 0;
    for (int i = 0; i != 10; ++i)
      nh += gemmi::find_tabulated_residue(names[i]).hydrogen_count;
    benchmark::DoNotOptimize(nh);
  }
}

static void sequence_weight(benchmark::State& state) {
  // example sequence from 1MRU
  std::vector<std::string> seq{
    "MET", "THR", "THR", "PRO", "SER", "HIS", "LEU", "SER", "ASP", "ARG",
    "TYR", "GLU", "LEU", "GLY", "GLU", "ILE", "LEU", "GLY", "PHE", "GLY",
    "GLY", "MET", "SER", "GLU", "VAL", "HIS", "LEU", "ALA", "ARG", "ASP",
    "LEU", "ARG", "LEU", "HIS", "ARG", "ASP", "VAL", "ALA", "VAL", "LYS",
    "VAL", "LEU", "ARG", "ALA", "ASP", "LEU", "ALA", "ARG", "ASP", "PRO",
    "SER", "PHE", "TYR", "LEU", "ARG", "PHE", "ARG", "ARG", "GLU", "ALA",
    "GLN", "ASN", "ALA", "ALA", "ALA", "LEU", "ASN", "HIS", "PRO", "ALA",
    "ILE", "VAL", "ALA", "VAL", "TYR", "ASP", "THR", "GLY", "GLU", "ALA",
    "GLU", "THR", "PRO", "ALA", "GLY", "PRO", "LEU", "PRO", "TYR", "ILE",
    "VAL", "MET", "GLU", "TYR", "VAL", "ASP", "GLY", "VAL", "THR", "LEU",
    "ARG", "ASP", "ILE", "VAL", "HIS", "THR", "GLU", "GLY", "PRO", "MET",
    "THR", "PRO", "LYS", "ARG", "ALA", "ILE", "GLU", "VAL", "ILE", "ALA",
    "ASP", "ALA", "CYS", "GLN", "ALA", "LEU", "ASN", "PHE", "SER", "HIS",
    "GLN", "ASN", "GLY", "ILE", "ILE", "HIS", "ARG", "ASP", "VAL", "LYS",
    "PRO", "ALA", "ASN", "ILE", "MET", "ILE", "SER", "ALA", "THR", "ASN",
    "ALA", "VAL", "LYS", "VAL", "MET", "ASP", "PHE", "GLY", "ILE", "ALA",
    "ARG", "ALA", "ILE", "THR", "ALA", "GLN", "TYR", "LEU", "SER", "PRO",
    "GLU", "GLN", "ALA", "ARG", "GLY", "ASP", "SER", "VAL", "ASP", "ALA",
    "ARG", "SER", "ASP", "VAL", "TYR", "SER", "LEU", "GLY", "CYS", "VAL",
    "LEU", "TYR", "GLU", "VAL", "LEU", "THR", "GLY", "GLU", "PRO", "PRO",
    "PHE", "THR", "GLY", "ASP", "SER", "PRO", "VAL", "SER", "VAL", "ALA",
    "TYR", "GLN", "HIS", "VAL", "ARG", "GLU", "ASP", "PRO", "ILE", "PRO",
    "PRO", "SER", "ALA", "ARG", "HIS", "GLU", "GLY", "LEU", "SER", "ALA",
    "ASP", "LEU", "ASP", "ALA", "VAL", "VAL", "LEU", "LYS", "ALA", "LEU",
    "ALA", "LYS", "ASN", "PRO", "GLU", "ASN", "ARG", "TYR", "GLN", "THR",
    "ALA", "ALA", "GLU", "MET", "ARG", "ALA", "ASP", "LEU", "VAL", "ARG",
    "VAL", "HIS", "ASN", "GLY", "GLU", "PRO", "PRO", "GLU", "ALA", "PRO",
    "LYS"};

  for (auto _ : state) {
    double weight = gemmi::calculate_sequence_weight(seq);
    benchmark::DoNotOptimize(weight);
  }
}

BENCHMARK(find_tabulated_residue_x10);
BENCHMARK(sequence_weight);
BENCHMARK_MAIN();

/* Output from my laptop:

------------------------------------------------------------------
Benchmark                           Time           CPU Iterations
------------------------------------------------------------------
find_tabulated_residue_x10         52 ns         52 ns   13458745

*/
