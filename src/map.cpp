// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "gemmi/util.hpp"  // for trim_str
#include "input.h"
#include <cstdlib>  // for strtod
#include <cstdio>  // for fprintf
#define EXE_NAME "gemmi-map"
#include "options.h"

namespace mol = gemmi::mol;

enum OptionIndex { Unknown, Verbose, OutputMode };

static const option::Descriptor Usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] CCP4_MAP[...]\n" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};


template<typename T>
void print_info(const gemmi::Grid<T>& grid) {
  std::printf("Map mode: %d\n", grid.header_u32(4));
  std::printf("Number of columns, rows, sections: %5d %5d %5d\n",
              grid.nu, grid.nv, grid.nw);
  std::printf("Number of data points: %d\n", grid.nu * grid.nv * grid.nw);
  int u0 = grid.header_u32(5);
  int v0 = grid.header_u32(6);
  int w0 = grid.header_u32(7);
  std::printf("                             from: %5d %5d %5d\n", u0, v0, w0);
  std::printf("                               to: %5d %5d %5d\n",
              u0 + grid.nu - 1, v0 + grid.nv - 1, w0 + grid.nw - 1);
  const char* xyz = "?XYZ";
  std::printf("Fast, medium, slow axes: %c %c %c\n",
              xyz[grid.header_u32(17)],
              xyz[grid.header_u32(18)],
              xyz[grid.header_u32(19)]);
  unsigned mx = grid.header_u32(8);
  unsigned my = grid.header_u32(9);
  unsigned mz = grid.header_u32(10);
  std::printf("Grid sampling on x, y, z: %5d %5d %5d\n", mx, my, mz);
  std::printf("Number of grid points: %d\n", mx * my * mz);
  const mol::UnitCell& cell = grid.unit_cell;
  std::printf("Space group number: %d\n", grid.header_u32(23));
  std::printf("Cell dimensions: %g %g %g  %g %g %g\n",
              cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);

  std::printf("\nSTATS from HEADER and DATA\n");
  gemmi::GridStats st = grid.calculate_statistics();
  std::printf("Minimum: %8.5f   %8.5f\n", grid.stats.dmin, st.dmin);
  std::printf("Maximum: %8.5f   %8.5f\n", grid.stats.dmax, st.dmax);
  std::printf("Mean:    %8.5f   %8.5f\n", grid.stats.dmean, st.dmean);
  std::printf("RMS:     %8.5f   %8.5f\n", grid.stats.rms, st.rms);
  uint32_t nlabl = grid.header_u32(56);
  for (uint32_t i = 0; i < nlabl && i < 10; ++i) {
    std::string label = gemmi::trim_str(grid.header_str(57 + i * 20));
    std::printf("Label #%u\n%s\n", i, label.c_str());
  }
}

int main(int argc, char **argv) {
  OptParser parse;
  auto options = parse.simple_parse(argc, argv, Usage);
  bool verbose = options[Verbose];

  if (parse.nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    return 0;
  }

  try {
    for (int i = 0; i < parse.nonOptionsCount(); ++i) {
      const char* input = parse.nonOption(i);
      gemmi::Grid<> grid;
      if (verbose)
        std::fprintf(stderr, "Reading %s ...\n", input);
      grid.read_ccp4(input);
      print_info(grid);

    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
