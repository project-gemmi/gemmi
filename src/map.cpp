// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "gemmi/util.hpp"  // for trim_str
#include "input.h"
#include <cstdlib>   // for strtod
#include <cstdio>    // for fprintf
#include <algorithm> // for nth_element
#define USE_UNICODE
#ifdef USE_UNICODE
#include <clocale>  // for setlocale
#include <cwchar>  // for wint_t
#endif
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
void print_histogram(const std::vector<T>& data, double min, double max) {
#ifdef USE_UNICODE
  std::setlocale(LC_ALL, "");
#endif
  int bins[81] = {0};
  double delta = max - min;
  for (T d : data) {
    int n = gemmi::iround((d - min) * (80 / delta));
    bins[n]++;
  }
  double max_h = *std::max_element(std::begin(bins), std::end(bins));
  const int rows = 24;
  for (int i = rows; i > 0; --i) {
    for (int j = 0; j < 80; ++j) {
      double h = bins[j] / max_h * rows;
#ifdef USE_UNICODE
      wint_t c = ' ';
      if (h > i) {
        c = 0x2588; // 0x2581 = one eighth block, ..., 0x2588 = full block
      } else if (h > i - 1) {
        c = 0x2581 + static_cast<int>((h - (i - 1)) * 7);
      }
      printf("%lc", c);
#else
      std::putchar(h > i + 0.5 ? '#' : ' ');
#endif
    }
    std::putchar('\n');
  }
}

template<typename T>
void print_info(const gemmi::Grid<T>& grid) {
  std::printf("Map mode: %d\n", grid.header_i32(4));
  std::printf("Number of columns, rows, sections: %5d %5d %5d %8s %d points\n",
              grid.nu, grid.nv, grid.nw, "->", grid.nu * grid.nv * grid.nw);
  int u0 = grid.header_i32(5);
  int v0 = grid.header_i32(6);
  int w0 = grid.header_i32(7);
  std::printf("                             from: %5d %5d %5d\n", u0, v0, w0);
  std::printf("                               to: %5d %5d %5d\n",
              u0 + grid.nu - 1, v0 + grid.nv - 1, w0 + grid.nw - 1);
  const char* xyz = "?XYZ";
  std::printf("Fast, medium, slow axes: %c %c %c\n",
              xyz[grid.header_i32(17)],
              xyz[grid.header_i32(18)],
              xyz[grid.header_i32(19)]);
  int mx = grid.header_i32(8);
  int my = grid.header_i32(9);
  int mz = grid.header_i32(10);
  std::printf("Grid sampling on x, y, z: %5d %5d %5d          %8s %d points\n",
              mx, my, mz, "->", mx * my * mz);
  const mol::UnitCell& cell = grid.unit_cell;
  std::printf("Space group number: %d\n", grid.header_i32(23));
  std::printf("Cell dimensions: %g %g %g  %g %g %g\n",
              cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
  int origin[3] = {
    grid.header_i32(50),
    grid.header_i32(51),
    grid.header_i32(52)
  };
  if (origin[0] != 0 || origin[1] != 0 || origin[2] != 0)
    std::printf("Non-zero origin: %d %d %d\n", origin[0], origin[1], origin[2]);

  std::printf("\nStatistics from HEADER and DATA\n");
  gemmi::GridStats st = grid.calculate_statistics();
  std::printf("Minimum: %12.5f  %12.5f\n", grid.stats.dmin, st.dmin);
  std::printf("Maximum: %12.5f  %12.5f\n", grid.stats.dmax, st.dmax);
  std::printf("Mean:    %12.5f  %12.5f\n", grid.stats.dmean, st.dmean);
  std::printf("RMS:     %12.5f  %12.5f\n", grid.stats.rms, st.rms);
  std::vector<T> data = grid.data;
  size_t mpos = data.size() / 2;
  std::nth_element(data.begin(), data.begin() + mpos, data.end());
  std::printf("Median:                %12.5f\n", data[mpos]);
  print_histogram(data, st.dmin, st.dmax);
  int nlabl = grid.header_i32(56);
  if (nlabl != 0)
    std::printf("\n");
  for (int i = 0; i < nlabl && i < 10; ++i) {
    std::string label = gemmi::trim_str(grid.header_str(57 + i * 20));
    std::printf("Label #%d\n%s\n", i, label.c_str());
  }
  int nsymbt = grid.header_i32(24);
  if (nsymbt != 0)
    std::printf("\n");
  for (int i = 0; i * 80 < nsymbt; i++) {
    std::string symop = grid.header_str(256 + i * 20 /*words not bytes*/, 80);
    std::printf("Sym op #%d: %s\n", i + 1, gemmi::trim_str(symop).c_str());
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
