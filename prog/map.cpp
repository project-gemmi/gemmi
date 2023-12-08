// Copyright 2017 Global Phasing Ltd.

#include <cmath>           // for floor
#include <cstdio>          // for fprintf
#include <algorithm>       // for nth_element, count_if
#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"    // for MaybeGzipped
#include "gemmi/util.hpp"  // for trim_str
#include "gemmi/symmetry.hpp"
#include "gemmi/floodfill.hpp"  // for mask_points_above_threshold
#include "histogram.h"     // for print_histogram
#define GEMMI_PROG map
#include "options.h"

namespace {

enum OptionIndex {
  Dump=4, Deltas, CheckSym, Reorder, Full, Mask, Threshold, Fraction
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] CCP4_MAP[...]\n" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Dump, 0, "d", "dump", Arg::None,
    "  -d, --dump  \tPrint a map summary (default action)." },
  { Deltas, 0, "", "deltas", Arg::None,
    "  --deltas  \tStatistics of dx, dy and dz." },
  { CheckSym, 0, "", "check-symmetry", Arg::None,
    "  --check-symmetry  \tCompare the values of symmetric points." },
  { Reorder, 0, "", "write-xyz", Arg::Required,
    "  --write-xyz=FILE  \tWrite transposed map with fast X axis and slow Z." },
  { Full, 0, "", "write-full", Arg::Required,
    "  --write-full=FILE  \tWrite map extended to cover whole unit cell." },
  { Mask, 0, "", "write-mask", Arg::Required,
    "  --write-mask=FILE  \tMake a mask by thresholding the map." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for making a mask:" },
  { Threshold, 0, "", "threshold", Arg::Float,
    "  --threshold  \tExplicit threshold value for 0/1 mask." },
  { Fraction, 0, "", "fraction", Arg::Float,
    "  --fraction  \tThreshold is selected to have this fraction of 1's." },
  { 0, 0, 0, 0, 0, 0 }
};

template<typename T>
void print_info(const gemmi::Ccp4<T>& map, const gemmi::DataStats& st) {
  const gemmi::Grid<T>& grid = map.grid;
  std::printf("Map mode: %d\n", map.header_i32(4));
  std::printf("Endianness: %snative\n", map.same_byte_order ? "" : "NOT ");
  std::printf("Number of columns, rows, sections: %5d %5d %5d %6s %d points\n",
              grid.nu, grid.nv, grid.nw, "->", grid.nu * grid.nv * grid.nw);
  int u0 = map.header_i32(5);
  int v0 = map.header_i32(6);
  int w0 = map.header_i32(7);
  std::printf("                             from: %5d %5d %5d\n", u0, v0, w0);
  std::printf("                               to: %5d %5d %5d\n",
              u0 + grid.nu - 1, v0 + grid.nv - 1, w0 + grid.nw - 1);
  std::printf("Fast, medium, slow axes: %c %c %c\n",
              'X' + map.header_i32(17) - 1,
              'X' + map.header_i32(18) - 1,
              'X' + map.header_i32(19) - 1);
  int mx = map.header_i32(8);
  int my = map.header_i32(9);
  int mz = map.header_i32(10);
  std::printf("Grid sampling on x, y, z: %5d %5d %5d    %12s %d points/cell\n",
              mx, my, mz, "->", mx * my * mz);
  const gemmi::UnitCell& cell = grid.unit_cell;
  const gemmi::SpaceGroup* sg = grid.spacegroup;
  std::printf("Space group: %d  (%s)\n",
              sg ? sg->ccp4 : 0, sg ? sg->hm : "unknown");
  int order = sg ? sg->operations().order() : 1;
  std::printf("SG order: %-3d      %40s %d points/ASU\n",
               order, "->", mx * my * mz / order);
  std::printf("Cell dimensions: %g %g %g  %g %g %g\n",
              cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
  gemmi::Position origin = map.get_origin();
  if (origin.x != 0 || origin.y != 0 || origin.z != 0)
    std::printf("Non-zero origin: %g %g %g\n", origin.x, origin.y, origin.z);
  if (map.has_skew_transformation())
    std::printf("Defines skew transformation with translation length %.5g A.\n",
                map.get_skew_transformation().vec.length());

  if (st.nan_count != 0)
    std::printf("\n*** Data includes NaNs: %zu of %zu points ***",
                st.nan_count, grid.data.size());
  std::printf("\nStatistics from HEADER and DATA\n");
  std::printf("Minimum: %12.5f  %12.5f\n", map.hstats.dmin, st.dmin);
  std::printf("Maximum: %12.5f  %12.5f\n", map.hstats.dmax, st.dmax);
  std::printf("Mean:    %12.5f  %12.5f\n", map.hstats.dmean, st.dmean);
  std::printf("RMS:     %12.5f  %12.5f\n", map.hstats.rms, st.rms);
  std::vector<T> data = grid.data;  // copy b/c nth_element() reorders data
  if (st.nan_count != 0)
    gemmi::vector_remove_if(data, [](T x) { return std::isnan(x); });
  if (!data.empty()) {
    size_t mpos = data.size() / 2;
    std::nth_element(data.begin(), data.begin() + mpos, data.end());
    std::printf("Median:                %12.5f\n", data[mpos]);
    bool mask = std::all_of(data.begin(), data.end(),
                            [&st](T x) { return x == st.dmin || x == st.dmax; });
    double margin = mask ? 7 * (st.dmax - st.dmin) : 0;
    print_histogram(data, st.dmin - margin, st.dmax + margin);
  }

  int nlabl = map.header_i32(56);
  if (nlabl != 0)
    std::printf("\n");
  for (int i = 0; i < nlabl && i < 10; ++i) {
    std::string label = gemmi::trim_str(map.header_str(57 + i * 20));
    std::printf("Label #%d\n%s\n", i, label.c_str());
  }
  int nsymbt = map.header_i32(24);
  if (nsymbt != 0) {
    std::printf("\n");
    std::vector<gemmi::Op> ops;
    int bad_counter = 0;
    for (int i = 0; i * 80 < nsymbt; i++) {
      std::string symop = map.header_str(256 + i * 20 /*words not bytes*/, 80);
      try {
        gemmi::Op op = gemmi::parse_triplet(symop);
        ops.push_back(op);
      } catch (std::exception&) {
        ++bad_counter;
      }
      std::printf("Sym op #%d: %s\n", i + 1, gemmi::trim_str(symop).c_str());
    }
    if (bad_counter == 0) {
      gemmi::GroupOps gops = gemmi::split_centering_vectors(ops);
      const gemmi::SpaceGroup* sg2 = gemmi::find_spacegroup_by_ops(gops);
      std::printf("Space group from the operators: ");
      if (sg2)
        std::printf("%d  (%s)\n", sg2->ccp4, sg2->xhm().c_str());
      else
        std::printf("unknown\n");
      if (sg2 && sg && sg != sg2)
        std::printf("NOTE: different than from the ISPG header.\n");
    } else {
      std::printf("NOTE: %d of symmetry lines do not parse as operators.\n", bad_counter);
    }
  }
}

template<typename T>
void print_deltas(const gemmi::Grid<T>& grid, double dmin, double dmax) {
  std::vector<double> deltas;
  deltas.reserve(grid.data.size());
  for (int i = 0; i < 3; ++i) {
    int f[3] = {0, 0, 0};
    f[i] = 1;
    for (int w = f[2]; w < grid.nw; ++w)
      for (int v = f[1]; v < grid.nv; ++v)
        for (int u = f[0]; u < grid.nu; ++u)
          deltas.push_back(grid.get_value_q(u, v, w) -
                           grid.get_value_q(u - f[0], v - f[1], w - f[2]));
    gemmi::DataStats st = gemmi::calculate_data_statistics(deltas);
    std::printf("\nd%c: min: %.5f  max: %.5f  mean: %.5f  std.dev: %.5f\n",
                "XYZ"[i], st.dmin, st.dmax, st.dmean, st.rms);
    print_histogram(deltas, dmin, dmax);
    deltas.clear();
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  p.check_exclusive_pair(Threshold, Fraction);
  //bool verbose = p.options[Verbose];

  if (p.nonOptionsCount() > 1 && (p.options[Reorder] || p.options[Full])) {
    std::fprintf(stderr, "Option --write-... can be only used "
                         "with a single input file.\n");
    return 1;
  }

  bool dump = (p.options[Dump] ||
               !(p.options[Deltas] || p.options[CheckSym] ||
                 p.options[Reorder] || p.options[Full] || p.options[Mask]));
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      const char* input = p.nonOption(i);
      gemmi::Ccp4<> map;
      if (i != 0)
        std::printf("\n\n");
      std::printf("Reading file: %s\n", input);
      map.read_ccp4(gemmi::MaybeGzipped(input));
      gemmi::DataStats stats = gemmi::calculate_data_statistics(map.grid.data);
      if (dump)
        print_info(map, stats);
      if (p.options[Deltas])
        print_deltas(map.grid, stats.dmin, stats.dmax);
      if (p.options[Reorder]) {
        map.setup(NAN, gemmi::MapSetup::ReorderOnly);
        map.write_ccp4_map(p.options[Reorder].arg);
      }
      double max_err = 0.;
      if (p.options[CheckSym]) {
        const double eps = 0.01;
        // this is equivalent to calling map.setup(Full), but additionally
        // it prints inconsistent points.
        map.setup(NAN, gemmi::MapSetup::NoSymmetry);
        map.grid.symmetrize([&](float a, float b) {
            if (a < b || a > b) {
              double diff = std::fabs(a - b);
              if (diff > eps)
                std::printf("Symmetry-equivalent values differ: "
                            "%g != %g  diff: %g\n", a, b, diff);
              max_err = std::max(max_err, diff);
            }
            return std::isnan(a) ? b : a;
        });
        map.grid.calculate_spacing();
      } else if (p.options[Full] || p.options[Mask]) {
        map.setup(NAN, gemmi::MapSetup::Full);
      }
      if (p.options[Full]) {
        size_t nn = std::count_if(map.grid.data.begin(), map.grid.data.end(),
                                  [](float x) { return std::isnan(x); });
        if (nn != 0)
          std::fprintf(stderr, "WARNING: %zu unknown values set to NAN\n", nn);
        map.write_ccp4_map(p.options[Full].arg);
      }
      if (p.options[Mask]) {
        double threshold;
        if (p.options[Threshold]) {
          threshold = std::atof(p.options[Threshold].arg);
        } else if (p.options[Fraction]) {
          double fraction = std::atof(p.options[Fraction].arg);
          if (fraction < 0. || fraction > 1.) {
            std::fprintf(stderr, "Fraction must be between 0 and 1.\n");
            return 2;
          }
          auto data = map.grid.data;  // making a copy for nth_element()
          size_t n = std::min(static_cast<size_t>(data.size() * (1.0 - fraction)),
                              data.size() - 1);
          std::nth_element(data.begin(), data.begin() + n, data.end());
          threshold = data[n];
        } else {
          std::fprintf(stderr, "ERROR: you need to use --threshold or --fraction.\n");
          return 2;
        }
        gemmi::Ccp4<std::int8_t> mask;
        gemmi::mask_nodes_above_threshold(mask.grid, map.grid, threshold);
        size_t ones = std::count(mask.grid.data.begin(), mask.grid.data.end(), 1);
        size_t all = mask.grid.data.size();
        std::fprintf(stderr, "Masked %zu of %zu points (%.1f%%) above %g\n",
                     ones, all, 100.0 * ones / all, threshold);
        mask.update_ccp4_header(0);
        mask.write_ccp4_map(p.options[Mask].arg);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
