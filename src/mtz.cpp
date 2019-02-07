// Copyright 2019 Global Phasing Ltd.
//
// MTZ info

#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#define GEMMI_PROG mtz
#include "options.h"
#include <stdio.h>

using gemmi::Mtz;

enum OptionIndex { Verbose=3, Headers, Dump, PrintTsv, PrintStats,
                   ToggleEndian };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] MTZ_FILE[...]"
    "\nPrint informations from an mtz file."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Headers, 0, "H", "headers", Arg::None,
    "  -H, --headers  \tPrint raw headers, until the END record." },
  { Dump, 0, "d", "dump", Arg::None,
    "  -d, --dump  \tPrint a subset of CCP4 mtzdmp informations." },
  { PrintTsv, 0, "", "tsv", Arg::None,
    "  --tsv  \tPrint all the data as tab-separated values." },
  { PrintStats, 0, "", "stats", Arg::None,
    "  --stats  \tPrint column statistics (completeness, mean, etc)." },
  { ToggleEndian, 0, "", "toggle-endian", Arg::None,
    "  --toggle-endian  \tToggle assumed endiannes (little <-> big)." },
  { 0, 0, 0, 0, 0, 0 }
};

static void dump(const Mtz& mtz) {
  printf("Title: %s\n", mtz.title.c_str());
  printf("Number of Datasets = %zu\n\n", mtz.datasets.size());
  for (const Mtz::Dataset& ds : mtz.datasets) {
    printf("Dataset %4d   %s > %s > %s:\n",
           ds.number, ds.project_name.c_str(),
           ds.crystal_name.c_str(), ds.dataset_name.c_str());
    printf("        cell  %g %7g %7g  %6g %6g %6g\n",
           ds.cell.a, ds.cell.b, ds.cell.c,
           ds.cell.alpha, ds.cell.beta, ds.cell.gamma);
    printf("  wavelength  %g\n", ds.wavelength);
  }
  printf("\nNumber of Columns = %d\n", mtz.ncol);
  printf("Number of Reflections = %d\n", mtz.nreflections);
  if (mtz.nbatches != 0)
    printf("Number of Batches = %d\n", mtz.nbatches);
  printf("Missing values marked as: %g\n", mtz.valm);
  printf("Global Cell (obsolete):  %7.3f %7.3f %7.3f  %g %g %g\n",
         mtz.cell.a, mtz.cell.b, mtz.cell.c,
         mtz.cell.alpha, mtz.cell.beta, mtz.cell.gamma);
  printf("Resolution: %.2f - %.2f A\n",
         mtz.resolution_high(), mtz.resolution_low());
  printf("Sort Order: %d %d %d %d %d\n",
         mtz.sort_order[0], mtz.sort_order[1], mtz.sort_order[2],
         mtz.sort_order[3], mtz.sort_order[4]);
  printf("Space Group: %s\n", mtz.spacegroup_name.c_str());
  printf("Space Group Number: %d\n", mtz.spacegroup_number);
  printf("\nColumn    Type Dataset    Min        Max\n");
  for (const Mtz::Column& col : mtz.columns)
    printf("%-12s %c %2d %12.6g %10.6g\n",
           col.label.c_str(), col.type, col.dataset_number,
           col.min_value, col.max_value);
  if (mtz.history.empty()) {
    printf("\nNo history in the file.\n");
  } else {
    printf("\nHistory (%zu lines):\n", mtz.history.size());
    for (const std::string& hline : mtz.history)
      printf("%s\n", hline.c_str());
  }
}

static void print_tsv(const Mtz& mtz) {
  int ncol = mtz.ncol;
  for (int i = 0; i < ncol; ++i)
    printf("%s%c", mtz.columns[i].label.c_str(), i + 1 != ncol ? '\t' : '\n');
  for (int i = 0; i < mtz.nreflections * ncol; ++i)
    printf("%g%c", mtz.data[i], (i + 1) % ncol != 0 ? '\t' : '\n');
}

struct ColumnStats {
  float min_value = INFINITY;
  float max_value = -INFINITY;
  gemmi::Variance var;
};

static void print_stats(const Mtz& mtz) {
  std::vector<ColumnStats> column_stats(mtz.columns.size());
  for (size_t i = 0; i != mtz.data.size(); ++i) {
    float v = mtz.data[i];
    if (!std::isnan(v)) {
      ColumnStats& stat = column_stats[i % column_stats.size()];
      if (v < stat.min_value)
        stat.min_value = v;
      if (v > stat.max_value)
        stat.max_value = v;
      stat.var.add_point(v);
    }
  }
  printf("column type @dataset  completeness        min       max"
         "       mean   stddev\n");
  for (size_t i = 0; i != column_stats.size(); ++i) {
    const Mtz::Column& col = mtz.columns[i];
    const ColumnStats& stat = column_stats[i];
    printf("%-14s %c @%d  %d (%6.2f%%) %9.5g %9.5g  %9.5g %8.4g\n",
           col.label.c_str(), col.type, col.dataset_number,
           stat.var.n, 100.0 * stat.var.n / mtz.nreflections,
           stat.min_value, stat.max_value,
           stat.var.mean_x, std::sqrt(stat.var.for_population()));
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  bool verbose = p.options[Verbose];


  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      const char* path = p.nonOption(i);
      if (i != 0)
        printf("\n\n");
      if (verbose)
        fprintf(stderr, "Reading %s ...\n", path);
      gemmi::fileptr_t f = gemmi::file_open(path, "rb");
      Mtz mtz;
      try {
        mtz.read_first_bytes(f.get());
        if (p.options[ToggleEndian])
          mtz.toggle_endiannes();
      } catch (std::runtime_error& e) {
        gemmi::fail(std::string(e.what()) + ": " + path);
      }
      if (p.options[Headers]) {
        char buf[81] = {0};
        mtz.seek_headers(f.get());
        while (fread(buf, 1, 80, f.get()) != 0) {
          printf("%s\n", gemmi::rtrim_str(buf).c_str());
          if (gemmi::ialpha3_id(buf) == gemmi::ialpha3_id("END"))
            break;
        }
      }
      if (verbose)
        mtz.warnings = stderr;
      mtz.read_main_headers(f.get());
      mtz.read_history_and_later_headers(f.get());
      if (p.options[Dump])
        dump(mtz);
      if (p.options[PrintTsv] || p.options[PrintStats])
        mtz.read_raw_data(f.get());
      if (p.options[PrintTsv])
        print_tsv(mtz);
      if (p.options[PrintStats])
        print_stats(mtz);
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
