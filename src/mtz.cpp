// Copyright 2019 Global Phasing Ltd.
//
// MTZ info

#include <cstdio>
#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/input.hpp>    // for FileStream, MemoryStream
#define GEMMI_PROG mtz
#include "options.h"

using gemmi::Mtz;
using std::printf;

namespace {

enum OptionIndex { Headers=4, Dump, PrintBatch, PrintBatches, PrintTsv,
                   PrintStats, CheckAsu, ToggleEndian, NoIsym, UpdateReso };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] MTZ_FILE[...]"
    "\nPrint informations from an mtz file."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Headers, 0, "H", "headers", Arg::None,
    "  -H, --headers  \tPrint raw headers, until the END record." },
  { Dump, 0, "d", "dump", Arg::None,
    "  -d, --dump  \tPrint a subset of CCP4 mtzdmp informations." },
  { PrintBatch, 0, "B", "batch", Arg::Int,
    "  -B N, --batch=N  \tPrint data from batch header N." },
  { PrintBatches, 0, "b", "batches", Arg::None,
    "  -b, --batches  \tPrint data from all batch headers." },
  { PrintTsv, 0, "", "tsv", Arg::None,
    "  --tsv  \tPrint all the data as tab-separated values." },
  { PrintStats, 0, "s", "stats", Arg::None,
    "  -s, --stats  \tPrint column statistics (completeness, mean, etc)." },
  { CheckAsu, 0, "", "check-asu", Arg::None,
    "  --check-asu  \tCheck if reflections are in conventional ASU." },
  { ToggleEndian, 0, "", "toggle-endian", Arg::None,
    "  --toggle-endian  \tToggle assumed endiannes (little <-> big)." },
  { NoIsym, 0, "", "no-isym", Arg::None,
    "  --no-isym  \tDo not apply symmetry from M/ISYM column." },
  { UpdateReso, 0, "", "update-reso", Arg::None,
    "  --update-reso  \tRecalculate resolution limits before printing." },
  { 0, 0, 0, 0, 0, 0 }
};

void dump(const Mtz& mtz) {
  printf("Title: %s\n", mtz.title.c_str());
  printf("Total Number of Datasets = %zu\n\n", mtz.datasets.size());
  for (const Mtz::Dataset& ds : mtz.datasets) {
    printf("Dataset %4d   %s > %s > %s:\n",
           ds.id, ds.project_name.c_str(),
           ds.crystal_name.c_str(), ds.dataset_name.c_str());
    printf("        cell  %g %7g %7g  %6g %6g %6g\n",
           ds.cell.a, ds.cell.b, ds.cell.c,
           ds.cell.alpha, ds.cell.beta, ds.cell.gamma);
    printf("  wavelength  %g\n", ds.wavelength);
  }
  printf("\nNumber of Columns = %zu\n", mtz.columns.size());
  printf("Number of Reflections = %d\n", mtz.nreflections);
  printf("Number of Batches = %zu\n", mtz.batches.size());
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
  printf("Space Group Number: %d\n\n", mtz.spacegroup_number);
  printf("Header info (run with option -s for recalculated statistics):\n");
  printf("Column    Type  Dataset    Min        Max\n");
  for (const Mtz::Column& col : mtz.columns)
    printf("%-12s %c %2d %12.6g %10.6g\n",
           col.label.c_str(), col.type, col.dataset_id,
           col.min_value, col.max_value);
  if (mtz.history.empty()) {
    printf("\nNo history in the file.\n");
  } else {
    printf("\nHistory (%zu lines):\n", mtz.history.size());
    for (const std::string& hline : mtz.history)
      printf("%s\n", hline.c_str());
  }
  if (!mtz.batches.empty()) {
    printf("\nBatch counts:\n");
    int max_dataset_id = 0;
    for (const Mtz::Dataset& ds : mtz.datasets)
      max_dataset_id = std::max(max_dataset_id, ds.id);
    if (max_dataset_id < 1000000) {
      std::vector<int> batch_count(max_dataset_id + 1, 0);
      std::vector<int> batch_min(batch_count.size(), INT_MAX);
      std::vector<int> batch_max(batch_count.size(), -1);
      for (const Mtz::Batch& batch : mtz.batches) {
        int ds_id = batch.dataset_id();
        batch_count[ds_id]++;
        batch_min[ds_id] = std::min(batch_min[ds_id], batch.number);
        batch_max[ds_id] = std::max(batch_max[ds_id], batch.number);
      }
      for (int i = 0; i <= max_dataset_id; ++i)
        if (batch_count[i] != 0)
          printf(" dataset %d: %d batches (%d - %d)\n",
                 i, batch_count[i], batch_min[i], batch_max[i]);
    }
  }
}

void print_batch(const Mtz::Batch& b) {
  printf("Batch %d - %s\n", b.number, b.title.c_str());
  printf("    %zu %s: %s\n", b.axes.size(),
         b.axes.size() == 1 ? "axis" : "axes",
         gemmi::join_str(b.axes, ", ").c_str());
  printf("  %4zu integers:", b.ints.size());
  for (size_t i = 0; i != b.ints.size(); ++i) {
    if (i != 0 && i % 10 == 0)
      printf("\n                ");
    printf(" %5d", b.ints[i]);
  }
  printf("\n  %4zu floats:", b.floats.size());
  size_t last_non_zero = b.floats.size();
  while (last_non_zero != 0 && b.floats[last_non_zero - 1] == 0.f)
    --last_non_zero;
  for (size_t i = 0; i != b.floats.size(); ++i) {
    if (i != 0 && i % 5 == 0) {
      printf("\n       %4zu|  ", i);
      if (i >= last_non_zero) {
        printf("          ...");
        break;
      }
    }
    printf(" %12.5g", b.floats[i]);
  }
  printf("\n    dataset: %d\n", b.dataset_id());
}

void print_batch_extra_info(const Mtz::Batch& b) {
  gemmi::UnitCell uc = b.get_cell();
  printf("    Unit cell parameters: %g %g %g   %g %g %g\n",
         uc.a, uc.b, uc.c, uc.alpha, uc.beta, uc.gamma);
  printf("    Phi start - end: %g - %g\n", b.phi_start(), b.phi_end());
  gemmi::Mat33 u = b.matrix_U();
  for (int i = 0; i != 3; ++i)
    printf("    %s % 10.6f % 10.6f % 10.6f\n",
           i == 0 ? "Orientation matrix U:" : "                     ",
           u.a[i][0], u.a[i][1], u.a[i][2]);
}

void print_tsv(const Mtz& mtz) {
  size_t ncol = mtz.columns.size();
  for (size_t i = 0; i < ncol; ++i)
    printf("%s%c", mtz.columns[i].label.c_str(), i + 1 != ncol ? '\t' : '\n');
  for (size_t i = 0; i < mtz.nreflections * ncol; ++i)
    printf("%g%c", mtz.data[i], (i + 1) % ncol != 0 ? '\t' : '\n');
}

struct ColumnStats {
  float min_value = INFINITY;
  float max_value = -INFINITY;
  gemmi::Variance var;
};

void print_stats(const Mtz& mtz) {
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
           col.label.c_str(), col.type, col.dataset_id,
           stat.var.n, 100.0 * stat.var.n / mtz.nreflections,
           stat.min_value, stat.max_value,
           stat.var.mean_x, std::sqrt(stat.var.for_population()));
  }
}

void check_asu(const Mtz& mtz) {
  size_t ncol = mtz.columns.size();
  const gemmi::SpaceGroup* sg = mtz.spacegroup;
  if (!sg)
    gemmi::fail("no spacegroup in the MTZ file.");
  int counter = 0;
  gemmi::ReciprocalAsu asu(sg);
  for (int i = 0; i < mtz.nreflections; ++i) {
    int h = (int) mtz.data[i * ncol + 0];
    int k = (int) mtz.data[i * ncol + 1];
    int l = (int) mtz.data[i * ncol + 2];
    if (asu.is_in({{h, k, l}}))
      ++counter;
  }
  printf("spacegroup: %s\n", sg->xhm().c_str());
  printf("ccp4 ASU convention wrt. standard setting: %s\n",
         asu.condition_str());
  printf("inside / outside of ASU: %d / %d\n",
         counter, mtz.nreflections - counter);
  double dmin = mtz.resolution_high() - 1e-6;
  printf("All unique reflections up to d=%g: %d\n",
         dmin, gemmi::count_reflections(mtz.cell, mtz.spacegroup, dmin));
}

template<typename Stream>
void print_mtz_info(Stream&& stream, const char* path,
                    const std::vector<option::Option>& options) {
  Mtz mtz;
  try {
    mtz.read_first_bytes(stream);
    if (options[ToggleEndian])
      mtz.toggle_endiannes();
  } catch (std::runtime_error& e) {
    gemmi::fail(std::string(e.what()) + ": " + path);
  }
  if (options[Headers]) {
    char buf[81] = {0};
    mtz.seek_headers(stream);
    while (stream.read(buf, 80)) {
      printf("%s\n", gemmi::rtrim_str(buf).c_str());
      if (gemmi::ialpha3_id(buf) == gemmi::ialpha3_id("END"))
        break;
    }
  }
  if (options[Verbose])
    mtz.warnings = stderr;
  mtz.read_main_headers(stream);
  mtz.read_history_and_batch_headers(stream);
  mtz.setup_spacegroup();
  if (options[PrintTsv] || options[PrintStats] || options[CheckAsu] ||
      options[UpdateReso])
    mtz.read_raw_data(stream);
  if (options[UpdateReso])
    mtz.update_reso();
  if (options[Dump] ||
      !(options[PrintBatch] || options[PrintBatches] || options[PrintTsv] ||
        options[PrintStats] || options[CheckAsu] || options[Headers]))
    dump(mtz);
  if (options[PrintBatch]) {
    for (const option::Option* o = options[PrintBatch]; o; o = o->next()) {
      int number = std::atoi(o->arg);
      for (const Mtz::Batch& b : mtz.batches)
        if (b.number == number) {
          print_batch(b);
          print_batch_extra_info(b);
        }
    }
  }
  if (options[PrintBatches])
    for (const Mtz::Batch& b : mtz.batches)
      print_batch(b);
  if (mtz.has_data() && !options[NoIsym])
    mtz.switch_to_original_hkl();
  if (options[PrintTsv])
    print_tsv(mtz);
  if (options[PrintStats])
    print_stats(mtz);
  if (options[CheckAsu])
    check_asu(mtz);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      const char* path = p.nonOption(i);
      if (i != 0)
        printf("\n\n");
      if (p.options[Verbose])
        std::fprintf(stderr, "Reading %s ...\n", path);
      gemmi::MaybeGzipped input(path);
      if (input.is_stdin()) {
        print_mtz_info(gemmi::FileStream{stdin}, path, p.options);
      } else if (std::unique_ptr<char[]> mem = input.memory()) {
        gemmi::MemoryStream stream(mem.get(), mem.get() + input.memory_size());
        print_mtz_info(std::move(stream), path, p.options);
      } else {
        gemmi::fileptr_t f = gemmi::file_open(input.path().c_str(), "rb");
        print_mtz_info(gemmi::FileStream{f.get()}, path, p.options);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
