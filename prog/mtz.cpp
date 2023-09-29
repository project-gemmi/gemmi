// Copyright 2019 Global Phasing Ltd.
//
// MTZ info

#include <cstdio>
#include <iostream>  // for cerr
#include <unordered_map>
#include <gemmi/mtz.hpp>
#include <gemmi/asudata.hpp>  // for AsuData
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/input.hpp>    // for FileStream, MemoryStream
#include <gemmi/reciproc.hpp> // for count_reflections
#include "histogram.h"        // for print_histogram
#define GEMMI_PROG mtz
#include "options.h"

using gemmi::Mtz;
using std::printf;

namespace {

struct MtzArg: public Arg {
  static option::ArgStatus AsuChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"ccp4", "tnt"});
  }
};

enum OptionIndex {
  Headers=4, Dump, PrintBatch, PrintBatches, ExpandedBatches, PrintAppendix,
  PrintTsv, PrintStats, PrintHistogram, PrintCells, CheckAsu,
  Compare, ToggleEndian, NoIsym, UpdateReso
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] MTZ_FILE[...]"
    "\nPrint information from an mtz file."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Headers, 0, "H", "headers", Arg::None,
    "  -H, --headers  \tPrint raw headers, until the END record." },
  { Dump, 0, "d", "dump", Arg::None,
    "  -d, --dump  \tPrint a subset of CCP4 mtzdmp information." },
  { PrintBatch, 0, "B", "batch", Arg::Int,
    "  -B N, --batch=N  \tPrint data from batch header N." },
  { PrintBatches, 0, "b", "batches", Arg::None,
    "  -b, --batches  \tPrint data from all batch headers." },
  { ExpandedBatches, 0, "e", "", Arg::None,
    "  -e  \t(with -B or -b) expanded info from batch headers." },
  { PrintAppendix, 0, "A", "appendix", Arg::None,
    "  -A, --appendix  \tPrint appended text." },
  { PrintTsv, 0, "", "tsv", Arg::None,
    "  --tsv  \tPrint all the data as tab-separated values." },
  { PrintStats, 0, "s", "stats", Arg::None,
    "  -s, --stats  \tPrint column statistics (completeness, mean, etc)." },
  { PrintHistogram, 0, "", "histogram", Arg::Required,
    "  --histogram=LABEL  \tPrint histogram of values in column LABEL." },
  { PrintCells, 0, "", "cells", Arg::None,
    "  --cells  \tPrint cell parameters only." },
  { CheckAsu, 0, "", "check-asu", MtzArg::AsuChoice,
    "  --check-asu=ccp4|tnt  \tCheck if reflections are in ASU." },
  { Compare, 0, "", "compare", Arg::Required,
    "  --compare=FILE  \tCompare two MTZ files." },
  { ToggleEndian, 0, "", "toggle-endian", Arg::None,
    "  --toggle-endian  \tToggle assumed endianness (little <-> big)." },
  { NoIsym, 0, "", "no-isym", Arg::None,
    "  --no-isym  \tDo not apply symmetry from M/ISYM column." },
  { UpdateReso, 0, "", "update-reso", Arg::None,
    "  --update-reso  \tRecalculate resolution limits before printing." },
  { 0, 0, 0, 0, 0, 0 }
};

void print_cell_parameters(const char* prefix, const gemmi::UnitCell& cell) {
  printf("%s %g %7g %7g  %6g %6g %6g\n", prefix,
         cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
}

void dump(const Mtz& mtz) {
  printf("Title: %s\n", mtz.title.c_str());
  printf("Total Number of Datasets = %zu\n\n", mtz.datasets.size());
  for (const Mtz::Dataset& ds : mtz.datasets) {
    printf("Dataset %4d   %s > %s > %s:\n",
           ds.id, ds.project_name.c_str(),
           ds.crystal_name.c_str(), ds.dataset_name.c_str());
    print_cell_parameters("        cell ", ds.cell);
    printf("  wavelength  %g\n", ds.wavelength);
  }
  printf("\nNumber of Columns = %zu\n", mtz.columns.size());
  printf("Number of Reflections = %d\n", mtz.nreflections);
  printf("Number of Batches = %zu\n", mtz.batches.size());
  printf("Missing values marked as: %g\n", mtz.valm);
  print_cell_parameters("Global Cell (obsolete): ", mtz.cell);
  printf("Resolution: %.2f - %.2f A\n",
         mtz.resolution_high(), mtz.resolution_low());
  printf("Sort Order: %d %d %d %d %d\n",
         mtz.sort_order[0], mtz.sort_order[1], mtz.sort_order[2],
         mtz.sort_order[3], mtz.sort_order[4]);
  printf("Space Group: %s\n", mtz.spacegroup_name.c_str());
  printf("Space Group Number: %d\n", mtz.spacegroup_number);
  if (mtz.symops.empty()) {
    printf("No SYMM records.\n");
  } else {
    gemmi::GroupOps gops = gemmi::split_centering_vectors(mtz.symops);
    const gemmi::SpaceGroup* symm_sg = find_spacegroup_by_ops(gops);
    if (symm_sg == nullptr) {
      printf("Space Group from SYMM Records: unknown, the operations are:\n");
      for (const gemmi::Op& op : mtz.symops)
        printf("    %s\n", op.triplet().c_str());
    } else {
      printf("Space Group from SYMM Records: %s\n", symm_sg->xhm().c_str());
      if (symm_sg != mtz.spacegroup)
        printf("  WARNING: the space group differs!\n");
    }
  }
  printf("\nHeader info (run with option -s for recalculated statistics):\n");
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
    int prev_ds_id = -INT_MAX;
    int bspan[2] = {-INT_MAX, -INT_MAX};
    printf("\nBatch numbers:");
    for (size_t i = 0; i < mtz.batches.size(); ++i) {
      const Mtz::Batch& batch = mtz.batches[i];
      int ds_id = batch.dataset_id();
      if (ds_id != prev_ds_id || batch.number != bspan[1] + 1) {
        if (i != 0)
          printf(" %d-%d", bspan[0], bspan[1]);
        bspan[0] = batch.number;
        if (ds_id != prev_ds_id) {
          printf("\n dataset %d:", ds_id);
          prev_ds_id = ds_id;
        }
      }
      bspan[1] = batch.number;
    }
    printf(" %d-%d\n", bspan[0], bspan[1]);
  }
  if (!mtz.appended_text.empty())
    printf("\nAppendix: %zu bytes.\n", mtz.appended_text.size());
}

const char* batch_int_desc[] = {
  "no. of words",
  "no. of integers",
  "no. of reals",
  "type of orientation block",  // 3
  "refinement flag for cell a",
  "refinement flag for cell b",
  "refinement flag for cell c",
  "refinement flag for cell alpha",
  "refinement flag for cell beta",
  "refinement flag for cell gamma",
  "no. of missetting angle sets (PhiXYZ)", // 10
  "reciprocal axis closest to rot. axis E1",
  "crystal number",
  "crystal mosaicity (0=iso, 1=anisotropic)",
  "type of data (1=2D, 2=3D, 3=Laue)",
  "goniostat scan axis number",  // 15
  "no. of batch scales & Bfactors + SD's",
  "no. of goniostat axes",
  "beam info (0=lab, 1=synchrotron)",
  "no. of detectors",
  "dataset id",  // 20
};

const char* batch_float_desc[] = {
  "unit cell a",  // 0
  "unit cell b",
  "unit cell c",
  "unit cell alpha",
  "unit cell beta",
  "unit cell gamma",
  "U(1,1)", "U(2,1)", "U(3,1)",  // 6-8
  "U(1,2)", "U(2,2)", "U(3,2)",
  "U(1,3)", "U(2,3)", "U(3,3)",
  "misseting angle PhiXYZ(1,1)",  // 15
  "misseting angle PhiXYZ(2,1)",
  "misseting angle PhiXYZ(3,1)",
  "misseting angle PhiXYZ(1,2)",
  "misseting angle PhiXYZ(2,2)",
  "misseting angle PhiXYZ(3,2)",
  "mosaicity(1) reflection width (deg)",  // 21
  "mosaicity(2) vertical width (deg)",
  nullptr, nullptr, nullptr, nullptr, nullptr, // mosaicity padding
  nullptr, nullptr, nullptr, nullptr, nullptr,
  "datum(1) value of goniostat axis",  // 33
  "datum(2) value of goniostat axis",
  "datum(3) value of goniostat axis",
  "initial phi relative to datum",  // 36
  "final phi relative to datum",
  "scanax(1) rotation axis in lab frame",  // 38
  "scanax(2) rotation axis in lab frame",
  "scanax(3) rotation axis in lab frame",
  "start time [minutes]",  // 41
  "stop time [minutes]",
  "batch scale",  // 43
  "batch temperature factor",
  "sd of batch scale",
  "sd of temperature factor",
  "range of phi values",  // 47
  nullptr, nullptr, nullptr, nullptr, nullptr,
  nullptr, nullptr, nullptr, nullptr, nullptr,
  nullptr,
  "E1(1) \"Cambridge\" lab axes...",  // 59
  "E1(2) ...defining goniostat axes",
  "E1(3)",
  "E2(1)",
  "E2(2)",
  "E2(3)",
  "E3(1)",
  "E3(2)",
  "E3(3)",
  nullptr, nullptr, nullptr, nullptr, nullptr,  // 68-79
  nullptr, nullptr, nullptr, nullptr, nullptr,
  nullptr, nullptr,
  "source(1) idealised (ie. excluding ...", // 80
  "source(2) ...tilts) source vector, ...",
  "source(3) ...in Cambridge lab frame",
  "S0(1) source vector (incl. tilts), ...", // 83
  "S0(2) ...antiparallel to beam, ...",
  "S0(3) ...in Cambridge lab frame",
  "wavelength [A]",  // 86
  "dispersion delta(lambda)/lambda",  // 87
  "correlated component of dispersion",
  "horizontal beam divergence (deg)",  // 89
  "vertical beam divergence (0=isotropic)",
};

const char* batch_det_desc[] = {
  "DX crystal to detector distance [mm]",  // 111+40n
  "THETA detector tilt angle [deg]",
  "minimum Y coordinate [pixel]",  // 113+40n
  "maximum Y coordinate [pixel]",
  "minimum Z coordinate [pixel]",
  "maximum Z coordinate [pixel]",
};

void print_batch(const Mtz::Batch& b, bool expanded) {
  printf("Batch %d - %s\n", b.number, b.title.c_str());
  printf("    %zu %s: %s\n", b.axes.size(),
         b.axes.size() == 1 ? "axis" : "axes",
         gemmi::join_str(b.axes, ", ").c_str());
  printf("  %4zu integers:", b.ints.size());
  for (size_t i = 0; i != b.ints.size(); ++i) {
    if (expanded) {
      constexpr size_t n = sizeof(batch_int_desc) / sizeof(batch_int_desc[0]);
      const char* desc = i < n ? batch_int_desc[i] : nullptr;
      if (desc || b.ints[i] != 0)
        printf("\n     %4zu %-40s %5d", i, desc ? desc : "", b.ints[i]);
    } else {
      if (i != 0 && i % 10 == 0)
        printf("\n                ");
      printf(" %5d", b.ints[i]);
    }
  }
  printf("\n  %4zu floats:", b.floats.size());
  size_t last_non_zero = b.floats.size();
  while (last_non_zero != 0 && b.floats[last_non_zero - 1] == 0.f)
    --last_non_zero;
  if (expanded) {
    size_t expected_end = 111;
    for (size_t i = 0; i < std::min(b.floats.size(), expected_end); ++i) {
      constexpr size_t n = sizeof(batch_float_desc) / sizeof(batch_float_desc[0]);
      const char* desc = i < n ? batch_float_desc[i] : nullptr;
      if (desc || b.floats[i] != 0)
        printf("\n     %4zu %-40s %-.5g", i, desc ? desc : "", b.floats[i]);
    }
    int ndet = b.ints[19];
    if (ndet < 0)
      gemmi::fail("NDET < 0");
    for (int det = 0; det < ndet; ++det) {
      printf("\n    detector #%d", det+1);
      size_t offset = expected_end;
      expected_end += 40;
      if (b.floats.size() < expected_end)
        gemmi::fail("header too short for no. detectors:", std::to_string(ndet));
      constexpr size_t n = sizeof(batch_det_desc) / sizeof(batch_det_desc[0]);
      for (size_t i = 0; i < 40; ++i) {
        const char* desc = i < n ? batch_det_desc[i] : nullptr;
        size_t oi = offset + i;
        if (desc || b.floats[oi] != 0)
          printf("\n     %4zu %-40s %-.5g", oi, desc ? desc : "", b.floats[oi]);
      }
    }
    for (size_t i = expected_end; i < b.floats.size(); ++i)
      if (b.floats[i] != 0)
        printf("\n     %4zu %-40s %-.5g", i, "??", b.floats[i]);
    printf("\n");
  } else {
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
}

void print_batch_extra_info(const Mtz::Batch& b) {
  gemmi::UnitCell uc = b.get_cell();
  print_cell_parameters("    Unit cell parameters:", uc);
  printf("    Phi start - end: %g - %g\n", b.phi_start(), b.phi_end());
  gemmi::Mat33 u = b.matrix_U();
  for (int i = 0; i != 3; ++i)
    printf("    %s % 10.6f % 10.6f % 10.6f\n",
           i == 0 ? "Orientation matrix U:" : "                     ",
           u.a[i][0], u.a[i][1], u.a[i][2]);
}

void print_cells(const Mtz& mtz) {
  print_cell_parameters("global cell param.:", mtz.cell);
  for (const Mtz::Dataset& ds : mtz.datasets) {
    printf("dataset %d %-8s:", ds.id, ds.dataset_name.c_str());
    print_cell_parameters("", ds.cell);
  }
}

void print_tsv(const Mtz& mtz) {
  size_t ncol = mtz.columns.size();
  for (size_t i = 0; i < ncol; ++i)
    printf("%s%c", mtz.columns[i].label.c_str(), i + 1 != ncol ? '\t' : '\n');
  for (size_t i = 0; i < mtz.nreflections * ncol; ++i)
    printf("%g%c", mtz.data[i], (i + 1) % ncol != 0 ? '\t' : '\n');
}

void print_stats(const Mtz& mtz) {
  struct ColumnStats {
    float min_value = INFINITY;
    float max_value = -INFINITY;
    gemmi::Variance var;
  };
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

  if (!mtz.batches.empty()) {
    const Mtz::Column* col = mtz.column_with_label("BATCH");
    if (!col) {
      printf("Missing column BATCH\n");
      return;
    }
    std::unordered_map<int,int> batch_stat;
    for (float b : *col)
      batch_stat[(int)b]++;
    int min_loc = 0;
    int min_val = INT_MAX;
    int max_loc = 0;
    int max_val = -INT_MAX;
    for (auto b : batch_stat) {
      if (b.second < min_val) {
        min_loc = b.first;
        min_val = b.second;
      }
      if (b.second > max_val) {
        max_loc = b.first;
        max_val = b.second;
      }
    }
    printf("\n%zu MTZ batches, including %d empty. Non-empty batches have\n",
           mtz.batches.size(), (int)mtz.batches.size() - (int)batch_stat.size());
    printf("from %d (in BATCH %d) to %d (in BATCH %d) reflections.\n",
           min_val, min_loc, max_val, max_loc);
  }
}

void print_column_statistics(const Mtz& mtz, const char* label) {
  const Mtz::Column& col = mtz.get_column_with_label(label);
  std::vector<float> data(col.size());
  for (size_t i = 0; i != data.size(); ++i)
    data[i] = col[i];
  gemmi::DataStats st = gemmi::calculate_data_statistics(data);
  std::printf("\nStatistics of column %s:\n", label);
  std::printf("NaN count:  %zu of %zu\n", st.nan_count, data.size());
  if (st.nan_count == data.size())
    return;
  std::printf("Minimum: %12.5f\n", st.dmin);
  std::printf("Maximum: %12.5f\n", st.dmax);
  std::printf("Mean:    %12.5f\n", st.dmean);
  std::printf("RMS:     %12.5f\n", st.rms);
  if (st.nan_count != 0)
    gemmi::vector_remove_if(data, [](float x) { return std::isnan(x); });
  size_t mpos = data.size() / 2;
  std::nth_element(data.begin(), data.begin() + mpos, data.end());
  std::printf("Median:  %12.5f\n", data[mpos]);
  double margin = 0;
  print_histogram(data, st.dmin - margin, st.dmax + margin);
}

void check_asu(const Mtz& mtz, bool tnt) {
  size_t ncol = mtz.columns.size();
  const gemmi::SpaceGroup* sg = mtz.spacegroup;
  if (!sg)
    gemmi::fail("no spacegroup in the MTZ file.");
  int counter = 0;
  gemmi::ReciprocalAsu asu(sg, tnt);
  for (int i = 0; i < mtz.nreflections; ++i) {
    int h = (int) mtz.data[i * ncol + 0];
    int k = (int) mtz.data[i * ncol + 1];
    int l = (int) mtz.data[i * ncol + 2];
    if (asu.is_in({{h, k, l}}))
      ++counter;
  }
  if (!mtz.is_merged())
    printf("NOTE: this is multirecord (unmerged) MTZ file\n");
  printf("spacegroup: %s\n", sg->xhm().c_str());
  printf("%s ASU convention wrt. standard setting: %s\n",
         tnt ? "TNT" : "CCP4", asu.condition_str());
  printf("inside / outside of ASU: %d / %d\n",
         counter, mtz.nreflections - counter);
  double dmin = mtz.resolution_high() - 1e-6;
  printf("All unique reflections up to d=%g: %d\n",
         dmin, gemmi::count_reflections(mtz.cell, mtz.spacegroup, dmin));
}

void compare_mtz(Mtz& mtz1, const char* path2, bool verbose) {
  Mtz mtz2;
  mtz2.read_input(gemmi::MaybeGzipped(path2), true);
  if (mtz1.spacegroup != mtz2.spacegroup)
    printf("Spacegroup differs:  %s  and  %s\n",
           mtz1.spacegroup_name.c_str(), mtz2.spacegroup_name.c_str());
  else if (verbose)
    printf("Spacegroup the same.\n");
  if (mtz1.cell != mtz2.cell)
    printf("Unit cell differs:\n    %g %g %g  %g %g %g\n    %g %g %g  %g %g %g\n",
           mtz1.cell.a, mtz1.cell.b, mtz1.cell.c,
           mtz1.cell.alpha, mtz1.cell.beta, mtz1.cell.gamma,
           mtz2.cell.a, mtz2.cell.b, mtz2.cell.c,
           mtz2.cell.alpha, mtz2.cell.beta, mtz2.cell.gamma);
  else if (verbose)
    printf("Unit cell the same.\n");
  mtz1.sort();
  mtz2.sort();
  // Check if indices are the same. "H" is a dummy value for AsuData.
  {
    gemmi::AsuData<int> ad1 = gemmi::make_asu_data<int>(mtz1, "H", true);
    gemmi::AsuData<int> ad2 = gemmi::make_asu_data<int>(mtz2, "H", true);
    int n = gemmi::count_equal_values(ad1.v, ad2.v);
    if (n != mtz1.nreflections || n != mtz2.nreflections)
      printf("Miller indices differ: %d common (all: %d and %d).\n",
              n, mtz1.nreflections, mtz2.nreflections);
    else
      printf("All Miller indices are the same. Count: %d\n", n);
  }
  for (auto col = mtz1.columns.begin() + 3; col < mtz1.columns.end(); ++col) {
    const Mtz::Column* col2 = mtz2.column_with_label(col->label);
    if (!col2) {
      printf("Missing column: %s\n", col->label.c_str());
      continue;
    }
    if (col->type != col2->type) {
      printf("Type of column %s differs: %c and %c\n",
             col->label.c_str(), col->type, col2->type);
      continue;
    }
    if (col->type == 'F' && col+1 != mtz1.columns.end() && (col+1)->type == 'P') {
      std::array<std::string,2> labels{{col->label, (col+1)->label}};
      auto ad1 = gemmi::make_asu_data<std::complex<float>,2>(mtz1, labels, true);
      auto ad2 = gemmi::make_asu_data<std::complex<float>,2>(mtz2, labels, true);
      gemmi::ComplexCorrelation cor = gemmi::calculate_hkl_complex_correlation(ad1.v, ad2.v);
      std::complex<double> cc = cor.coefficient();
      printf("Column %s/%s: |CC|=%.8g  phase(CC)=%g deg  ratio=%g  n=%d\n",
             col->label.c_str(), (col+1)->label.c_str(),
             std::abs(cc), gemmi::deg(std::arg(cc)), cor.mean_ratio(), cor.n);
      ++col;
    } else if (col->type == 'I' || col->type == 'B' || col->type == 'Y') {
      auto ad1 = gemmi::make_asu_data<float>(mtz1, col->label, true);
      auto ad2 = gemmi::make_asu_data<float>(mtz2, col->label, true);
      int nid = gemmi::count_equal_values(ad1.v, ad2.v);
      printf("Column %s: identical: %d  (all: %zu and %zu)\n",
             col->label.c_str(), nid, ad1.size(), ad2.size());
    } else { // J, D, Q, G, L, K, M, E, P, A, Y
      auto ad1 = gemmi::make_asu_data<float>(mtz1, col->label, true);
      auto ad2 = gemmi::make_asu_data<float>(mtz2, col->label, true);
      int nid = gemmi::count_equal_values(ad1.v, ad2.v);
      printf("Column %s: identical: %d ", col->label.c_str(), nid);
      if ((size_t)nid == ad1.size() && ad1.size() == ad2.size()) {
        printf("(all)\n");
      } else {
        gemmi::Correlation cor = gemmi::calculate_hkl_value_correlation(ad1.v, ad2.v);
        printf("CC=%.8g  ratio=%.8g  n=%d\n",
               cor.coefficient(), cor.mean_ratio(), cor.n);
      }
    }
  }
}

template<typename Stream>
void print_mtz_info(Stream&& stream, const char* path,
                    const std::vector<option::Option>& options) {
  Mtz mtz;
  try {
    mtz.read_first_bytes(stream);
    if (options[ToggleEndian])
      mtz.toggle_endianness();
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
    mtz.warnings = &std::cerr;
  mtz.read_main_headers(stream);
  mtz.read_history_and_batch_headers(stream);
  mtz.setup_spacegroup();
  if (options[PrintTsv] || options[PrintStats] || options[PrintHistogram] ||
      options[CheckAsu] || options[Compare] || options[UpdateReso])
    mtz.read_raw_data(stream);
  if (options[UpdateReso])
    mtz.update_reso();
  if (options[Dump] ||
      !(options[PrintBatch] || options[PrintBatches] || options[PrintTsv] ||
        options[PrintStats] || options[PrintHistogram] ||
        options[PrintAppendix] || options[PrintCells] ||
        options[CheckAsu] || options[Compare] || options[Headers]))
    dump(mtz);
  if (options[PrintBatch]) {
    for (const option::Option* o = options[PrintBatch]; o; o = o->next()) {
      int number = std::atoi(o->arg);
      for (const Mtz::Batch& b : mtz.batches)
        if (b.number == number) {
          bool expanded = options[ExpandedBatches];
          print_batch(b, expanded);
          if (!expanded)
            print_batch_extra_info(b);
        }
    }
  }
  if (options[PrintBatches])
    for (const Mtz::Batch& b : mtz.batches)
      print_batch(b, options[ExpandedBatches]);
  if (options[PrintAppendix])
    printf("%s", mtz.appended_text.c_str());
  if (mtz.has_data() && !options[NoIsym])
    mtz.switch_to_original_hkl();
  if (options[PrintCells])
    print_cells(mtz);
  for (const option::Option* opt = options[PrintHistogram]; opt; opt = opt->next())
    print_column_statistics(mtz, opt->arg);
  if (options[PrintTsv])
    print_tsv(mtz);
  if (options[PrintStats])
    print_stats(mtz);
  if (options[CheckAsu])
    check_asu(mtz, options[CheckAsu].arg[0] == 't');
  if (options[Compare])
    // here mtz gets sorted, so this option must be at the end
    compare_mtz(mtz, options[Compare].arg, options[Verbose]);
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
      if (p.options[Verbose]) {
        std::fflush(stdout);
        std::fprintf(stderr, "Reading %s ...\n", path);
        std::fflush(stderr);
      }
      gemmi::MaybeGzipped input(path);
      if (input.is_stdin()) {
        print_mtz_info(gemmi::FileStream{stdin}, path, p.options);
      } else if (gemmi::CharArray mem = input.uncompress_into_buffer()) {
        print_mtz_info(mem.stream(), path, p.options);
      } else {
        gemmi::fileptr_t f = gemmi::file_open(input.path().c_str(), "rb");
        print_mtz_info(gemmi::FileStream{f.get()}, path, p.options);
      }
    }
  } catch (std::runtime_error& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    std::fflush(stderr);
    return 1;
  }
  return 0;
}
