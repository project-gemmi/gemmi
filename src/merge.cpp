// Copyright 2020 Global Phasing Ltd.
//
// Merges multi-record (unmerged) data - I and sigma(I)
// Uses Inverse-variance weighting.

#include <cmath>              // for sqrt
#include <cstdio>
#include <algorithm>          // for sort
#include <iostream>
#include <tuple>              // for tie
#ifndef GEMMI_ALL_IN_ONE
# define GEMMI_WRITE_IMPLEMENTATION 1
#endif
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/gzread.hpp>
#include <gemmi/mtz2cif.hpp>  // for MtzToCif
#include <gemmi/ofstream.hpp> // for Ofstream
#include <gemmi/refln.hpp>
#define GEMMI_PROG merge
#include "options.h"

enum OptionIndex { BlockName=4, Compare };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n  " EXE_NAME " --compare [options] UNMERGED_FILE MERGED_FILE"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \toutput mmCIF block name: data_NAME (default: merged)." },
  { Compare, 0, "", "compare", Arg::None, "" },
  { NoOp, 0, "", "", Arg::None,
    "The input file can be either SF-mmCIF with _diffrn_refln or MTZ.\n"
    "The output file can be either SF-mmCIF or MTZ."
  },
  { 0, 0, 0, 0, 0, 0 }
};

struct Intensity {
  gemmi::Miller hkl;
  bool plus;  // I(+) if true, I(-) otherwise
  double value;
  double sigma;

  bool operator<(const Intensity& o) const {
    return std::tie(hkl[0], hkl[1], hkl[2]) < std::tie(o.hkl[0], o.hkl[1], o.hkl[2]);
  }
};

struct Intensities {
  std::vector<Intensity> data;
  const gemmi::SpaceGroup* spacegroup = nullptr;
  gemmi::UnitCell unit_cell;
  double wavelength;

  template<typename DataProxy>
  void read(const DataProxy& proxy, size_t value_idx, size_t sigma_idx) {
    unit_cell = proxy.unit_cell();
    spacegroup = proxy.spacegroup();
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      Intensity intensity;
      intensity.hkl = proxy.get_hkl(i);
      intensity.value = proxy.get_num(i + value_idx);
      intensity.sigma = proxy.get_num(i + sigma_idx);
      // XDS marks rejected reflections with negative sigma
      if (!std::isnan(intensity.value) && intensity.sigma >= 0)
        data.push_back(intensity);
    }
  }

  std::string spacegroup_str() const { return spacegroup ? spacegroup->xhm() : "none"; }

  void remove_systematic_absences() {
    if (!spacegroup)
      return;
    gemmi::GroupOps gops = spacegroup->operations();
    gemmi::vector_remove_if(data, [&](Intensity& x) {
        return gops.is_systematically_absent(x.hkl);
    });
  }

  void sort() { std::sort(data.begin(), data.end()); }
};

[[noreturn]]
inline void exit_message(const char* msg, const char* path=nullptr) {
  if (path)
    std::fprintf(stderr, "%s: %s\n", msg, path);
  else
    std::fprintf(stderr, "%s\n", msg);
  std::exit(1);
}

static
Intensities read_intensities(bool merged, const char* input_path,
                             const char* block_name, bool verbose) {
  Intensities intensities;
  try {
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.warnings = stderr;
      mtz.read_input(gemmi::MaybeGzipped(input_path), /*with_data=*/true);
      if (merged && !mtz.batches.empty())
        exit_message("The second file should be merged.");
      if (!merged && mtz.batches.empty() && !mtz.column_with_label("M/ISYM"))
        exit_message("The input mtz file should be unmerged.");
      const std::string value_label = merged ? "IMEAN" : "I";
      const std::string sigma_label = "SIG" + value_label;
      const gemmi::Mtz::Column& col = mtz.get_column_with_label(value_label);
      size_t sigma_idx = mtz.get_column_with_label(sigma_label).idx;
      size_t value_idx = col.idx;
      intensities.wavelength = mtz.dataset(col.dataset_id).wavelength;
      intensities.read(gemmi::MtzDataProxy{mtz}, value_idx, sigma_idx);
    } else {
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      const std::string value_label = merged ? "intensity_meas" : "intensity_net";
      const std::string sigma_label = "intensity_sigma";
      for (gemmi::ReflnBlock& rb : rblocks) {
        if (block_name && rb.block.name != block_name)
          continue;
        rb.use_unmerged(!merged);
        if (rb.default_loop) {
          size_t value_idx = rb.get_column_index(value_label);
          size_t sigma_idx = rb.get_column_index(sigma_label);
          if (verbose)
            std::fprintf(stderr, "Taking data from %s ...\n", rb.block.name.c_str());
          intensities.wavelength = rb.wavelength;
          intensities.read(gemmi::ReflnDataProxy(rb), value_idx, sigma_idx);
          if (!intensities.spacegroup)
            exit_message("Unknown space group", input_path);
          if (!merged) {
            gemmi::GroupOps gops = intensities.spacegroup->operations();
            gemmi::ReciprocalAsu asu(intensities.spacegroup);
            for (Intensity& intensity : intensities.data)
              intensity.hkl = asu.to_asu(intensity.hkl, gops);
          }
          break;
        }
      }
      if (intensities.data.empty())
        exit_message("Unmerged data not found in", input_path);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while reading %s: %s\n", input_path, e.what());
    std::exit(1);
  }
  return intensities;
}

static
void merge_intensities_in_place(Intensities& intensities) {
  intensities.sort();
  std::vector<Intensity>::iterator out = intensities.data.begin();
  double sum_wI = 0.;
  double sum_w = 0.;
  for (auto in = intensities.data.begin(); in != intensities.data.end(); ++in) {
    if (out->hkl != in->hkl) {
      out->value = sum_wI / sum_w;
      out->sigma = 1.0 / std::sqrt(sum_w);
      sum_wI = sum_w = 0.;
      ++out;
      if (out != in)
        out->hkl = in->hkl;
    }
    double w = 1. / (in->sigma * in->sigma);
    sum_wI += w * in->value;
    sum_w += w;
  }
  out->value = sum_wI / sum_w;
  out->sigma = 1.0 / std::sqrt(sum_w);
  intensities.data.erase(++out, intensities.data.end());
}

static
void write_merged_intensities(const Intensities& intensities, const char* output_path) {
  gemmi::Mtz mtz;
  mtz.cell = intensities.unit_cell;
  mtz.spacegroup = intensities.spacegroup;
  mtz.add_dataset("HKL_base");
  mtz.add_column("H", 'H');
  mtz.add_column("K", 'H');
  mtz.add_column("L", 'H');
  mtz.add_dataset("unknown").wavelength = intensities.wavelength;
  mtz.add_column("IMEAN", 'J');
  mtz.add_column("SIGIMEAN", 'Q');
  mtz.nreflections = (int) intensities.data.size();
  mtz.data.resize(intensities.data.size() * mtz.columns.size());
  for (size_t i = 0; i != intensities.data.size(); ++i) {
    size_t offset = i * mtz.columns.size();
    mtz.set_hkl(offset, intensities.data[i].hkl);
    mtz.data[offset + 3] = (float) intensities.data[i].value;
    mtz.data[offset + 4] = (float) intensities.data[i].sigma;
  }
  try {
    if (gemmi::giends_with(output_path, ".mtz")) {
      mtz.write_to_file(output_path);
    } else {
      gemmi::MtzToCif mtz_to_cif;
      mtz_to_cif.with_comments = false;
      mtz_to_cif.block_name = "merged";
      gemmi::Ofstream os(output_path, /*dash=*/&std::cout);
      mtz_to_cif.write_cif(mtz, os.ref());
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while writing %s: %s\n", output_path, e.what());
    std::exit(1);
  }
}

static
void compare_intensities(Intensities& intensities, Intensities& ref) {
  if (intensities.spacegroup == ref.spacegroup)
    printf("Space group: %s\n", intensities.spacegroup_str().c_str());
  else
    printf("Space groups: %s   %s\n", intensities.spacegroup_str().c_str(),
           ref.spacegroup_str().c_str());
  printf("Reflections: %zu  %zu\n", intensities.data.size(), ref.data.size());
  intensities.remove_systematic_absences();
  ref.remove_systematic_absences();
  printf("Excluding sys. absences: %zu  %zu\n",
         intensities.data.size(), ref.data.size());
  const gemmi::UnitCell& cell = intensities.unit_cell;
  ref.sort();
  auto a = intensities.data.begin();
  gemmi::Correlation vc, sc;
  for (const Intensity& r : ref.data) {
    if (r.hkl != a->hkl) {
      while (*a < r && ++a != intensities.data.end()) {
      }
      if (a == intensities.data.end())
        break;
      if (r.hkl != a->hkl)
        continue;
    }
    vc.add_point(r.value, a->value);
    sc.add_point(r.sigma, a->sigma);
    ++a;
    if (a == intensities.data.end())
      break;
  }
  printf("Common reflections: %d\n", vc.n);
  double min_a_1_d2 = INFINITY;
  double max_a_1_d2 = 0;
  double min_r_1_d2 = INFINITY;
  double max_r_1_d2 = 0;
  for (const Intensity& x : intensities.data) {
    double a_1_d2 = cell.calculate_1_d2(x.hkl);
    if (a_1_d2 < min_a_1_d2)
      min_a_1_d2 = a_1_d2;
    if (a_1_d2 > max_a_1_d2)
      max_a_1_d2 = a_1_d2;
  }
  for (const Intensity& x : ref.data) {
    double r_1_d2 = cell.calculate_1_d2(x.hkl);
    if (r_1_d2 < min_r_1_d2)
      min_r_1_d2 = r_1_d2;
    if (r_1_d2 > max_r_1_d2)
      max_r_1_d2 = r_1_d2;
  }
  printf("Resolution: %g-%g    %g-%g\n",
         1/std::sqrt(min_a_1_d2), 1/std::sqrt(max_a_1_d2),
         1/std::sqrt(min_r_1_d2), 1/std::sqrt(max_r_1_d2));
  printf("Value CC: %.9g%% (slope: %g)\n", 100 * vc.coefficient(), vc.slope());
  printf("Sigma CC: %.9g%% (slope: %g)\n", 100 * sc.coefficient(), sc.slope());
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);
  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path);
  const char* block_name = nullptr;
  if (p.options[BlockName])
    block_name = p.options[BlockName].arg;
  Intensities intensities = read_intensities(false, input_path, block_name, verbose);
  if (verbose)
    std::fprintf(stderr, "Merging observations (%zu total) ...\n", intensities.data.size());
  merge_intensities_in_place(intensities);
  if (p.options[Compare]) {
    if (verbose)
      std::fprintf(stderr, "Reading merged reflections from %s ...\n", output_path);
    Intensities ref = read_intensities(true, input_path, nullptr, verbose);
    compare_intensities(intensities, ref);
  } else {
    if (verbose)
      std::fprintf(stderr, "Writing %zu reflections to %s ...\n",
                   intensities.data.size(), output_path);
    write_merged_intensities(intensities, output_path);
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
