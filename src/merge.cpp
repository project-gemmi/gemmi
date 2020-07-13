// Copyright 2020 Global Phasing Ltd.
//
// merge multi-record (unmerged) data - I and sigma(I)

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

enum OptionIndex { BlockName=4, };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \toutput mmCIF block name: data_NAME (default: merged)." },
  { NoOp, 0, "", "", Arg::None,
    "\nThe input file can be either SF-mmCIF with _diffrn_refln or MTZ."
    "\nThe output file can be either SF-mmCIF or MTZ."
  },
  { 0, 0, 0, 0, 0, 0 }
};

struct Intensity {
  gemmi::Miller hkl;
  bool plus;  // I(+) if true, I(-) otherwise
  double value;
  double sigma;
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
      if (intensity.sigma >= 0)
        data.push_back(intensity);
    }
  }
};

static
Intensities read_unmerged_intensities(const char* input_path,
                                      const char* block_name, bool verbose) {
  Intensities intensities;
  try {
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.warnings = stderr;
      mtz.read_input(gemmi::MaybeGzipped(input_path), /*with_data=*/true);
      if (mtz.batches.empty() && !mtz.column_with_label("M/ISYM")) {
        std::fprintf(stderr, "The input mtz file should be unmerged.\n");
        std::exit(1);
      }
      const gemmi::Mtz::Column& col = mtz.get_column_with_label("I");
      size_t value_idx = col.idx;
      size_t sigma_idx = mtz.get_column_with_label("SIGI").idx;
      intensities.wavelength = mtz.dataset(col.dataset_id).wavelength;
      intensities.read(gemmi::MtzDataProxy{mtz}, value_idx, sigma_idx);
    } else {
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      for (gemmi::ReflnBlock& rb : rblocks) {
        if (block_name && rb.block.name != block_name)
          continue;
        if (rb.diffrn_refln_loop) {
          size_t value_idx = rb.get_column_index("intensity_net");
          size_t sigma_idx = rb.get_column_index("intensity_sigma");
          rb.use_unmerged(true);
          if (verbose)
            std::fprintf(stderr, "Taking data from %s ...\n", rb.block.name.c_str());
          intensities.wavelength = rb.wavelength;
          intensities.read(gemmi::ReflnDataProxy(rb), value_idx, sigma_idx);
          if (!intensities.spacegroup) {
            std::fprintf(stderr, "Unknown space group\n");
            std::exit(1);
          }
          gemmi::GroupOps gops = intensities.spacegroup->operations();
          gemmi::ReciprocalAsu asu(intensities.spacegroup);
          for (Intensity& intensity : intensities.data)
            intensity.hkl = asu.to_asu(intensity.hkl, gops);
          break;
        }
      }
      if (intensities.data.empty()) {
        std::fprintf(stderr, "Unmerged data not found in: %s\n", input_path);
        std::exit(1);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while reading %s: %s\n", input_path, e.what());
    std::exit(1);
  }
  return intensities;
}

static
void merge_intensities_in_place(Intensities& intensities) {
  std::sort(intensities.data.begin(), intensities.data.end(),
            [&](const Intensity& a, const Intensity& b) {
    return std::tie(a.hkl[0], a.hkl[1], a.hkl[2]) < std::tie(b.hkl[0], b.hkl[1], b.hkl[2]);
  });
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
  Intensities intensities = read_unmerged_intensities(input_path, block_name, verbose);
  if (verbose)
    std::fprintf(stderr, "Merging observations (%zu total) ...\n", intensities.data.size());
  merge_intensities_in_place(intensities);
  if (verbose)
    std::fprintf(stderr, "Merged reflections: %zu\n", intensities.data.size());
  if (verbose)
    std::fprintf(stderr, "Writing %s ...\n", output_path);
  write_merged_intensities(intensities, output_path);
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
