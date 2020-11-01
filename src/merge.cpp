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
#include <gemmi/fstream.hpp>  // for Ofstream
#include <gemmi/refln.hpp>
#include <gemmi/merge.hpp>
#define GEMMI_PROG merge
#include "options.h"

namespace {

enum OptionIndex { WriteAnom=4, BlockName, Compare, PrintAll };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n  " EXE_NAME " --compare [options] UNMERGED_FILE MERGED_FILE"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { WriteAnom, 0, "a", "write-anom", Arg::None,
    "  --write-anom  \toutput I(+) and I(-) instead of IMEAN." },
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \toutput mmCIF block name: data_NAME (default: merged)." },
  { Compare, 0, "", "compare", Arg::None,
    "  --compare  \tcompare unmerged and merged data (no output file)." },
  { PrintAll, 0, "", "print-all", Arg::None,
    "  --print-all  \tprint all compared reflections." },
  { NoOp, 0, "", "", Arg::None,
    "\nThe input file can be either SF-mmCIF with _diffrn_refln or MTZ."
    "\nThe output file can be either SF-mmCIF or MTZ."
  },
  { 0, 0, 0, 0, 0, 0 }
};

using gemmi::Intensities;

Intensities read_unmerged_intensities_from_mmcif(const gemmi::ReflnBlock& rb) {
  size_t value_idx = rb.get_column_index("intensity_net");
  size_t sigma_idx = rb.get_column_index("intensity_sigma");
  Intensities intensities;
  intensities.copy_metadata(rb);
  intensities.wavelength = rb.wavelength;
  intensities.read_data(gemmi::ReflnDataProxy(rb), value_idx, sigma_idx);
  // switch to ASU indices
  gemmi::GroupOps gops = intensities.spacegroup->operations();
  gemmi::ReciprocalAsu asu(intensities.spacegroup);
  for (Intensities::Refl& refl : intensities.data) {
    auto hkl_isym = asu.to_asu(refl.hkl, gops);
    refl.hkl = hkl_isym.first;
    refl.isign = (hkl_isym.second % 2 == 0 ? -1 : 1);
  }
  return intensities;
}

Intensities read_mean_intensities_from_mmcif(const gemmi::ReflnBlock& rb) {
  size_t value_idx = rb.get_column_index("intensity_meas");
  size_t sigma_idx = rb.get_column_index("intensity_sigma");
  Intensities intensities;
  intensities.copy_metadata(rb);
  intensities.wavelength = rb.wavelength;
  intensities.read_data(gemmi::ReflnDataProxy(rb), value_idx, sigma_idx);
  return intensities;
}

Intensities read_anomalous_intensities_from_mmcif(const gemmi::ReflnBlock& rb) {
  size_t value_idx[2] = {rb.get_column_index("pdbx_I_plus"),
                         rb.get_column_index("pdbx_I_minus")};
  size_t sigma_idx[2] = {rb.get_column_index("pdbx_I_plus_sigma"),
                         rb.get_column_index("pdbx_I_minus_sigma")};
  Intensities intensities;
  intensities.copy_metadata(rb);
  intensities.wavelength = rb.wavelength;
  gemmi::ReflnDataProxy proxy(rb);
  for (size_t i = 0; i < proxy.size(); i += proxy.stride())
    for (int j = 0; j < 2; ++j) {
      Intensities::Refl refl;
      refl.hkl = proxy.get_hkl(i);
      refl.isign = (j == 0 ? 1 : -1);
      refl.value = proxy.get_num(i + value_idx[j]);
      refl.sigma = proxy.get_num(i + sigma_idx[j]);
      intensities.add_if_valid(refl);
    }
  return intensities;
}

enum class ReadType { Unmerged, Mean, Anomalous };

Intensities read_intensities(ReadType itype, const char* input_path,
                             const char* block_name, bool verbose) {
  try {
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.warnings = stderr;
      mtz.read_input(gemmi::MaybeGzipped(input_path), /*with_data=*/true);
      if (itype == ReadType::Mean && !mtz.column_with_label("IMEAN")) {
        std::fprintf(stderr, "No IMEAN, using I(+) and I(-) ...\n");
        if (!mtz.column_with_label("I(+)"))
          gemmi::fail("I(+) not found");
        itype = ReadType::Anomalous;
      }
      switch (itype) {
        case ReadType::Unmerged:
          return gemmi::read_unmerged_intensities_from_mtz(mtz);
        case ReadType::Mean:
          return gemmi::read_mean_intensities_from_mtz(mtz);
        case ReadType::Anomalous:
          return gemmi::read_anomalous_intensities_from_mtz(mtz);
      }
      gemmi::unreachable();
    } else {
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      for (gemmi::ReflnBlock& rb : rblocks) {
        if (!block_name || rb.block.name == block_name) {
          rb.use_unmerged(itype == ReadType::Unmerged);
          if (rb.default_loop) {
            if (itype == ReadType::Mean && rb.find_column_index("intensity_meas") < 0) {
              if (rb.find_column_index("pdbx_I_plus") < 0)
                gemmi::fail("merged intensities not found");
              std::fprintf(stderr, "No _refln.intensity_meas, using pdbx_I_plus/minus ...\n");
              itype = ReadType::Anomalous;
            }
            if (verbose)
              std::fprintf(stderr, "Taking %s from block %s ...\n",
                           itype == ReadType::Unmerged ? "unmerged I" :
                           itype == ReadType::Mean ? "<I>" : "I+/I-",
                           rb.block.name.c_str());
            switch (itype) {
              case ReadType::Unmerged:
                return read_unmerged_intensities_from_mmcif(rb);
              case ReadType::Mean:
                return read_mean_intensities_from_mmcif(rb);
              case ReadType::Anomalous:
                return read_anomalous_intensities_from_mmcif(rb);
            }
          }
        }
      }
      gemmi::fail("data not found");
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while reading %s: %s\n", input_path, e.what());
    std::exit(1);
  }
}

void write_merged_intensities(const Intensities& intensities, const char* output_path) {
  gemmi::Mtz mtz;
  mtz.cell = intensities.unit_cell;
  mtz.spacegroup = intensities.spacegroup;
  mtz.add_dataset("HKL_base");
  mtz.add_column("H", 'H');
  mtz.add_column("K", 'H');
  mtz.add_column("L", 'H');
  mtz.add_dataset("unknown").wavelength = intensities.wavelength;
  if (!intensities.have_sign()) {
    mtz.add_column("IMEAN", 'J');
    mtz.add_column("SIGIMEAN", 'Q');
  } else {
    mtz.add_column("I(+)", 'K');
    mtz.add_column("SIGI(+)", 'M');
    mtz.add_column("I(-)", 'K');
    mtz.add_column("SIGI(-)", 'M');
  }
  mtz.data.resize(intensities.data.size() * mtz.columns.size(), NAN);
  gemmi::Miller prev_hkl = intensities.data[0].hkl;
  mtz.set_hkl(0, prev_hkl);
  size_t offset = 0;
  for (const Intensities::Refl& refl : intensities.data) {
    if (refl.hkl != prev_hkl) {
      offset += mtz.columns.size();
      mtz.set_hkl(offset, refl.hkl);
      prev_hkl = refl.hkl;
    }
    size_t value_offset = offset + (refl.isign >= 0 ? 3 : 5);
    mtz.data[value_offset] = (float) refl.value;
    mtz.data[value_offset + 1] = (float) refl.sigma;
  }
  mtz.data.resize(offset + mtz.columns.size());
  mtz.nreflections = int(mtz.data.size() / mtz.columns.size());
  try {
    if (gemmi::giends_with(output_path, ".mtz")) {
      mtz.write_to_file(output_path);
    } else {
      gemmi::MtzToCif mtz_to_cif;
      mtz_to_cif.with_comments = false;
      mtz_to_cif.block_name = "merged";
      gemmi::Ofstream os(output_path, /*dash=*/&std::cout);
      mtz_to_cif.write_cif(mtz, nullptr, os.ref());
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while writing %s: %s\n", output_path, e.what());
    std::exit(1);
  }
}

void print_reflection(const Intensities::Refl* a, const Intensities::Refl* r) {
  const Intensities::Refl& h = a ? *a : *r;
  printf("   %3d %3d %3d  %c ",
         h.hkl[0], h.hkl[1], h.hkl[2],
         h.isign == 0 ? 'm' : h.isign > 0 ? '+' : '-');
  if (a && r)
    printf("%8.2f vs %8.2f\n", a->value, r->value);
  else if (a)
    printf("%8.2f vs   N/A\n", a->value);
  else
    printf("  N/A    vs %8.2f\n", r->value);
}

void compare_intensities(Intensities& intensities, Intensities& ref, bool print_all) {
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
  if (intensities.data.empty())
    return;
  ref.sort();
  auto a = intensities.data.begin();
  gemmi::Correlation ci, cs; // correlation of <I>, Isigma, I+, I-
  for (const Intensities::Refl& r : ref.data) {
    if (r.hkl != a->hkl || r.isign != a->isign) {
      while (*a < r) {
        if (print_all)
          print_reflection(&*a, nullptr);
        ++a;
        if (a == intensities.data.end())
          break;
      }
      if (a == intensities.data.end())
        break;
      if (r.hkl != a->hkl || r.isign != a->isign) {
        if (print_all)
          print_reflection(nullptr, &r);
        continue;
      }
    }
    ci.add_point(r.value, a->value);
    cs.add_point(r.sigma, a->sigma);
    if (print_all)
      print_reflection(&*a, &r);
    ++a;
    if (a == intensities.data.end())
      break;
  }
  printf("Common reflections: %d\n", ci.n);
  auto a_resol = intensities.resolution_range();
  ref.unit_cell = intensities.unit_cell;
  auto r_resol = ref.resolution_range();
  printf("Resolution: %g-%g    %g-%g\n", a_resol[0], a_resol[1], r_resol[0], r_resol[1]);
  const char* what = intensities.have_sign() ? "I+/I-" : "Imean";
  printf("%s CC: %.9g%% (mean ratio: %g)\n", what, 100 * ci.coefficient(), ci.mean_ratio());
  printf("Sigma CC: %.9g%% (mean ratio: %g)\n", 100 * cs.coefficient(), cs.mean_ratio());
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);

  Intensities ref;
  if (p.options[Compare]) {
    if (verbose)
      std::fprintf(stderr, "Reading merged reflections from %s ...\n", output_path);
    ReadType itype = p.options[WriteAnom] ? ReadType::Anomalous : ReadType::Mean;
    ref = read_intensities(itype, output_path, nullptr, verbose);
  }

  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path);
  const char* block_name = nullptr;
  if (p.options[BlockName])
    block_name = p.options[BlockName].arg;
  Intensities intensities = read_intensities(ReadType::Unmerged,
                                             input_path, block_name, verbose);
  if (verbose) {
    size_t plus_count = 0;
    size_t minus_count = 0;
    for (const Intensities::Refl& refl : intensities.data) {
      if (refl.isign > 0)
        ++plus_count;
      else if (refl.isign < 0)
        ++minus_count;
    }
    std::fprintf(stderr, "Merging observations (%zu total, %zu for I+, %zu for I-) ...\n",
                 intensities.data.size(), plus_count, minus_count);
  }
  bool output_plus_minus = (p.options[Compare] ? ref.have_sign() : p.options[WriteAnom]);
  intensities.merge_in_place(output_plus_minus);
  if (p.options[Compare]) {
    compare_intensities(intensities, ref, p.options[PrintAll]);
  } else {
    if (verbose)
      std::fprintf(stderr, "Writing %zu reflections to %s ...\n",
                   intensities.data.size(), output_path);
    write_merged_intensities(intensities, output_path);
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
