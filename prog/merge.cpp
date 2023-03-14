// Copyright 2020 Global Phasing Ltd.
//
// Merges multi-record (unmerged) data - I and sigma(I)
// Uses Inverse-variance weighting.

#include <cmath>              // for sqrt
#include <cstdio>             // for fprintf
#include <algorithm>          // for sort
#include <iostream>           // for cout, cerr
#include <gemmi/asudata.hpp>  // for calculate_hkl_value_correlation
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/mtz2cif.hpp>  // for MtzToCif
#include <gemmi/fstream.hpp>  // for Ofstream
#include <gemmi/merge.hpp>
#include <gemmi/read_cif.hpp> // for read_cif_gz
#define GEMMI_PROG merge
#include "options.h"

namespace {

enum OptionIndex {
  WriteAnom=4, NoSysAbs, NumObs, BlockName, Compare, PrintAll
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n  " EXE_NAME " --compare [options] UNMERGED_FILE MERGED_FILE"
    "\n  " EXE_NAME " --compare [options] MMCIF_FILE_WITH_BOTH"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { WriteAnom, 0, "a", "anom", Arg::None,
    "  --anom  \toutput/compare I(+) and I(-) instead of IMEAN." },
  { NoSysAbs, 0, "", "no-sysabs", Arg::None,
    "  --no-sysabs  \tDo not output systematic absences." },
  { NumObs, 0, "", "nobs", Arg::None,
    "  --nobs  \tAdd MTZ column NOBS with the number of merged reflections." },
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \toutput mmCIF block name: data_NAME (default: merged)." },
  { Compare, 0, "", "compare", Arg::None,
    "  --compare  \tcompare unmerged and merged data (no output file)." },
  { PrintAll, 0, "", "print-all", Arg::None,
    "  --print-all  \tprint all compared reflections." },
  { NoOp, 0, "", "", Arg::None,
    "\nThe input file can be SF-mmCIF with _diffrn_refln, MTZ or XDS_ASCII.HKL."
    "\nThe output file can be either SF-mmCIF or MTZ."
  },
  { 0, 0, 0, 0, 0, 0 }
};

using gemmi::Intensities;
using gemmi::DataType;

void output_intensity_statistics(const Intensities& intensities) {
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

void read_intensities_from_rblocks(Intensities& intensities,
                                   DataType data_type,
                                   std::vector<gemmi::ReflnBlock>& rblocks,
                                   const char* block_name, bool verbose) {
  for (gemmi::ReflnBlock& rb : rblocks) {
    if (block_name && rb.block.name != block_name)
      continue;
    rb.use_unmerged(data_type == DataType::Unmerged);
    if (!rb.default_loop)
      continue;
    if (data_type == DataType::Mean && rb.find_column_index("intensity_meas") < 0) {
      if (rb.find_column_index("pdbx_I_plus") < 0)
        gemmi::fail("merged intensities not found");
      std::fprintf(stderr, "No _refln.intensity_meas, using pdbx_I_plus/minus ...\n");
      data_type = DataType::Anomalous;
    }
    if (verbose)
      std::fprintf(stderr, "Reading %s from block %s ...\n",
                   Intensities::type_str(data_type), rb.block.name.c_str());
    intensities.read_mmcif(rb, data_type);
    if (data_type != DataType::Unmerged)
      intensities.take_staraniso_b_from_mmcif(rb.block);
    break;
  }
}

Intensities read_intensities(DataType data_type, const char* input_path,
                             const char* block_name, bool verbose) {
  try {
    Intensities intensities;
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.warnings = &std::cerr;
      mtz.read_input(gemmi::MaybeGzipped(input_path), /*with_data=*/true);
      if (data_type == DataType::Unknown)
        data_type = mtz.batches.empty() ? DataType::Mean : DataType::Unmerged;
      if (data_type == DataType::Mean && !mtz.imean_column()) {
        std::fprintf(stderr, "No IMEAN, using I(+) and I(-) ...\n");
        if (!mtz.iplus_column())
          gemmi::fail("I(+) not found");
        data_type = DataType::Anomalous;
      }
      intensities.read_mtz(mtz, data_type);
      if (data_type != DataType::Unmerged)
        intensities.take_staraniso_b_from_mtz(mtz);
    } else if (gemmi::giends_with(input_path, ".hkl")) {
      gemmi::XdsAscii xds_ascii;
      xds_ascii.read_input(gemmi::MaybeGzipped(input_path));
      intensities.read_unmerged_intensities_from_xds(xds_ascii);
    } else {
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      read_intensities_from_rblocks(intensities, data_type, rblocks, block_name, verbose);
    }
    if (intensities.data.empty())
      gemmi::fail("data not found");
    return intensities;
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR while reading %s: %s\n", input_path, e.what());
    std::exit(1);
  }
}

void write_merged_intensities(const Intensities& intensities, bool write_nobs,
                              const char* output_path) {
  gemmi::Mtz mtz(/*with_base=*/true);
  mtz.spacegroup = intensities.spacegroup;
  mtz.set_cell_for_all(intensities.unit_cell);
  mtz.add_dataset("unknown").wavelength = intensities.wavelength;
  if (intensities.type == DataType::Mean) {
    mtz.add_column("IMEAN", 'J', -1, -1, false);
    mtz.add_column("SIGIMEAN", 'Q', -1, -1, false);
    if (write_nobs)
      mtz.add_column("NOBS", 'I', -1, -1, false);
  } else if (intensities.type == DataType::Anomalous) {
    mtz.add_column("I(+)", 'K', -1, -1, false);
    mtz.add_column("SIGI(+)", 'M', -1, -1, false);
    mtz.add_column("I(-)", 'K', -1, -1, false);
    mtz.add_column("SIGI(-)", 'M', -1, -1, false);
    if (write_nobs) {
      mtz.add_column("NOBS(+)", 'I', -1, -1, false);
      mtz.add_column("NOBS(-)", 'I', -1, -1, false);
    }
  } else {
    return;
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
    if (write_nobs) {
      size_t nobs_offset = offset + 5;  // for "NOBS"
      if (intensities.type == DataType::Anomalous)
        nobs_offset += (refl.isign >= 0 ? 2 : 3);
      mtz.data[nobs_offset] = (float) refl.nobs;
    }
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
      mtz_to_cif.write_cif(mtz, nullptr, nullptr, os.ref());
    }
  } catch (std::exception& e) {
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
  printf("Comparing unmerged and merged reflections...\n");
  if (intensities.spacegroup == ref.spacegroup)
    printf("Space group: %s\n", intensities.spacegroup_str().c_str());
  else
    printf("Space groups: %s   %s\n", intensities.spacegroup_str().c_str(),
           ref.spacegroup_str().c_str());
  printf("Reflections after merging: %zu  %zu\n", intensities.data.size(), ref.data.size());
  intensities.remove_systematic_absences();
  ref.remove_systematic_absences();
  printf("Excluding sys. absences: %zu  %zu\n",
         intensities.data.size(), ref.data.size());
  if (intensities.data.empty())
    return;
  ref.sort();
  if (ref.staraniso_b.ok()) {
    printf("Applying the same anisotropic scaling.\n");
    for (Intensities::Refl& refl : intensities.data) {
      double scale = ref.staraniso_b.scale(refl.hkl, intensities.unit_cell);
      refl.value *= scale;
      refl.sigma *= scale;
    }
  }

  gemmi::Correlation ci = calculate_hkl_value_correlation(intensities.data, ref.data);
  double intensity_ratio = ci.mean_ratio();
  printf("Ratio of compared intensities (merged : unmerged): %g\n", intensity_ratio);
  for (Intensities::Refl& refl : intensities.data) {
    refl.value *= intensity_ratio;
    refl.sigma *= intensity_ratio;
  }

  auto a = intensities.data.begin();
  gemmi::Correlation cs; // correlation of sigma
  int different_intensity_count = 0;
  int different_sigma_count = 0;
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
    //ci.add_point(r.value, a->value);
    cs.add_point(r.sigma, a->sigma);

    using gemmi::sq;
    double sq_max = std::max(sq(r.value), sq(a->value));
    double sq_diff = sq(r.value - a->value);
    if (sq_diff > 1e-4 && sq_diff > sq(0.005) * sq_max)
      different_intensity_count++;
    sq_max = std::max(sq(r.sigma), sq(a->sigma));
    sq_diff = sq(r.sigma - a->sigma);
    if (sq_diff > 1e-4 && sq_diff > sq(0.005) * sq_max) {
      if (different_sigma_count == 0)
        printf("First difference: %s %g vs %g\n",
               r.hkl_label().c_str(), a->value, r.value);
      different_sigma_count++;
    }

    if (print_all)
      print_reflection(&*a, &r);
    ++a;
    if (a == intensities.data.end())
      break;
  }
  printf("Common reflections: %d\n", ci.n);
  printf("%d of intensities and %d of sigmas differ by >0.5%%.\n",
         different_intensity_count, different_sigma_count);
  auto a_resol = intensities.resolution_range();
  ref.unit_cell = intensities.unit_cell;
  auto r_resol = ref.resolution_range();
  printf("Resolution: %g-%g    %g-%g\n", a_resol[0], a_resol[1], r_resol[0], r_resol[1]);
  printf("%s CC: %.9g%%\n", intensities.type_str(), 100 * ci.coefficient());
  printf("Sigma CC: %.9g%% (mean ratio: %g)\n", 100 * cs.coefficient(), 1./cs.mean_ratio());
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.nonOptionsCount() != 2 &&
      !(p.nonOptionsCount() == 1 && p.options[Compare])) {
    fprintf(stderr, "%s requires 2 arguments (or single arg with --compare), got %d.",
                    p.program_name, p.nonOptionsCount());
    p.print_try_help_and_exit("");
  }
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = nullptr;
  if (p.nonOptionsCount() == 2)
    output_path = p.nonOption(1);
  DataType otype = p.options[WriteAnom] ? DataType::Anomalous : DataType::Mean;
  const char* block_name = nullptr;
  if (p.options[BlockName])
    block_name = p.options[BlockName].arg;

  Intensities ref;
  if (p.options[Compare] && output_path) {
    if (verbose)
      std::fprintf(stderr, "Reading merged reflections from %s ...\n", output_path);
    ref = read_intensities(otype, output_path, nullptr, verbose);
  }

  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path);
  try {
    Intensities intensities;
    if (output_path) {
      DataType data_type = DataType::Unmerged;
      if (p.options[Compare] && gemmi::giends_with(input_path, ".mtz"))
        // it's OK to compare also two merged files
        data_type = DataType::Unknown;
      intensities = read_intensities(data_type, input_path, block_name, verbose);
    } else { // special case of --compare with one mmCIF file
      if (gemmi::giends_with(input_path, ".mtz") ||
          gemmi::giends_with(input_path, ".hkl"))
        gemmi::fail("`--compare ONE_FILE' make sense only with mmCIF.");
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      read_intensities_from_rblocks(ref, otype, rblocks, nullptr, verbose);
      read_intensities_from_rblocks(intensities, DataType::Unmerged,
                                    rblocks, block_name, verbose);
      if (intensities.data.empty())
        gemmi::fail("unmerged data not found");
    }
    if (verbose)
      output_intensity_statistics(intensities);
    if (p.options[Compare]) {
      if (intensities.type != ref.type)
        intensities.merge_in_place(ref.type);
      compare_intensities(intensities, ref, p.options[PrintAll]);
    } else {
      intensities.merge_in_place(otype);
      if (p.options[NoSysAbs])
        intensities.remove_systematic_absences();
      if (verbose)
        std::fprintf(stderr, "Writing %zu reflections to %s ...\n",
                     intensities.data.size(), output_path);
      write_merged_intensities(intensities, p.options[NumObs], output_path);
    }
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
