// Copyright 2020 Global Phasing Ltd.
//
// Merges multi-record (unmerged) data - I and sigma(I)
// Uses Inverse-variance weighting.

#include <cstdio>             // for fprintf
#include <algorithm>          // for sort
#include <iostream>           // for cout
#include <gemmi/asudata.hpp>  // for calculate_hkl_value_correlation
#include <gemmi/binner.hpp>   // for Binner
#include <gemmi/mtz2cif.hpp>  // for MtzToCif
#include <gemmi/fstream.hpp>  // for Ofstream
#include <gemmi/intensit.hpp> // for Intensities
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/read_cif.hpp> // for read_cif_gz
#define GEMMI_PROG merge
#include "options.h"

namespace {

enum OptionIndex {
  WriteAnom=4, NoSysAbs, BatchCell, NumObs, InputBlock, OutputBlock, Stats, Compare, PrintAll
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n  " EXE_NAME " --stats [options] UNMERGED_FILE"
    "\n  " EXE_NAME " --compare [options] UNMERGED_FILE MERGED_FILE"
    "\n  " EXE_NAME " --compare [options] MMCIF_FILE_WITH_BOTH"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { WriteAnom, 0, "a", "anom", Arg::None,
    "  --anom  \tOutput/compare I(+) and I(-) instead of IMEAN." },
  { NoSysAbs, 0, "", "no-sysabs", Arg::None,
    "  --no-sysabs  \tRemove systematic absences." },
  { BatchCell, 0, "", "batch-cell", Arg::None,
    "  --batch-cell  \tUse cell parameters from MTZ batch headers." },
  { NumObs, 0, "", "nobs", Arg::None,
    "  --nobs  \tAdd MTZ column NOBS with the number of merged reflections." },
  { InputBlock, 0, "", "input-block", Arg::Required,
    "  --input-block=NAME  \tInput mmCIF block name (default: first unmerged)." },
  { OutputBlock, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \tOutput mmCIF block name: data_NAME (default: merged)." },
  { Stats, 0, "", "stats", Arg::Optional,
    "  --stats[=N]  \tPrint data metrics in N resol. shells (default: 10)."
    "\n\tAdd U (e.g. =10U or =U) for unweighted metrics."
    "\n\tAdd X for XDS-like weighted metrics (details in docs)."
    "\n\tAdd s or e for different binning (more in docs)." },
  { Compare, 0, "", "compare", Arg::None,
    "  --compare  \tCompare unmerged and merged data (no output file)." },
  { PrintAll, 0, "", "print-all", Arg::None,
    "  --print-all  \tPrint all compared reflections." },
  { NoOp, 0, "", "", Arg::None,
    "\nThe input file can be SF-mmCIF with _diffrn_refln, MTZ or XDS_ASCII.HKL."
    "\nThe output file can be either SF-mmCIF or MTZ."
  },
  { 0, 0, 0, 0, 0, 0 }
};

using gemmi::Intensities;
using gemmi::DataType;
using gemmi::MergingStats;

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

void read_intensities(Intensities& intensities, DataType data_type,
                      const std::string& input_path, const char* block_name,
                      bool verbose, bool batch_cell) {
  try {
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.logger.callback = gemmi::Logger::to_stderr;
      mtz.read_file_gz(input_path);
      intensities.import_mtz(mtz, data_type);
      if (data_type != DataType::UAM)
        intensities.take_staraniso_b_from_mtz(mtz);
      if (!mtz.is_merged()) {
        if (intensities.unit_cell != mtz.cell) {
          auto print_cell = [](const char* label, gemmi::UnitCellParameters& p) {
            std::fprintf(stderr, "%s:  %-8g %-8g %-8g   %-8g %-8g %-8g\n",
                         label, p.a, p.b, p.c, p.alpha, p.beta, p.gamma);
          };
          print_cell("Global unit cell from MTZ header", mtz.cell);
          print_cell("Mean of batch headers unit cells", intensities.unit_cell);
          if (!batch_cell) {
            intensities.unit_cell = mtz.cell;
            std::fprintf(stderr, "Using the former. Use option --batch-cell for the latter.\n");
          }
          std::fprintf(stderr, "\n");
        }
      }
    } else if (gemmi::giends_with(input_path, "hkl")) {  // .hkl or .ahkl
      gemmi::XdsAscii xds_ascii = gemmi::read_xds_ascii(input_path);
      intensities.import_xds(xds_ascii);

    } else {  // mmCIF
      auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
      gemmi::ReflnBlock* rblock = nullptr;
      if (block_name) {
        for (gemmi::ReflnBlock& rb : rblocks)
          if (rb.block.name == block_name) {
            rblock = &rb;
            break;
          }
        if (!rblock) {
          std::fprintf(stderr, "Error. No block data_%s in mmCIF file!\n", block_name);
          std::exit(1);
        }
      } else if (data_type == DataType::UAM) {
        for (gemmi::ReflnBlock& rb : rblocks)
          if (rb.diffrn_refln_loop) {
            rblock = &rb;
            break;
          }
      }
      if (!rblock)
        rblock = &rblocks.at(0);
      intensities.import_refln_block(*rblock, data_type);
      if (verbose)
        std::fprintf(stderr, "Got %s from block %s ...\n",
                     Intensities::type_str(intensities.type), rblock->block.name.c_str());
      if (data_type != DataType::UAM)
        intensities.take_staraniso_b_from_mmcif(rblock->block);
    }

    if (intensities.data.empty())
      gemmi::fail("data not found");
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR while reading %s: %s\n", input_path.c_str(), e.what());
    std::exit(1);
  }
}

void write_merged_intensities(const gemmi::Mtz& mtz, const std::string& output_path) {
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
    std::fprintf(stderr, "ERROR while writing %s: %s\n", output_path.c_str(), e.what());
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

// cf. validate_merged_intensities()
void compare_intensities(Intensities& intensities, Intensities& ref, bool print_all) {
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

double take_average(const std::vector<MergingStats>& stats,
                    double (MergingStats::*method)() const) {
  double sum = 0;
  int count = 0;
  for (const MergingStats& ms : stats) {
    sum += ms.stats_refl * (ms.*method)();
    count += ms.stats_refl;
  }
  return sum / count;
}

void print_merging_statistics(const std::vector<MergingStats>& stats,
                              const gemmi::Binner* binner) {
  if (binner) {
    printf("In resolution shells:\n"
           "  d_max  d_min   #obs  #uniq  #used  Rmerge   Rmeas   Rpim    CC1/2\n");
    for (size_t i = 0; i < binner->size(); ++i) {
      const MergingStats& ms = stats[i];
      printf("%7.3f %5.3f %7d %6d %6d %7.3f %7.3f %7.3f %8.4f\n",
             binner->dmax_of_bin(i), binner->dmin_of_bin(i),
             ms.all_refl, ms.unique_refl, ms.stats_refl,
             ms.r_merge(), ms.r_meas(), ms.r_pim(), ms.cc_half());
    }
  }

  MergingStats total;
  for (const MergingStats& ms : stats)
    total.add_other(ms);
  printf("\nObservations (all reflections): %d\n", total.all_refl);
  printf("Unique reflections: %d\n", total.unique_refl);
  printf("Used refl. (those with multiplicity 2+): %d\n", total.stats_refl);
  printf("          Overall");
  if (binner)
    printf("    Avg of %zu shells weighted by #used", binner->size());
  printf("\nR-merge: %7.4f", total.r_merge());
  if (binner)
    printf("         %7.4f", take_average(stats, &MergingStats::r_merge));
  printf("\nR-meas:  %7.4f", total.r_meas());
  if (binner)
    printf("         %7.4f", take_average(stats, &MergingStats::r_meas));
  printf("\nR-pim:   %7.4f", total.r_pim());
  if (binner)
    printf("         %7.4f", take_average(stats, &MergingStats::r_pim));
  // xdscc12 prints CC1/2 with 5 significant digits, do the same for easy comparison
  printf("\nCC1/2:   %8.5f", total.cc_half());
  if (binner)
    printf("        %8.5f", take_average(stats, &MergingStats::cc_half));
  printf("\n");
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.check_exclusive_pair(Stats, Compare);
  int n_args = p.nonOptionsCount();
  if (n_args != 2 && !(n_args == 1 && (p.options[Stats] || p.options[Compare]))) {
    fprintf(stderr, "%s requires 2 arguments or 1 arg with --stats/--compare, got %d.\n",
                    p.program_name, p.nonOptionsCount());
    p.print_try_help_and_exit("");
  }
  bool verbose = p.options[Verbose];
  bool batch_cell = p.options[BatchCell];
  bool two_files = p.nonOptionsCount() == 2;
  std::string input_path = p.coordinate_input_file(0, 'S');
  const char* output_path = two_files ? p.nonOption(1) : input_path.c_str();
  const char* input_block = p.options[InputBlock].arg;  // nullptr if option not given
  const char* output_block = p.options[OutputBlock].arg;
  bool to_anom = p.options[WriteAnom];

  if (p.options[Compare] && !two_files && gemmi::giends_with(input_path, "hkl")) {
    fprintf(stderr, "ERROR. Option --compare doesn't work with one XDS file.\n");
    return 1;
  }

  Intensities intensities;
  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path.c_str());
  read_intensities(intensities, DataType::UAM, input_path, input_block, verbose, batch_cell);
  if (intensities.type != DataType::Unmerged)
    std::fprintf(stderr, "NOTE: Got merged %s instead of unmerged data.\n",
                 intensities.type_str());
  if (p.options[NoSysAbs]) {
    size_t before = intensities.data.size();
    intensities.remove_systematic_absences();
    size_t diff = intensities.data.size() - before;
    if (verbose && diff != 0)
      std::fprintf(stderr, "Removed systematic absences: %zu observations.\n", diff);
  }
  if (verbose)
    output_intensity_statistics(intensities);

  auto binning_method = gemmi::Binner::Method::Dstar3;
  if (p.options[Stats]) {
    int nbins = 10;
    char use_weights = 'Y';
    if (const char* arg = p.options[Stats].arg) {
      char* endptr;
      long n = std::strtol(arg, &endptr, 10);
      if (endptr != arg) {
        nbins = (int) n;
        if (nbins <= 0) {
          fprintf(stderr, "Wrong argument for option --stats (expected N>0): %s\n", arg);
          return 1;
        }
      }
      for (; *endptr != '\0'; ++endptr) {
        char c = *endptr;
        if (c == 'U' || c == 'X') {
          use_weights = c;
        } else if (c == 'e') {
          binning_method = gemmi::Binner::Method::EqualCount;
        } else if (c == 's') {
          binning_method = gemmi::Binner::Method::Dstar2;
        } else {
          fprintf(stderr, "Unexpected character '%c' in option --stats: %s\n", c, arg);
          return 1;
        }
      }
    }
    if (verbose) {
      if (use_weights == 'Y')
        printf("Using weights for R-metrics like Aimless\n");
      else if (use_weights == 'X')
        printf("Using weights for R-metrics like XDS\n");
      else if (use_weights == 'U')
        printf("NOT using weights for R-metrics.\n");
    }
    try {
      if (to_anom)
        intensities.prepare_for_merging(DataType::Anomalous);
      else
        intensities.sort();
      gemmi::Binner binner;
      const gemmi::Binner* binner_ptr = nullptr;
      if (nbins > 1) {
        if (verbose) {
          const char* binning_desc = "equal volumes (equispaced in d*^3)";
          if (binning_method == gemmi::Binner::Method::Dstar2)
            binning_desc = "increasing volumes (equispaced in d*^2)";
          else if (binning_method == gemmi::Binner::Method::EqualCount)
            binning_desc = "equal reflection count";
          printf("Setting up resolution shells with %s.\n", binning_desc);
        }
        binner.setup(nbins, binning_method, gemmi::IntensitiesDataProxy{intensities});
        binner_ptr = &binner;
      }
      auto stats = intensities.calculate_merging_stats(binner_ptr, use_weights);
      print_merging_statistics(stats, binner_ptr);
    } catch (std::exception& e) {
      std::fprintf(stderr, "ERROR: %s\n", e.what());
      return 1;
    }
  }

  if (p.options[Compare]) {
    if (verbose)
      std::fprintf(stderr, "Reading merged reflections from %s ...\n", output_path);
    auto type_we_want = to_anom ? DataType::Anomalous : DataType::MergedMA;
    Intensities ref;  // for --compare
    read_intensities(ref, type_we_want, output_path, output_block, verbose, false);
    if (!to_anom && ref.type == DataType::Anomalous)
      std::fprintf(stderr, "Using I(+)/I(-) because <I> is absent.\n");
    if (intensities.type != ref.type)
      intensities.merge_in_place(ref.type);
    printf("Comparing %s ...%s\n", intensities.type_str(),
           verbose ? "" : "   (use -v for more details)");
    compare_intensities(intensities, ref, p.options[PrintAll]);
  }

  if (two_files) {
    intensities.merge_in_place(to_anom ? DataType::Anomalous : DataType::Mean);
    if (verbose)
      std::fprintf(stderr, "Writing %zu reflections to %s ...\n",
                   intensities.data.size(), output_path);
    bool with_nobs = p.options[NumObs];
    write_merged_intensities(intensities.prepare_merged_mtz(with_nobs), output_path);
  }
  return 0;
}
