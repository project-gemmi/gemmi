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

struct Intensity {
  gemmi::Miller hkl;
  int isign = 0;  // 1 for I(+), -1 for I(-)
  double value;
  double sigma;

  bool operator<(const Intensity& o) const {
    return std::tie(hkl[0], hkl[1], hkl[2], isign) <
           std::tie(o.hkl[0], o.hkl[1], o.hkl[2], o.isign);
  }
};

struct Intensities {
  std::vector<Intensity> data;
  const gemmi::SpaceGroup* spacegroup = nullptr;
  gemmi::UnitCell unit_cell;
  double wavelength;

  bool have_sign() const { return !data.empty() && data[0].isign != 0; }

  std::array<double,2> resolution_range() const {
    double min_1_d2 = INFINITY;
    double max_1_d2 = 0;
    for (const Intensity& x : data) {
      double a_1_d2 = unit_cell.calculate_1_d2(x.hkl);
      if (a_1_d2 < min_1_d2)
        min_1_d2 = a_1_d2;
      if (a_1_d2 > max_1_d2)
        max_1_d2 = a_1_d2;
    }
    return { 1 / std::sqrt(min_1_d2), 1 / std::sqrt(max_1_d2) };
  }
  void copy_metadata(const gemmi::Mtz& mtz, int dataset_id) {
    unit_cell = mtz.cell;
    spacegroup = mtz.spacegroup;
    if (!spacegroup)
      gemmi::fail("unknown space group");
    wavelength = mtz.dataset(dataset_id).wavelength;
  }

  void copy_metadata(const gemmi::ReflnBlock& rb) {
    unit_cell = rb.cell;
    spacegroup = rb.spacegroup;
    if (!spacegroup)
      gemmi::fail("unknown space group");
    wavelength = rb.wavelength;
  }

  void add_if_valid(const Intensity& intensity) {
      // XDS marks rejected reflections with negative sigma.
      // Sigma 0.0 is also problematic - it rarely happens (e.g. 5tkn).
    if (!std::isnan(intensity.value) && intensity.sigma > 0)
      data.push_back(intensity);
  }

  template<typename DataProxy>
  void read_data(const DataProxy& proxy, size_t value_idx, size_t sigma_idx) {
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      Intensity intensity;
      intensity.hkl = proxy.get_hkl(i);
      intensity.value = proxy.get_num(i + value_idx);
      intensity.sigma = proxy.get_num(i + sigma_idx);
      add_if_valid(intensity);
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

enum class ReadType { Unmerged, Mean, Anomalous };

Intensities read_unmerged_intensities_from_mtz(const gemmi::Mtz& mtz) {
  if (mtz.batches.empty())
    gemmi::fail("expected unmerged file");
  const gemmi::Mtz::Column* isym_col = mtz.column_with_label("M/ISYM");
  if (!isym_col || isym_col->idx != 3)
    gemmi::fail("unmerged file should have M/ISYM as 4th column");
  const gemmi::Mtz::Column& col = mtz.get_column_with_label("I");
  size_t value_idx = col.idx;
  size_t sigma_idx = mtz.get_column_with_label("SIGI").idx;
  Intensities intensities;
  intensities.copy_metadata(mtz, col.dataset_id);
  for (size_t i = 0; i < mtz.data.size(); i += mtz.columns.size()) {
    Intensity intensity;
    intensity.hkl = mtz.get_hkl(i);
    intensity.isign = ((int)mtz.data[i + 3] % 2 == 0 ? -1 : 1);
    intensity.value = mtz.data[i + value_idx];
    intensity.sigma = mtz.data[i + sigma_idx];
    intensities.add_if_valid(intensity);
  }
  return intensities;
}

Intensities read_mean_intensities_from_mtz(const gemmi::Mtz& mtz) {
  if (!mtz.batches.empty())
    gemmi::fail("expected merged file");
  const gemmi::Mtz::Column& col = mtz.get_column_with_label("IMEAN");
  size_t value_idx = col.idx;
  size_t sigma_idx = mtz.get_column_with_label("SIGIMEAN").idx;
  Intensities intensities;
  intensities.copy_metadata(mtz, col.dataset_id);
  intensities.read_data(gemmi::MtzDataProxy{mtz}, value_idx, sigma_idx);
  return intensities;
}

Intensities read_anomalous_intensities_from_mtz(const gemmi::Mtz& mtz) {
  if (!mtz.batches.empty())
    gemmi::fail("expected merged file");
  const gemmi::Mtz::Column& col = mtz.get_column_with_label("I(+)");
  size_t value_idx[2] = {col.idx, mtz.get_column_with_label("I(-)").idx};
  size_t sigma_idx[2] = {mtz.get_column_with_label("SIGI(+)").idx,
                         mtz.get_column_with_label("SIGI(-)").idx};
  Intensities intensities;
  intensities.copy_metadata(mtz, col.dataset_id);
  for (size_t i = 0; i < mtz.data.size(); i += mtz.columns.size())
    for (int j = 0; j < 2; ++j) {
      Intensity intensity;
      intensity.hkl = mtz.get_hkl(i);
      intensity.isign = (j == 0 ? 1 : -1);
      intensity.value = mtz.data[i + value_idx[j]];
      intensity.sigma = mtz.data[i + sigma_idx[j]];
      intensities.add_if_valid(intensity);
    }
  return intensities;
}

Intensities read_unmerged_intensities_from_mmcif(const gemmi::ReflnBlock& rb) {
  size_t value_idx = rb.get_column_index("intensity_net");
  size_t sigma_idx = rb.get_column_index("intensity_sigma");
  Intensities intensities;
  intensities.copy_metadata(rb);
  intensities.read_data(gemmi::ReflnDataProxy(rb), value_idx, sigma_idx);
  // switch to ASU indices
  gemmi::GroupOps gops = intensities.spacegroup->operations();
  gemmi::ReciprocalAsu asu(intensities.spacegroup);
  for (Intensity& intensity : intensities.data) {
    auto hkl_isym = asu.to_asu(intensity.hkl, gops);
    intensity.hkl = hkl_isym.first;
    intensity.isign = (hkl_isym.second % 2 == 0 ? -1 : 1);
  }
  return intensities;
}

Intensities read_mean_intensities_from_mmcif(const gemmi::ReflnBlock& rb) {
  size_t value_idx = rb.get_column_index("intensity_meas");
  size_t sigma_idx = rb.get_column_index("intensity_sigma");
  Intensities intensities;
  intensities.copy_metadata(rb);
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
  gemmi::ReflnDataProxy proxy(rb);
  for (size_t i = 0; i < proxy.size(); i += proxy.stride())
    for (int j = 0; j < 2; ++j) {
      Intensity intensity;
      intensity.hkl = proxy.get_hkl(i);
      intensity.isign = (j == 0 ? 1 : -1);
      intensity.value = proxy.get_num(i + value_idx[j]);
      intensity.sigma = proxy.get_num(i + sigma_idx[j]);
      intensities.add_if_valid(intensity);
    }
  return intensities;
}

Intensities read_intensities(ReadType itype, const char* input_path,
                             const char* block_name, bool verbose) {
  try {
    if (gemmi::giends_with(input_path, ".mtz")) {
      gemmi::Mtz mtz;
      if (verbose)
        mtz.warnings = stderr;
      mtz.read_input(gemmi::MaybeGzipped(input_path), /*with_data=*/true);
      if (itype == ReadType::Mean && !mtz.column_with_label("IMEAN")) {
        if (!mtz.column_with_label("I(+)"))
          gemmi::fail("merged intensities not found");
        std::fprintf(stderr, "No IMEAN, using I(+) and I(-) ...\n");
        itype = ReadType::Anomalous;
      }
      switch (itype) {
        case ReadType::Unmerged:
          return read_unmerged_intensities_from_mtz(mtz);
        case ReadType::Mean:
          return read_mean_intensities_from_mtz(mtz);
        case ReadType::Anomalous:
          return read_anomalous_intensities_from_mtz(mtz);
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

void merge_intensities_in_place(Intensities& intensities) {
  intensities.sort();
  std::vector<Intensity>::iterator out = intensities.data.begin();
  double sum_wI = 0.;
  double sum_w = 0.;
  for (auto in = intensities.data.begin(); in != intensities.data.end(); ++in) {
    if (out->hkl != in->hkl || out->isign != in->isign) {
      out->value = sum_wI / sum_w;
      out->sigma = 1.0 / std::sqrt(sum_w);
      sum_wI = sum_w = 0.;
      ++out;
      if (out != in) {
        out->hkl = in->hkl;
        out->isign = in->isign;
      }
    }
    double w = 1. / (in->sigma * in->sigma);
    sum_wI += w * in->value;
    sum_w += w;
  }
  out->value = sum_wI / sum_w;
  out->sigma = 1.0 / std::sqrt(sum_w);
  intensities.data.erase(++out, intensities.data.end());
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
  for (const Intensity& inten : intensities.data) {
    if (inten.hkl != prev_hkl) {
      offset += mtz.columns.size();
      mtz.set_hkl(offset, inten.hkl);
      prev_hkl = inten.hkl;
    }
    size_t value_offset = offset + (inten.isign >= 0 ? 3 : 5);
    mtz.data[value_offset] = (float) inten.value;
    mtz.data[value_offset + 1] = (float) inten.sigma;
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
      mtz_to_cif.write_cif(mtz, os.ref());
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR while writing %s: %s\n", output_path, e.what());
    std::exit(1);
  }
}

void print_reflection(const Intensity* a, const Intensity* r) {
  const Intensity& h = a ? *a : *r;
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
  for (const Intensity& r : ref.data) {
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
    ref = read_intensities(itype, input_path, nullptr, verbose);
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
    for (const Intensity& inten : intensities.data) {
      if (inten.isign > 0)
        ++plus_count;
      else if (inten.isign < 0)
        ++minus_count;
    }
    std::fprintf(stderr, "Merging observations (%zu total, %zu for I+, %zu for I-) ...\n",
                 intensities.data.size(), plus_count, minus_count);
  }
  if (p.options[Compare] ? !ref.have_sign() : !p.options[WriteAnom])
    // discard signs so that merging produces Imean
    for (Intensity& inten : intensities.data)
      inten.isign = 0;
  merge_intensities_in_place(intensities);
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
