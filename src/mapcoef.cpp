// Copyright 2019 Global Phasing Ltd.
//
// Common options for programs that read MTZ or SF-mmCIF map coefficients
// and transform them to a map.

#include <stdio.h>
#include <cstring>            // for strcmp
#include <cstdlib>            // for strtod, exit
#include <array>
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for get_f_phi_on_grid, transform_f_phi_..
#include <gemmi/gzread.hpp>   // for read_cif_gz
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/util.hpp>     // for fail, giends_with
#include "mapcoef.h"

#define GEMMI_PROG n/a
#include "options.h"

using gemmi::Mtz;

const option::Descriptor MapUsage[] = {
  { 0, 0, 0, 0, 0, 0 }, // The first 3 entries are empty to make MapUsage[enum]
  { 0, 0, 0, 0, 0, 0 }, // work as expected.
  { 0, 0, 0, 0, 0, 0 }, // NoOp=0, Help=1, Version=2.
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Diff, 0, "d", "diff", Arg::None,
    "  -d, --diff  \tUse difference map coefficients." },
  { Section, 0, "", "section", Arg::Required,
    "  --section=NAME  \tMTZ dataset name or CIF block name" },
  { FLabel, 0, "f", "", Arg::Required,
    "  -f COLUMN  \tF column (MTZ label or mmCIF tag)." },
  { PhLabel, 0, "p", "", Arg::Required,
    "  -p COLUMN  \tPhase column (MTZ label or mmCIF tag)." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid size (user-specified minimum)." },
  { Sample, 0, "s", "sample", Arg::Float,
    "  -s, --sample=NUMBER  \tSet spacing to d_min/NUMBER (3 is usual)." },
  { GridQuery, 0, "G", "", Arg::None,
    "  -G  \tPrint size of the grid that would be used and exit." },
};


static std::array<const Mtz::Column*, 2>
get_mtz_map_columns(const Mtz& mtz, const char* section, bool diff_map,
                    const char* f_label, const char* phi_label) {
  static const char* default_labels[] = {
    "FWT", "PHWT", "DELFWT", "PHDELWT",
    "2FOFCWT", "PH2FOFCWT", "FOFCWT", "PHFOFCWT",
    nullptr
  };
  const Mtz::Column* f_col = nullptr;
  const Mtz::Column* phi_col = nullptr;
  const Mtz::Dataset* ds = nullptr;
  if (section) {
    ds = mtz.dataset_with_name(section);
    if (!ds)
      gemmi::fail(std::string("No such dataset in the MTZ file: ") + section);
  }
  if (f_label) {
    f_col = mtz.column_with_label(f_label, ds);
    if (!f_col)
      gemmi::fail(std::string("Column not found: ") + f_label);
    if (phi_label) {
      phi_col = mtz.column_with_label(phi_label, ds);
    } else {
      for (int i = 0; ; i += 2) {
        if (default_labels[i] == nullptr)
          gemmi::fail("Unknown phase column label.\n");
        if (std::strcmp(default_labels[i], f_label) == 0) {
          phi_col = mtz.column_with_label(default_labels[i + 1], ds);
          break;
        }
      }
    }
  } else {
    for (int i = (diff_map ? 2 : 0); default_labels[i]; i += 4)
      if ((f_col = mtz.column_with_label(default_labels[i], ds)) &&
          (phi_col = mtz.column_with_label(default_labels[i+1], ds))) {
        break;
      }
  }
  if (!f_col || !phi_col)
    gemmi::fail("Default map coefficient labels not found.\n");
  return {{f_col, phi_col}};
}

template<typename DataProxy>
void print_grid_size(const DataProxy& data, std::array<int, 3> min_size,
                     double rate) {
  std::array<int, 3> size = gemmi::get_size_for_hkl(data, min_size, rate);
  printf("Grid size: %d x %d x %d\n", size[0], size[1], size[2]);
}


gemmi::Grid<float>
read_sf_and_fft_to_map(const char* input_path,
                       const std::vector<option::Option>& options,
                       FILE* output,
                       bool oversample_by_default) {
  if (options[PhLabel] && !options[FLabel])
    gemmi::fail("Option -p can be given only together with -f");
  if (options[FLabel] && options[Diff])
    gemmi::fail("Option -d has no effect together with -f");
  if (output)
    fprintf(output, "Reading reflections from %s ...\n", input_path);
  std::vector<int> vsize{0, 0, 0};
  if (options[GridDims])
    vsize = parse_comma_separated_ints(options[GridDims].arg);
  std::array<int,3> min_size = {{vsize[0], vsize[1], vsize[2]}};
  double sample_rate = 0.;
  if (options[Sample])
    sample_rate = std::strtod(options[Sample].arg, nullptr);
  else if (oversample_by_default && !options[GridDims])
    sample_rate = 3.;
  const char* section = options[Section] ? options[Section].arg : nullptr;
  const char* f_label = options[FLabel] ? options[FLabel].arg : nullptr;
  const char* ph_label = options[PhLabel] ? options[PhLabel].arg : nullptr;
  bool diff_map = options[Diff];
  gemmi::Grid<std::complex<float>> grid;
  if (gemmi::giends_with(input_path, ".cif") ||
      gemmi::giends_with(input_path, ".ent")) {
    if (!f_label)
      f_label = diff_map ? "pdbx_DELFWT" : "pdbx_FWT";
    if (!ph_label)
      ph_label = diff_map ?  "pdbx_DELPHWT" : "pdbx_PHWT";
    if (output)
      fprintf(output, "Looking for tags _refln.%s and _refln.%s...\n",
              f_label, ph_label);
    gemmi::ReflnBlock rblock = gemmi::get_refln_block(
                                   gemmi::read_cif_gz(input_path).blocks,
                                   {f_label, ph_label}, section);
    if (options[GridQuery]) {
      print_grid_size(gemmi::ReflnDataProxy{rblock}, min_size, sample_rate);
      std::exit(0);
    }
    if (output)
      fprintf(output, "Putting data from block %s into matrix...\n",
              rblock.block.name.c_str());
    grid = gemmi::get_f_phi_on_grid<float>(gemmi::ReflnDataProxy{rblock},
                                           rblock.find_column_index(f_label),
                                           rblock.find_column_index(ph_label),
                                           /*half_l=*/true,
                                           min_size, sample_rate);
  } else {
    Mtz mtz = gemmi::read_mtz(gemmi::MaybeGzipped(input_path), true);
    auto cols = get_mtz_map_columns(mtz, section, diff_map, f_label, ph_label);
    if (output)
      fprintf(output, "Putting data from columns %s and %s into matrix...\n",
              cols[0]->label.c_str(), cols[1]->label.c_str());
    if (options[GridQuery]) {
      print_grid_size(gemmi::MtzDataProxy{mtz}, min_size, sample_rate);
      std::exit(0);
    }
    grid = gemmi::get_f_phi_on_grid<float>(gemmi::MtzDataProxy{mtz},
                                           cols[0]->idx, cols[1]->idx, true,
                                           min_size, sample_rate);
  }
  if (output)
    fprintf(output, "Fourier transform -> grid %d x %d x %d...\n",
            grid.nu, grid.nv, (grid.nw - 1) * 2);
  return gemmi::transform_f_phi_grid_to_map(std::move(grid));
}

