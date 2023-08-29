// Copyright 2019 Global Phasing Ltd.
//
// Common options for programs that read MTZ or SF-mmCIF map coefficients
// and transform them to a map.

#include "mapcoef.h"
#include <stdio.h>
#include <cstring>            // for strcmp
#include <cstdlib>            // for strtod, exit
#include <array>
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for get_f_phi_on_grid, transform_f_phi_..
#include <gemmi/recgrid.hpp>  // for ReciprocalGrid
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/util.hpp>     // for fail, giends_with
#include <gemmi/read_cif.hpp> // for read_cif_gz
#include "timer.h"

#define GEMMI_PROG n/a
#include "options.h"

using gemmi::Mtz;

const option::Descriptor MapUsage[] = {
  { 0, 0, 0, 0, 0, 0 }, // The first 3 entries are empty to make MapUsage[enum]
  { 0, 0, 0, 0, 0, 0 }, // work as expected.
  { 0, 0, 0, 0, 0, 0 }, // NoOp=0, Help=1, Version=2, Verbose=3.
  { 0, 0, 0, 0, 0, 0 },
  { Diff, 0, "d", "diff", Arg::None,
    "  -d, --diff  \tUse difference map coefficients." },
  { Section, 0, "", "section", Arg::Required,
    "  --section=NAME  \tMTZ dataset name or CIF block name" },
  { FLabel, 0, "f", "", Arg::Required,
    "  -f COLUMN  \tF column (MTZ label or mmCIF tag)." },
  { PhLabel, 0, "p", "", Arg::Required,
    "  -p COLUMN  \tPhase column (MTZ label or mmCIF tag)." },
  { WeightLabel, 0, "", "weight", Arg::Required,
    "  --weight=COLUMN  \t(normally not needed) weighting for F." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid size (user-specified minimum)." },
  { ExactDims, 0, "", "exact", Arg::None,
    "  --exact  \tUse the exact grid size specified by --grid." },
  { Sample, 0, "s", "sample", Arg::Float,
    "  -s, --sample=NUMBER  \tSet spacing to d_min/NUMBER (3 is usual)." },
  { AxesZyx, 0, "", "zyx", Arg::None,
    "  --zyx  \tOutput axes Z Y X as fast, medium, slow (default is X Y Z)." },
  { GridQuery, 0, "G", "", Arg::None,
    "  -G  \tPrint size of the grid that would be used and exit." },
  { TimingFft, 0, "", "timing", Arg::None,
    "  --timing  \tPrint calculation times." },
};


static std::array<const Mtz::Column*, 2>
get_mtz_map_columns(const Mtz& mtz, const char* section, bool diff_map,
                    const char* f_label, const char* phi_label) {
  const Mtz::Column* f_col = nullptr;
  const Mtz::Column* phi_col = nullptr;
  const Mtz::Dataset* ds = nullptr;
  if (section) {
    ds = mtz.dataset_with_name(section);
    if (!ds)
      gemmi::fail("No such dataset in the MTZ file: ", section);
  }
  if (f_label) {
    f_col = mtz.column_with_label(f_label, ds);
    if (!f_col)
      gemmi::fail("Column not found: ", f_label);
    if (phi_label)
      phi_col = mtz.column_with_label(phi_label, ds);
    else
      phi_col = f_col->get_next_column_if_type('P');
    if (!f_col || !phi_col)
      gemmi::fail("Specified map coefficient labels not found.\n");
  } else if (diff_map) {
    for (const char* label : {"DELFWT", "FOFCWT"})
      if ((f_col = mtz.column_with_label(label, ds)))
        if ((phi_col = f_col->get_next_column_if_type('P')))
          break;
    if (!f_col || !phi_col)
      gemmi::fail("Default difference map labels (DELFWT/FOFCWT + phase) not found.\n");
  } else {
    for (const char* label : {"FWT", "2FOFCWT"})
      if ((f_col = mtz.column_with_label(label, ds)))
        if ((phi_col = f_col->get_next_column_if_type('P')))
          break;
    if (!f_col || !phi_col)
      gemmi::fail("Default map labels (FWT/2FOFCWT + phase) not found.\n");
  }
  return {{f_col, phi_col}};
}

static const Mtz::Column&
get_mtz_column(const Mtz& mtz, const char* section, const char* label) {
  const Mtz::Dataset* ds = nullptr;
  if (section)
    ds = mtz.dataset_with_name(section);
  const Mtz::Column* col = mtz.column_with_label(label, ds);
  if (!col)
    gemmi::fail("Column not found: ", label);
  return *col;
}

template<typename DataProxy>
void adjust_size(const DataProxy& data, std::array<int, 3>& size,
                 double sample_rate, bool exact_dims, bool grid_query) {
  if (exact_dims) {
    if (!gemmi::data_fits_into(data, size))
      gemmi::fail("grid size is too small for hkl data");
    gemmi::check_grid_factors(data.spacegroup(), size);
  } else {
    size = gemmi::get_size_for_hkl(data, size, sample_rate);
  }
  if (grid_query) {
    printf("Grid size: %d x %d x %d\n", size[0], size[1], size[2]);
    std::exit(0);
  }
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
  if (options[ExactDims] && !options[GridDims])
    gemmi::fail("Option --exact requires option --grid");
  if (options[ExactDims] && options[Sample])
    gemmi::fail("Option --sample has not effect together with --exact");
  if (output)
    fprintf(output, "Reading reflections from %s ...\n", input_path);
  std::vector<int> vsize{0, 0, 0};
  if (options[GridDims])
    vsize = parse_comma_separated_ints(options[GridDims].arg);
  Timer timer(options[TimingFft]);
  std::array<int,3> size = {{vsize[0], vsize[1], vsize[2]}};
  double sample_rate = 0.;
  if (options[Sample])
    sample_rate = std::strtod(options[Sample].arg, nullptr);
  else if (oversample_by_default && !options[GridDims])
    sample_rate = 3.;
  const char* section = options[Section] ? options[Section].arg : nullptr;
  const char* f_label = options[FLabel] ? options[FLabel].arg : nullptr;
  const char* ph_label = options[PhLabel] ? options[PhLabel].arg : nullptr;
  const char* weight_label = options[WeightLabel] ? options[WeightLabel].arg
                                                  : nullptr;
  bool diff_map = options[Diff];
  bool half_l = true;
  gemmi::AxisOrder axis_order = options[AxesZyx] ? gemmi::AxisOrder::ZYX
                                                 : gemmi::AxisOrder::XYZ;
  gemmi::FPhiGrid<float> grid;
  gemmi::ReciprocalGrid<float> weight_grid;
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
    gemmi::ReflnDataProxy data(rblock);
    adjust_size(data, size, sample_rate,
                options[ExactDims], options[GridQuery]);
    if (output)
      fprintf(output, "Putting data from block %s into matrix...\n",
              rblock.block.name.c_str());
    gemmi::ReflnDataProxy data_proxy(rblock);
    int f_col = rblock.find_column_index(f_label);
    int phi_col = rblock.find_column_index(ph_label);
    gemmi::FPhiProxy<gemmi::ReflnDataProxy> fphi(data_proxy, f_col, phi_col);
    grid = gemmi::get_f_phi_on_grid<float>(fphi, size, half_l, axis_order);
    if (weight_label)
      weight_grid = gemmi::get_value_on_grid<float>(
          data_proxy, rblock.find_column_index(weight_label),
          size, half_l, axis_order);
  } else {
    timer.start();
    Mtz mtz;
    mtz.read_file_gz(input_path);
    timer.print("MTZ read in");
    auto cols = get_mtz_map_columns(mtz, section, diff_map, f_label, ph_label);
    gemmi::MtzDataProxy data_proxy{mtz};
    adjust_size(data_proxy, size, sample_rate,
                options[ExactDims], options[GridQuery]);
    if (output)
      fprintf(output, "Putting data from columns %s and %s into matrix...\n",
              cols[0]->label.c_str(), cols[1]->label.c_str());
    timer.start();
    gemmi::FPhiProxy<gemmi::MtzDataProxy> fphi(data_proxy, cols[0]->idx, cols[1]->idx);
    grid = gemmi::get_f_phi_on_grid<float>(fphi, size, half_l, axis_order);
    timer.print("F/Phi grid prepared in");
    if (weight_label) {
      const Mtz::Column& col = get_mtz_column(mtz, section, weight_label);
      weight_grid = gemmi::get_value_on_grid<float>(data_proxy, col.idx,
                                                    size, half_l, axis_order);
    }
  }
  if (weight_grid.data.size() == grid.data.size())
    for (size_t i = 0; i != grid.data.size(); ++i)
      grid.data[i] *= weight_grid.data[i];
  if (output)
    fprintf(output, "Fourier transform...\n");
  timer.start();
  gemmi::Grid<float> map = gemmi::transform_f_phi_grid_to_map(std::move(grid));
  timer.print("FFT in");
  assert(map.axis_order == axis_order);
  if (output)
    fprintf(output, "Map size: %d x %d x %d\n", map.nu, map.nv, map.nw);
  return map;
}
