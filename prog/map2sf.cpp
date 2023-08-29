// Copyright 2019 Global Phasing Ltd.
//
// Transform CCP4 map to map coefficient columns in MTZ file.

#include <stdio.h>
#include <cctype>             // for toupper
#include <cstdlib>            // for strtod
#include <gemmi/fail.hpp>     // for fail
#include <gemmi/grid.hpp>     // for Grid, ReciprocalGrid, ReciprocalGrid<>...
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/ccp4.hpp>     // for Ccp4
#include <gemmi/fourier.hpp>  // for transform_map_to_f_phi
#include <gemmi/util.hpp>     // for iends_with
#include <gemmi/gz.hpp>       // for MaybeGzipped

#define GEMMI_PROG map2sf
#include "options.h"

namespace {

using gemmi::Mtz;

enum OptionIndex { Base=4, Section, DMin, FType, PhiType, Spacegroup };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] MAP_FILE OUTPUT_FILE COL_F COL_PH\n\n"
    "Writes map coefficients (amplitude and phase) of a map to OUTPUT_FILE.\n"
    "The output is MTZ if it has mtz extension, otherwise it is mmCIF.\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Base, 0, "b", "base", Arg::Required,
    "  -b, --base=PATH  \tAdd new columns to the data from this file." },
  { Section, 0, "", "section", Arg::Required,
    "  --section=NAME  \tAdd new columns to this MTZ dataset or CIF block." },
  { DMin, 0, "", "dmin", Arg::Float,
    "  --dmin=D_MIN  \tResolution limit." },
  { FType, 0, "", "ftype", Arg::Char,
    "  --ftype=TYPE   \tMTZ amplitude column type (default: F)." },
  { PhiType, 0, "", "phitype", Arg::Char,
    "  --phitype=TYPE  \tMTZ phase column type (default: P)." },
  { Spacegroup, 0, "", "spacegroup", Arg::Required,
    "  --spacegroup=SG  \tOverwrite space group from map header." },
  { 0, 0, 0, 0, 0, 0 }
};

void transform_map_to_sf(OptParser& p) {
  bool verbose = p.options[Verbose];
  const char* map_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);
  const char* f_col = p.nonOption(2);
  const char* phi_col = p.nonOption(3);
  char f_type = p.options[FType] ? std::toupper(p.options[FType].arg[0]) : 'F';
  char phi_type = p.options[PhiType] ? std::toupper(p.options[PhiType].arg[0]) : 'P';
  if (verbose)
    fprintf(stderr, "Reading %s ...\n", map_path);
  gemmi::Ccp4<float> map;
  map.read_ccp4(gemmi::MaybeGzipped(map_path));
  if (p.options[Spacegroup]) {
    map.grid.spacegroup = gemmi::find_spacegroup_by_name(p.options[Spacegroup].arg);
    if (map.grid.spacegroup == nullptr)
      gemmi::fail("unknown space group: ", p.options[Spacegroup].arg);
  }
  map.setup(NAN);
  if (std::any_of(map.grid.data.begin(), map.grid.data.end(),
                  [](float x) { return std::isnan(x); }))
    gemmi::fail("Map does not cover all the ASU");
  if (verbose)
    fprintf(stderr, "Fourier transform of grid %d x %d x %d...\n",
            map.grid.nu, map.grid.nv, map.grid.nw);
  gemmi::FPhiGrid<float> hkl = gemmi::transform_map_to_f_phi(map.grid, /*half_l=*/true);
  if (gemmi::giends_with(output_path, ".mtz")) {
    gemmi::Mtz mtz;
    if (p.options[Base]) {
      if (verbose)
        fprintf(stderr, "Reading %s ...\n", p.options[Base].arg);
      mtz.read_file_gz(p.options[Base].arg);
      int dataset_id = -1;
      if (p.options[Section]) {
        const char* ds_name = p.options[Section].arg;
        if (Mtz::Dataset* ds = mtz.dataset_with_name(ds_name))
          dataset_id = ds->id;
        else
          mtz.add_dataset(ds_name);
      }
      if (verbose)
        fprintf(stderr, "Copied columns: %s\n",
                gemmi::join_str(mtz.columns, ", ",
                                [](const Mtz::Column& c) { return c.label; }
                               ).c_str());
      mtz.add_column(f_col, f_type, dataset_id, -1, false);
      size_t f_idx = mtz.columns.back().idx;
      mtz.add_column(phi_col, phi_type, dataset_id, -1, false);
      mtz.expand_data_rows(2);
      for (int i = 0; i != mtz.nreflections; ++i) {
        size_t offset = i * mtz.columns.size();
        int h = (int) mtz.data[offset + 0];
        int k = (int) mtz.data[offset + 1];
        int l = (int) mtz.data[offset + 2];
        std::complex<float> v = hkl.get_value(h, k, l);
        mtz.data[offset + f_idx] = (float) std::abs(v);
        mtz.data[offset + f_idx + 1] = (float) gemmi::phase_in_angles(v);
      }
    } else {
      double dmin = 0;
      if (p.options[DMin])
        dmin = std::strtod(p.options[DMin].arg, nullptr);
      mtz.cell = map.grid.unit_cell;
      mtz.spacegroup = map.grid.spacegroup;
      mtz.sort_order = {{1, 2, 3, 0, 0}};
      mtz.add_base();
      mtz.add_dataset(p.options[Section] ? p.options[Section].arg : "unknown");
      mtz.add_column(f_col, f_type, -1, -1, false);
      mtz.add_column(phi_col, phi_type, -1, -1, false);
      gemmi::AsuData<std::complex<float>> data = hkl.prepare_asu_data<>(dmin);
      mtz.nreflections = (int) data.v.size();
      add_asu_f_phi_to_float_vector(mtz.data, data);
    }
    if (verbose)
      fprintf(stderr, "Writing %s ...\n", output_path);
    mtz.write_to_file(output_path);
  } else {
    gemmi::fail("mmCIF support not implemented yet");
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(4);
  try {
    transform_map_to_sf(p);
  } catch (std::exception& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
