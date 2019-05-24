// Copyright 2019 Global Phasing Ltd.
//
// Transform MTZ or SF-mmCIF map coefficients to CCP4 map.

#include <stdio.h>
#include <cstring>            // for strcmp
#include <cstdlib>            // for strtod
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/ccp4.hpp>     // for Ccp4
#include <gemmi/fourier.hpp>  // for get_f_phi_on_grid, transform_f_phi_...
#include <gemmi/gzread.hpp>   // for read_cif_gz
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/util.hpp>     // for fail, giends_with
//#include <gemmi/version.hpp>  // for GEMMI_VERSION

#define GEMMI_PROG sf2map
#include "options.h"

using gemmi::Mtz;
using options_type = std::vector<option::Option>;

enum OptionIndex { Verbose=3, Diff, Section, FLabel, PhiLabel, GridDims,
                   Sample, GridQuery, Normalize };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE MAP_FILE\n\n"
    "INPUT_FILE must be either MTZ or mmCIF with map coefficients.\n\n"
    "By default, the program searches for 2mFo-DFc map coefficients in:\n"
    "  - MTZ columns FWT/PHWT or 2FOFCWT/PH2FOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_FWT/pdbx_PHWT.\n"
    "If option \"-d\" is given, mFo-DFc map coefficients are searched in:\n"
    "  - MTZ columns DELFWT/PHDELWT or FOFCWT/PHFOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_DELFWT/pdbx_DELPHWT\n\n"
    "\nOptions:"},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Diff, 0, "d", "diff", Arg::None,
    "  -d, --diff  \tUse difference map coefficients." },
  { Section, 0, "", "section", Arg::Required,
    "  --section=NAME  \tMTZ dataset name or CIF block name" },
  { FLabel, 0, "f", "", Arg::Required,
    "  -f COLUMN  \tF column (MTZ label or mmCIF tag)." },
  { PhiLabel, 0, "p", "", Arg::Required,
    "  -p COLUMN  \tPhase column (MTZ label or mmCIF tag)." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid size (user-specified minimum)." },
  { Sample, 0, "s", "sample", Arg::Float,
    "  -s, --sample=NUMBER  \tSet spacing to d_min/NUMBER (3 is usual)." },
  { GridQuery, 0, "G", "", Arg::None,
    "  -G  \tPrint size of the grid that would be used and exit." },
  { Normalize, 0, "", "normalize", Arg::None,
    "  --normalize  \tScale the map to standard deviation 1 and mean 0." },
  { 0, 0, 0, 0, 0, 0 }
};


static std::array<const Mtz::Column*, 2>
get_mtz_columns(const Mtz& mtz, const options_type& options) {
  static const char* default_labels[] = {
    "FWT", "PHWT", "DELFWT", "PHDELWT",
    "2FOFCWT", "PH2FOFCWT", "FOFCWT", "PHFOFCWT",
    nullptr
  };
  const Mtz::Column* f_col = nullptr;
  const Mtz::Column* phi_col = nullptr;
  const Mtz::Dataset* ds = nullptr;
  if (options[Section]) {
    std::string dname = options[Section].arg;
    ds = mtz.dataset_with_name(dname);
    if (!ds)
      gemmi::fail("No such dataset in the MTZ file: " + dname);
  }
  if (options[FLabel]) {
    f_col = mtz.column_with_label(options[FLabel].arg, ds);
    if (!f_col)
      gemmi::fail("Column not found: " + std::string(options[FLabel].arg));
    if (options[PhiLabel]) {
      phi_col = mtz.column_with_label(options[PhiLabel].arg, ds);
    } else {
      for (int i = 0; ; i += 2) {
        if (default_labels[i] == nullptr)
          gemmi::fail("Unknown phase column label.\n");
        if (std::strcmp(default_labels[i], options[FLabel].arg) == 0) {
          phi_col = mtz.column_with_label(default_labels[i + 1], ds);
          break;
        }
      }
    }
  } else {
    for (int i = (options[Diff] ? 2 : 0); default_labels[i]; i += 4)
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

static void transform_sf_to_map(OptParser& p) {
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* map_path = p.nonOption(1);
  if (p.options[PhiLabel] && !p.options[FLabel])
    gemmi::fail("Option -p can be given only together with -f");
  if (p.options[FLabel] && p.options[Diff])
    gemmi::fail("Option -d has no effect together with -f");
  if (verbose)
    fprintf(stderr, "Reading %s ...\n", input_path);
  std::vector<int> vsize{0, 0, 0};
  if (p.options[GridDims])
    vsize = parse_comma_separated_ints(p.options[GridDims].arg);
  std::array<int,3> min_size = {{vsize[0], vsize[1], vsize[2]}};
  double sample_rate = 0;
  if (p.options[Sample])
    sample_rate = std::strtod(p.options[Sample].arg, nullptr);
  gemmi::Grid<std::complex<float>> grid;
  if (gemmi::giends_with(input_path, ".cif") ||
      gemmi::giends_with(input_path, ".ent")) {
    const char* f_label = p.options[Diff] ? "pdbx_DELFWT" : "pdbx_FWT";
    const char* ph_label = p.options[Diff] ?  "pdbx_DELPHWT" : "pdbx_PHWT";
    const char* block_name = p.options[Section] ? p.options[Section].arg
                                                : nullptr;
    if (p.options[FLabel])
      f_label = p.options[FLabel].arg;
    if (p.options[PhiLabel])
      ph_label = p.options[PhiLabel].arg;
    std::vector<gemmi::ReflnBlock> rblocks = gemmi::as_refln_blocks(
                                        gemmi::read_cif_gz(input_path).blocks);
    if (verbose)
      fprintf(stderr, "Looking for tags _refln.%s and _refln.%s...\n",
              f_label, ph_label);
    bool has_section = false;
    for (gemmi::ReflnBlock& rb : rblocks) {
      if (block_name && rb.block.name != block_name)
        continue;
      has_section = true;
      int f_idx = rb.find_column_index(f_label);
      int ph_idx = rb.find_column_index(ph_label);
      if (f_idx != -1 && ph_idx != -1) {
        if (verbose)
          fprintf(stderr, "Putting data from block %s into matrix...\n",
                  rb.block.name.c_str());
        if (p.options[GridQuery]) {
          print_grid_size(gemmi::ReflnDataProxy{rb}, min_size, sample_rate);
          return;
        }
        grid = gemmi::get_f_phi_on_grid<float>(gemmi::ReflnDataProxy{rb},
                                               f_idx, ph_idx, /*half_l=*/true,
                                               min_size, sample_rate);
        break;
      }
    }
    if (!has_section)
      gemmi::fail("Section not found in the mmCIF file.");
    if (grid.data.empty())
      gemmi::fail("Required tags not found in the mmCIF file.");
  } else {
    Mtz mtz = gemmi::read_mtz_file(input_path);
    auto cols = get_mtz_columns(mtz, p.options);
    if (p.options[Verbose])
      fprintf(stderr, "Putting data from columns %s and %s into matrix...\n",
              cols[0]->label.c_str(), cols[1]->label.c_str());
    if (p.options[GridQuery]) {
      print_grid_size(gemmi::MtzDataProxy{mtz}, min_size, sample_rate);
      return;
    }
    grid = gemmi::get_f_phi_on_grid<float>(gemmi::MtzDataProxy{mtz},
                                           cols[0]->idx, cols[1]->idx, true,
                                           min_size, sample_rate);
  }

  if (verbose)
    fprintf(stderr, "Fourier transform -> grid %d x %d x %d...\n",
            grid.nu, grid.nv, (grid.nw - 1) * 2);
  gemmi::Ccp4<float> ccp4;
  ccp4.grid = gemmi::transform_f_phi_grid_to_map(std::move(grid));
  if (verbose)
    fprintf(stderr, "Writing %s ...\n", map_path);
  ccp4.update_ccp4_header(2, true);
  if (p.options[Normalize]) {
    double mult = 1.0 / ccp4.hstats.rms;
    for (float& x : ccp4.grid.data)
      x = float((x - ccp4.hstats.dmean) * mult);
  }
  ccp4.write_ccp4_map(map_path);
}


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  try {
    transform_sf_to_map(p);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
