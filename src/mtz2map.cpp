// Copyright 2019 Global Phasing Ltd.
//
// MTZ (or SF-mmCIF) map coefficients -> CCP4 map.

#include <stdio.h>
#include <cstring>            // for strcmp
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/ccp4.hpp>     // for Ccp4
#include <gemmi/fourier.hpp>  // for get_f_phi_on_grid, transform_f_phi_...
#include <gemmi/gzread.hpp>   // for read_cif_gz
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/util.hpp>     // for fail, giends_with
//#include <gemmi/version.hpp>  // for GEMMI_VERSION

#define GEMMI_PROG mtz2map
#include "options.h"

using gemmi::Mtz;
using options_type = std::vector<option::Option>;

enum OptionIndex { Verbose=3, Diff, FLabel, PhiLabel, GridDims };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE MAP_FILE\n"
    "INPUT_FILE must be either MTZ or mmCIF with map coefficients."
    "\nOptions:"},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Diff, 0, "d", "diff", Arg::None,
    "  -d, --diff  \tUse data from columns (PH)DELFWT or (PH)FOFCWT." },
  { FLabel, 0, "f", "", Arg::Required,
    "  -fCOLUMN  \tF label in the MTZ file." },
  { PhiLabel, 0, "p", "", Arg::Required,
    "  -pCOLUMN  \tPhase label in the MTZ file." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid sampling." },
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
  if (options[FLabel]) {
    f_col = mtz.column_with_label(options[FLabel].arg);
    if (!f_col)
      gemmi::fail("Column not found: " + std::string(options[FLabel].arg));
    if (options[PhiLabel]) {
      phi_col = mtz.column_with_label(options[PhiLabel].arg);
    } else {
      for (int i = 0; ; i += 2) {
        if (default_labels[i] == nullptr)
          gemmi::fail("Unknown phase column label.\n");
        if (std::strcmp(default_labels[i], options[FLabel].arg) == 0) {
          phi_col = mtz.column_with_label(default_labels[i + 1]);
          break;
        }
      }
    }
  } else {
    for (int i = (options[Diff] ? 2 : 0); default_labels[i]; i += 4)
      if ((f_col = mtz.column_with_label(default_labels[i])) &&
          (phi_col = mtz.column_with_label(default_labels[i+1]))) {
        break;
      }
  }
  if (!f_col || !phi_col)
    gemmi::fail("Default map coefficient labels not found.\n");
  return {{f_col, phi_col}};
}

static gemmi::Grid<std::complex<float>>
grid_from_mmcif(const char* input_path, std::array<int,3> size,
                const char* f_label, const char* ph_label, bool verbose) {
  auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(input_path).blocks);
  if (verbose)
    fprintf(stderr, "Looking for tags _refln.%s and _refln.%s...\n",
            f_label, ph_label);
  for (gemmi::ReflnBlock& rb : rblocks) {
    int f_idx = rb.find_column_index(f_label);
    int ph_idx = rb.find_column_index(ph_label);
    if (f_idx != -1 && ph_idx != -1) {
      if (verbose)
        fprintf(stderr, "Putting data from block %s into matrix...\n",
                rb.block.name.c_str());
      return gemmi::get_f_phi_on_grid<float>(gemmi::ReflnDataProxy{rb},
                                             f_idx, ph_idx,
                                             /*half_l*/true, size);
    }
  }
  gemmi::fail("Corresponding tags not found in the mmCIF file.");
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* map_path = p.nonOption(1);
  if (p.options[PhiLabel] && !p.options[FLabel]) {
    fprintf(stderr, "Option -p can be given only together with -f\n");
    return 1;
  }
  if (p.options[FLabel] && p.options[Diff])
    fprintf(stderr, "Option -d has no effect together with -f\n");
  if (verbose)
    fprintf(stderr, "Reading %s ...\n", input_path);
  std::vector<int> vsize{0, 0, 0};
  if (p.options[GridDims])
    vsize = parse_comma_separated_ints(p.options[GridDims].arg);
  std::array<int,3> size = {{vsize[0], vsize[1], vsize[2]}};
  try {
    gemmi::Grid<std::complex<float>> grid;
    if (gemmi::giends_with(input_path, ".cif") ||
        gemmi::giends_with(input_path, ".ent")) {
      const char* f_label = p.options[Diff] ? "pdbx_DELFWT" : "pdbx_FWT";
      const char* ph_label = p.options[Diff] ?  "pdbx_DELPHWT" : "pdbx_PHWT";
      if (p.options[FLabel])
        f_label = p.options[FLabel].arg;
      if (p.options[PhiLabel])
        ph_label = p.options[PhiLabel].arg;
      grid = grid_from_mmcif(input_path, size, f_label, ph_label, verbose);
    } else {
      Mtz mtz = gemmi::read_mtz_file(input_path);
      auto cols = get_mtz_columns(mtz, p.options);
      if (p.options[Verbose])
        fprintf(stderr, "Putting data from columns %s and %s into matrix...\n",
                cols[0]->label.c_str(), cols[1]->label.c_str());
      grid = gemmi::get_f_phi_on_grid<float>(gemmi::MtzDataProxy{mtz},
                                             cols[0]->idx, cols[1]->idx,
                                             /*half_l*/true, size);
    }
    if (verbose)
      fprintf(stderr, "Fourier transform -> grid %d x %d x %d...\n",
              grid.nu, grid.nv, (grid.nw - 1) * 2);
    gemmi::Ccp4<float> ccp4;
    ccp4.grid = gemmi::transform_f_phi_half_to_map(std::move(grid));

    if (verbose)
      fprintf(stderr, "Writing %s ...\n", map_path);
    ccp4.update_ccp4_header(2, true);
    ccp4.write_ccp4_map(map_path);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
