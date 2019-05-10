// Copyright 2019 Global Phasing Ltd.
//
// MTZ map coefficients -> CCP4 map.

#include <stdio.h>
#include <cstring>              // for strcmp
#include <gemmi/mtz.hpp>
#include <gemmi/ccp4.hpp>
#include <gemmi/fourier.hpp>
//#include <gemmi/fileutil.hpp> // for file_open
//#include <gemmi/atox.hpp>     // for read_word
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#include <pocketfft/pocketfft_hdronly.h>

#define GEMMI_PROG mtz2map
#include "options.h"

using gemmi::Mtz;

enum OptionIndex { Verbose=3, Diff, FLabel, PhiLabel, GridDims };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] MTZ_FILE MAP_FILE"
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
get_mtz_columns(const Mtz& mtz, const std::vector<option::Option>& options) {
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

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* mtz_path = p.nonOption(0);
  const char* map_path = p.nonOption(1);
  if (p.options[PhiLabel] && !p.options[FLabel]) {
    fprintf(stderr, "Option -p can be given only together with -f\n");
    return 1;
  }
  if (p.options[FLabel] && p.options[Diff])
    fprintf(stderr, "Option -d has no effect together with -f\n");
  if (verbose)
    fprintf(stderr, "Reading %s ...\n", mtz_path);
  try {
    Mtz mtz = gemmi::read_mtz_file(mtz_path);
    auto cols = get_mtz_columns(mtz, p.options);
    if (verbose)
      fprintf(stderr, "Putting data from columns %s and %s into matrix...\n",
              cols[0]->label.c_str(), cols[1]->label.c_str());
    std::vector<int> size{0, 0, 0};
    if (p.options[GridDims])
      size = parse_comma_separated_ints(p.options[GridDims].arg);
    bool half_l = true;
    auto grid = gemmi::get_f_phi_on_grid<>(mtz, cols[0]->idx, cols[1]->idx,
                                           half_l, {size[0], size[1], size[2]});
    if (verbose)
      fprintf(stderr, "Fourier transform...\n");
    // x -> conj(x) is equivalent to changing axis direction before FFT
    for (std::complex<float>& x : grid.data)
      x.imag(-x.imag());
    gemmi::Ccp4<float> ccp4;
    ccp4.grid.space_group = mtz.spacegroup;
    ccp4.grid.unit_cell = mtz.cell;
    int full_nw = 2 * (grid.nw - 1);
    ccp4.grid.set_size(grid.nu, grid.nv, full_nw);
    float norm = float(1.0 / ccp4.grid.unit_cell.volume);
    pocketfft::shape_t shape{(size_t)grid.nw, (size_t)grid.nv, (size_t)grid.nu};
    ptrdiff_t s = sizeof(float);
    pocketfft::stride_t stride{grid.nv * grid.nu * 2*s, grid.nu * 2*s, 2*s};
    pocketfft::c2c<float>(shape, stride, stride, {1, 2}, pocketfft::BACKWARD,
                          &grid.data[0], &grid.data[0], 1.0f);
    pocketfft::stride_t stride_out{grid.nv * grid.nu * s, grid.nu * s, s};
    shape[0] = full_nw;
    pocketfft::c2r<float>(shape, stride, stride_out, 0,
                          &grid.data[0], &ccp4.grid.data[0], norm);
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
