// Copyright 2019 Global Phasing Ltd.
//
// Transform MTZ or SF-mmCIF map coefficients to CCP4 map.

#include <stdio.h>
#include <gemmi/ccp4.hpp>      // for Ccp4
#include <gemmi/calculate.hpp> // for calculate_fractional_box
#include <gemmi/select.hpp>    // for Selection
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include "mapcoef.h"

#define GEMMI_PROG sf2map
#include "options.h"

namespace {

enum OptionIndex { Normalize=AfterMapOptions, MapMask, Margin, Select };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE MAP_FILE\n\n"
    "INPUT_FILE must be either MTZ or mmCIF with map coefficients.\n\n"
    "By default, the program searches for 2mFo-DFc map coefficients in:\n"
    "  - MTZ columns FWT/PHWT or 2FOFCWT/PH2FOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_FWT/pdbx_PHWT.\n"
    "If option \"-d\" is given, mFo-DFc map coefficients are searched in:\n"
    "  - MTZ columns DELFWT/PHDELWT or FOFCWT/PHFOFCWT,\n"
    "  - mmCIF tags _refln.pdbx_DELFWT/pdbx_DELPHWT.\n\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  MapUsage[Diff],
  MapUsage[Section],
  MapUsage[FLabel],
  MapUsage[PhLabel],
  MapUsage[WeightLabel],
  MapUsage[GridDims],
  MapUsage[ExactDims],
  MapUsage[Sample],
  MapUsage[AxesZyx],
  MapUsage[GridQuery],
  MapUsage[TimingFft],
  { Normalize, 0, "", "normalize", Arg::None,
    "  --normalize  \tScale the map to standard deviation 1 and mean 0." },
  { MapMask, 0, "", "mapmask", Arg::Required,
    "  --mapmask=FILE  \tOutput only map covering the structure from FILE,"
    " similarly to CCP4 MAPMASK with XYZIN." },
  { Margin, 0, "", "margin", Arg::Float,
    "  --margin=N  \t(w/ --mapmask) Border in Angstrom (default: 5)." },
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL  \t(w/ --mapmask) Selection of atoms in FILE, MMDB syntax." },
  { 0, 0, 0, 0, 0, 0 }
};


void transform_sf_to_map(OptParser& p) {
  const char* input_path = p.nonOption(0);
  const char* map_path = p.options[GridQuery] ? nullptr : p.nonOption(1);
  gemmi::Ccp4<float> ccp4;
  ccp4.grid = read_sf_and_fft_to_map(input_path, p.options,
                                     p.options[Verbose] ? stderr : nullptr);
  if (p.options[Verbose])
    fprintf(stderr, "Writing %s ...\n", map_path);
  ccp4.update_ccp4_header(2);
  if (p.options[Normalize]) {
    double mult = 1.0 / ccp4.hstats.rms;
    for (float& x : ccp4.grid.data)
      x = float((x - ccp4.hstats.dmean) * mult);
    ccp4.update_ccp4_header(2);
  }
  if (p.options[MapMask]) {
    double margin = 5;
    if (p.options[Margin])
      margin = std::atof(p.options[Margin].arg);
    gemmi::Structure st = gemmi::read_structure_gz(p.options[MapMask].arg);
    st.cell = ccp4.grid.unit_cell; // needed in case the unit cells differ
    gemmi::Box<gemmi::Fractional> box;
    if (p.options[Select]) {
      gemmi::Selection sel(p.options[Select].arg);
      for (gemmi::Model& model : sel.models(st))
        for (gemmi::Chain& chain : sel.chains(model))
          for (gemmi::Residue& res : sel.residues(chain))
            for (gemmi::Atom& atom : sel.atoms(res))
              box.extend(st.cell.fractionalize(atom.pos));
      box.add_margins({margin * st.cell.ar, margin * st.cell.br, margin * st.cell.cr});
    } else {
      box = gemmi::calculate_fractional_box(st, margin);
    }
    ccp4.set_extent(box);
  }
  ccp4.write_ccp4_map(map_path);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[GridQuery])
    p.require_input_files_as_args(0);
  else
    p.require_positional_args(2);
  try {
    transform_sf_to_map(p);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
