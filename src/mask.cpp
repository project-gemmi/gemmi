// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "gemmi/pdb.hpp"
#include <cstdlib>  // for strtod

#define EXE_NAME "gemmi-mask"
#include "options.h"

namespace mol = gemmi::mol;

enum OptionIndex { Unknown, Verbose, FormatIn, Threshold, GridDims, Radius};

struct MaskArg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"ccp4", "pdb", "cif", "none"});
  }
};

static const option::Descriptor Usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT output.msk"
    "\n\nINPUT is either or a CCP4 map or (not yet) a coordinate file." },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { FormatIn, 0, "", "from", MaskArg::FileFormat,
    "  --from=ccp4|pdb|cif  \tInput format (default: from file extension)." },
  { Unknown, 0, "", "", Arg::None, "\nOptions for map to mask conversions:" },
  { Threshold, 0, "t", "threshold", Arg::Float,
    "  -t, --threshold  \tMask map below the threshold." },
  { Unknown, 0, "", "", Arg::None, "\nOptions for model masking:" },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid sampling (default: ~1A spacing)." },
  { Radius, 0, "r", "radius", Arg::Float,
    "  -r, --radius  \tRadius of atom spheres (default: 3.0A)." },
  { 0, 0, 0, 0, 0, 0 }
};


int main(int argc, char **argv) {
  OptParser parse;
  auto options = parse.simple_parse(argc, argv, Usage);
  parse.require_positional_args(2);
  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  if (options[Verbose])
    fprintf(stderr, "Converting %s ...\n", input);

  gemmi::Grid<> grid;
  double threshold = 0.5;
  try {
    if (gemmi::iends_with(input, ".pdb")) {
      double radius = (options[Radius]
                       ? std::strtod(options[Radius].arg, nullptr)
                       : 3.0);
      mol::Structure st = mol::read_pdb(input);
      grid.unit_cell = st.cell;
      if (options[GridDims]) {
        auto dims = parse_comma_separated_ints(options[GridDims].arg);
        grid.set_size(dims[0], dims[1], dims[2]);
      } else {
        grid.set_spacing(1);
      }
      if (st.models.size() > 1)
        fprintf(stderr, "Note: only the first model is used.\n");
      for (const mol::Chain& chain : st.models[0].chains)
        for (const mol::Residue& res : chain.residues)
          for (const mol::Atom& atom : res.atoms)
            grid.set_points_around(atom.pos, radius, 1.0);
      grid.calculate_statistics();
    } else {
      if (!options[Threshold]) {
        fprintf(stderr, "You need to specify threshold (-t).\n");
        return 2;
      }
      threshold = std::strtod(options[Threshold].arg, nullptr);
      grid.read_ccp4(input);
    }
    //grid.write_ccp4_mask(output, threshold);
    grid.write_ccp4_map(output);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
