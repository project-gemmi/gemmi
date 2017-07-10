// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "gemmi/pdb.hpp"
#include <cstdlib>  // for strtod

#define EXE_NAME "gemmi-mask"
#include "options.h"

enum OptionIndex { Unknown, Verbose, FormatIn, Threshold, Radius};

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
  { Radius, 0, "r", "radius", Arg::Float,
    "  -r, --radius  \tRadius of atom spheres (default: 3.0A)." },
  { 0, 0, 0, 0, 0, 0 }
};

void map_to_mask(const char* input, const char* output, double threshold) {
  gemmi::Grid<> grid;
  grid.read_ccp4(input);
  grid.write_ccp4_mask(output, threshold);
}

void mol_to_mask(const char* input, const char* output, double radius) {
  gemmi::mol::Structure model = gemmi::mol::read_pdb(input);
  gemmi::Grid<> grid;
  grid.set_size(128, 128, 128);
  grid.unit_cell = model.cell;
  // TODO
  grid.write_ccp4_mask(output, 0.5);
}


int main(int argc, char **argv) {
  OptParser parse;
  auto options = parse.simple_parse(argc, argv, Usage);
  parse.require_positional_args(2);
  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  if (options[Verbose])
    fprintf(stderr, "Converting %s ...\n", input);

  try {
    if (gemmi::iends_with(input, ".pdb")) {
      double radius = 3.0;
      if (options[Radius])
        radius = std::strtod(options[Radius].arg, nullptr);
      mol_to_mask(input, output, radius);
    } else {
      if (!options[Threshold]) {
        fprintf(stderr, "You need to specify threshold (-t).\n");
        return 2;
      }
      double threshold = std::strtod(options[Threshold].arg, nullptr);
      map_to_mask(input, output, threshold);
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=include,third_party
