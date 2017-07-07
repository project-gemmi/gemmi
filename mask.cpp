// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include <cstdlib>  // for strtod

#define EXE_NAME "gemmi-mask"
#include "options.h"

enum OptionIndex { Unknown, Verbose, Threshold, FormatIn };

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
  { Threshold, 0, "t", "threshold", Arg::Float,
    "  -t, --threshold  \tMask map below the threshold." },
  { FormatIn, 0, "", "from", MaskArg::FileFormat,
    "  --from=ccp4|pdb|cif  \tInput format (default: from file extension)." },
  { 0, 0, 0, 0, 0, 0 }
};

void map_to_mask(const char* input, const char* output, double threshold) {
  gemmi::Grid grid;
  grid.read_ccp4(input);
  grid.write_ccp4_mask(output, threshold);
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
    double threshold = std::strtod(options[Threshold].arg, nullptr);
    map_to_mask(input, output, threshold);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=include,third_party
