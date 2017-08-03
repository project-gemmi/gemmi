// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "input.h"
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
    "\n\nINPUT is either a CCP4 map or a coordinate file." },
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

enum class InputType : char { Pdb, Mmcif, Ccp4, Unknown };


int main(int argc, char **argv) {
  OptParser parse;
  auto options = parse.simple_parse(argc, argv, Usage);
  parse.require_positional_args(2);
  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  if (options[Verbose])
    std::fprintf(stderr, "Converting %s ...\n", input);

  InputType in_type = InputType::Unknown;
  if (options[FormatIn]) {
    const char* arg = options[FormatIn].arg;
    if (strcmp(arg, "pdb") == 0)
      in_type = InputType::Pdb;
    else if (strcmp(arg, "cif") == 0)
      in_type = InputType::Mmcif;
    else if (strcmp(arg, "ccp4") == 0)
      in_type = InputType::Ccp4;
  }
  if (in_type == InputType::Unknown) {
    using gemmi::iends_with;
    if (iends_with(input, ".pdb") || iends_with(input, ".ent") ||
        iends_with(input, ".pdb.gz") || iends_with(input, ".ent.gz"))
      in_type = InputType::Pdb;
    else if (iends_with(input, ".cif") || iends_with(input, ".cif.gz"))
      in_type = InputType::Mmcif;
    else if (iends_with(input, ".ccp4") || iends_with(input, ".map"))
      in_type = InputType::Ccp4;
  }
  if (in_type == InputType::Unknown) {
    std::fprintf(stderr, "Cannot determine input type for extension."
                         " Use --from=...");
    return 1;
  }

  try {
    // map -> mask
    if (in_type == InputType::Ccp4) {
      if (!options[Threshold]) {
        std::fprintf(stderr, "You need to specify threshold (-t).\n");
        return 2;
      }
      double threshold = std::strtod(options[Threshold].arg, nullptr);
      gemmi::Grid<> grid;
      grid.read_ccp4(input);
      grid.write_ccp4_mask(output, threshold);

    // model -> mask
    } else {
      double radius = (options[Radius]
                       ? std::strtod(options[Radius].arg, nullptr)
                       : 3.0);

      mol::Structure st;
      if (in_type == InputType::Pdb)
        st = pdb_read_any(input);
      else if (in_type == InputType::Mmcif)
        st = mmcif_read_atoms(cif_read_any(input));
      gemmi::Grid<> grid;
      grid.unit_cell = st.cell;
      if (options[GridDims]) {
        auto dims = parse_comma_separated_ints(options[GridDims].arg);
        grid.set_size(dims[0], dims[1], dims[2]);
      } else {
        grid.set_spacing(1);
      }
      if (st.models.size() > 1)
        std::fprintf(stderr, "Note: only the first model is used.\n");
      for (const mol::Chain& chain : st.models[0].chains)
        for (const mol::Residue& res : chain.residues)
          for (const mol::Atom& atom : res.atoms)
            grid.set_points_around(atom.pos, radius, 1.0);
      grid.calculate_statistics();
      //grid.write_ccp4_mask(output, 0.5);
      grid.write_ccp4_map(output);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
