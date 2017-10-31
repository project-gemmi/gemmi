// Copyright 2017 Global Phasing Ltd.

#include "gemmi/grid.hpp"
#include "gemmi/crystal.hpp"
#include "input.h"
#include <cstdlib>  // for strtod
#define EXE_NAME "gemmi-mask"
#include "options.h"

enum OptionIndex { Verbose=3, FormatIn, Threshold, Fraction, GridDims, Radius };

struct MaskArg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"ccp4", "pdb", "cif", "none"});
  }
};

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT output.msk"
    "\n\nMakes a mask in the CCP4 format."
    "\nIf INPUT is a CCP4 map the mask is created by thresholding the map."
    "\nIf INPUT is a coordinate file (mmCIF, PDB, etc) the atoms are masked." },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { FormatIn, 0, "", "from", MaskArg::FileFormat,
    "  --from=ccp4|pdb|cif  \tInput format (default: from file extension)." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for making a mask from a map:" },
  { Threshold, 0, "t", "threshold", Arg::Float,
    "  -t, --threshold  \tThe density cutoff value." },
  { Fraction, 0, "f", "fraction", Arg::Float,
    "  -f, --fraction  \tThe volume fraction to be above the threshold." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for masking a model:" },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid sampling (default: ~1A spacing)." },
  { Radius, 0, "r", "radius", Arg::Float,
    "  -r, --radius  \tRadius of atom spheres (default: 3.0A)." },
  { 0, 0, 0, 0, 0, 0 }
};

enum class InputType : char { Coordinates, Ccp4, Unknown };


int main(int argc, char **argv) {
  OptParser p;
  p.exclusive_groups.push_back({Threshold, Fraction});
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);

  if (p.options[Verbose])
    std::fprintf(stderr, "Converting %s ...\n", input);

  InputType in_type = InputType::Unknown;
  gemmi::CoorFormat in_subtype = gemmi::CoorFormat::Unknown;
  if (p.options[FormatIn]) {
    const char* arg = p.options[FormatIn].arg;
    if (strcmp(arg, "pdb") == 0) {
      in_type = InputType::Coordinates;
      in_subtype = gemmi::CoorFormat::Pdb;
    } else if (strcmp(arg, "cif") == 0) {
      in_type = InputType::Coordinates;
      in_subtype = gemmi::CoorFormat::Cif;
    } else if (strcmp(arg, "ccp4") == 0) {
      in_type = InputType::Ccp4;
    }
  }
  if (in_type == InputType::Unknown) {
    in_subtype = coordinate_format_from_extension(input);
    if (in_subtype != gemmi::CoorFormat::Unknown) {
      in_type = InputType::Coordinates;
    } else if (gemmi::iends_with(input, ".ccp4") ||
               gemmi::iends_with(input, ".map")) {
      in_type = InputType::Ccp4;
    }
  }
  if (in_type == InputType::Unknown) {
    std::fprintf(stderr, "Cannot determine input type for extension."
                         " Use --from=...");
    return 1;
  }

  try {
    // map -> mask
    if (in_type == InputType::Ccp4) {
      double threshold;
      gemmi::Grid<> grid;
      grid.read_ccp4(input);
      if (p.options[Threshold]) {
        threshold = std::strtod(p.options[Threshold].arg, nullptr);
      } else if (p.options[Fraction]) {
        double fraction = std::strtod(p.options[Fraction].arg, nullptr);
        if (fraction < 0) {
          std::fprintf(stderr, "Cannot use negative fraction.\n");
          return 2;
        }
        auto data = grid.data;  // making a copy for nth_element()
        size_t n = std::min(static_cast<size_t>(data.size() * fraction),
                            data.size() - 1);
        std::nth_element(data.begin(), data.begin() + n, data.end());
        threshold = data[n];
      } else {
        std::fprintf(stderr, "You need to specify threshold (-t or -f).\n");
        return 2;
      }
      size_t count = grid.write_ccp4_mask(output, threshold);
      std::fprintf(stderr, "Masked %zu of %zu points (%.1f%%) above %g\n",
                   count, grid.data.size(), 100.0 * count / grid.data.size(),
                   threshold);

    // model -> mask
    } else {
      double radius = (p.options[Radius]
                       ? std::strtod(p.options[Radius].arg, nullptr)
                       : 3.0);

      gemmi::Structure st = read_structure(input);
      gemmi::Grid<> grid;
      grid.unit_cell = st.cell;
      if (p.options[GridDims]) {
        auto dims = parse_comma_separated_ints(p.options[GridDims].arg);
        grid.set_size(dims[0], dims[1], dims[2]);
      } else {
        grid.set_spacing(1);
      }
      if (st.models.size() > 1)
        std::fprintf(stderr, "Note: only the first model is used.\n");
      for (const gemmi::Chain& chain : st.models[0].chains)
        for (const gemmi::Residue& res : chain.residues)
          for (const gemmi::Atom& atom : res.atoms)
            grid.set_points_around(atom.pos, radius, 1.0);
      grid.stats = grid.calculate_statistics();
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
