// Copyright 2017 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/symmetry.hpp"
#include "input.h"
#include <cstdlib>  // for strtod

#define GEMMI_PROG mask
#include "options.h"

enum OptionIndex { Verbose=3, FormatIn, Threshold, Fraction, GridSpac,
                   GridDims, Radius };

struct MaskArg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"coor", "map", "none"});
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
    "  --from=coor|map  \tInput type (default: from file extension)." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for making a mask from a map:" },
  { Threshold, 0, "t", "threshold", Arg::Float,
    "  -t, --threshold  \tThe density cutoff value." },
  { Fraction, 0, "f", "fraction", Arg::Float,
    "  -f, --fraction  \tThe volume fraction to be above the threshold." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for masking a model:" },
  { GridSpac, 0, "s", "spacing", Arg::Float,
    "  -s, --spacing=D  \tMax. sampling for the grid (default: 1A)." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid sampling." },
  { Radius, 0, "r", "radius", Arg::Float,
    "  -r, --radius  \tRadius of atom spheres (default: 3.0A)." },
  { 0, 0, 0, 0, 0, 0 }
};

enum class InputType : char { Coordinates, Ccp4, Unknown };


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.exclusive_groups.push_back({Threshold, Fraction});
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);

  if (p.options[Verbose])
    std::fprintf(stderr, "Converting %s ...\n", input);

  InputType in_type = InputType::Unknown;
  if (p.options[FormatIn]) {
    const char* arg = p.options[FormatIn].arg;
    if (strcmp(arg, "coor") == 0)
      in_type = InputType::Coordinates;
    else if (strcmp(arg, "map") == 0)
      in_type = InputType::Ccp4;
  }
  if (in_type == InputType::Unknown) {
    if (coordinate_format_from_extension(input) != gemmi::CoorFormat::Unknown) {
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
  if (p.options[GridDims] && p.options[GridSpac]) {
    std::fprintf(stderr, "Options --grid and --spacing are exclusive.");
    return 1;
  }

  try {
    // map -> mask
    if (in_type == InputType::Ccp4) {
      double threshold;
      gemmi::Ccp4<signed char> mask;
      mask.read_ccp4_map(input);
      if (p.options[Threshold]) {
        threshold = std::strtod(p.options[Threshold].arg, nullptr);
      } else if (p.options[Fraction]) {
        double fraction = std::strtod(p.options[Fraction].arg, nullptr);
        if (fraction < 0) {
          std::fprintf(stderr, "Cannot use negative fraction.\n");
          return 2;
        }
        auto data = mask.grid.data;  // making a copy for nth_element()
        size_t n = std::min(static_cast<size_t>(data.size() * fraction),
                            data.size() - 1);
        std::nth_element(data.begin(), data.begin() + n, data.end());
        threshold = data[n];
      } else {
        std::fprintf(stderr, "You need to specify threshold (-t or -f).\n");
        return 2;
      }
      mask.grid.make_zeros_and_ones(threshold);
      size_t ones = std::count(mask.grid.data.begin(), mask.grid.data.end(), 1);
      size_t all = mask.grid.data.size();
      std::fprintf(stderr, "Masked %zu of %zu points (%.1f%%) above %g\n",
                   ones, all, 100.0 * ones / all, threshold);
      mask.update_ccp4_header(0);
      mask.write_ccp4_map(output);

    // model -> mask
    } else {
      double radius = (p.options[Radius]
                       ? std::strtod(p.options[Radius].arg, nullptr)
                       : 3.0);
      gemmi::Structure st = read_structure(input);
      gemmi::Ccp4<signed char> mask;
      mask.grid.unit_cell = st.cell;
      mask.grid.space_group = gemmi::find_spacegroup_by_name(st.sg_hm);
      if (p.options[GridDims]) {
        auto dims = parse_comma_separated_ints(p.options[GridDims].arg);
        mask.grid.set_size(dims[0], dims[1], dims[2]);
      } else {
        double spac = 1;
        if (p.options[GridSpac])
          spac = std::strtod(p.options[GridSpac].arg, nullptr);
        mask.grid.set_size_from_spacing(spac, true);
      }
      if (p.options[Verbose]) {
        const auto& g = mask.grid;
        std::fprintf(stderr, "Grid: %d x %d x %d\n", g.nu, g.nv, g.nw);
        std::fprintf(stderr, "Spacing: %.3f, %.3f, %.3f\n",
                             g.spacing[0], g.spacing[1], g.spacing[2]);
        int np = g.data.size();
        double vol = st.cell.volume;
        std::fprintf(stderr, "Total points: %d\n", np);
        std::fprintf(stderr, "Unit cell volume: %.1f A^3\n", vol);
        std::fprintf(stderr, "Volume per point: %.3f A^3\n", vol / np);
        if (g.space_group) {
          std::fprintf(stderr, "Spacegroup: %s\n", g.space_group->hm);
          int na = g.space_group->operations().order();
          std::fprintf(stderr, "ASU volume: %.1f A^3\n", vol / na);
          std::fprintf(stderr, "Points per ASU: %d\n", np / na);
        } else {
          std::fprintf(stderr, "No spacegroup\n");
        }
      }
      if (st.models.size() > 1)
        std::fprintf(stderr, "Note: only the first model is used.\n");
      for (const gemmi::Chain& chain : st.models[0].chains)
        for (const gemmi::Residue& res : chain.residues)
          for (const gemmi::Atom& atom : res.atoms)
            mask.grid.set_points_around(atom.pos, radius, 1);
      if (p.options[Verbose]) {
        int n = std::count(mask.grid.data.begin(), mask.grid.data.end(), 1);
        std::fprintf(stderr, "Points masked by model: %d\n", n);
      }
      mask.grid.symmetrize(
          [](signed char a, signed char b) { return std::max(a,b); });
      if (p.options[Verbose]) {
        int n = std::count(mask.grid.data.begin(), mask.grid.data.end(), 1);
        std::fprintf(stderr, "After symmetrizing: %d\n", n);
      }
      mask.update_ccp4_header(0, true);
      mask.write_ccp4_map(output);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
