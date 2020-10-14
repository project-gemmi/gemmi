// Copyright 2017 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <cstdlib>  // for atof
#include "gemmi/ccp4.hpp"
#include "gemmi/floodfill.hpp" // for remove_islands_in_mask
#include "gemmi/gz.hpp"        // for MaybeGzipped
#include "gemmi/gzread.hpp"
#include "gemmi/rhogrid.hpp"   // for mask_points_in_constant_radius, etc
#include "gemmi/symmetry.hpp"
#include "timer.h"

#define GEMMI_PROG mask
#include "options.h"

namespace {

enum OptionIndex { Timing=4, GridSpac, GridDims, Radius, RProbe, RShrink,
                   IslandLimit, CctbxCompat, RefmacCompat, Invert };

struct MaskArg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"coor", "map", "none"});
  }
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT output.msk"
    "\n\nMakes a bulk-solvent mask in the CCP4 format."
    "\nINPUT is a coordinate file (mmCIF, PDB, etc)." },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Timing, 0, "", "timing", Arg::None,
    "  --timing  \tPrint how long individual steps take." },
  { GridSpac, 0, "s", "spacing", Arg::Float,
    "  -s, --spacing=D  \tMax. sampling for the grid (default: 1A)." },
  { GridDims, 0, "g", "grid", Arg::Int3,
    "  -g, --grid=NX,NY,NZ  \tGrid sampling." },
  { Radius, 0, "r", "radius", Arg::Float,
    "  -r, --radius=R  \tUse constant radius of atom spheres." },
  { RProbe, 0, "", "r-probe", Arg::Float,
    "  --r-probe=Rp  \tUse VdW radius + Rp (default: 1.0A)." },
  { RShrink, 0, "", "r-shrink", Arg::Float,
    "  --r-shrink=Rs  \tFinally, remove a shell of thickness Rs (default: 1.1A)." },
  { IslandLimit, 0, "", "island-limit", Arg::Float,
    "  --island-limit=VOL  \tRemove \"islands\" up to VOL A^3." },
  { CctbxCompat, 0, "", "cctbx-compat", Arg::None,
    "  --cctbx-compat  \tUse vdW, Rprobe, Rshrink radii from cctbx." },
  { RefmacCompat, 0, "", "refmac-compat", Arg::None,
    "  --refmac-compat  \tUse radii compatible with Refmac." },
  { Invert, 0, "I", "invert", Arg::None,
    "  -I, --invert  \t0 for solvent, 1 for molecule." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);

  if (p.options[GridDims] && p.options[GridSpac]) {
    std::fprintf(stderr, "Options --grid and --spacing are exclusive.");
    return 1;
  }

  if (p.options[Verbose])
    std::fprintf(stderr, "Converting %s ...\n", input);

  try {
    using std::int8_t;
    Timer timer(p.options[Timing]);
    gemmi::Structure st = gemmi::read_structure_gz(input);
    gemmi::Ccp4<int8_t> mask;
    mask.grid.unit_cell = st.cell;
    mask.grid.spacegroup = st.find_spacegroup();
    if (p.options[GridDims]) {
      auto dims = parse_comma_separated_ints(p.options[GridDims].arg);
      mask.grid.set_size(dims[0], dims[1], dims[2]);
    } else {
      double spac = 1;
      if (p.options[GridSpac])
        spac = std::atof(p.options[GridSpac].arg);
      mask.grid.set_size_from_spacing(spac, true);
    }
    if (p.options[Verbose]) {
      const auto& g = mask.grid;
      std::fprintf(stderr, "Grid: %d x %d x %d\n", g.nu, g.nv, g.nw);
      std::fprintf(stderr, "Spacing: %.3f, %.3f, %.3f\n",
                           g.spacing[0], g.spacing[1], g.spacing[2]);
      size_t np = g.data.size();
      double vol = st.cell.volume;
      std::fprintf(stderr, "Total points: %zu\n", np);
      std::fprintf(stderr, "Unit cell volume: %.1f A^3\n", vol);
      std::fprintf(stderr, "Volume per point: %.3f A^3\n", vol / np);
      if (g.spacegroup) {
        std::fprintf(stderr, "Spacegroup: %s\n", g.spacegroup->hm);
        int na = g.spacegroup->operations().order();
        std::fprintf(stderr, "ASU volume: %.1f A^3\n", vol / na);
        std::fprintf(stderr, "Points per ASU: %.1f\n", double(np) / na);
      } else {
        std::fprintf(stderr, "No spacegroup\n");
      }
    }
    if (st.models.size() > 1)
      std::fprintf(stderr, "Note: only the first model is used.\n");

    const int8_t vmol = 1;
    const int8_t vsol = 0;
    double rshrink = 1.1;
    double island_vol = 0.;
    if (p.options[Radius]) {
      if (p.options[RProbe])
        p.exit_exclusive(Radius, RProbe);
      if (p.options[CctbxCompat])
        p.exit_exclusive(Radius, CctbxCompat);
      if (p.options[RefmacCompat])
        p.exit_exclusive(Radius, RefmacCompat);
      double radius = std::atof(p.options[Radius].arg);
      timer.start();
      gemmi::mask_points_in_constant_radius(mask.grid, st.models[0], radius, vmol);
      timer.print("Points in constant radius masked in");
    } else {
      double rprobe = 1.0;
      gemmi::AtomicRadiiSet radii_set = gemmi::AtomicRadiiSet::VanDerWaals;
      if (p.options[CctbxCompat]) {
        if (p.options[RefmacCompat])
          p.exit_exclusive(CctbxCompat, RefmacCompat);
        radii_set = gemmi::AtomicRadiiSet::Cctbx;
        rprobe = 1.11;
        rshrink = 0.9;
      } else if (p.options[RefmacCompat]) {
        radii_set = gemmi::AtomicRadiiSet::Refmac;
        rprobe = 1.0;
        rshrink = 0.8;
        island_vol = 50;  // the exact value used in Refmac is yet to be found
      }
      if (p.options[RProbe])
        rprobe = std::atof(p.options[RProbe].arg);
      timer.start();
      gemmi::mask_points_in_varied_radius(mask.grid, st.models[0], radii_set, rprobe, vmol);
      timer.print("Points masked in");
    }
    timer.start();
    mask.grid.symmetrize(
        [&](int8_t a, int8_t b) { return a == vmol || b == vmol ? vmol : vsol; });
    timer.print("Mask symmetrized in");
    if (p.options[RShrink])
      rshrink = std::atof(p.options[RShrink].arg);
    if (rshrink > 0) {
      timer.start();
      gemmi::set_margin_around(mask.grid, rshrink, vsol, (int8_t)-1);
      mask.grid.change_values(-1, vsol);
      timer.print("Mask shrunken in");
    }

    if (p.options[IslandLimit])
      island_vol = std::atof(p.options[IslandLimit].arg);
    if (island_vol > 0.) {
      size_t limit = static_cast<size_t>(island_vol * mask.grid.point_count()
                                         / mask.grid.unit_cell.volume);
      timer.start();
      int n = gemmi::remove_islands_in_mask(mask.grid, limit);
      timer.print("Islands removed in");
      if (p.options[Verbose])
        std::fprintf(stderr, "Islands removed: %d\n", n);
    }
    if (p.options[Verbose]) {
      size_t n = std::count(mask.grid.data.begin(), mask.grid.data.end(), vmol);
      std::fprintf(stderr, "Points masked by model: %zu\n", n);
    }
    if (!p.options[Invert])
      for (int8_t& v : mask.grid.data)
        v = (int8_t)1 - v;
    mask.update_ccp4_header(0, true);
    mask.write_ccp4_map(output);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
