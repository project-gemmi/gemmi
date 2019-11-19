// Copyright 2019 Global Phasing Ltd.
//
// Calculate structure factors from a molecular model.

#include <stdio.h>
#include <complex>
#include <gemmi/fourier.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/it92.hpp>

#define GEMMI_PROG sfcalc
#include "options.h"

enum OptionIndex { Hkl=4, Dmin, Rate, Smear, RadiusMult };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE\n\n"
    "Calculates structure factors of a model (PDB or mmCIF file).\n\n"
    "Either directly calculates reflections specified by option --hkl\n"
    "or uses FFT to calculate all reflections up to requested resolution.\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Hkl, 0, "", "hkl", Arg::Int3,
    "  --hkl=H,K,L  \tCalculate structure factor F_hkl." },
  { Dmin, 0, "", "dmin", Arg::Float,
    "  --dmin=D  \tCalculate structure factors up to given resolution." },
  { NoOp, 0, "", "", Arg::None, "\nOptions affecting FFT-based calculations:" },
  { Rate, 0, "", "rate", Arg::Float,
    "  --rate=R  \tShannon rate used for grid spacing (default: 1.5)." },
  { Smear, 0, "", "smear", Arg::Float,
    "  --smear=X  \tB added for Gaussian smearing (default: 0.0)." },
  { RadiusMult, 0, "", "radius-mult", Arg::Float,
    "  --radius-mult=X  \tMultiply atomic radius by X (default: 1.0)." },
  { 0, 0, 0, 0, 0, 0 }
};

static
std::complex<double> calculate_structure_factor(const gemmi::Structure& st,
                                                int h, int k, int l) {
  using namespace gemmi;
  const Model& model = st.models[0];
  std::complex<double> sf = 0.;
  std::vector<double> scattering_factors((int) El::END, 0.);
  double stol2 = 0.25 * st.cell.calculate_1_d2(h, k, l);
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        El el = atom.element.elem;
        double& sfactor = scattering_factors[(int)el];
        if (sfactor == 0.) {
          using table = IT92<double>;
          if (!table::has(el))
            fail("Missing scattering factor for " +
                 std::string(atom.element.name()));
          sfactor = table::get(el).calculate_sf(stol2);
        }
        Fractional fpos = st.cell.fractionalize(atom.pos);
        std::complex<double> part = structure_factor_part(fpos, h, k, l);
        for (const FTransform& image : st.cell.images)
          part += structure_factor_part(image.apply(fpos), h, k, l);
        sf += atom.occ * sfactor * std::exp(-atom.b_iso * stol2) * part;
      }
  return sf;
}

static
void print_sf(std::complex<double> sf, int h, int k, int l) {
  printf(" (%d %d %d)\t%8.2f\t%6.2f\n",
         h, k, l, std::abs(sf), gemmi::phase_in_angles(sf));
}

namespace {
using namespace gemmi;

struct RhoGridOptions {
  double d_min;
  double rate = 1.5;
  double smear = 0.0;
  double radius_mult = 1.0;
};

template <typename T>
void add_atom_density_to_grid(const Atom& atom, Grid<T>& grid,
                              const RhoGridOptions& opt) {
  const double radius = 6 * opt.radius_mult; // FIXME
  Fractional fpos = grid.unit_cell.fractionalize(atom.pos).wrap_to_unit();
  auto& scat = IT92<T>::get(atom.element);
  grid.use_points_around(fpos, radius, [&](T& point, double r2) {
      double b = atom.b_iso + opt.smear;
      point += (T)atom.occ * scat.calculate_density((T)r2, (T)b);
  });
}

void add_model_density_to_grid(const Model& model, Grid<float>& grid,
                               const RhoGridOptions& opt) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        add_atom_density_to_grid(atom, grid, opt);
}

template <typename T>
void set_grid_cell_and_spacegroup(Grid<T>& grid, const Structure& st) {
  grid.unit_cell = st.cell;
  grid.spacegroup = find_spacegroup_by_name(st.spacegroup_hm);
}

void put_first_model_density_on_grid(const Structure& st, Grid<float>& grid,
                                     const RhoGridOptions& opt) {
  grid.data.clear();
  set_grid_cell_and_spacegroup(grid, st);
  grid.set_size_from_spacing(opt.d_min / (2 * opt.rate), true);
  add_model_density_to_grid(st.models.at(0), grid, opt);
  grid.symmetrize([](double a, double b) { return a + b; });
}

void print_structure_factors(const Structure& st, const RhoGridOptions& opt,
                             bool verbose) {
  Grid<float> grid;
  if (verbose) {
    fprintf(stderr, "Preparing electron density on a grid...\n");
    fflush(stderr);
  }
  put_first_model_density_on_grid(st, grid, opt);
  if (verbose) {
    fprintf(stderr, "FFT of grid %d x %d x %d\n", grid.nu, grid.nv, grid.nw);
    fflush(stderr);
  }
  Grid<std::complex<float>> sf = transform_map_to_f_phi(grid, /*half_l=*/true);
  if (verbose) {
    fprintf(stderr, "Printing results...\n");
    fflush(stderr);
  }
  double max_1_d2 = 1. / std::sqrt(opt.d_min);
  for (int h = -sf.nu / 2; h < sf.nu / 2; ++h)
    for (int k = -sf.nv / 2; k < sf.nv / 2; ++k)
      for (int l = 0; l < sf.nw / 2; ++l) {
        double hkl_1_d2 = sf.unit_cell.calculate_1_d2(h, k, l);
        if (hkl_1_d2 < max_1_d2) {
          std::complex<double> value = sf.data[sf.index_n(h, k, l)];
          double unsmear = std::exp(opt.smear * 0.25 * hkl_1_d2);
          print_sf(unsmear * value, h, k, l);
        }
      }
}

} // namespace


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (p.options[Verbose]) {
        fprintf(stderr, "Reading file %s...\n", input.c_str());
        fflush(stderr);
      }
      gemmi::Structure st = gemmi::read_structure_gz(input);
      for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
        std::vector<int> hkl = parse_comma_separated_ints(opt->arg);
        if (p.options[Verbose])
          fprintf(stderr, "hkl=(%d %d %d) -> d=%g\n", hkl[0], hkl[1], hkl[2],
                  st.cell.calculate_d(hkl[0], hkl[1], hkl[2]));
        std::complex<double> sf =
          calculate_structure_factor(st, hkl[0], hkl[1], hkl[2]);
        print_sf(sf, hkl[0], hkl[1], hkl[2]);
      }
      if (p.options[Dmin]) {
        RhoGridOptions opt;
        opt.d_min = std::strtod(p.options[Dmin].arg, nullptr);
        if (p.options[Rate])
          opt.rate = std::strtod(p.options[Rate].arg, nullptr);
        if (p.options[Smear])
          opt.smear = std::strtod(p.options[Smear].arg, nullptr);
        if (p.options[RadiusMult])
          opt.radius_mult = std::strtod(p.options[RadiusMult].arg, nullptr);
        print_structure_factors(st, opt, p.options[Verbose]);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
