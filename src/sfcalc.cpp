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

enum OptionIndex { Hkl=4, Dmin };

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

template <typename T>
void add_atom_density_to_grid(const Atom& atom, Grid<T>& grid) {
  const double radius = 6; // FIXME
  Fractional fpos = grid.unit_cell.fractionalize(atom.pos).wrap_to_unit();
  auto& scat = IT92<T>::get(atom.element);
  grid.use_points_around(fpos, radius, [&](T& point, double r2) {
      point += (T) atom.occ * scat.calculate_density((T)r2, atom.b_iso);
  });
}

void add_model_density_to_grid(const Model& model, Grid<float>& grid) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        add_atom_density_to_grid(atom, grid);
}

template <typename T>
void set_grid_cell_and_spacegroup(Grid<T>& grid, const Structure& st) {
  grid.unit_cell = st.cell;
  grid.spacegroup = find_spacegroup_by_name(st.spacegroup_hm);
}

void put_first_model_density_on_grid(const Structure& st, Grid<float>& grid,
                                     double grid_spacing) {
  grid.data.clear();
  set_grid_cell_and_spacegroup(grid, st);
  grid.set_size_from_spacing(grid_spacing, true);
  add_model_density_to_grid(st.models.at(0), grid);
  grid.symmetrize([](double a, double b) { return a + b; });
}

void print_structure_factors(const Structure& st, double d_min) {
  Grid<float> grid;
  put_first_model_density_on_grid(st, grid, d_min / 3.);
  Grid<std::complex<float>> sf = transform_map_to_f_phi(grid, /*half_l=*/true);
  double max_1_d2 = 1. / std::sqrt(d_min);
  for (int h = -sf.nu / 2; h < sf.nu / 2; ++h)
    for (int k = -sf.nv / 2; k < sf.nv / 2; ++k)
      for (int l = 0; l < sf.nw / 2; ++l)
        if (sf.unit_cell.calculate_1_d2(h, k, l) < max_1_d2)
          print_sf(sf.data[sf.index_n(h, k, l)], h, k, l);
}

} // namespace


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      gemmi::Structure st = gemmi::read_structure_gz(input);
      for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
        std::vector<int> hkl = parse_comma_separated_ints(opt->arg);
        std::complex<double> sf =
          calculate_structure_factor(st, hkl[0], hkl[1], hkl[2]);
        print_sf(sf, hkl[0], hkl[1], hkl[2]);
      }
      if (p.options[Dmin]) {
        double dmin = std::strtod(p.options[Dmin].arg, nullptr);
        print_structure_factors(st, dmin);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
