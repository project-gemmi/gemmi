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

enum OptionIndex { Hkl=4, Dmin, Rate, Smear, RCutoff, Check };

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
  { RCutoff, 0, "", "r-cutoff", Arg::Float,
    "  --r-cutoff=Y  \tUse atomic radius R such that rho(R) < Y (default: 1e-4)." },
  { Check, 0, "", "check", Arg::None,
    "  --check  \tCalculate exact values and report differences (slow)." },
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
double determine_effective_radius(const T& coef, float b, float min_value) {
  float x1 = 3.5f;
  float y1 = coef.calculate_density(x1*x1, b);
  float x2 = x1;
  float y2 = y1;
  if (y1 < min_value)
    while (y1 < min_value) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5f;
      y1 = coef.calculate_density(x1*x1, b);
    }
  else
    while (y2 > min_value) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5f;
      y2 = coef.calculate_density(x2*x2, b);
    }
  while (x2 - x1 > 0.02f) {
    float new_x = 0.5f * (x2 + x1);
    float new_y = coef.calculate_density(new_x*new_x, b);
    if (new_y < min_value) {
      x2 = new_x;
      y2 = new_y;
    } else {
      x1 = new_x;
      y1 = new_y;
    }
  }
  return x2;
}

struct RhoGridOptions {
  double d_min;
  double rate = 1.5;
  double smear = 0.0;
  float r_cutoff = 1e-4f;
};

template <typename T>
void add_atom_density_to_grid(const Atom& atom, Grid<T>& grid,
                              const RhoGridOptions& opt) {
  auto& scat = IT92<T>::get(atom.element);
  double b = atom.b_iso + opt.smear;
  double radius = determine_effective_radius(scat, (float) b, opt.r_cutoff);
  Fractional fpos = grid.unit_cell.fractionalize(atom.pos).wrap_to_unit();
  grid.use_points_around(fpos, radius, [&](T& point, double r2) {
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
                             bool verbose, bool check) {
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
  double sum_sq_diff = 0.;
  double sum_abs = 0.;
  double max_abs_df = 0.;
  int count = 0;
  double max_1_d2 = 1. / (opt.d_min * opt.d_min);
  for (int h = -sf.nu / 2; h < sf.nu / 2; ++h)
    for (int k = -sf.nv / 2; k < sf.nv / 2; ++k)
      for (int l = 0; l < sf.nw / 2; ++l) {
        double hkl_1_d2 = sf.unit_cell.calculate_1_d2(h, k, l);
        if (hkl_1_d2 < max_1_d2) {
          std::complex<double> value = sf.data[sf.index_n(h, k, l)];
          value *= std::exp(opt.smear * 0.25 * hkl_1_d2);
          if (check) {
            std::complex<double> exact = calculate_structure_factor(st, h,k,l);
            double abs_df = std::abs(value - exact);
            sum_sq_diff += sq(abs_df);
            sum_abs += std::abs(exact);
            if (abs_df > max_abs_df)
              max_abs_df = abs_df;
            ++count;
            printf(" (%d %d %d)\t%7.2f\t%7.2f  \t%6.2f\t%6.2f\td=%5.2f\n",
                   h, k, l, std::abs(value), std::abs(exact),
                   gemmi::phase_in_angles(value), gemmi::phase_in_angles(exact),
                   1. / std::sqrt(hkl_1_d2));
          } else {
            print_sf(value, h, k, l);
          }
        }
      }
  if (check) {
    double rmse = std::sqrt(sum_sq_diff / count);
    double abs_avg = sum_abs / count;
    fprintf(stderr, "RMSE: %g\tNormalized RMSE: %g%%\tMax |dF|: %g\n",
            rmse, 100. * rmse / abs_avg, max_abs_df);
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
        if (p.options[RCutoff])
          opt.r_cutoff = (float) std::strtod(p.options[RCutoff].arg, nullptr);
        print_structure_factors(st, opt, p.options[Verbose], p.options[Check]);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
