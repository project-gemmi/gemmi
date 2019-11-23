// Copyright 2019 Global Phasing Ltd.
//
// Calculate structure factors from a molecular model.

#include <stdio.h>
#include <chrono>
#include <complex>
#include <gemmi/fourier.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/it92.hpp>
#include <gemmi/fileutil.hpp>  // for file_open

#define GEMMI_PROG sfcalc
#include "options.h"

enum OptionIndex { Hkl=4, Dmin, Rate, Smear, RCut, Check };

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
    "  --smear=X  \tB added for Gaussian smearing (default: auto)." },
  { RCut, 0, "", "rcut", Arg::Float,
    "  --rcut=Y  \tUse atomic radius r such that rho(r) < Y (default: 5e-5)." },
  { Check, 0, "", "check", Arg::Optional,
    "  --check[=CACHE]  \tCalculate exact values and report differences (slow)." },
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
  double smear = 0.;
  float r_cut = 5e-5f;
};

template <typename T>
void add_atom_density_to_grid(const Atom& atom, Grid<T>& grid,
                              const RhoGridOptions& opt) {
  auto& scat = IT92<T>::get(atom.element);
  double b = atom.b_iso + opt.smear;
  double radius = determine_effective_radius(scat, (float) b, opt.r_cut);
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
                             bool verbose, bool check, const char* cache_file) {
  using Clock = std::chrono::steady_clock;
  Grid<float> grid;
  if (verbose) {
    fprintf(stderr, "Preparing electron density on a grid...\n");
    fflush(stderr);
  }
  auto start = Clock::now();
  put_first_model_density_on_grid(st, grid, opt);
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "FFT of grid %d x %d x %d\n", grid.nu, grid.nv, grid.nw);
    fflush(stderr);
    start = Clock::now();
  }
  Grid<std::complex<float>> sf = transform_map_to_f_phi(grid, /*half_l=*/true);
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "Printing results...\n");
    fflush(stderr);
  }
  gemmi::fileptr_t cache(nullptr, nullptr);
  if (cache_file)
    cache = gemmi::file_open(cache_file, "r");
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
            std::complex<double> exact;
            if (cache_file) {
              char buf[100];
              if (fgets(buf, 99, cache.get()) == nullptr)
                gemmi::fail("cannot read line from file");
              int h_, k_, l_;
              double f_abs, f_deg;
              sscanf(buf, " (%d %d %d) %*f %lf %*f %lf",
                     &h_, &k_, &l_, &f_abs, &f_deg);
              if (h_ != h || k_ != k || l_ != l)
                gemmi::fail("Different h k l order in file.");
              exact = std::polar(f_abs, gemmi::rad(f_deg));
            } else {
              exact = calculate_structure_factor(st, h, k, l);
            }
            double abs_df = std::abs(value - exact);
            sum_sq_diff += sq(abs_df);
            sum_abs += std::abs(exact);
            if (abs_df > max_abs_df)
              max_abs_df = abs_df;
            ++count;
            printf(" (%d %d %d)\t%7.2f\t%8.3f \t%6.2f\t%7.3f\td=%5.2f\n",
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
    fprintf(stderr, "RMSE: %#.5g\t%#.5g%%\tMax |dF|: %#.5g",
            rmse, 100. * rmse / abs_avg, max_abs_df);
    if (!verbose) {
      std::chrono::duration<double> elapsed = Clock::now() - start;
      fprintf(stderr, "\t%#.5gs", elapsed.count());
    }
    fprintf(stderr, "\n");
  }
}

double get_minimum_b_iso(const Model& model) {
  double b_min = 1000.;
  for (const Chain& chain : model.chains)
    for (const Residue& residue : chain.residues)
      for (const Atom& atom : residue.atoms)
        if (atom.b_iso < b_min)
          b_min = atom.b_iso;
  return b_min;
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
        if (p.options[RCut])
          opt.r_cut = (float) std::strtod(p.options[RCut].arg, nullptr);

        if (p.options[Smear]) {
          opt.smear = std::strtod(p.options[Smear].arg, nullptr);
        } else if (opt.rate < 3) {
          // ITfC vol B section 1.3.4.4.5 has formula
          // B = log Q / (sigma * (sigma - 1) * d^*_max^2)
          // But it is more complex. If we take into account numerical
          // accuracy, optimal B depends also on atomic cutoff radius
          // and on distribution of B-factors. On the other hand increasing
          // B increases the radius. He we use an ad-hoc rule.
          double sqrtB = 4 * opt.d_min * (1./opt.rate - 0.2);
          double b_min = get_minimum_b_iso(st.models.at(0));
          opt.smear = sqrtB * sqrtB - b_min;
          if (p.options[Verbose])
            fprintf(stderr, "B_min=%g, B_add=%g\n", b_min, opt.smear);
        }

        print_structure_factors(st, opt, p.options[Verbose],
                                p.options[Check], p.options[Check].arg);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
