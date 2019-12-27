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
#include <gemmi/math.hpp>      // for sq
#include <gemmi/rhogrid.hpp>   // for put_first_model_density_on_grid
#include <gemmi/sfcalc.hpp>    // for calculate_structure_factor
#include <gemmi/smcif.hpp>     // for make_atomic_structure_from_block
#include <gemmi/fprime.hpp>    // for cromer_libermann

#define GEMMI_PROG sfcalc
#include "options.h"

enum OptionIndex { Hkl=4, Dmin, Rate, Smear, RCut, Check,
                   NoFp, NoFileFp, Measured, Scale };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE\n\n"
    "Calculates structure factors of a model (PDB, mmCIF or SMX CIF file).\n\n"
    "Uses FFT to calculate all reflections up to requested resolution for MX\n"
    "files. Otherwise, for SMX and --hkl, F's are calculated directly.\n"
    "This program can also compare F's calculated directly with values\n"
    "calculated through FFT or with values read from a reflection file.\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Hkl, 0, "", "hkl", Arg::Int3,
    "  --hkl=H,K,L  \tCalculate structure factor F_hkl." },
  { Dmin, 0, "", "dmin", Arg::Float,
    "  --dmin=D  \tCalculate structure factors up to given resolution." },
  { NoFp, 0, "", "nofp", Arg::None,
    "  --nofp  \tIgnore f' (anomalous dispersion scattering)." },
  { NoFileFp, 0, "", "no-file-fp", Arg::None,
    "  --no-file-fp  \tIgnore _atom_type_scat_dispersion_real." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for FFT-based calculations:" },
  { Rate, 0, "", "rate", Arg::Float,
    "  --rate=R  \tShannon rate used for grid spacing (default: 1.5)." },
  { Smear, 0, "", "smear", Arg::Float,
    "  --smear=X  \tB added for Gaussian smearing (default: auto)." },
  { RCut, 0, "", "rcut", Arg::Float,
    "  --rcut=Y  \tUse atomic radius r such that rho(r) < Y (default: 5e-5)." },
  { Check, 0, "", "check", Arg::Optional,
    "  --check[=CACHE]  \tCalculate exact values and report differences (slow)." },
  { NoOp, 0, "", "", Arg::None,
    "\nOptions for comparing calculated values with values from a file:" },
  { NoOp, 0, "", "", Arg::None,
    "  --check[=FILE]  \tRe-calculate Fcalc and report differences." },
  { Measured, 0, "", "meas", Arg::None,
    "  --meas  \tCompare against measured values _refln_F_(squared_)meas." },
  { Scale, 0, "", "scale", Arg::Float,
    "  --scale=S  \tMultiply calculated F by sqrt(S) (default: 1)." },
  { 0, 0, 0, 0, 0, 0 }
};

static
void print_sf(std::complex<double> sf, const gemmi::Miller& hkl) {
  printf(" (%d %d %d)\t%8.2f\t%6.2f\n",
         hkl[0], hkl[1], hkl[2], std::abs(sf), gemmi::phase_in_angles(sf));
}

struct Comparator {
  double sum_sq_diff = 0.;
  double sum_sq1 = 0.;
  double sum_sq2 = 0.;
  double sum_abs = 0.;
  double max_abs_df = 0.;
  double sum_abs_diff = 0.;
  int count = 0;

  template<typename T> void add(T value, T exact) {
    double abs_df = std::abs(value - exact);
    sum_sq_diff += abs_df * abs_df;
    sum_sq1 += gemmi::sq(std::abs(value));
    sum_sq2 += gemmi::sq(std::abs(exact));
    sum_abs += std::abs(exact);
    sum_abs_diff += std::abs(std::abs(value) - std::abs(exact));
    if (abs_df > max_abs_df)
      max_abs_df = abs_df;
    ++count;
  }

  double rmse() const { return std::sqrt(sum_sq_diff / count); }
  double abs_avg() const { return sum_abs / count; }
  double weighted_rmse() const { return rmse() / abs_avg(); }
  double rfactor() const { return sum_abs_diff / sum_abs; }
  double scale() const { return std::sqrt(sum_sq1 / sum_sq2); }
};

using namespace gemmi;

static
void print_structure_factors(const Structure& st, const RhoGridOptions& opt,
                             bool verbose, bool check, const char* cache_file) {
  using Clock = std::chrono::steady_clock;
  Grid<float> grid;
  if (verbose) {
    fprintf(stderr, "Preparing electron density on a grid...\n");
    fflush(stderr);
  }
  auto start = Clock::now();
  using Table = IT92<double>;
  put_first_model_density_on_grid<Table>(st, grid, opt);
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "FFT of grid %d x %d x %d\n", grid.nu, grid.nv, grid.nw);
    fflush(stderr);
    start = Clock::now();
  }
  Grid<std::complex<float>> sf = transform_map_to_f_phi(grid, /*half_l=*/true);
  StructureFactorCalculator<Table> calc(st.cell);
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "Printing results...\n");
    fflush(stderr);
  }
  gemmi::fileptr_t cache(nullptr, nullptr);
  if (cache_file)
    cache = gemmi::file_open(cache_file, "r");
  Comparator comparator;
  double max_1_d2 = 1. / (opt.d_min * opt.d_min);
  for (int h = -sf.nu / 2; h < sf.nu / 2; ++h)
    for (int k = -sf.nv / 2; k < sf.nv / 2; ++k)
      for (int l = 0; l < sf.nw / 2; ++l) {
        Miller hkl{{h, k, l}};
        double hkl_1_d2 = sf.unit_cell.calculate_1_d2(hkl);
        if (hkl_1_d2 < max_1_d2) {
          std::complex<double> value = sf.data[sf.index_n(h, k, l)];
          value *= std::exp(opt.smear * 0.25 * hkl_1_d2);
          if (check) {
            std::complex<double> exact;
            if (cache_file) {
              char cache_line[100];
              if (fgets(cache_line, 99, cache.get()) == nullptr)
                gemmi::fail("cannot read line from file");
              int cache_h, cache_k, cache_l;
              double f_abs, f_deg;
              sscanf(cache_line, " (%d %d %d) %*f %lf %*f %lf",
                     &cache_h, &cache_k, &cache_l, &f_abs, &f_deg);
              if (cache_h != h || cache_k != k || cache_l != l)
                gemmi::fail("Different h k l order than in cache file.");
              exact = std::polar(f_abs, gemmi::rad(f_deg));
            } else {
              exact = calc.calculate_sf_from_model(st.models[0], hkl);
            }
            comparator.add(value, exact);
            printf(" (%d %d %d)\t%7.2f\t%8.3f \t%6.2f\t%7.3f\td=%5.2f\n",
                   h, k, l, std::abs(value), std::abs(exact),
                   gemmi::phase_in_angles(value), gemmi::phase_in_angles(exact),
                   1. / std::sqrt(hkl_1_d2));
          } else {
            print_sf(value, hkl);
          }
        }
      }
  if (check) {
    fprintf(stderr, "RMSE: %#.5g\t%#.4g%%\tmax|dF|=%#.4g  \tR-factor: %#.3g%%",
            comparator.rmse(), 100. * comparator.weighted_rmse(),
            comparator.max_abs_df, 100 * comparator.rfactor());
    if (!verbose) {
      std::chrono::duration<double> elapsed = Clock::now() - start;
      fprintf(stderr, "\t%#.5gs", elapsed.count());
    }
    fprintf(stderr, "\n");
  }
}

static
double get_minimum_b_iso(const Model& model) {
  double b_min = 1000.;
  for (const Chain& chain : model.chains)
    for (const Residue& residue : chain.residues)
      for (const Atom& atom : residue.atoms)
        if (atom.b_iso < b_min)
          b_min = atom.b_iso;
  return b_min;
}

void fill_atom_types(gemmi::AtomicStructure& ast) {
  if (ast.wavelength > 0.) {
    double energy = hc() / ast.wavelength;
    for (gemmi::AtomicStructure::Site& site : ast.sites) {
      if (!ast.get_atom_type(site.type_symbol)) {
        ast.atom_types.emplace_back();
        gemmi::AtomicStructure::AtomType& at = ast.atom_types.back();
        at.symbol = site.type_symbol;
        gemmi::split_element_and_charge(at.symbol, &at);
        int z = site.element.atomic_number();
        at.dispersion_real = gemmi::cromer_libermann(z, energy,
                                                     &at.dispersion_imag);
      }
    }
  }
}

static
void verify_f_calc(const AtomicStructure& ast, const std::string& suffix,
                   double scale, bool verbose, const char* path) {
  StructureFactorCalculator<IT92<double>> calc(ast.cell);
  cif::Document hkl_doc = read_cif_gz(path);
  cif::Block& block = hkl_doc.blocks.at(0);
  std::string f_cols[2] = {"?F_calc", "?F_squared_calc"};
  cif::Table table = block.find("_refln_",
      {"index_h", "index_k", "index_l", "?F_"+suffix, "?F_squared_"+suffix});
  if (!table.ok())
    fail("_refln_index_ category not found in ", path);
  if (table.has_column(3)) {
    if (verbose)
      fprintf(stderr, "Checking _refln_F_calc from %s\n", path);
  } else if (table.has_column(4)) {
    if (verbose)
      fprintf(stderr, "Checking sqrt of _refln_F_squared_calc from %s\n", path);
  } else {
    fail("neither _refln_F_calc nor _refln_F_squared_calc not found");
  }
  if (verbose && !ast.atom_types.empty())
    fprintf(stderr, "Using f' read from cif file (%u atom types)\n",
            (unsigned) ast.atom_types.size());
  Comparator comparator;
  Miller hkl;
  for (auto row : table) {
    double fcalc_from_file = NAN;
    try {
      for (int i = 0; i != 3; ++i)
        hkl[i] = cif::as_int(row[i]);
      if (row.has(3))
        fcalc_from_file = cif::as_number(row[3]);
      else if (row.has(4))
        fcalc_from_file = std::sqrt(cif::as_number(row[4]));
    } catch(std::runtime_error& e) {
      fprintf(stderr, "Error in _refln_[] in %s: %s\n", path, e.what());
      continue;
    } catch(std::invalid_argument& e) {
      fprintf(stderr, "Error in _refln_[] in %s: %s\n", path, e.what());
      continue;
    }
    double f = std::abs(calc.calculate_sf_from_atomic_structure(ast, hkl));
    f *= scale;
    comparator.add(fcalc_from_file, f);
    if (verbose)
      printf(" (%d %d %d)\t%7.2f\t%8.3f \td=%5.2f\n",
             hkl[0], hkl[1], hkl[2], fcalc_from_file, f,
             ast.cell.calculate_d(hkl));
  }
  fprintf(stderr,
          "RMSE=%#.5g  %#.4g%%  max|dF|=%#.4g  R=%#.3g%%  sum(F^2)_ratio=%g\n",
          comparator.rmse(), 100. * comparator.weighted_rmse(),
          comparator.max_abs_df, 100*comparator.rfactor(), comparator.scale());
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (p.options[Verbose]) {
        fprintf(stderr, "Reading file %s ...\n", input.c_str());
        fflush(stderr);
      }
      gemmi::Structure st = gemmi::read_structure_gz(input);
      gemmi::AtomicStructure ast;
      bool use_st = !st.models.empty();
      if (!use_st) {
        if (giends_with(input, ".cif"))
          ast = gemmi::make_atomic_structure_from_block(
                                    gemmi::read_cif_gz(input).sole_block());
        if (ast.sites.empty() ||
            // COD can have a row of nulls as a placeholder (e.g. 2211708)
            (ast.sites.size() == 1 && ast.sites[0].element == El::X))
          gemmi::fail("no atoms in the file");
        // SM CIF files specify full occupancy for atoms on special positions.
        // We need to adjust it for symmetry calculations.
        for (AtomicStructure::Site& site : ast.sites) {
          int n_mates = ast.cell.is_special_position(site.fract, 0.4);
          if (n_mates != 0)
            site.occ /= (n_mates + 1);
        }
        if (p.options[NoFp] || p.options[NoFileFp])
          ast.atom_types.clear();
        if (!p.options[NoFp] && ast.atom_types.empty())
          fill_atom_types(ast);
      }
      const UnitCell& cell = use_st ? st.cell : ast.cell;
      StructureFactorCalculator<IT92<double>> calc(cell);
      for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
        std::vector<int> hkl_ = parse_comma_separated_ints(opt->arg);
        gemmi::Miller hkl{{hkl_[0], hkl_[1], hkl_[2]}};
        if (p.options[Verbose])
          fprintf(stderr, "hkl=(%d %d %d) -> d=%g\n", hkl[0], hkl[1], hkl[2],
                  cell.calculate_d(hkl));
        if (use_st)
          print_sf(calc.calculate_sf_from_model(st.models[0], hkl), hkl);
        else
          print_sf(calc.calculate_sf_from_atomic_structure(ast, hkl), hkl);
      }
      if (p.options[Dmin]) {
        if (!use_st)
          gemmi::fail("for now small-molecule CIF files work only with --hkl");
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
          double b_min = get_minimum_b_iso(st.models[0]);
          opt.smear = sqrtB * sqrtB - b_min;
          if (p.options[Verbose])
            fprintf(stderr, "B_min=%g, B_add=%g\n", b_min, opt.smear);
        }

        print_structure_factors(st, opt, p.options[Verbose],
                                p.options[Check], p.options[Check].arg);
      } else if (p.options[Check]) {
        if (use_st)
          gemmi::fail("not a SM CIF");
        double scale = 1.0;
        if (p.options[Scale])
          scale = std::strtod(p.options[Scale].arg, nullptr);
        std::string suffix = p.options[Measured] ? "meas" : "calc";
        const char* path = p.options[Check].arg;
        verify_f_calc(ast, suffix, scale, p.options[Verbose], path);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
