// Copyright 2019 Global Phasing Ltd.
//
// Calculate structure factors from a molecular model.

#include <stdio.h>
#include <chrono>
#include <complex>
#include <gemmi/bulksol.hpp>
#include <gemmi/fileutil.hpp>  // for file_open
#include <gemmi/fourier.hpp>
#include <gemmi/gz.hpp>        // for MaybeGzipped
#include <gemmi/gzread.hpp>
#include <gemmi/it92.hpp>      // for IT92
#include <gemmi/itc4322.hpp>   // for ITC4322
#include <gemmi/math.hpp>      // for sq
#include <gemmi/mtz.hpp>       // for read_mtz_file
#include <gemmi/rhogrid.hpp>   // for put_model_density_on_grid
#include <gemmi/sfcalc.hpp>    // for calculate_structure_factor
#include <gemmi/smcif.hpp>     // for make_small_structure_from_block

#define GEMMI_PROG sfcalc
#include "options.h"

namespace {

enum OptionIndex { Hkl=4, Dmin, For, Rate, Blur, RCut, Test, Compare,
                   CifFp, Wavelength, Unknown, NoAniso, ScaleTo, FLabel,
                   PhiLabel, Ksolv, Bsolv, Rprobe, Rshrink };

struct SfCalcArg: public Arg {
  static option::ArgStatus FormFactors(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"xray", "electron"});
  }
};


const option::Descriptor Usage[] = {
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
    "  --dmin=NUM  \tCalculate structure factors up to given resolution." },
  { For, 0, "", "for", SfCalcArg::FormFactors,
    "  --for=TYPE  \tTYPE is xray (default) or electron." },
  { CifFp, 0, "", "ciffp", Arg::None,
    "  --ciffp  \tRead f' from _atom_type_scat_dispersion_real in CIF." },
  { Wavelength, 0, "w", "wavelength", Arg::Float,
    "  --wavelength=NUM  \tWavelength [A] for calculation of f' "
    "(use --wavelength=0 or -w0 to ignore anomalous scattering)." },
  { Unknown, 0, "", "unknown", Arg::Required,
    "  --unknown=SYMBOL  \tUse form factor of SYMBOL for unknown atoms." },
  { NoAniso, 0, "", "noaniso", Arg::None,
    "  --noaniso  \tIgnore anisotropic ADPs." },

  { NoOp, 0, "", "", Arg::None, "\nOptions for FFT-based calculations:" },
  { Rate, 0, "", "rate", Arg::Float,
    "  --rate=NUM  \tShannon rate used for grid spacing (default: 1.5)." },
  { Blur, 0, "", "blur", Arg::Float,
    "  --blur=NUM  \tB added for Gaussian blurring (default: auto)." },
  { RCut, 0, "", "rcut", Arg::Float,
    "  --rcut=Y  \tUse atomic radius r such that rho(r) < Y (default: 5e-5)." },
  { Test, 0, "", "test", Arg::Optional,
    "  --test[=CACHE]  \tCalculate exact values and report differences (slow)." },
  { ScaleTo, 0, "", "scale-to", Arg::ColonPair,
    "  --scale-to=FILE:COL  \tAnisotropic scaling to F from MTZ file." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for bulk solvent correction (only w/ FFT):" },
  { Rprobe, 0, "", "r-probe", Arg::Float,
    "  --r-probe=NUM  \tValue added to VdW radius (default: 1.0A)." },
  { Rshrink, 0, "", "r-shrink", Arg::Float,
    "  --r-shrink=NUM  \tValue for shrinking the solvent area (default: 1.1A)." },
  { Ksolv, 0, "", "ksolv", Arg::Float,
    "  --ksolv=NUM  \tValue (if optimizing: initial value) of k_solv." },
  { Bsolv, 0, "", "bsolv", Arg::Float,
    "  --bsolv=NUM  \tValue (if optimizing: initial value) of B_solv." },

  { NoOp, 0, "", "", Arg::None,
    "\nOptions for comparing calculated values with values from a file:" },
  { Compare, 0, "", "compare", Arg::Required,
    "  --compare=FILE  \tRe-calculate Fcalc and report differences." },
  { FLabel, 0, "", "f", Arg::Required,
    "  --f=LABEL  \tMTZ column label (default: FC) or small molecule cif"
    " tag (default: F_calc or F_squared_calc)." },
  { PhiLabel, 0, "", "phi", Arg::Required,
    "  --phi=LABEL  \tMTZ column label (default: PHIC)" },
  { 0, 0, 0, 0, 0, 0 }
};

struct RefFile {
  enum class Mode { None, Test, Compare };
  Mode mode = Mode::None;
  const char* path = nullptr;
  std::string f_label;
  std::string phi_label;
};

struct SolventParam {
  bool use_solvent = false;
  // initialize with average values (Fokine & Urzhumtsev, 2002)
  double k = 0.35;
  double B = 46.0;

  double calculate_scale(double stol2) const {
    return k * std::exp(-B * stol2);
  }
};

void print_sf(std::complex<double> sf, const gemmi::Miller& hkl) {
  printf(" (%d %d %d)\t%.8f\t%.6f\n",
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
  double phi_diff_weighted_sum = 0.;

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

  void add_complex(std::complex<double> value, std::complex<double> exact) {
    add(value, exact);
    double d = gemmi::angle_abs_diff(gemmi::deg(std::arg(value)),
                                     gemmi::deg(std::arg(exact)));
    phi_diff_weighted_sum += std::abs(exact) * d;
  }

  double rmse() const { return std::sqrt(sum_sq_diff / count); }
  double abs_avg() const { return sum_abs / count; }
  double weighted_rmse() const { return rmse() / abs_avg(); }
  double rfactor() const { return sum_abs_diff / sum_abs; }
  double scale() const { return std::sqrt(sum_sq1 / sum_sq2); }
  double mean_dphi() const { return phi_diff_weighted_sum / sum_abs; }
};

void print_to_stderr(const Comparator& c) {
  fflush(stdout);
  fprintf(stderr, "RMSE=%#.5g  %#.4g%%  max|dF|=%#.4g  R=%.3f%%",
          c.rmse(), 100 * c.weighted_rmse(), c.max_abs_df, 100 * c.rfactor());
}

template<typename Table, typename Real>
void print_structure_factors(const gemmi::Structure& st,
                             gemmi::DensityCalculator<Table, Real>& dencalc,
                             const SolventParam& solvent,
                             bool verbose, const RefFile& file,
                             const gemmi::AsuData<std::array<Real,2>>& scale_to) {
  using Clock = std::chrono::steady_clock;
  if (verbose) {
    fprintf(stderr, "Preparing electron density on a grid...\n");
    fflush(stderr);
  }
  auto start = Clock::now();
  dencalc.set_grid_cell_and_spacegroup(st);
  dencalc.put_model_density_on_grid(st.models[0]);
  const gemmi::Grid<Real>& grid = dencalc.grid;
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "FFT of grid %d x %d x %d\n", grid.nu, grid.nv, grid.nw);
    fflush(stderr);
    start = Clock::now();
  }
  gemmi::FPhiGrid<Real> sf = transform_map_to_f_phi(grid, /*half_l=*/true);
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fprintf(stderr, "...took %g s.\n", elapsed.count());
    fprintf(stderr, "Printing results...\n");
    fflush(stderr);
  }
  gemmi::StructureFactorCalculator<Table> calc(st.cell);
  for (int i = 0; i != (int)gemmi::El::END; ++i)
    if (dencalc.fprimes[i] != 0.f)
      calc.set_fprime((gemmi::El)i, dencalc.fprimes[i]);
  gemmi::fileptr_t cache(nullptr, nullptr);
  gemmi::AsuData<std::complex<double>> compared_data;
  if (file.path) {
    if (file.mode == RefFile::Mode::Test) {
      cache = gemmi::file_open(file.path, "r");
    } else if (file.mode == RefFile::Mode::Compare) {
      gemmi::Mtz mtz;
      mtz.read_input(gemmi::MaybeGzipped(file.path), true);
      compared_data = mtz.get_f_phi<double>(file.f_label, file.phi_label);
      compared_data.ensure_sorted();
    }
  }
  auto asu_data = sf.prepare_asu_data(dencalc.d_min, dencalc.blur);

  if (solvent.use_solvent) {
    dencalc.put_solvent_mask_on_grid(st.models[0]);
    auto asu_mask = transform_map_to_f_phi(dencalc.grid, /*half_l=*/true)
                    .prepare_asu_data(dencalc.d_min, 0);
    assert(asu_mask.v.size() == asu_data.v.size());
    for (size_t i = 0; i != asu_data.v.size(); ++i) {
      const gemmi::HklValue<std::complex<Real>>& m = asu_mask.v[i];
      assert(asu_data.v[i].hkl == m.hkl);
      double stol2 = asu_data.unit_cell().calculate_stol_sq(m.hkl);
      asu_data.v[i].value += (Real) solvent.calculate_scale(stol2) * m.value;
    }
  }

  Comparator comparator;
  if (scale_to.size() != 0) {
    std::vector<gemmi::FcFo<Real>> data = prepare_fc_fo(asu_data, scale_to);
    gemmi::BulkSolvent<Real> bulk(asu_data.unit_cell(), asu_data.spacegroup());
    printf("Calculating scale factors using %zu points...\n", data.size());
    bulk.quick_iso_fit(data);
    //fprintf(stderr, "k_ov=%g B_ov=%g\n", bulk.k_overall, bulk.B_aniso.u11);
    bulk.aniso_fit(data);
    fprintf(stderr, "k_ov=%g B11=%g B22=%g B33=%g B12=%g B13=%g B23=%g\n",
            bulk.k_overall, bulk.B_aniso.u11, bulk.B_aniso.u22, bulk.B_aniso.u33,
                            bulk.B_aniso.u12, bulk.B_aniso.u13, bulk.B_aniso.u23);
    for (typename gemmi::FPhiGrid<Real>::HklValue& hv : asu_data.v)
      hv.value *= bulk.get_scale_factor_aniso(hv.hkl);
  }
  for (typename gemmi::FPhiGrid<Real>::HklValue& hv : asu_data.v) {
    if (file.mode == RefFile::Mode::None) {
      print_sf(hv.value, hv.hkl);
    } else {
      std::complex<double> exact;
      if (file.path) {
        if (file.mode == RefFile::Mode::Test) {
          char cache_line[100];
          if (fgets(cache_line, 99, cache.get()) == nullptr)
            gemmi::fail("cannot read line from file");
          gemmi::Miller cache_hkl;
          double f_abs, f_deg;
          sscanf(cache_line, " (%d %d %d) %*f %lf %*f %lf",
                 &cache_hkl[0], &cache_hkl[1], &cache_hkl[2], &f_abs, &f_deg);
          if (cache_hkl != hv.hkl)
            gemmi::fail("Different h k l order than in cache file.");
          exact = std::polar(f_abs, gemmi::rad(f_deg));
        } else if (file.mode == RefFile::Mode::Compare) {
          auto it = std::lower_bound(compared_data.v.begin(), compared_data.v.end(), hv.hkl);
          if (it == compared_data.v.end() || it->hkl != hv.hkl)
            continue;
          exact = it->value;
        }
      } else {
        exact = calc.calculate_sf_from_model(st.models[0], hv.hkl);
      }
      comparator.add_complex(hv.value, exact);
      printf(" (%d %d %d)\t%7.2f\t%8.3f \t%6.2f\t%7.3f\td=%5.2f\n",
             hv.hkl[0], hv.hkl[1], hv.hkl[2], std::abs(hv.value), std::abs(exact),
             gemmi::phase_in_angles(hv.value), gemmi::phase_in_angles(exact),
             sf.unit_cell.calculate_d(hv.hkl));
    }
  }
  if (file.mode != RefFile::Mode::None) {
    print_to_stderr(comparator);
    fprintf(stderr, "  <dPhi>=%#.4g", comparator.mean_dphi());
    if (!verbose) {
      std::chrono::duration<double> elapsed = Clock::now() - start;
      fprintf(stderr, "   %#.5gs", elapsed.count());
    }
    fprintf(stderr, "\n");
  }
}

template<typename Table>
void print_structure_factors_sm(const gemmi::SmallStructure& small,
                                gemmi::StructureFactorCalculator<Table>& calc,
                                double d_min, bool verbose) {
  using Clock = std::chrono::steady_clock;
  auto start = Clock::now();
  int counter = 0;
  double max_1_d = 1. / d_min;
  int max_h = int(max_1_d / small.cell.ar);
  int max_k = int(max_1_d / small.cell.br);
  int max_l = int(max_1_d / small.cell.cr);
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name(small.spacegroup_hm,
                                                       small.cell.alpha, small.cell.gamma);
  gemmi::ReciprocalAsu asu(sg ? sg : &gemmi::get_spacegroup_p1());
  for (int h = -max_h; h <= max_h; ++h)
    for (int k = -max_k; k <= max_k; ++k)
      for (int l = 0; l <= max_l; ++l) {
        gemmi::Miller hkl{{h, k, l}};
        if (!asu.is_in(hkl))
          continue;
        double hkl_1_d2 = small.cell.calculate_1_d2(hkl);
        if (hkl_1_d2 < max_1_d * max_1_d) {
          auto value = calc.calculate_sf_from_small_structure(small, hkl);
          print_sf(value, hkl);
          ++counter;
        }
      }
  if (verbose) {
    std::chrono::duration<double> elapsed = Clock::now() - start;
    fflush(stdout);
    fprintf(stderr, "Calculated %d SFs in %g s.\n", counter, elapsed.count());
    fflush(stderr);
  }
}

double get_minimum_b_iso(const gemmi::Model& model) {
  double b_min = 1000.;
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& residue : chain.residues)
      for (const gemmi::Atom& atom : residue.atoms)
        if (atom.b_iso < b_min)
          b_min = atom.b_iso;
  return b_min;
}

template<typename Table>
void compare_with_hkl(const gemmi::SmallStructure& small,
                      gemmi::StructureFactorCalculator<Table>& calc,
                      const RefFile& file,
                      bool verbose,
                      Comparator& comparator) {
  namespace cif = gemmi::cif;
  cif::Document hkl_doc = gemmi::read_cif_gz(file.path);
  cif::Block& block = hkl_doc.blocks.at(0);
  std::vector<std::string> tags =
    {"index_h", "index_k", "index_l", "?F_calc", "?F_squared_calc"};
  if (!file.f_label.empty()) {
    tags.pop_back();
    tags.back().replace(1, std::string::npos, file.f_label);
  }
  cif::Table table = block.find("_refln_", tags);
  if (!table.ok())
    gemmi::fail("_refln_index_ category not found in ", file.path);
  int col = 0;
  if (table.has_column(3))
    col = 3;
  else if (tags.size() > 4 && table.has_column(4))
    col = 4;
  if (col == 0) {
    std::string msg;
    if (file.f_label.empty())
      msg = "Neither _refln_F_calc nor _refln_F_squared_calc";
    else
      msg = "_refln_" + file.f_label;
    gemmi::fail(msg + " not found in: ", file.path);
  }
  bool use_sqrt = (col == 4 ||
                   file.f_label == "F_squared_calc" ||
                   file.f_label == "F_squared_meas");
  if (verbose)
    fprintf(stderr, "Checking %s_refln_%s from %s\n",
            use_sqrt ? "sqrt of " : "", tags[col].c_str()+1, file.path);
  gemmi::Miller hkl;
  int missing = 0;
  int negative = 0;
  for (auto row : table) {
    if (!row.has2(col)) {
      missing++;
      continue;
    }
    double f_from_file = NAN;
    try {
      for (int i = 0; i != 3; ++i)
        hkl[i] = cif::as_int(row[i]);
      f_from_file = cif::as_number(row[col]);
      if (use_sqrt) {
        if (f_from_file >= 0) {
          f_from_file = std::sqrt(f_from_file);
        } else {
          negative++;
          f_from_file = 0;
        }
      }
    } catch(std::runtime_error& e) {
      fprintf(stderr, "Error in _refln_[] in %s: %s\n", file.path, e.what());
      continue;
    } catch(std::invalid_argument& e) {
      fprintf(stderr, "Error in _refln_[] in %s: %s\n", file.path, e.what());
      continue;
    }
    double f = std::abs(calc.calculate_sf_from_small_structure(small, hkl));
    comparator.add(f_from_file, f);
    if (verbose)
      printf(" (%d %d %d)\t%7.2f\t%8.3f \td=%5.2f\n",
             hkl[0], hkl[1], hkl[2], f_from_file, f,
             small.cell.calculate_d(hkl));
  }
  if (missing)
    fprintf(stderr, "missing value in %d rows\n", missing);
  if (negative)
    fprintf(stderr, "negative value in %d rows\n", negative);
}

template<typename Table>
void compare_with_mtz(const gemmi::Model& model, const gemmi::UnitCell& cell,
                      gemmi::StructureFactorCalculator<Table>& calc,
                      const RefFile& file, bool verbose, Comparator& comparator) {
  gemmi::Mtz mtz;
  mtz.read_input(gemmi::MaybeGzipped(file.path), true);
  gemmi::Mtz::Column* col = mtz.column_with_label(file.f_label);
  if (!col)
    gemmi::fail("MTZ file has no column with label: " + file.f_label);
  gemmi::MtzDataProxy data_proxy{mtz};
  for (size_t i = 0; i < data_proxy.size(); i += data_proxy.stride()) {
    gemmi::Miller hkl = data_proxy.get_hkl(i);
    double f_from_file = data_proxy.get_num(i + col->idx);
    double f = std::abs(calc.calculate_sf_from_model(model, hkl));
    comparator.add(f_from_file, f);
    if (verbose)
      printf(" (%d %d %d)\t%7.2f\t%8.3f \td=%5.2f\n",
             hkl[0], hkl[1], hkl[2], f_from_file, f, cell.calculate_d(hkl));
  }
}

template<typename Table>
void process_with_table(bool use_st, gemmi::Structure& st, const gemmi::SmallStructure& small,
                        double wavelength, const OptParser& p) {
  const gemmi::UnitCell& cell = use_st ? st.cell : small.cell;
  gemmi::StructureFactorCalculator<Table> calc(cell);

  // assign f'
  if (p.options[CifFp]) {
    if (use_st) {
      // _atom_type.scat_dispersion_real is almost never used,
      // so for now we ignore it.
    } else { // small molecule
      if (p.options[Verbose])
        fprintf(stderr, "Using f' read from cif file (%u atom types)\n",
                (unsigned) small.atom_types.size());
      for (const gemmi::SmallStructure::AtomType& atom_type : small.atom_types)
        calc.set_fprime(atom_type.element, atom_type.dispersion_real);
    }
  }

  auto present_elems = use_st ? st.models[0].present_elements()
                              : small.present_elements();
  if (present_elems[(int)gemmi::El::X])
    gemmi::fail("unknown element. Add --unknown=O to treat unknown atoms as oxygen.");
  for (size_t i = 1; i != present_elems.size(); ++i)
    if (present_elems[i] && !Table::has((gemmi::El)i))
      gemmi::fail("Missing form factor for element ", element_name((gemmi::El)i));
  if (wavelength > 0) {
    double energy = gemmi::hc() / wavelength;
    for (int z = 1; z <= 92; ++z)
      if (present_elems[z]) {
        double fprime = gemmi::cromer_libermann(z, energy, nullptr);
        calc.set_fprime_if_not_set((gemmi::El)z, fprime);
      }
  }

  // handle option --hkl
  for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
    std::vector<int> hkl_ = parse_comma_separated_ints(opt->arg);
    gemmi::Miller hkl{{hkl_[0], hkl_[1], hkl_[2]}};
    if (p.options[Verbose])
      fprintf(stderr, "hkl=(%d %d %d) -> d=%g\n", hkl[0], hkl[1], hkl[2],
              cell.calculate_d(hkl));
    if (use_st)
      print_sf(calc.calculate_sf_from_model(st.models[0], hkl), hkl);
    else
      print_sf(calc.calculate_sf_from_small_structure(small, hkl), hkl);
  }

  RefFile file;
  if (p.options[Test]) {
    file.mode = RefFile::Mode::Test;
    file.path = p.options[Test].arg;
  } else if (p.options[Compare]) {
    file.mode = RefFile::Mode::Compare;
    file.path = p.options[Compare].arg;
  }
  if (p.options[FLabel])
    file.f_label = p.options[FLabel].arg;
  else if (use_st)
    file.f_label = "FC";
  if (p.options[PhiLabel])
   file.phi_label = p.options[PhiLabel].arg;
  else if (use_st)
    file.phi_label = "PHIC";

  using Real = float;
  gemmi::AsuData<std::array<Real,2>> scale_to;
  if (p.options[ScaleTo]) {
    const char* arg = p.options[ScaleTo].arg;
    const char* sep = std::strchr(arg, ':');
    std::string path(arg, sep);
    std::string label(sep+1);
    gemmi::Mtz mtz;
    mtz.read_input(gemmi::MaybeGzipped(path), true);
    scale_to = mtz.get_values<Real,2>({label, "SIG"+label});
    scale_to.ensure_sorted();
  }

  // handle option --dmin
  if (p.options[Dmin]) {
    double d_min = std::atof(p.options[Dmin].arg);
    if (use_st) {
      gemmi::DensityCalculator<Table, Real> dencalc;
      dencalc.d_min = d_min;
      if (p.options[Rate])
        dencalc.rate = std::atof(p.options[Rate].arg);
      if (p.options[RCut])
        dencalc.r_cut = (float) std::atof(p.options[RCut].arg);
      for (auto& it : calc.fprimes())
        dencalc.fprimes[(int)it.first] = (float) it.second;
      if (p.options[Blur]) {
        dencalc.blur = std::atof(p.options[Blur].arg);
      } else if (dencalc.rate < 3) {
        // ITfC vol B section 1.3.4.4.5 has formula
        // B = log Q / (sigma * (sigma - 1) * d*_max ^2)
        // where Q is quality factor, sigma is the oversampling rate.
        // This value is not optimal.
        // The optimal value would depend on the distribution of B-factors
        // and on the atomic cutoff radius, and probably it would be too
        // hard to estimate.
        // Here we use a simple ad-hoc rule:
        double sqrtB = 4 * dencalc.d_min * (1./dencalc.rate - 0.2);
        double b_min = get_minimum_b_iso(st.models[0]);
        dencalc.blur = sqrtB * sqrtB - b_min;
        if (p.options[Verbose])
          fprintf(stderr, "B_min=%g, B_add=%g\n", b_min, dencalc.blur);
      }
      if (p.options[Rprobe])
        dencalc.rprobe = std::atof(p.options[Rprobe].arg);
      if (p.options[Rshrink])
        dencalc.rshrink = std::atof(p.options[Rshrink].arg);

      SolventParam solvent;
      if (p.options[Ksolv] || p.options[Bsolv]) {
        solvent.use_solvent = true;
        if (p.options[Ksolv])
          solvent.k = std::atof(p.options[Ksolv].arg);
        if (p.options[Bsolv])
          solvent.B = std::atof(p.options[Bsolv].arg);
      }

      print_structure_factors(st, dencalc, solvent, p.options[Verbose], file, scale_to);
    } else {
      if (p.options[Rate] || p.options[RCut] || p.options[Blur] ||
          p.options[Test])
        gemmi::fail("Small molecule SFs are calculated directly. Do not use any\n"
                    "of the FFT-related options: --rate, --blur, --rcut, --test.");
      print_structure_factors_sm(small, calc, d_min, p.options[Verbose]);
    }

  // handle option --compare
  } else if (file.mode == RefFile::Mode::Compare) {
    Comparator comparator;
    if (use_st)
      compare_with_mtz(st.models[0], st.cell, calc, file, p.options[Verbose], comparator);
    else
      compare_with_hkl(small, calc, file, p.options[Verbose], comparator);
    print_to_stderr(comparator);
    fprintf(stderr, "  sum(F^2)_ratio=%g\n", comparator.scale());
  }
}

void process(const std::string& input, const OptParser& p) {
  // read (Small)Structure
  gemmi::Structure st = gemmi::read_structure_gz(input);
  gemmi::SmallStructure small;
  bool use_st = !st.models.empty();
  if (!use_st) {
    if (gemmi::giends_with(input, ".cif"))
      small = gemmi::make_small_structure_from_block(
                                gemmi::read_cif_gz(input).sole_block());
    if (small.sites.empty() ||
        // COD can have a row of nulls as a placeholder (e.g. 2211708)
        (small.sites.size() == 1 && small.sites[0].element == gemmi::El::X))
      gemmi::fail("no atoms in the file");
    // SM CIF files specify full occupancy for atoms on special positions.
    // We need to adjust it for symmetry calculations.
    for (gemmi::SmallStructure::Site& site : small.sites) {
      int n_mates = small.cell.is_special_position(site.fract, 0.4);
      if (n_mates != 0)
        site.occ /= (n_mates + 1);
    }
  }

  if (p.options[ScaleTo] || p.options[Ksolv] || p.options[Bsolv] ||
      p.options[Rprobe] || p.options[Rshrink]) {
    static const char* msg =
      "Options --scale-to, --ksolv, --bsolv, --r-probe, --r-shrink works only with ";
    if (!p.options[Dmin])
      gemmi::fail(msg, "--dmin");
    if (!use_st)
      gemmi::fail(msg, "macromolecular structures");
  }

  if (p.options[NoAniso]) {
    if (use_st) {
      for (gemmi::CRA& cra : st.models[0].all())
        cra.atom->aniso.u11 = cra.atom->aniso.u22 = cra.atom->aniso.u33 = 0;
    } else {
      for (gemmi::SmallStructure::Site& site : small.sites)
        site.aniso.u11 = site.aniso.u22 = site.aniso.u33 = 0;
    }
  }

  double wavelength = 0;
  if (p.options[Wavelength]) {
    wavelength = std::atof(p.options[Wavelength].arg);
  } else {
    if (use_st) {
      // reading wavelength from PDB and mmCIF files needs to be revisited
      //if (!st.crystals.empty() && !st.crystals[0].diffractions.empty())
      //  wavelength_list = st.crystals[0].diffractions[0].wavelengths;
    } else {
      wavelength = small.wavelength;
    }
  }

  if (p.options[Unknown]) {
    gemmi::El new_el = gemmi::find_element(p.options[Unknown].arg);
    if (new_el == gemmi::El::X)
      gemmi::fail("--unknown must specify chemical element symbol.");
    if (use_st) {
      for (gemmi::Chain& chain : st.models[0].chains)
        for (gemmi::Residue& residue : chain.residues)
          for (gemmi::Atom& atom : residue.atoms)
            if (atom.element == gemmi::El::X)
              atom.element = new_el;
    } else {
      for (gemmi::SmallStructure::Site& atom : small.sites)
        if (atom.element == gemmi::El::X)
          atom.element = new_el;
    }
  }

  char table = p.options[For] ? p.options[For].arg[0] : 'x';
  if (table == 'x') {
    process_with_table<gemmi::IT92<double>>(use_st, st, small, wavelength, p);
  } else if (table == 'e') {
    if (p.options[CifFp])
      gemmi::fail("Electron scattering has no dispersive part (--ciffp)");
    process_with_table<gemmi::ITC4322<double>>(use_st, st, small, 0., p);
  }
}

} // anonymous namespace

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
      process(input, p);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
