// Copyright 2019 Global Phasing Ltd.
//
// Calculate structure factors from a molecular model.

#include <stdio.h>
#include <complex>
#include <gemmi/assembly.hpp>  // for expand_ncs
#include <gemmi/calculate.hpp> // for expand_box
#include <gemmi/ccp4.hpp>      // for Ccp4
#include <gemmi/fileutil.hpp>  // for file_open
#include <gemmi/fourier.hpp>
#include <gemmi/fprime.hpp>    // for cromer_liberman
#include <gemmi/gz.hpp>        // for MaybeGzipped
#include <gemmi/it92.hpp>      // for IT92
#include <gemmi/c4322.hpp>     // for C4322
#include <gemmi/neutron92.hpp> // for Neutron92
#include <gemmi/math.hpp>      // for sq
#include <gemmi/mtz.hpp>       // for Mtz
#include <gemmi/dencalc.hpp>   // for DensityCalculator
#include <gemmi/scaling.hpp>   // for Scaling
#include <gemmi/sfcalc.hpp>    // for calculate_structure_factor
#include <gemmi/smcif.hpp>     // for make_small_structure_from_block
#include <gemmi/solmask.hpp>   // for SolventMasker
#include <gemmi/read_cif.hpp>  // for read_cif_gz
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include "timer.h"             // for Timer

#define GEMMI_PROG sfcalc
#include "options.h"

namespace {

enum OptionIndex {
  Hkl=4, Dmin, For, NormalizeIt92, Rate, Blur, RCut, Test, ToMtz, Compare,
  CifFp, Wavelength, Unknown, NoAniso, Margin, ScaleTo, SigmaCutoff, FLabel,
  PhiLabel, Ksolv, Bsolv, Baniso, RadiiSet, Rprobe, Rshrink, WriteMap
};

struct SfCalcArg: public Arg {
  static option::ArgStatus FormFactors(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"xray", "electron", "neutron", "mott-bethe"});
  }
  static option::ArgStatus Radii(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"vdw", "cctbx", "refmac"});
  }

  static option::ArgStatus Float6(const option::Option& option, bool msg) {
    if (option.arg) {
      char* endptr = nullptr;
      int counter = 0;
      do {
        (void) std::strtod(endptr ? endptr + 1 : option.arg, &endptr);
        ++counter;
      } while (*endptr == ':');
      if (counter == 6 && *endptr == '\0')
        return option::ARG_OK;
    }
    if (msg)
      fprintf(stderr, "Option '%.*s' requires six colon-separated numbers "
                      "as an argument,\n for example: %.*s=2.1:3:4:0:0:0\n",
                      option.namelen, option.name, option.namelen, option.name);
    return option::ARG_ILLEGAL;
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
    "  --for=TYPE  \tTYPE is xray (default), electron, neutron or mott-bethe." },
  { NormalizeIt92, 0, "", "normalize-it92", Arg::None,
    "  --normalize-it92  \tNormalize X-ray form factors (a tiny change)." },
  { CifFp, 0, "", "ciffp", Arg::None,
    "  --ciffp  \tRead f' from _atom_type_scat_dispersion_real in CIF." },
  { Wavelength, 0, "w", "wavelength", Arg::Float,
    "  --wavelength=NUM  \tWavelength [A] for calculation of f' "
    "(use --wavelength=0 or -w0 to ignore anomalous scattering)." },
  { Unknown, 0, "", "unknown", Arg::Required,
    "  --unknown=SYMBOL  \tUse form factor of SYMBOL for unknown atoms." },
  { NoAniso, 0, "", "noaniso", Arg::None,
    "  --noaniso  \tIgnore anisotropic ADPs." },
  { Margin, 0, "", "margin", Arg::Float,
    "  --margin=NUM  \tFor non-crystal use bounding box w/ margin (default: 10)." },

  { NoOp, 0, "", "", Arg::None,
    "\nOptions for density and FFT calculations (with --dmin):" },
  { Rate, 0, "", "rate", Arg::Float,
    "  --rate=NUM  \tShannon rate used for grid spacing (default: 1.5)." },
  { Blur, 0, "", "blur", Arg::Float,
    "  --blur=NUM  \tB added for Gaussian blurring (default: auto)." },
  { RCut, 0, "", "rcut", Arg::Float,
    "  --rcut=Y  \tUse atomic radius r such that rho(r) < Y (default: 1e-5)." },
  { Test, 0, "", "test", Arg::Optional,
    "  --test[=CACHE]  \tCalculate exact values and report differences (slow)." },
  { WriteMap, 0, "", "write-map", Arg::Required,
    "  --write-map=FILE  \tWrite density (excl. bulk solvent) as CCP4 map." },
  { ToMtz, 0, "", "to-mtz", Arg::Required,
    "  --to-mtz=FILE  \tWrite Fcalc to a new MTZ file." },

  { NoOp, 0, "", "", Arg::None, "\nOptions for anisotropic scaling (only w/ FFT):" },
  { ScaleTo, 0, "", "scale-to", Arg::Required,
    "  --scale-to=FILE:COL  \tAnisotropic scaling to F from MTZ file."
    "\n\tArgument: FILE[:FCOL[:SIGFCOL]] (defaults: F and SIGF)." },
  { SigmaCutoff, 0, "", "sigma-cutoff", Arg::Float,
    "  --sigma-cutoff=NUM  \tUse only data with F/SIGF > NUM (default: 0)." },
  // TODO: solvent option: mask, babinet, none

  { NoOp, 0, "", "", Arg::None, "\nOptions for bulk solvent correction (only w/ FFT):" },
  { RadiiSet, 0, "", "radii-set", SfCalcArg::Radii,
    "  --radii-set=SET  \tSet of per-element radii, one of: vdw, cctbx, refmac." },
  { Rprobe, 0, "", "r-probe", Arg::Float,
    "  --r-probe=NUM  \tValue added to VdW radius (default: 1.0A)." },
  { Rshrink, 0, "", "r-shrink", Arg::Float,
    "  --r-shrink=NUM  \tValue for shrinking the solvent area (default: 1.1A)." },
  { Ksolv, 0, "", "ksolv", Arg::Float,
    "  --ksolv=NUM  \tValue (if optimizing: initial value) of k_solv." },
  { Bsolv, 0, "", "bsolv", Arg::Float,
    "  --bsolv=NUM  \tValue (if optimizing: initial value) of B_solv." },
  { Baniso, 0, "", "baniso", SfCalcArg::Float6,
    "  --baniso=B11:...:B23 \tAnisotropic scale matrix (6 colon-separated numbers: "
      "B11, B22, B33, B12, B13, B23)." },

  { NoOp, 0, "", "", Arg::None,
    "\nOptions for comparing calculated values with values from a file:" },
  { Compare, 0, "", "compare", Arg::Required,
    "  --compare=FILE  \tRe-calculate Fcalc and report differences." },
  { FLabel, 0, "", "f", Arg::Required,
    "  --f=LABEL  \tMTZ column label (default: FC) or small molecule cif"
    " tag (default: F_calc or F_squared_calc)." },
  { PhiLabel, 0, "", "phi", Arg::Required,
    "  --phi=LABEL  \tMTZ column label (default: PHIC)." },
  { 0, 0, 0, 0, 0, 0 }
};

struct RefFile {
  enum class Mode { None, Test, Compare, WriteMtz };
  Mode mode = Mode::None;
  const char* path = nullptr;
  std::string f_label;
  std::string phi_label;
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

template<typename Real>
void write_asudata_to_mtz(const gemmi::AsuData<std::complex<Real>>& asu_data,
                          const RefFile& file) {
  std::unique_ptr<gemmi::Mtz> output_mtz(new gemmi::Mtz(true));
  output_mtz->set_cell_for_all(asu_data.unit_cell());
  output_mtz->spacegroup = asu_data.spacegroup();
  output_mtz->add_dataset("calculated");
  output_mtz->add_column(file.f_label, 'F', -1, -1, false);
  output_mtz->add_column(file.phi_label, 'P', -1, -1, false);
  output_mtz->title = "Fcalc generated by gemmi";
  output_mtz->nreflections = (int) asu_data.v.size();
  output_mtz->data.reserve(5 * asu_data.v.size());
  for (const gemmi::HklValue<std::complex<Real>>& hv : asu_data.v) {
    for (int i = 0; i != 3; ++i)
      output_mtz->data.push_back((float) hv.hkl[i]);
    output_mtz->data.push_back((float) std::abs(hv.value));
    output_mtz->data.push_back((float) gemmi::phase_in_angles(hv.value));
  }
  output_mtz->write_to_file(file.path);
}

template<typename Table, typename Real>
void process_with_fft(const gemmi::Structure& st,
                      gemmi::DensityCalculator<Table, Real>& dencalc,
                      bool mott_bethe,
                      const gemmi::SolventMasker& masker,
                      gemmi::Scaling<Real>& scaling,
                      bool verbose, const RefFile& file,
                      const gemmi::AsuData<gemmi::ValueSigma<Real>>& scale_to,
                      const char* map_file) {
  // prepare electron density map
  if (verbose) {
    fprintf(stderr, "Preparing electron density on a grid...\n");
    fflush(stderr);
  }
  Timer timer(verbose);
  timer.start();
  dencalc.set_grid_cell_and_spacegroup(st);
  dencalc.put_model_density_on_grid(st.models[0]);
  timer.print("...took");
  if (map_file) {
    gemmi::Ccp4<Real> ccp4;
    ccp4.grid = dencalc.grid;
    ccp4.update_ccp4_header(2);
    ccp4.write_ccp4_map(map_file);
  }

  if (verbose) {
#if GEMMI_COUNT_DC
    fprintf(stderr, "Density-points calculated: %zu (avg per atom: %g)\n",
            dencalc.density_computations,
            double(dencalc.density_computations) / dencalc.atoms_added);
#endif
    fprintf(stderr, "FFT of grid %d x %d x %d\n",
            dencalc.grid.nu, dencalc.grid.nv, dencalc.grid.nw);
    fflush(stderr);
    timer.start();
  }
  gemmi::FPhiGrid<Real> sf = transform_map_to_f_phi(dencalc.grid, /*half_l=*/true);
  if (verbose) {
    timer.print("...took");
    fprintf(stderr, "Preparing results...\n");
    fflush(stderr);
  }
  gemmi::StructureFactorCalculator<Table> calc(st.cell);
  calc.addends = dencalc.addends;
  gemmi::fileptr_t cache(nullptr, nullptr);
  gemmi::AsuData<std::complex<double>> compared_data;
  if (file.path) {
    if (file.mode == RefFile::Mode::Test) {
      cache = gemmi::file_open(file.path, "r");
    } else if (file.mode == RefFile::Mode::Compare) {
      gemmi::Mtz mtz;
      mtz.read_input(gemmi::MaybeGzipped(file.path), true);
      compared_data.load_values<2>(gemmi::MtzDataProxy{mtz}, {file.f_label, file.phi_label});
    }
  }
  auto asu_data = sf.prepare_asu_data(dencalc.d_min, dencalc.blur, false, false, mott_bethe);

  gemmi::AsuData<std::complex<Real>> mask_data;
  if (scaling.use_solvent) {
    // uses scaling.grid as a temporary array
    masker.put_mask_on_grid(dencalc.grid, st.models[0]);
    mask_data = transform_map_to_f_phi(dencalc.grid, /*half_l=*/true)
                .prepare_asu_data(dencalc.d_min, 0);
  }

  if (scale_to.size() != 0) {
    scaling.prepare_points(asu_data, scale_to, &mask_data);
    printf("Calculating scale factors using %zu points...\n", scaling.points.size());
    scaling.fit_isotropic_b_approximately();
    //fprintf(stderr, "k_ov=%g B_ov=%g\n", scaling.k_overall, scaling.get_b_overall().u11);
    scaling.fit_parameters();
    gemmi::SMat33<double> b_aniso = scaling.get_b_overall();
    if (scaling.use_solvent)
      fprintf(stderr, "Bulk solvent parameters: k_sol=%g B_sol=%g\n",
              scaling.k_sol, scaling.b_sol);
    fprintf(stderr, "k_ov=%g B11=%g B22=%g B33=%g B12=%g B13=%g B23=%g\n",
            scaling.k_overall, b_aniso.u11, b_aniso.u22, b_aniso.u33,
                               b_aniso.u12, b_aniso.u13, b_aniso.u23);
    if (verbose) {
      std::vector<double> computed = scaling.compute_values();
      Comparator comparator;
      for (size_t i = 0; i != scaling.points.size(); ++i)
        comparator.add(computed[i], (double)scaling.points[i].fobs);
      fprintf(stderr, "After scaling: ");
      print_to_stderr(comparator);
      fprintf(stderr, "\n");
    }
  }
  scaling.scale_data(asu_data, &mask_data);

  if (file.mode == RefFile::Mode::WriteMtz) {
    write_asudata_to_mtz(asu_data, file);
  } else if (file.mode == RefFile::Mode::None) {
    for (gemmi::HklValue<std::complex<Real>>& hv : asu_data.v)
      print_sf(hv.value, hv.hkl);
  } else {
    Comparator comparator;
    for (gemmi::HklValue<std::complex<Real>>& hv : asu_data.v) {
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
        if (mott_bethe)
          exact *= calc.mott_bethe_factor();
      }
      comparator.add_complex(hv.value, exact);
      printf(" (%d %d %d)\t%7.2f\t%8.3f \t%6.2f\t%7.3f\td=%5.2f\n",
             hv.hkl[0], hv.hkl[1], hv.hkl[2], std::abs(hv.value), std::abs(exact),
             gemmi::phase_in_angles(hv.value), gemmi::phase_in_angles(exact),
             sf.unit_cell.calculate_d(hv.hkl));
    }
    print_to_stderr(comparator);
    fprintf(stderr, "  <dPhi>=%#.4g", comparator.mean_dphi());
    if (!verbose)
      fprintf(stderr, "   %#.5gs", timer.count());
    fprintf(stderr, "\n");
  }
}

template<typename Table>
void print_structure_factors_sm(const gemmi::SmallStructure& small,
                                gemmi::StructureFactorCalculator<Table>& calc,
                                bool mott_bethe, double d_min, bool verbose,
                                const RefFile& file) {
  Timer timer(verbose);
  timer.start();
  int counter = 0;
  // cf. prepare_asu_data()
  double max_1_d = 1. / d_min;
  int max_h = int(max_1_d / small.cell.ar);
  int max_k = int(max_1_d / small.cell.br);
  int max_l = int(max_1_d / small.cell.cr);
  const gemmi::SpaceGroup* sg = small.find_spacegroup();
  if (!sg)
    sg = &gemmi::get_spacegroup_p1();
  gemmi::ReciprocalAsu asu(sg);
  gemmi::AsuData<std::complex<double>> asu_data;
  gemmi::GroupOps gops = sg->operations();
  for (int h = -max_h; h <= max_h; ++h)
    for (int k = -max_k; k <= max_k; ++k)
      for (int l = 0; l <= max_l; ++l) {
        gemmi::Miller hkl{{h, k, l}};
        if (!asu.is_in(hkl) || (hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0))
          continue;
        if (gops.is_systematically_absent(hkl))
          continue;
        double hkl_1_d2 = small.cell.calculate_1_d2(hkl);
        if (hkl_1_d2 < max_1_d * max_1_d) {
          auto value = calc.calculate_sf_from_small_structure(small, hkl);
          if (mott_bethe)
            value *= calc.mott_bethe_factor();
          if (file.mode == RefFile::Mode::WriteMtz)
            asu_data.v.push_back({hkl, value});
          else
            print_sf(value, hkl);
          ++counter;
        }
      }
  if (verbose) {
    fflush(stdout);
    fprintf(stderr, "Calculated %d SFs in %g s.\n", counter, timer.count());
    fflush(stderr);
  }
  if (file.mode == RefFile::Mode::WriteMtz) {
    asu_data.unit_cell_ = small.cell;
    asu_data.spacegroup_ = sg;
    write_asudata_to_mtz(asu_data, file);
  }
}

template<typename Table>
void compare_with_hkl(const gemmi::SmallStructure& small,
                      gemmi::StructureFactorCalculator<Table>& calc,
                      const RefFile& file,
                      bool verbose,
                      Comparator& comparator,
                      bool mott_bethe) {
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
    } catch(std::exception& e) {
      fprintf(stderr, "Error in _refln_[] in %s: %s\n", file.path, e.what());
      continue;
    }
    double f = std::abs(calc.calculate_sf_from_small_structure(small, hkl));
    if (mott_bethe)
      f *= calc.mott_bethe_factor();
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
                      const RefFile& file, bool verbose, Comparator& comparator,
                      bool mott_bethe) {
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
    if (mott_bethe)
      f *= calc.mott_bethe_factor();
    comparator.add(f_from_file, f);
    if (verbose)
      printf(" (%d %d %d)\t%7.2f\t%8.3f \td=%5.2f\n",
             hkl[0], hkl[1], hkl[2], f_from_file, f, cell.calculate_d(hkl));
  }
}

template<typename Table>
void process_with_table(bool use_st, gemmi::Structure& st, const gemmi::SmallStructure& small,
                        double wavelength, bool mott_bethe, const OptParser& p) {
  const gemmi::UnitCell& cell = use_st ? st.cell : small.cell;
  gemmi::StructureFactorCalculator<Table> calc(cell);

  // assign f' given explicitly in a file
  if (p.options[CifFp]) {
    if (use_st) {
      // _atom_type.scat_dispersion_real is almost never used,
      // so for now we ignore it.
    } else { // small molecule
      if (p.options[Verbose])
        fprintf(stderr, "Using f' read from cif file (%u atom types)\n",
                (unsigned) small.atom_types.size());
      for (const gemmi::SmallStructure::AtomType& atom_type : small.atom_types)
        calc.addends.set(atom_type.element, (float)atom_type.dispersion_real);
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
      if (present_elems[z] && calc.addends.values[z] == 0) {
        calc.addends.values[z] = (float) gemmi::cromer_liberman(z, energy, nullptr);
      }
  }
  if (mott_bethe)
    calc.addends.subtract_z();

  // handle option --hkl
  for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
    std::vector<int> hkl_ = parse_comma_separated_ints(opt->arg);
    gemmi::Miller hkl{{hkl_[0], hkl_[1], hkl_[2]}};
    if (p.options[Verbose])
      fprintf(stderr, "hkl=(%d %d %d) -> d=%g\n", hkl[0], hkl[1], hkl[2],
              cell.calculate_d(hkl));
    std::complex<double> sf;
    if (use_st)
      sf = calc.calculate_sf_from_model(st.models[0], hkl);
    else
      sf = calc.calculate_sf_from_small_structure(small, hkl);
    if (mott_bethe)
      sf *= calc.mott_bethe_factor();
    print_sf(sf, hkl);
  }

  RefFile file;
  if (p.options[Test]) {
    file.mode = RefFile::Mode::Test;
    file.path = p.options[Test].arg;
  } else if (p.options[Compare]) {
    file.mode = RefFile::Mode::Compare;
    file.path = p.options[Compare].arg;
  } else if (p.options[ToMtz]) {
    file.mode = RefFile::Mode::WriteMtz;
    file.path = p.options[ToMtz].arg;
  }
  if (p.options[FLabel])
    file.f_label = p.options[FLabel].arg;
  // we leave f_label empty for use in compare_with_hkl()
  else if (use_st || file.mode == RefFile::Mode::WriteMtz)
    file.f_label = "FC";
  if (p.options[PhiLabel])
   file.phi_label = p.options[PhiLabel].arg;
  else
    file.phi_label = "PHIC";

  using Real = float;
  gemmi::AsuData<gemmi::ValueSigma<Real>> scale_to;
  if (p.options[ScaleTo]) {
    std::string path = p.options[ScaleTo].arg;
    std::string flabel = "F";
    std::string siglabel = "SIGF";
    size_t sep2 = path.rfind(':');
    if (sep2 != std::string::npos && sep2 != 0) {
      size_t sep = path.rfind(':', sep2 - 1);
      if (sep == std::string::npos)
        std::swap(sep, sep2);
      flabel = path.substr(sep+1, sep2 - (sep+1));
      if (sep2 != std::string::npos)
        siglabel = path.substr(sep2+1);
      path.resize(sep);
    }
    double sigma_cutoff = 0;
    if (p.options[SigmaCutoff])
      sigma_cutoff = std::atof(p.options[SigmaCutoff].arg);
    gemmi::Mtz mtz;
    mtz.read_input(gemmi::MaybeGzipped(path), true);
    if (siglabel.empty()) {
      scale_to.load_values<2>(gemmi::MtzDataProxy{mtz}, {flabel, flabel});
      for (auto& hkl_value : scale_to.v)
        hkl_value.value.sigma = std::sqrt(hkl_value.value.sigma);
    } else {
      scale_to.load_values<2>(gemmi::MtzDataProxy{mtz}, {flabel, siglabel});
      size_t size_before = scale_to.size();
      vector_remove_if(scale_to.v, [=](const gemmi::HklValue<gemmi::ValueSigma<Real>>& x) {
          return x.value.value <= sigma_cutoff * x.value.sigma;
      });
      if (p.options[Verbose])
        fprintf(stderr, "Sigma cutoff (F/sigF > %g) excluded %zu out of %zu points.\n",
                sigma_cutoff, size_before - scale_to.size(), size_before);
    }
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
        dencalc.cutoff = (float) std::atof(p.options[RCut].arg);
      dencalc.addends = calc.addends;
      if (p.options[Blur]) {
        dencalc.blur = std::atof(p.options[Blur].arg);
      } else if (dencalc.rate < 3) {
        // ITfC vol B section 1.3.4.4.5 has formula
        // B = log Q / (sigma * (sigma - 1) * d*_max ^2)
        // where Q is quality factor, sigma is the oversampling rate.
        // This value is not optimal.
        // The optimal value would depend on the distribution of B-factors
        // and on the atomic cutoff radius, and probably it would be too
        // hard to estimate. Here we use the same formula as in Refmac.
        dencalc.set_refmac_compatible_blur(st.models[0]);
        if (p.options[Verbose])
          fprintf(stderr, "B_min=%g, B_add=%g\n",
                  gemmi::get_minimum_b(st.models[0]), dencalc.blur);
      }
      gemmi::AtomicRadiiSet radii_choice = gemmi::AtomicRadiiSet::VanDerWaals;
      if (p.options[RadiiSet]) {
        char c = p.options[RadiiSet].arg[0];
        if (c == 'v')
          radii_choice = gemmi::AtomicRadiiSet::VanDerWaals;
        else if (c == 'c')
          radii_choice = gemmi::AtomicRadiiSet::Cctbx;
        else if (c == 'r')
          radii_choice = gemmi::AtomicRadiiSet::Refmac;
      }
      gemmi::SolventMasker masker(radii_choice);
      if (p.options[Rprobe])
        masker.rprobe = std::atof(p.options[Rprobe].arg);
      if (p.options[Rshrink])
        masker.rshrink = std::atof(p.options[Rshrink].arg);

      gemmi::Scaling<Real> scaling(cell, st.find_spacegroup());
      if (p.options[Ksolv] || p.options[Bsolv] || scale_to.size() != 0) {
        scaling.use_solvent = true;
        if (p.options[Ksolv])
          scaling.k_sol = std::atof(p.options[Ksolv].arg);
        if (p.options[Bsolv])
          scaling.b_sol = std::atof(p.options[Bsolv].arg);
      }
      if (p.options[Baniso]) {
        char* endptr = nullptr;
        gemmi::SMat33<double> b_aniso;
        b_aniso.u11 = std::strtod(p.options[Baniso].arg, &endptr);
        b_aniso.u22 = std::strtod(endptr + 1, &endptr);
        b_aniso.u33 = std::strtod(endptr + 1, &endptr);
        b_aniso.u12 = std::strtod(endptr + 1, &endptr);
        b_aniso.u13 = std::strtod(endptr + 1, &endptr);
        b_aniso.u23 = std::strtod(endptr + 1, &endptr);
        scaling.set_b_overall(b_aniso);
      }
      const char* map_file = p.options[WriteMap] ? p.options[WriteMap].arg : nullptr;
      process_with_fft(st, dencalc, mott_bethe, masker, scaling,
                       p.options[Verbose], file, scale_to, map_file);
    } else {
      if (p.options[Rate] || p.options[RCut] || p.options[Blur] ||
          p.options[Test])
        gemmi::fail("Small molecule SFs are calculated directly. Do not use any\n"
                    "of the FFT-related options: --rate, --blur, --rcut, --test.");
      print_structure_factors_sm(small, calc, mott_bethe, d_min, p.options[Verbose], file);
    }

  // handle option --compare
  } else if (file.mode == RefFile::Mode::Compare) {
    Comparator comparator;
    if (use_st)
      compare_with_mtz(st.models[0], st.cell, calc, file, p.options[Verbose],
                       comparator, mott_bethe);
    else
      compare_with_hkl(small, calc, file, p.options[Verbose], comparator, mott_bethe);
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
    small.change_occupancies_to_crystallographic();
  }

  if (!p.options[Dmin]) {
    for (OptionIndex opt : {ToMtz, WriteMap, ScaleTo, Ksolv, Bsolv,
                            RadiiSet, Rprobe, Rshrink})
      if (p.options[opt])
        gemmi::fail("Option ", p.options[opt].name, " works only with --dmin");
  }
  if (!use_st) {
    for (OptionIndex opt : {WriteMap, ScaleTo, Ksolv, Bsolv, RadiiSet, Rprobe, Rshrink})
      if (p.options[opt])
        gemmi::fail("Option ", p.options[opt].name, " is only used in density-FFT "
                    "route which is only used for macromolecular structures");
  }

  if (p.options[NoAniso]) {
    if (use_st) {
      for (gemmi::CRA cra : st.models[0].all())
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
  if (wavelength < 0)
    gemmi::fail("wavelength should not be negative");


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

  if (use_st && st.ncs_not_expanded())
    gemmi::expand_ncs(st, gemmi::HowToNameCopiedChain::Dup);

  if (use_st && !st.cell.is_crystal()) {
    double margin = 10.;
    if (p.options[Margin])
      margin = std::atof(p.options[Margin].arg);
    gemmi::Box<gemmi::Position> box;
    expand_box(st.models[0], box);
    gemmi::Position size = box.get_size();
    st.cell.set(size.x + margin, size.y + margin, size.z + margin, 90, 90, 90);
    if (p.options[Verbose]) {
      fprintf(stderr, "Unit cell set to %g x %g x %g\n", st.cell.a, st.cell.b, st.cell.c);
      fflush(stderr);
    }
  }

  char table = p.options[For] ? p.options[For].arg[0] : 'x';
  if (p.options[CifFp] && table != 'x')
    gemmi::fail("Electron scattering has no dispersive part (--ciffp)");
  if (table == 'x' || table == 'm') {
    if (p.options[NormalizeIt92])
      gemmi::IT92<float>::normalize();
    if (table == 'm')
      gemmi::IT92<float>::ignore_charge = true;
    process_with_table<gemmi::IT92<float>>(use_st, st, small, wavelength,
                                           table == 'm', p);
  } else if (table == 'e') {
    process_with_table<gemmi::C4322<float>>(use_st, st, small, 0., false, p);
  } else if (table == 'n') {
    process_with_table<gemmi::Neutron92<double>>(use_st, st, small, 0., false, p);
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.check_exclusive_group({ToMtz, Test, Compare});
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
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
