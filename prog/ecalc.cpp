// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <cstdlib>  // for atof
#include "gemmi/binner.hpp"  // for Binner
#include "gemmi/mtz.hpp"     // for Mtz
#include "gemmi/stats.hpp"   // for Mean

#define GEMMI_PROG ecalc
#include "options.h"

namespace {

enum OptionIndex {
  LabelF=4, LabelE, NoSigma, Method, BinSize
};

struct EcalcArg: public Arg {
  static option::ArgStatus MethodChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"ma", "2", "3", "ec"});
  }
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT.mtz OUTPUT.mtz"
    "\n\nCalculates E (normalised structure amplitudes) from F." },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { LabelF, 0, "", "F", Arg::Required,
    "  --F=LABEL  \tLabel of the input column (default: F)." },
  { LabelE, 0, "", "E", Arg::Required,
    "  --E=LABEL  \tLabel of the output column (default: E)." },
  { NoSigma, 0, "", "no-sigma", Arg::Required,
    "  --no-sigma  \tDo not use sigma columns (SIGF->SIGE)." },
  { Method, 0, "", "method", EcalcArg::MethodChoice,
    "  --method=METHOD  \tMethod of calculating <F^2>, one of:\n\t"
    "ma - moving average,\n\t"
    "2 - resolution bins equispaced in 1/d^2,\n\t"
    "3 - resolution bins in 1/d^3 (default),\n\t"
    "ec - bins with equal count of reflections." },
  { BinSize, 0, "n", "binsize", Arg::Int,
    "  -n, --binsize=N  \tNumber of reflections per bin or in moving window\n\t"
    "(default: 200)." },
  { NoOp, 0, "", "", Arg::None,
    "\nColumn for SIGF is always the next column after F (the type must be Q)."
    "\nIf INPUT.mtz has column E, it is replaced in OUTPUT.mtz."
    "\nOtherwise, a new column is appended."
    "\nThe name for SIGE, if it is added as a column, is 'SIG' + label of E."
    "\nTo print the list of resolution shells add -vv." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  using gemmi::Mtz;
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);
  const char* f_label = p.options[LabelF] ? p.options[LabelF].arg : "F";
  const char* e_label = p.options[LabelE] ? p.options[LabelE].arg : "E";
  bool use_sigma = !p.options[NoSigma];
  int verbose = p.options[Verbose].count();
  if (verbose > 0)
    std::fprintf(stderr, "Reading %s ...\n", input);
  try {
    Mtz mtz;
    mtz.read_file_gz(input);
    if (mtz.spacegroup == nullptr)
      gemmi::fail("Unknown space group of the input data: ", input);
    const Mtz::Column* fcol = &mtz.get_column_with_label(f_label);
    if (use_sigma && !fcol->get_next_column_if_type('Q'))
        gemmi::fail("Column ", f_label, " not followed by sigma column. Use --no-sigma.\n");
    int e_idx = -1;
    if (Mtz::Column* e = mtz.column_with_label(e_label)) {
      if (e->type != 'E')
        gemmi::fail("Column for E exists, but is not of type E");
      if (use_sigma && !e->get_next_column_if_type('Q'))
        gemmi::fail("Column ", e_label, " not followed by sigma column. Use --no-sigma.\n");
      e_idx = (int)e->idx;
    }
    if (verbose > 0)
      std::fprintf(stderr, "%s column %s ...\n",
                   (e_idx >= 0 ? "Replacing existing" : "Adding"), e_label);
    std::vector<std::string> trailing_cols(use_sigma ? 1 : 0);
    int fcol_idx = fcol->idx; // fcol gets invalidated in the next line
    Mtz::Column& ecol = mtz.copy_column(e_idx, *fcol, trailing_cols);
    fcol = &mtz.columns[fcol_idx];
    ecol.label = e_label;
    ecol.type = 'E';
    if (use_sigma) {
      assert(ecol.get_next_column_if_type('Q') != nullptr);
      (&ecol + 1)->label = "SIG" + ecol.label;
    }
    if (!mtz.has_data())  // shouldn't happen
      gemmi::fail("no data");
    gemmi::Binner binner;
    int f_count = 0;
    for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
      if (!std::isnan(mtz.data[n + fcol_idx]))
        ++f_count;
    int binsize = 200;
    if (p.options[BinSize])
      binsize = std::atoi(p.options[BinSize].arg);
    int nbins = gemmi::iround((double)f_count / binsize);
    if (verbose > 0)
      std::fprintf(stderr, "%d reflections, %d Fs, %d shells\n",
                   mtz.nreflections, f_count, nbins);
    bool use_moving_average = false;
    auto method = gemmi::Binner::Method::Dstar3;
    if (p.options[Method])
      switch (p.options[Method].arg[0]) {
        case '2': method = gemmi::Binner::Method::Dstar3; break;
        case 'e': method = gemmi::Binner::Method::EqualCount; break;
        case 'm': use_moving_average = true; break;
      }
    if (verbose > 0)
      std::fprintf(stderr, "Calculating E ...\n");
    gemmi::GroupOps gops = mtz.spacegroup->operations();
    std::vector<double> computed_values(mtz.nreflections, NAN);
    if (use_moving_average) {
      //TODO
    } else {
      std::vector<double> mids = binner.setup_mid(nbins, method, gemmi::MtzDataProxy{mtz});
      std::vector<double> inv_d2(mtz.nreflections);
      for (size_t i = 0, n = 0; n < mtz.data.size(); n += mtz.columns.size(), ++i)
        inv_d2[i] = mtz.cell.calculate_1_d2(mtz.get_hkl(n));
      std::vector<int> bin_index = binner.get_bins_from_1_d2(inv_d2);
      std::vector<gemmi::Mean> denom(binner.size());
      for (size_t i = 0, n = 0; n < mtz.data.size(); n += mtz.columns.size(), i++) {
        gemmi::Miller hkl = mtz.get_hkl(n);
        double f = mtz.data[n + fcol_idx];
        if (!std::isnan(f)) {
          double inv_epsilon = 1.0 / gops.epsilon_factor(hkl);
          double f2 = f * f * inv_epsilon;
          computed_values[i] = std::sqrt(inv_epsilon);
          denom[bin_index[i]].add_point(f2);
        }
      }
      if (verbose > 1) {
        // print shell statistics
        std::vector<int> refl_counts(binner.size());
        printf(" shell\t    #F\t    d\t <F^2>\t  #refl\t mid d\n");
        for (int idx : bin_index)
          ++refl_counts[idx];
        for (size_t i = 0; i < binner.size(); ++i) {
          double d = 1 / std::sqrt(binner.limits[i]);
          double mid_d = 1 / std::sqrt(mids[i]);
          printf("%6zu\t%6d\t%6.2f\t%7.1f\t%6d\t%6.2f\n",
                 i+1, denom[i].n, d, denom[i].get_mean(), refl_counts[i], mid_d);
        }
        printf("\n");
      }
      // re-use Mean::sum for RMS
      for (gemmi::Mean& m : denom)
        m.sum = std::sqrt(m.get_mean());
      for (size_t i = 0; i < computed_values.size(); ++i) {
        double x = inv_d2[i];
        int bin = bin_index[i];
        double rms = denom[bin].sum;
        if (x > mids.front() && x < mids.back()) {
          // linear interpolation in 1/d^2
          if (x > mids[bin])
            ++bin;
          double x0 = mids[bin - 1];
          double x1 = mids[bin];
          double y0 = denom[bin - 1].sum;
          double y1 = denom[bin].sum;
          assert(x0 <= x && x <= x1);
          rms = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
        }
        computed_values[i] /= rms;
      }
    }
    if (verbose > 0)
      std::fprintf(stderr, "Writing %s ...\n", output);
    for (size_t i = 0, n = 0; n < mtz.data.size(); n += mtz.columns.size(), ++i) {
      mtz.data[n + ecol.idx] *= float(computed_values[i]);
      if (use_sigma)
        mtz.data[n + ecol.idx + 1] *= float(computed_values[i]);
    }
    mtz.write_to_file(output);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

