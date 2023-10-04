// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include "gemmi/binner.hpp"  // for Binner
#include "gemmi/mtz.hpp"     // for Mtz
#include "gemmi/ecalc.hpp"   // for Mean

#define GEMMI_PROG ecalc
#include "options.h"

namespace {

enum OptionIndex {
  LabelF=4, LabelE, NoSigma, Method, BinSize, MaxBins
};

struct EcalcArg: public Arg {
  static option::ArgStatus MethodChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {/*"ma",*/ "2", "3", "ec"});
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
    "ma - moving average (not implemented yet),\n\t"
    "2 - resolution bins equispaced in 1/d^2,\n\t"
    "3 - resolution bins in 1/d^3 (default),\n\t"
    "ec - bins with equal count of reflections." },
  { BinSize, 0, "", "binsize", Arg::Int,
    "  --binsize=N  \tNumber of reflections per bin or in moving window\n\t"
    "(default: 200)." },
  { MaxBins, 0, "", "maxbins", Arg::Int,
    "  --maxbins=N  \tMaximum number of bins (default: 100)." },
  { NoOp, 0, "", "", Arg::None,
    "\nColumn for SIGF is always the next column after F (the type must be Q)."
    "\nIf INPUT.mtz has column E, it is replaced in OUTPUT.mtz."
    "\nOtherwise, a new column is appended."
    "\nThe name for SIGE, if it is added as a column, is 'SIG' + label of E." },
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
    const Mtz::Column& fcol_before = mtz.get_column_with_label(f_label);
    if (use_sigma && !fcol_before.get_next_column_if_type('Q'))
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
    int fcol_idx = fcol_before.idx; // fcol_before gets invalidated in the next line
    Mtz::Column& ecol = mtz.copy_column(e_idx, fcol_before, trailing_cols);
    //const Mtz::Column* fcol = &mtz.columns[fcol_idx];
    ecol.label = e_label;
    ecol.type = 'E';
    if (use_sigma) {
      assert(ecol.get_next_column_if_type('Q') != nullptr);
      (&ecol + 1)->label = "SIG" + ecol.label;
    }
    if (!mtz.has_data())  // shouldn't happen
      gemmi::fail("no data");
    int f_count = 0;
    for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
      if (!std::isnan(mtz.data[n + fcol_idx]))
        ++f_count;
    int binsize = p.integer_or(BinSize, 200);
    int nbins = gemmi::iround((double)f_count / binsize);
    nbins = std::min(nbins, p.integer_or(MaxBins, 100));

    if (verbose > 0 || nbins < 3)
      std::fprintf(stderr, "%d reflections, %d Fs, %d shells\n",
                   mtz.nreflections, f_count, nbins);
    if (nbins < 3)
      gemmi::fail("not enough resolution bins");
    //bool use_moving_average = false;
    auto method = gemmi::Binner::Method::Dstar3;
    if (p.options[Method])
      switch (p.options[Method].arg[0]) {
        case '2': method = gemmi::Binner::Method::Dstar2; break;
        case 'e': method = gemmi::Binner::Method::EqualCount; break;
        //case 'm': use_moving_average = true; break;
      }
    if (verbose > 0)
      std::fprintf(stderr, "Calculating E ...\n");
    gemmi::Binner binner;
    gemmi::MtzDataProxy data_proxy{mtz};
    binner.setup(nbins, method, data_proxy, nullptr, /*with_mids=*/true, fcol_idx);
    std::vector<double> multipliers
      = gemmi::calculate_amplitude_normalizers(data_proxy, fcol_idx, binner);

    if (verbose > 0)
      std::fprintf(stderr, "Writing %s ...\n", output);
    for (size_t i = 0, n = 0; n < mtz.data.size(); n += mtz.columns.size(), ++i) {
      mtz.data[n + ecol.idx] *= float(multipliers[i]);
      if (use_sigma)
        mtz.data[n + ecol.idx + 1] *= float(multipliers[i]);
    }
    mtz.write_to_file(output);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

