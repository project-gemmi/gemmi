// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <cstdlib>  // for atof
#include "gemmi/binner.hpp"
#include "gemmi/mtz.hpp"

#define GEMMI_PROG ecalc
#include "options.h"

namespace {

enum OptionIndex {
  LabelF=4, LabelE, NoSigma
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
  bool verbose = p.options[Verbose];
  if (verbose)
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
    if (verbose)
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
    int nbins = gemmi::iround(f_count / 200.);
    if (verbose)
      std::fprintf(stderr, "%d reflections, %d Fs, %d shells\n",
                   mtz.nreflections, f_count, nbins);
    binner.setup(nbins, gemmi::Binner::Method::Dstar2, gemmi::MtzDataProxy{mtz});
    std::vector<int> bin_index = binner.get_bins(gemmi::MtzDataProxy{mtz});
    gemmi::GroupOps gops = mtz.spacegroup->operations();
    std::vector<double> denom(binner.size(), 0.);
    std::vector<int> count(denom.size(), 0);
    if (verbose)
      std::fprintf(stderr, "Calculating E ...\n");
    auto bin_num = bin_index.begin();
    for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size(), bin_num++) {
      gemmi::Miller hkl = mtz.get_hkl(n);
      double f = mtz.data[n + fcol_idx];
      if (std::isnan(f))
        continue;
      double f2 = f * f / gops.epsilon_factor(hkl);
      mtz.data[n + ecol.idx] = std::sqrt(float(f2));
      denom[*bin_num] += f2;
      count[*bin_num]++;
    }
    for (size_t i = 0; i < denom.size(); ++i)
      denom[i] = std::sqrt(denom[i] / count[i]);
    for (size_t i = 0, n = 0; n < mtz.data.size(); n += mtz.columns.size(), ++i)
      mtz.data[n + ecol.idx] /= (float) denom[bin_index[i]];
    if (verbose)
      std::fprintf(stderr, "Writing %s ...\n", output);
    mtz.write_to_file(output);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

