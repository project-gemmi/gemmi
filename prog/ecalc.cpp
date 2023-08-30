// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <cstdlib>  // for atof
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
    const Mtz::Column& fcol = mtz.get_column_with_label(f_label);
    if (use_sigma && !fcol.get_next_column_if_type('Q'))
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
    Mtz::Column& ecol = mtz.copy_column(e_idx, fcol, trailing_cols);
    ecol.label = e_label;
    ecol.type = 'E';
    if (use_sigma) {
      assert(ecol.get_next_column_if_type('Q') != nullptr);
      (&ecol + 1)->label = "SIG" + ecol.label;
    }
    gemmi::GroupOps gops = mtz.spacegroup->operations();
    for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size()) {
      gemmi::Miller hkl = mtz.get_hkl(n);
      int epsilon = gops.epsilon_factor(hkl);
      double f = mtz.data[n + fcol.idx];
      mtz.data[n + ecol.idx] = float(f / std::sqrt(epsilon));
    }
    // TODO
    if (verbose)
      std::fprintf(stderr, "Writing %s ...\n", output);
    mtz.write_to_file(output);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

