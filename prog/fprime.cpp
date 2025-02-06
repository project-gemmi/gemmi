// Copyright 2019 Global Phasing Ltd.

#include "gemmi/elem.hpp"    // for Element, find_element
#include "gemmi/fprime.hpp"  // for cromer_liberman_for_array
#include "gemmi/math.hpp"    // for hc
#include <cstdlib>           // for atof
#include <stdio.h>

#define GEMMI_PROG fprime
#include "options.h"

namespace {

enum OptionIndex { Energy=4, Wavelen, Step, Nvalues };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] ELEMENT[...]"
    "\n\nPrints anomalous scattering factors f' and f\"."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  //CommonUsage[Verbose],
  { Energy, 0, "e", "energy", Arg::NumberOrRange,
    "  -e, --energy=ENERGY  \tEnergy [eV] or range of energies (e.g. 8000:14000)." },
  { Wavelen, 0, "w", "wavelength", Arg::NumberOrRange,
    "  -w, --wavelength=LAMBDA  \tWavelength [A] or range (e.g. 0.5:0.9)." },
  { Step, 0, "s", "step", Arg::Float,
    "  -s, --step=STEP  \tStep size for a range." },
  { Nvalues, 0, "n", "", Arg::Int,
    "  -n N  \tNumber of values in a range." },
  { NoOp, 0, "", "", Arg::None,
    "\nOptions -e/-w can be given multiple times:"
    "\n  " EXE_NAME " -e 12400 -e 11500 -e 9800 Se"
    "\nIf -e or -w specifies a range, -n or -s must be provided, e.g.:"
    "\n  " EXE_NAME " -w 0.5:0.9 -s 0.01 Os"
  },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (!p.options[Energy] && !p.options[Wavelen]) {
    fprintf(stderr, "Neither energy nor wavelength was specified.\n");
    return -1;
  }
  if (p.nonOptionsCount() == 0) {
    fprintf(stderr, "No elements given.\n");
    return -1;
  }
  p.check_exclusive_pair(Step, Nvalues);
  bool range_error = false;
  auto expand_range = [&](double xmin, double xmax) {
    std::vector<double> result;
    int nvalues = 0;
    double step = 0.;
    if (p.options[Step]) {
      step = std::atof(p.options[Step].arg);
      if ((xmax - xmin) * step < 0)
        step = -step;
      nvalues = int((xmax - xmin) / step + 1.01);
    } else if (p.options[Nvalues]) {
      nvalues = std::atoi(p.options[Nvalues].arg);
      step = (xmax - xmin) / (nvalues - 1);
    } else {
      fprintf(stderr, "Specifying range requires either -s or -n.\n");
      range_error = true;
    }
    if (nvalues > 0) {
      result.reserve(nvalues);
      for (int i = 0; i < nvalues; ++i)
        result.push_back(xmin + i * step);
    }
    return result;
  };

  std::vector<double> energies;
  for (const option::Option* opt = p.options[Energy]; opt; opt = opt->next()) {
    double xmin, xmax;
    parse_number_or_range(opt->arg, &xmin, &xmax);
    if (xmin == xmax) {
      energies.push_back(xmin);
    } else {
      for (double energy : expand_range(xmin, xmax))
        energies.push_back(energy);
      if (range_error)
        return -1;
    }
  }
  double hc = gemmi::hc();
  for (const option::Option* opt = p.options[Wavelen]; opt; opt = opt->next()) {
    double xmin, xmax;
    parse_number_or_range(opt->arg, &xmin, &xmax);
    if (xmin == xmax) {
      energies.push_back(hc / xmin);
    } else {
      for (double wavelength : expand_range(xmin, xmax))
        energies.push_back(hc / wavelength);
      if (range_error)
        return -1;
    }
  }
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* name = p.nonOption(i);
    gemmi::Element elem = gemmi::find_element(name);
    if (elem == gemmi::El::X) {
      fprintf(stderr, "Error: element name not recognized: '%s'\n", name);
      return -1;
    }
    std::vector<double> fp(energies.size(), 0);
    std::vector<double> fpp(energies.size(), 0);
    if (i == 0)
      printf("Element\t E[eV]\tWavelength[A]\t   f'   \t  f\"\n");
    gemmi::cromer_liberman_for_array(elem.atomic_number(),
                                     (int) energies.size(), energies.data(),
                                     &fp[0], &fpp[0]);
    for (size_t j = 0; j != energies.size(); ++j) {
      printf("%s\t%#g\t% -8g\t%8.5g\t%.5g\n",
             elem.name(), energies[j], hc / energies[j], fp[j], fpp[j]);
    }
  }
  return 0;
}
