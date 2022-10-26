// Copyright 2019 Global Phasing Ltd.

#include "gemmi/elem.hpp"    // for Element, find_element
#include "gemmi/fprime.hpp"  // for cromer_liberman_for_array
#include "gemmi/math.hpp"    // for hc
#include <cstdlib>           // for atof
#include <stdio.h>

#define GEMMI_PROG fprime
#include "options.h"

namespace {

enum OptionIndex { Energy=4, Wavelen };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] ELEMENT[...]"
    "\nPrints anomalous scattering factors f' and f\"."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  //CommonUsage[Verbose],
  { Energy, 0, "e", "energy", Arg::Float,
    "  -e, --energy=ENERGY  \tEnergy [eV]" },
  { Wavelen, 0, "w", "wavelength", Arg::Float,
    "  -w, --wavelength=LAMBDA  \tWavelength [A]" },
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
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* name = p.nonOption(i);
    gemmi::Element elem = gemmi::find_element(name);
    if (elem == gemmi::El::X) {
      fprintf(stderr, "Error: element name not recognized: '%s'\n", name);
      return -1;
    }
    std::vector<double> energies;
    for (const option::Option* opt = p.options[Energy]; opt; opt = opt->next())
      energies.push_back(std::atof(opt->arg));
    double hc = gemmi::hc();
    for (const option::Option* opt = p.options[Wavelen]; opt; opt = opt->next())
      energies.push_back(hc / std::atof(opt->arg));
    std::vector<double> fp(energies.size(), 0);
    std::vector<double> fpp(energies.size(), 0);
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
