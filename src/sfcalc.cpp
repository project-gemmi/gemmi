// Copyright 2019 Global Phasing Ltd.
//
// Calculate structure factors from a molecular model.

#include <stdio.h>
#include <complex>
#include <gemmi/fourier.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/it92.hpp>

#define GEMMI_PROG sfcalc
#include "options.h"

enum OptionIndex { Verbose=3, Hkl };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_FILE\n\n"
    "Calculates structure factors of a model (PDB or mmCIF file).\n\n"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  { Hkl, 0, "", "hkl", Arg::Int3,
    "  --hkl=H,K,L  \tCalculate structure factor F_hkl." },
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

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      gemmi::Structure st = gemmi::read_structure_gz(input);
      for (const option::Option* opt = p.options[Hkl]; opt; opt = opt->next()) {
        std::vector<int> hkl = parse_comma_separated_ints(opt->arg);
        std::complex<double> sf =
          calculate_structure_factor(st, hkl[0], hkl[1], hkl[2]);
        printf(" (%d %d %d)\t%8.2f\t%6.2f\n", hkl[0], hkl[1], hkl[2],
               std::abs(sf), gemmi::phase_in_angles(sf));

      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
