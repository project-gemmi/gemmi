// Copyright 2017 Global Phasing Ltd.
//
// This program analyses PDB or mmCIF files, printing similar things
// as CCP4 RWCONTENTS. For example: Matthews Coefficient.

#include <gemmi/symmetry.hpp>
#include "input.h"
#define EXE_NAME "gemmi-contents"
#include "options.h"
#include <stdio.h>

void analyse(const gemmi::Structure& st, bool verbose) {
  using namespace gemmi;
  printf(" Spacegroup   %s\n", st.sg_hm.c_str());
  int order = 1;
  const SpaceGroup* sg = find_spacegroup_by_name(st.sg_hm);
  if (sg) {
    order = sg->operations().order();
    printf("   Group no. %d with %d operations.\n", sg->number, order);
  } else {
    std::fprintf(stderr, "Unrecognized space group name! Assuming P1.\n");
  }
  printf(" Cell volume: %20.3f\n", st.cell.volume);
  printf(" ASU volume:  %20.3f\n", st.cell.volume / order);
  // TODO: strict NCS
  if (st.models.size() > 1)
    std::fprintf(stderr, "Warning: using only the first model out of %zu.\n",
                 st.models.size());
  double weight = 0;
  double water_count = 0;
  const Model& model = st.models.at(0);
  for (const Chain& chain : model.chains) {
    for (const Residue& res : chain.residues) {
      for (const Atom& atom : res.atoms) {
        // skip hydrogens
        if (atom.element == El::H || atom.element == El::D)
          continue;
        if (atom.element == El::O && res.is_water()) {
          water_count += atom.occ;
          break; // move to next residue
        }
        bool rwcontents_compat = true;
        double w = rwcontents_compat ? 2 * atom.element.atomic_number()
                                     : atom.element.weight();
        weight += atom.occ * w;
      }
    }
  }
  weight += 18 * water_count;
  printf(" Water count: %5.1g\n", water_count);
  printf(" Molecular Weight of all atoms: %20.3f\n", weight);
}

enum OptionIndex { Verbose=3 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nAnalyses content of a PDB or mmCIF."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

int main(int argc, char **argv) {
  OptParser p;
  p.simple_parse(argc, argv, Usage);
  bool verbose = p.options[Verbose];
  if (p.nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    return 0;
  }
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      const char* input = p.nonOption(i);
      if (verbose)
        std::fprintf(stderr, "Reading %s ...\n", input);
      gemmi::Structure st = read_structure(input);
      analyse(st, verbose);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
