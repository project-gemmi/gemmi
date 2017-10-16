// Copyright 2017 Global Phasing Ltd.
//
// This program analyses PDB or mmCIF files, printing similar things
// as CCP4 RWCONTENTS: weight, Matthews coefficient, etc.

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
  int ncs = 1;
  // TODO: strict NCS
  printf(" Number of molecules: %8d\n", order * ncs);
  printf(" Cell volume: %20.3f\n", st.cell.volume);
  printf(" ASU volume:  %20.3f\n", st.cell.volume / order);
  if (st.models.size() > 1)
    std::fprintf(stderr, "Warning: using only the first model out of %zu.\n",
                 st.models.size());
  double water_count = 0;
  int h_count = 0;
  double weight = 0;
  double protein_weight = 0;
  double atom_count = 0;
  double protein_atom_count = 0;
  const Model& model = st.models.at(0);
  for (const Chain& chain : model.chains) {
    for (const Residue& res : chain.residues) {
      ResidueInfo res_info = res.get_info();
      if (res_info.is_water())
        if (const Atom* oxygen = res.find_by_element(El::O))
          water_count += oxygen->occ;
      bool is_protein = false;
      if (res_info.is_amino() || res_info.is_nucleic() ||
          res.name == "HEM" || res.name == "SO4" || res.name == "SUL") {
        is_protein = true;
        h_count += res_info.hydrogen_count;
      }
      for (const Atom& atom : res.atoms) {
        // skip hydrogens
        if (atom.element == El::H || atom.element == El::D)
          continue;
        if (is_protein) {
          protein_atom_count += atom.occ;
          protein_weight += atom.occ * atom.element.weight();
        }
        atom_count += atom.occ;
        weight += atom.occ * atom.element.weight();
      }
    }
  }
  double h_weight = Element(El::H).weight();
  weight += (2 * water_count + h_count) * h_weight;
  protein_weight += h_count * h_weight;
  printf(" Heavy (not H) atom count: %25.3f\n", + atom_count + water_count);
  printf(" Estimate of the protein hydrogens: %12d\n", h_count);
  printf(" Estimated total atom count (incl. H): %13.3f\n",
                                atom_count + 3 * water_count + h_count);
  printf(" Estimated protein atom count (incl. H): %11.3f\n",
                                          protein_atom_count + h_count);
  printf(" Water count: %38.3f\n", water_count);
  printf(" Molecular weight of all atoms: %20.3f\n", weight);
  printf(" Molecular weight of protein atoms: %16.3f\n", protein_weight);
  double total_protein_weight = protein_weight * order * ncs;
  double Vm = st.cell.volume / total_protein_weight;
  printf(" Matthews coefficient: %29.3f\n", Vm);
  double Na = 0.602214;
  double ro = 1.34;
  for (double ro : { 1.34, 1.33 })
    printf(" Solvent %% (for protein density %g): %12.3f\n",
           ro, 100. * (1 - 1 / (ro * Vm * Na)));
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
