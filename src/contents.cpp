// Copyright 2017 Global Phasing Ltd.
//
// This program analyses PDB or mmCIF files, printing similar things
// as CCP4 RWCONTENTS: weight, Matthews coefficient, etc.

#include <gemmi/symmetry.hpp>
#include <gemmi/resinfo.hpp>
#include <gemmi/calculate.hpp>
#include "input.h"
#define GEMMI_PROG contents
#include "options.h"
#include <stdio.h>

using namespace gemmi;

enum OptionIndex { Verbose=3, Dihedrals };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nAnalyses content of a PDB or mmCIF."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Dihedrals, 0, "", "dihedrals", Arg::None,
    "  --dihedrals  \tPrint peptide dihedral angles." },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_content_info(const Structure& st, bool /*verbose*/) {
  printf(" Spacegroup   %s\n", st.sg_hm.c_str());
  int order = 1;
  const SpaceGroup* sg = find_spacegroup_by_name(st.sg_hm);
  if (sg) {
    order = sg->operations().order();
    printf("   Group no. %d with %d operations.\n", sg->number, order);
  } else {
    std::fprintf(stderr, "%s space group name! Assuming P1.\n",
                 st.sg_hm.empty() ? "No" : "Unrecognized");
  }
  double n_molecules = order * st.get_ncs_multiplier();
  printf(" Number of images (symmetry * strict NCS): %5g\n", n_molecules);
  printf(" Cell volume [A^3]: %30.1f\n", st.cell.volume);
  printf(" ASU volume [A^3]:  %30.1f\n", st.cell.volume / order);
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
      ResidueInfo res_info = find_tabulated_residue(res.name);
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
  double total_protein_weight = protein_weight * n_molecules;
  double Vm = st.cell.volume / total_protein_weight;
  printf(" Matthews coefficient: %29.3f\n", Vm);
  double Na = 0.602214;  // Avogadro number x 10^-24 (cm^3->A^3)
  // rwcontents uses 1.34, Rupp's papers 1.35
  for (double ro : { 1.35, 1.34 })
    printf(" Solvent %% (for protein density %g): %13.3f\n",
           ro, 100. * (1. - 1. / (ro * Vm * Na)));
}

static void print_dihedrals(const Structure& st) {
  printf(" Chain Residue      Psi      Phi    Omega\n");
  const Model& model = st.models.at(0);
  for (const Chain& chain : model.chains) {
    const char* cname = chain.name_for_pdb().c_str();
    for (const Residue& res : chain.residues) {
      printf("%3s %4d%c %5s", cname, res.seq_num_for_pdb(),
             res.printable_icode(), res.name.c_str());
      const Residue* prev = chain.prev_bonded_aa(res);
      const Residue* next = chain.next_bonded_aa(res);
      double omega = next ? gemmi::calculate_omega(res, *next) : NAN;
      auto phi_psi = gemmi::calculate_phi_psi(prev, res, next);
      if (prev || next)
        printf(" % 8.2f % 8.2f % 8.2f\n",
               gemmi::deg(phi_psi[0]), gemmi::deg(phi_psi[1]),
               gemmi::deg(omega));
      else
        printf("\n");
    }
  }
  printf("\n");
}

static void print_atoms_on_special_positions(const Structure& st) {
  printf(" Atoms on special positions:");
  bool found = false;
  for (const Chain& chain : st.models.at(0).chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (int n = st.cell.is_special_position(atom.pos)) {
          found = true;
          NearbyImage im = st.cell.find_nearest_image(atom.pos, atom.pos,
                                                      SymmetryImage::Different);
          printf("\n    %s %4d %3s %-3s %c fold=%d  occ=%.2f  d_image=%.4f",
                 chain.name_for_pdb().c_str(),
                 res.seq_num_for_pdb(), res.name.c_str(),
                 atom.name.c_str(), (atom.altloc | 0x20),
                 n+1, atom.occ, im.dist());
        }
  if (!found)
    printf(" none");
  printf("\n");
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  bool verbose = p.options[Verbose];
  if (p.nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    return 0;
  }
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.nonOption(i);
      if (is_pdb_code(input))
        input = expand_pdb_code_to_path_or_fail(input);
      if (i > 0)
        std::printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        std::printf("File: %s\n", input.c_str());
      Structure st = read_structure(input);
      if (p.options[Dihedrals])
        print_dihedrals(st);
      print_atoms_on_special_positions(st);
      print_content_info(st, verbose);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
