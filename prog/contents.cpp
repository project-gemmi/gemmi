// Copyright 2017 Global Phasing Ltd.
//
// This program analyses PDB or mmCIF files, printing similar things
// as CCP4 RWCONTENTS: weight, Matthews coefficient, etc.

#include <cassert>
#include <cstdio>
#include <gemmi/symmetry.hpp>
#include <gemmi/resinfo.hpp>
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/seqtools.hpp>  // for calculate_sequence_weight
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/select.hpp>    // for Selection
#include <gemmi/stats.hpp>     // for DataStats
#include <gemmi/calculate.hpp> // for expand_box, calculate_omega
#include "histogram.h"         // for print_histogram
#define GEMMI_PROG contents
#include "options.h"

using namespace gemmi;
using std::printf;

namespace {

enum OptionIndex { Select=4, Bfactors, Dihedrals, NoContentInfo };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nAnalyses content of a PDB or mmCIF."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL  \tUse only the selection." },
  { Bfactors, 0, "b", "", Arg::None,
    "  -b  \tPrint statistics of isotropic ADPs (B-factors)." },
  { Dihedrals, 0, "", "dihedrals", Arg::None,
    "  --dihedrals  \tPrint peptide dihedral angles." },
  { NoContentInfo, 0, "n", "", Arg::None,
    "  -n  \tDo not print content (for use with other options)." },
  { 0, 0, 0, 0, 0, 0 }
};

void print_atoms_on_special_positions(const Structure& st) {
  printf(" Atoms on special positions:");
  bool found = false;
  for (const Chain& chain : st.first_model().chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (int n = st.cell.is_special_position(atom.pos)) {
          found = true;
          NearestImage im = st.cell.find_nearest_image(atom.pos, atom.pos,
                                                       Asu::Different);
          printf("\n    %s %4d%c %3s %-3s %c fold=%d  occ=%.2f  d_image=%.4f",
                 chain.name.c_str(), *res.seqid.num, res.seqid.icode,
                 res.name.c_str(), atom.name.c_str(), (atom.altloc | 0x20),
                 n+1, atom.occ, im.dist());
        }
  if (!found)
    printf(" none");
  printf("\n");
}

void print_solvent_content(const UnitCell& cell, double mol_weight) {
  if (cell.is_crystal()) {
    double Vm = cell.volume_per_image() / mol_weight;
    printf(" Matthews coefficient: %29.3f\n", Vm);
    double Na = 0.602214;  // Avogadro number x 10^-24 (cm^3->A^3)
    // rwcontents uses 1.34, Rupp's papers 1.35
    for (double ro : { 1.35, 1.34 })
      printf(" Solvent %% (for protein density %g): %13.3f\n",
             ro, 100. * (1. - 1. / (ro * Vm * Na)));
  } else {
    printf(" Not a crystal / unit cell not known.\n");
  }
}

void print_content_info(const Structure& st, bool /*verbose*/) {
  printf(" Spacegroup   %s\n", st.spacegroup_hm.c_str());
  const Model& model = st.first_model();
  int order = 1;
  if (st.cell.is_crystal()) {
    if (const SpaceGroup* sg = st.find_spacegroup()) {
      order = sg->operations().order();
      printf("   Group no. %d with %d operations.\n", sg->number, order);
    } else {
      std::fprintf(stderr, "%s space group name! Assuming P1.\n",
                   st.spacegroup_hm.empty() ? "No" : "Unrecognized");
    }
  } else {
    printf("   Not a crystal.\n");
    Box<Position> box;
    expand_box(model, box);
    printf("   Atoms in: x [%g, %g]  y [%g, %g]  z [%g, %g]\n",
           box.minimum.x, box.maximum.x,
           box.minimum.y, box.maximum.y,
           box.minimum.z, box.maximum.z);
    if (st.ncs_not_expanded()) {
      for (const NcsOp& ncs_op : st.ncs) {
        if (!ncs_op.given)
          for (const_CRA cra : model.all())
            box.extend(ncs_op.apply(cra.atom->pos));
      }
      printf("   With NCS: x [%g, %g]  y [%g, %g]  z [%g, %g]\n",
             box.minimum.x, box.maximum.x,
             box.minimum.y, box.maximum.y,
             box.minimum.z, box.maximum.z);
    }
  }
  if (!st.origx.is_identity())
    printf("   The ORIGX matrix is not identity.\n");
  if (st.cell.explicit_matrices)
    printf("   Non-standard fractionalization matrix is given.\n");
  if (st.cell.is_crystal())
    print_atoms_on_special_positions(st);
  double n_molecules = order * st.get_ncs_multiplier();
  printf(" Number of images (symmetry * strict NCS): %5g\n", n_molecules);
  assert(n_molecules == st.cell.images.size() + 1);
  if (st.cell.is_crystal()) {
    printf(" Cell volume [A^3]: %30.1f\n", st.cell.volume);
    printf(" ASU volume [A^3]:  %30.1f\n", st.cell.volume / order);
  }
  double water_count = 0;
  int residue_count = 0;
  int mol_h_count = 0;
  double mol_weight = 0;
  double mol_atom_count = 0;
  double buffer_atom_count = 0;
  double file_h_count = 0;
  for (const Chain& chain : model.chains) {
    for (const Residue& res : chain.residues) {
      ResidueInfo res_info = find_tabulated_residue(res.name);
      bool is_buffer = res_info.is_buffer_or_water();
      if (!is_buffer && chain.is_first_in_group(res)) {
        residue_count++;
        mol_h_count += std::max(res_info.hydrogen_count - 2, 0);
      }
      for (const Atom& atom : res.atoms) {
        if (atom.is_hydrogen()) {
          file_h_count += atom.occ;
        } else if (is_buffer) {
          if (res_info.is_water())
            water_count += atom.occ;
          buffer_atom_count += atom.occ;
        } else {
          mol_atom_count += atom.occ;
          mol_weight += atom.occ * atom.element.weight();
        }

        // sanity check: occupancies
        if (atom.occ > 1.0f || atom.occ < 0.f)
          printf("WARNING: Occupancy of %s: %g\n",
                 atom_str(chain, res, atom).c_str(), atom.occ);
        if (atom.altloc && (&atom == &res.atoms[0] || (&atom - 1)->name != atom.name)) {
          float occ_sum = atom.occ;
          for (const Atom* a = &atom + 1; a < res.atoms.data() + res.atoms.size(); ++a)
            if (a->name == atom.name)
              occ_sum += a->occ;
          if (occ_sum > 1.0f)
            printf("WARNING: Sum of altloc occupancies of %s/%s %s/%s: %g\n",
                   chain.name.c_str(), res.name.c_str(), res.seqid.str().c_str(),
                   atom.name.c_str(), occ_sum);
        }
      }
    }
  }
  // add weight of hydrogens
  mol_weight += mol_h_count * Element(El::H).weight();

  printf(" Residue count excl. solvent and buffer: %7d\n", residue_count);
  printf(" Water count: %38.3f\n", water_count);
  printf(" Heavy (not H) atom count: %25.3f\n",
         mol_atom_count + buffer_atom_count);
  printf("     in macromolecules and ligands: %16.3f\n", mol_atom_count);
  printf("     in solvent and buffer: %24.3f\n", buffer_atom_count);
  printf(" Hydrogens in the file: %28.3f\n", file_h_count);
  printf("Solvent content based on the model (excl. solvent and buffer)\n");
  printf(" Estimated hydrogen count: %21d\n", mol_h_count);
  printf(" Estimated molecular weight: %23.3f\n", mol_weight);
  print_solvent_content(st.cell, mol_weight);
  printf("Solvent content based on SEQRES\n");
  mol_weight = 0.;
  bool missing = false;
  for (const Chain& chain : model.chains)
    if (ConstResidueSpan polymer = chain.get_polymer()) {
      const Entity* entity = st.get_entity_of(polymer);
      if (entity && !entity->full_sequence.empty()) {
        mol_weight += calculate_sequence_weight(entity->full_sequence, 100.);
      } else {
        printf(" Missing sequence for chain %s.\n", chain.name.c_str());
        missing = true;
      }
    }
  if (missing)
    return;
  printf(" Molecular weight from sequence: %19.3f\n", mol_weight);
  print_solvent_content(st.cell, mol_weight);
}

void print_dihedrals(const Structure& st) {
  printf(" Chain Residue      Psi      Phi    Omega\n");
  const Model& model = st.first_model();
  for (const Chain& chain : model.chains) {
    for (const Residue& res : chain.residues) {
      printf("%3s %4d%c %5s", chain.name.c_str(), *res.seqid.num,
                              res.seqid.icode, res.name.c_str());
      const Residue* prev = chain.previous_residue(res);
      if (!are_connected(*prev, res, PolymerType::PeptideL))
        prev = nullptr;
      const Residue* next = chain.next_residue(res);
      if (!are_connected(res, *next, PolymerType::PeptideL))
        next = nullptr;
      double omega = next ? calculate_omega(res, *next) : NAN;
      auto phi_psi = calculate_phi_psi(prev, res, next);
      if (prev || next)
        printf(" % 8.2f % 8.2f % 8.2f\n",
               deg(phi_psi[0]), deg(phi_psi[1]), deg(omega));
      else
        printf("\n");
    }
  }
  printf("\n");
}

void print_bfactor_info(const gemmi::Model& model) {
  std::vector<double> bfactors;
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.occ > 0)
          bfactors.push_back(atom.b_iso);
  gemmi::DataStats stats = gemmi::calculate_data_statistics(bfactors);
  printf("\nIsotropic ADPs: %zu values\n", bfactors.size());
  printf("  min: %.2f  max: %.2f  mean: %.2f  std.dev: %.2f\n",
         stats.dmin, stats.dmax, stats.dmean, stats.rms);
  if (stats.dmin < stats.dmax)
    print_histogram(bfactors, stats.dmin, stats.dmax);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  bool verbose = p.options[Verbose];
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (i > 0)
        std::printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        std::printf("File: %s\n", input.c_str());
      Structure st = read_structure_gz(input);
      setup_entities(st);
      if (p.options[Select])
        gemmi::Selection(p.options[Select].arg).remove_not_selected(st);
      if (st.models.size() > 1)
        std::fprintf(stderr,
                     "Warning: using only the first model out of %zu.\n",
                     st.models.size());
      if (!p.options[NoContentInfo])
        print_content_info(st, verbose);
      if (p.options[Bfactors])
        print_bfactor_info(st.first_model());
      if (p.options[Dihedrals])
        print_dihedrals(st);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
