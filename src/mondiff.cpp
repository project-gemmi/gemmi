// Copyright 2019 Global Phasing Ltd.
//
// Compare two CIF files with monomer restraints.

#include <stdio.h>
#include <gemmi/gzread.hpp>    // for read_cif_gz
#include <gemmi/chemcomp.hpp>  // for make_chemcomp_from_block
#define GEMMI_PROG mondiff
#include "options.h"

using namespace gemmi;

enum OptionIndex { Verbose=3 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n  " EXE_NAME " [options] FILE1 FILE2"
    "\nOptions:"},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

static
void compare_chemcomps(const ChemComp& cc1, const ChemComp& cc2) {
  // atoms
  for (const ChemComp::Atom& a : cc1.atoms)
    if (cc2.find_atom(a.id) == cc2.atoms.end())
      printf("- atom %s\n", a.id.c_str());
  for (const ChemComp::Atom& a : cc2.atoms)
    if (cc1.find_atom(a.id) == cc1.atoms.end())
      printf("+ atom %s\n", a.id.c_str());

  // bonds
  for (const Restraints::Bond& b1 : cc1.rt.bonds) {
    auto b2 = cc2.rt.find_bond(b1.id1, b1.id2);
    if (b2 == cc2.rt.bonds.end()) {
      printf("- bond %s\n", b1.str().c_str());
      continue;
    }
    if (b1.type != b2->type)
      printf("! bond %-8s %s : %s\n", b1.str().c_str(),
             bond_type_to_string(b1.type).c_str(),
             bond_type_to_string(b2->type).c_str());
    if (std::fabs(b1.value - b2->value) > 0.01 ||
        std::fabs(b1.esd - b2->esd) > 0.1)
      printf("! bond %-8s  value %.3f : %.3f   esd %.3f : %.3f\n",
             b1.str().c_str(), b1.value, b2->value, b1.esd, b2->esd);
  }
  for (const Restraints::Bond& b2 : cc2.rt.bonds)
    if (cc1.rt.find_bond(b2.id1, b2.id2) == cc1.rt.bonds.end())
      printf("+ bond %s\n", b2.str().c_str());

}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* path1 = p.nonOption(0);
  const char* path2 = p.nonOption(1);
  try {
    if (verbose)
      fprintf(stderr, "Reading %s ...\n", path1);
    cif::Document doc1 = read_cif_gz(path1);
    const cif::Block* block1 = &doc1.blocks.at(0);
    if (block1->name == "comp_list")
      block1 = &doc1.blocks.at(1);
    if (verbose)
      fprintf(stderr, "Reading %s ...\n", path2);
    cif::Document doc2 = read_cif_gz(path2);
    const cif::Block* block2 = doc2.find_block(block1->name);
    if (!block2)
      fail("Block " + block1->name + " not found in " + path2);
    ChemComp cc1 = make_chemcomp_from_block(*block1);
    ChemComp cc2 = make_chemcomp_from_block(*block2);
    compare_chemcomps(cc1, cc2);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
