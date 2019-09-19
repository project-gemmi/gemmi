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
std::string str(const ChemComp& cc, const Restraints::Bond& b) {
  return gemmi::tostr("bond ", b.id1.atom, '-', b.id2.atom,
                      " (", cc.get_atom(b.id1.atom).chem_type,
                      '-', cc.get_atom(b.id2.atom).chem_type, ')');
}

static
std::string str(const ChemComp& cc, const Restraints::Angle& a) {
  return gemmi::tostr("angle ", a.id1.atom, '-', a.id2.atom, '-', a.id3.atom,
                      " (", cc.get_atom(a.id1.atom).chem_type,
                      '-', cc.get_atom(a.id2.atom).chem_type,
                      '-', cc.get_atom(a.id3.atom).chem_type, ')');
}

static
std::string str(const ChemComp&, const Restraints::Torsion& a) {
  return "torsion " + a.str();
}

static
std::string str(const ChemComp&, const Restraints::Chirality& a) {
  return "chirality " + a.str();
}

static
const char* mark(double delta, double eps) {
  if (delta < eps) return "";
  if (delta < 2*eps) return "*";
  if (delta < 4*eps) return "**";
  if (delta < 8*eps) return "***";
  return "****";
}

static
void compare_chemcomps(const ChemComp& cc1, const ChemComp& cc2) {
  // atoms
  for (const ChemComp::Atom& a : cc1.atoms) {
    auto b = cc2.find_atom(a.id);
    if (b == cc2.atoms.end())
      printf("- atom %s (%s)\n", a.id.c_str(), a.chem_type.c_str());
    else if (a.chem_type != b->chem_type)
      printf("! atom %s (%s : %s)\n",
             a.id.c_str(), a.chem_type.c_str(), b->chem_type.c_str());
  }
  for (const ChemComp::Atom& a : cc2.atoms)
    if (cc1.find_atom(a.id) == cc1.atoms.end())
      printf("+ atom %s (%s)\n", a.id.c_str(), a.chem_type.c_str());

  // bonds
  for (const Restraints::Bond& a : cc1.rt.bonds) {
    auto b = cc2.rt.find_bond(a.id1, a.id2);
    if (b == cc2.rt.bonds.end()) {
      printf("- bond %s\n", str(cc1, a).c_str());
    } else {
      if (a.type != b->type)
        printf("! %-30s %s : %s\n", str(cc1, a).c_str(),
               bond_type_to_string(a.type),
               bond_type_to_string(b->type));
      double delta = std::fabs(a.value - b->value);
      if (delta > 0.01 || std::fabs(a.esd - b->esd) > 0.1)
        printf("! %-30s %4s %.3f : %.3f   esd %.3f : %.3f\n",
               str(cc1, a).c_str(), mark(delta, a.esd),
               a.value, b->value, a.esd, b->esd);
    }
  }
  for (const Restraints::Bond& a : cc2.rt.bonds)
    if (cc1.rt.find_bond(a.id1, a.id2) == cc1.rt.bonds.end())
      printf("+ %s\n", str(cc2, a).c_str());

  // angles
  for (const Restraints::Angle& a : cc1.rt.angles) {
    auto b = cc2.rt.find_angle(a.id1, a.id2, a.id3);
    if (b == cc2.rt.angles.end()) {
      printf("- %s\n", str(cc1, a).c_str());
    } else {
      double delta = std::fabs(a.value - b->value);
      if (delta > 0.1 || std::fabs(a.esd - b->esd) > 1.0)
        printf("! %-30s %4s %6.2f : %6.2f   esd %.2f : %.2f\n",
               str(cc1, a).c_str(), mark(delta, a.esd),
               a.value, b->value, a.esd, b->esd);
    }
  }
  for (const Restraints::Angle& a : cc2.rt.angles)
    if (cc1.rt.find_angle(a.id1, a.id2, a.id3) == cc1.rt.angles.end())
      printf("+ %s\n", str(cc2, a).c_str());

  // torsion angles
  for (const Restraints::Torsion& a : cc1.rt.torsions) {
    auto b = cc2.rt.find_torsion(a.id1, a.id2, a.id3, a.id4);
    if (b == cc2.rt.torsions.end()) {
      printf("- %s\n", str(cc1, a).c_str());
    } else {
      double delta = std::fabs(a.value - b->value);
      if (delta > 0.1 || std::fabs(a.esd - b->esd) > 1.0)
        printf("! %-30s %4s %6.2f : %6.2f   esd %.2f : %.2f\n",
               str(cc1, a).c_str(), mark(delta, a.esd),
               a.value, b->value, a.esd, b->esd);
    }
  }
  for (const Restraints::Torsion& a : cc2.rt.torsions)
    if (cc1.rt.find_torsion(a.id1, a.id2, a.id3, a.id4) == cc1.rt.torsions.end())
      printf("+ %s\n", str(cc2, a).c_str());

  // chiralities
  for (const Restraints::Chirality& a : cc1.rt.chirs) {
    auto b = cc2.rt.find_chir(a.id_ctr, a.id1, a.id2, a.id3);
    if (b == cc2.rt.chirs.end())
      printf("- %s\n", str(cc1, a).c_str());
    else if (a.sign != b->sign)
      printf("! %-30s %s : %s\n", str(cc1, a).c_str(),
             chirality_to_string(a.sign),
             chirality_to_string(b->sign));
  }
  for (const Restraints::Chirality& a : cc2.rt.chirs)
    if (cc1.rt.find_chir(a.id_ctr, a.id1, a.id2, a.id3) == cc1.rt.chirs.end())
      printf("+ %s\n", str(cc2, a).c_str());
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
