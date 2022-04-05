// Copyright 2019 Global Phasing Ltd.
//
// Compare two CIF files with monomer restraints.

#include <stdio.h>
#include <cstdlib>
#include <set>
#include <gemmi/chemcomp.hpp>  // for make_chemcomp_from_block
#include <gemmi/read_cif.hpp>  // for read_cif_gz
#define GEMMI_PROG mondiff
#include "options.h"

using namespace gemmi;

namespace {

enum OptionIndex { Bond=4, BondEsd, Angle, AngleEsd, Relative };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n  " EXE_NAME " [options] FILE1 FILE2"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { NoOp, 0, "", "", Arg::None, "\nMinimal reported differences:" },
  { Bond, 0, "", "bond", Arg::Float,
    "  --bond=DELTA  \tdifference in distance value (default: 0.01)." },
  { BondEsd, 0, "", "bond-esd", Arg::Float,
    "  --bond-esd=DELTA  \tdifference in distance esd (default: 0.1)." },
  { Angle, 0, "", "angle", Arg::Float,
    "  --angle=DELTA  \tdifference in angle value (default: 0.1)." },
  { AngleEsd, 0, "", "angle-esd", Arg::Float,
    "  --angle-esd=DELTA  \tdifference in angle esd (default: 1.0)." },
  { Relative, 0, "", "rel", Arg::Float,
    "  --rel=SIGMA  \t abs(value difference) / esd > SIGMA (default: 0.0)." },
  { 0, 0, 0, 0, 0, 0 }
};

struct MinDelta {
  double bond = 0.01;
  double bond_esd = 0.1;
  double angle = 0.1;
  double angle_esd = 1.0;
  double rel = 0.0;
};

std::string str(const ChemComp& cc, const Restraints::Bond& b) {
  return gemmi::cat("bond ", b.id1.atom, '-', b.id2.atom,
                    " (", cc.get_atom(b.id1.atom).chem_type,
                    '-', cc.get_atom(b.id2.atom).chem_type, ')');
}

std::string str(const ChemComp& cc, const Restraints::Angle& a) {
  return gemmi::cat("angle ", a.id1.atom, '-', a.id2.atom, '-', a.id3.atom,
                    " (", cc.get_atom(a.id1.atom).chem_type,
                    '-', cc.get_atom(a.id2.atom).chem_type,
                    '-', cc.get_atom(a.id3.atom).chem_type, ')');
}

std::string str(const ChemComp&, const Restraints::Torsion& a) {
  return "torsion " + a.str();
}

std::string str(const ChemComp&, const Restraints::Chirality& a) {
  return "chirality " + a.str();
}

const char* mark(double delta, double eps) {
  if (delta < eps) return "";
  if (delta < 2*eps) return "*";
  if (delta < 4*eps) return "**";
  if (delta < 8*eps) return "***";
  return "****";
}

void compare_chemcomps(const ChemComp& cc1, const ChemComp& cc2,
                       const MinDelta& delta) {
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
      printf("- %s\n", str(cc1, a).c_str());
    } else {
      if (a.type != b->type)
        printf("! %-30s %s : %s\n", str(cc1, a).c_str(),
               bond_type_to_string(a.type),
               bond_type_to_string(b->type));
      double d = std::fabs(a.value - b->value);
      if ((d > delta.bond || std::fabs(a.esd - b->esd) > delta.bond_esd) &&
          d > delta.rel * std::min(a.esd, b->esd))
        printf("! %-30s %4s %.3f : %.3f   esd %.3f : %.3f\n",
               str(cc1, a).c_str(), mark(d, a.esd),
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
      double d = std::fabs(a.value - b->value);
      if ((d > delta.angle || std::fabs(a.esd - b->esd) > delta.angle_esd) &&
          d > delta.rel * std::min(a.esd, b->esd))
        printf("! %-30s %4s %6.2f : %6.2f   esd %.2f : %.2f\n",
               str(cc1, a).c_str(), mark(d, a.esd),
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
      double d = std::fabs(a.value - b->value);
      if ((d > delta.angle || std::fabs(a.esd - b->esd) > delta.angle_esd) &&
          d > delta.rel * std::min(a.esd, b->esd))
        printf("! %-30s %4s %6.2f : %6.2f   esd %.2f : %.2f\n",
               str(cc1, a).c_str(), mark(d, a.esd),
               a.value, b->value, a.esd, b->esd);
    }
  }
  for (const Restraints::Torsion& a : cc2.rt.torsions)
    if (cc1.rt.find_torsion(a.id1, a.id2, a.id3,a.id4) == cc1.rt.torsions.end())
      printf("+ %s\n", str(cc2, a).c_str());

  // chiralities
  std::vector<bool> matched_chir(cc2.rt.chirs.size(), false);
  for (const Restraints::Chirality& a : cc1.rt.chirs) {
    if (a.sign == ChiralityType::Both)  // not interesting
      continue;
    auto b = cc2.rt.find_chir(a.id_ctr, a.id1, a.id2, a.id3);
    if (b != cc2.rt.chirs.end()) {
      if (a.sign != b->sign)
        printf("! %-30s %s : %s\n", str(cc1, a).c_str(),
               chirality_to_string(a.sign), chirality_to_string(b->sign));
      matched_chir[b - cc2.rt.chirs.begin()] = true;
      continue;
    }
    b = cc2.rt.find_chir(a.id_ctr, a.id1, a.id3, a.id2);
    if (b == cc2.rt.chirs.end()) {
      printf("- %s  %s\n", str(cc1, a).c_str(), chirality_to_string(a.sign));
    } else {
      matched_chir[b - cc2.rt.chirs.begin()] = true;
      if (b->sign == ChiralityType::Both || b->sign == a.sign)
        printf("! %-30s %s : -%s\n", str(cc1, a).c_str(),
               chirality_to_string(a.sign), chirality_to_string(b->sign));
    }
  }
  for (size_t i = 0; i != cc2.rt.chirs.size(); ++i)
    if (!matched_chir[i]) {
      const Restraints::Chirality& a = cc2.rt.chirs[i];
      printf("+ %s  %s\n", str(cc2, a).c_str(), chirality_to_string(a.sign));
    }

  // planes
  std::vector<bool> matched_planes(cc2.rt.planes.size(), false);
  std::vector<std::set<Restraints::AtomId>> planes2;
  planes2.reserve(cc2.rt.planes.size());
  for (const Restraints::Plane& b : cc2.rt.planes)
    planes2.emplace_back(b.ids.begin(), b.ids.end());
  for (const Restraints::Plane& a : cc1.rt.planes) {
    std::set<Restraints::AtomId> plane1(a.ids.begin(), a.ids.end());
    auto b = std::find(planes2.begin(), planes2.end(), plane1);
    if (b == planes2.end()) {
      printf("- plane %s\n", a.str().c_str());
      continue;
    }
    double b_esd = cc2.rt.planes[b - planes2.begin()].esd;
    if (std::fabs(a.esd - b_esd) > 0.02)
      printf("! plane %-53s esd %.2f : %.2f\n", a.str().c_str(), a.esd, b_esd);
    b->clear();
  }
  for (size_t i = 0; i != planes2.size(); ++i)
    if (!planes2[i].empty())
      printf("+ plane %s\n", cc2.rt.planes[i].str().c_str());
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* path1 = p.nonOption(0);
  const char* path2 = p.nonOption(1);

  MinDelta delta;
  if (p.options[Bond])
    delta.bond = std::atof(p.options[Bond].arg);
  if (p.options[BondEsd])
    delta.bond_esd = std::atof(p.options[BondEsd].arg);
  if (p.options[Angle])
    delta.angle = std::atof(p.options[Angle].arg);
  if (p.options[AngleEsd])
    delta.angle_esd = std::atof(p.options[AngleEsd].arg);
  if (p.options[Relative])
    delta.rel = std::atof(p.options[Relative].arg);

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
      fail("Block ", block1->name, " not found in ", path2);
    ChemComp cc1 = make_chemcomp_from_block(*block1);
    ChemComp cc2 = make_chemcomp_from_block(*block2);
    compare_chemcomps(cc1, cc2, delta);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (std::out_of_range& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }
  return 0;
}
