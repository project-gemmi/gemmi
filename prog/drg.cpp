// Copyright 2025 Global Phasing Ltd.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include "gemmi/read_cif.hpp"     // for read_cif_gz
#include "gemmi/chemcomp.hpp"     // for ChemComp, make_chemcomp_from_block
#include "gemmi/acedrg_tables.hpp" // for AcedrgTables
#include "gemmi/to_cif.hpp"       // for write_cif_to_stream
#include "gemmi/to_chemcomp.hpp"  // for add_chemcomp_to_block
#include "gemmi/fstream.hpp"      // for Ofstream

#define GEMMI_PROG drg
#include "options.h"
#include "timer.h"

using namespace gemmi;

namespace {

template <typename Range>
int count_missing_values(const Range& range) {
  int missing = 0;
  for (const auto& item : range)
    if (std::isnan(item.value)) missing++;
  return missing;
}

void remove_atom_by_id(ChemComp& cc, const std::string& atom_id) {
  auto is_id = [&](const Restraints::AtomId& id) { return id.atom == atom_id; };
  cc.rt.bonds.erase(std::remove_if(cc.rt.bonds.begin(), cc.rt.bonds.end(),
                                  [&](const Restraints::Bond& b) {
                                    return is_id(b.id1) || is_id(b.id2);
                                  }),
                    cc.rt.bonds.end());
  cc.rt.angles.erase(std::remove_if(cc.rt.angles.begin(), cc.rt.angles.end(),
                                   [&](const Restraints::Angle& a) {
                                     return is_id(a.id1) || is_id(a.id2) || is_id(a.id3);
                                   }),
                     cc.rt.angles.end());
  cc.rt.torsions.erase(std::remove_if(cc.rt.torsions.begin(), cc.rt.torsions.end(),
                                     [&](const Restraints::Torsion& t) {
                                       return is_id(t.id1) || is_id(t.id2) ||
                                              is_id(t.id3) || is_id(t.id4);
                                     }),
                       cc.rt.torsions.end());
  cc.rt.chirs.erase(std::remove_if(cc.rt.chirs.begin(), cc.rt.chirs.end(),
                                  [&](const Restraints::Chirality& c) {
                                    return is_id(c.id_ctr) || is_id(c.id1) ||
                                           is_id(c.id2) || is_id(c.id3);
                                  }),
                    cc.rt.chirs.end());
  for (auto it = cc.rt.planes.begin(); it != cc.rt.planes.end(); ) {
    auto& ids = it->ids;
    ids.erase(std::remove_if(ids.begin(), ids.end(),
                             [&](const Restraints::AtomId& id) { return is_id(id); }),
              ids.end());
    if (ids.empty())
      it = cc.rt.planes.erase(it);
    else
      ++it;
  }
  cc.atoms.erase(std::remove_if(cc.atoms.begin(), cc.atoms.end(),
                               [&](const ChemComp::Atom& a) { return a.id == atom_id; }),
                 cc.atoms.end());
}

void adjust_terminal_carboxylate(ChemComp& cc) {
  bool has_oxt = false;
  bool has_hxt = false;
  for (const auto& atom : cc.atoms) {
    if (atom.id == "OXT" && atom.el == El::O)
      has_oxt = true;
    else if (atom.id == "HXT" && atom.el == El::H)
      has_hxt = true;
  }
  if (!has_oxt || !has_hxt)
    return;
  for (auto& atom : cc.atoms)
    if (atom.id == "OXT")
      atom.charge = -1.0f;
  remove_atom_by_id(cc, "HXT");
}

void add_n_terminal_h3(ChemComp& cc) {
  if (cc.find_atom("N") == cc.atoms.end())
    return;
  if (cc.find_atom("H3") != cc.atoms.end())
    return;
  if (cc.find_atom("H") == cc.atoms.end() || cc.find_atom("H2") == cc.atoms.end())
    return;

  bool n_has_h = false;
  bool n_has_h2 = false;
  bool n_has_h3 = false;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == "N" || bond.id2.atom == "N") {
      const std::string& other = (bond.id1.atom == "N") ? bond.id2.atom : bond.id1.atom;
      if (other == "H")
        n_has_h = true;
      else if (other == "H2")
        n_has_h2 = true;
      else if (other == "H3")
        n_has_h3 = true;
    }
  }
  if (!n_has_h || !n_has_h2 || n_has_h3)
    return;

  cc.atoms.push_back(ChemComp::Atom{"H3", "", El::H, 0.0f, "H", Position()});
  // Bond/angle values set to NAN - will be filled by fill_restraints()
  cc.rt.bonds.push_back({{1, "N"}, {1, "H3"}, BondType::Single, false,
                        NAN, NAN, NAN, NAN});

  auto add_angle_if_present = [&](const std::string& a1, const std::string& a2,
                                  const std::string& a3) {
    if (cc.find_atom(a1) == cc.atoms.end() ||
        cc.find_atom(a2) == cc.atoms.end() ||
        cc.find_atom(a3) == cc.atoms.end())
      return;
    cc.rt.angles.push_back({{1, a1}, {1, a2}, {1, a3}, NAN, NAN});  // filled later
  };

  add_angle_if_present("CA", "N", "H3");
  add_angle_if_present("H", "N", "H3");
  add_angle_if_present("H2", "N", "H3");
}

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n\nFill missing restraint values (bonds, angles) in a monomer CIF file"
    "\nusing COD/CSD statistical data from AceDRG tables."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Tables, 0, "t", "tables", Arg::Required,
    "  -t, --tables=DIR  \tDirectory with AceDRG tables." },
  { Sigma, 0, "", "sigma", Arg::Float,
    "  --sigma=NUM  \tMaximum sigma for bond restraints (default: 0.02)." },
  { Timing, 0, "", "timing", Arg::None,
    "  --timing  \tPrint timing information." },
  { CifStyle, 0, "", "style", Arg::CifStyle,
    "  --style=STYLE  \tOutput style: default, pdbx, aligned." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  std::string input = p.nonOption(0);
  std::string output = p.nonOption(1);
  int verbose = p.options[Verbose].count();

  // Get tables directory
  std::string tables_dir;
  if (p.options[Tables]) {
    tables_dir = p.options[Tables].arg;
  } else {
    // Try environment variable
    const char* env = std::getenv("ACEDRG_TABLES");
    if (env)
      tables_dir = env;
  }

  if (tables_dir.empty()) {
    std::fprintf(stderr, "ERROR: No tables directory specified.\n"
                         "Use --tables=DIR or set ACEDRG_TABLES environment variable.\n");
    return 1;
  }

  Timer timer(p.options[Timing]);

  try {
    // Load tables
    if (verbose)
      std::fprintf(stderr, "Loading tables from %s ...\n", tables_dir.c_str());
    timer.start();
    AcedrgTables tables;
    tables.load_tables(tables_dir);
    timer.print("Tables loaded in");

    if (p.options[Sigma])
      tables.lower_bond_sigma = std::strtod(p.options[Sigma].arg, nullptr);

    // Read input CIF
    if (verbose)
      std::fprintf(stderr, "Reading %s ...\n", input.c_str());
    timer.start();
    cif::Document doc = read_cif_gz(input);
    timer.print("Input CIF read in");

    int filled_count = 0;
    for (cif::Block& block : doc.blocks) {
      // Skip blocks that don't look like monomer definitions
      if (!block.find_values("_chem_comp_atom.atom_id"))
        continue;

      if (verbose)
        std::fprintf(stderr, "Processing block %s ...\n", block.name.c_str());

      ChemComp cc = make_chemcomp_from_block(block);
      adjust_terminal_carboxylate(cc);
      add_n_terminal_h3(cc);

      // Count missing values before
      int missing_bonds = count_missing_values(cc.rt.bonds);
      int missing_angles = count_missing_values(cc.rt.angles);

      if (missing_bonds == 0 && missing_angles == 0) {
        if (verbose)
          std::fprintf(stderr, "  No missing values.\n");
        continue;
      }

      // Fill restraints
      timer.start();
      tables.fill_restraints(cc);
      tables.assign_ccp4_types(cc);
      timer.print("Restraints filled in");

      // Count filled values
      int filled_bonds = missing_bonds - count_missing_values(cc.rt.bonds);
      int filled_angles = missing_angles - count_missing_values(cc.rt.angles);

      if (verbose)
        std::fprintf(stderr, "  Filled %d/%d bonds, %d/%d angles.\n",
                     filled_bonds, missing_bonds, filled_angles, missing_angles);

      filled_count += filled_bonds + filled_angles;

      // Update the block with new values
      add_chemcomp_to_block(cc, block);
    }

    if (verbose)
      std::fprintf(stderr, "Writing %s ...\n", output.c_str());

    timer.start();
    Ofstream os(output, &std::cout);
    write_cif_to_stream(os.ref(), doc, cif_write_options(p.options[CifStyle]));
    timer.print("Output written in");

    if (verbose)
      std::fprintf(stderr, "Done. Filled %d restraint values.\n", filled_count);

  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
