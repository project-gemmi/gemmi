// Copyright 2025 Global Phasing Ltd.

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

template <typename Range>
int count_missing_values(const Range& range) {
  int missing = 0;
  for (const auto& item : range)
    if (std::isnan(item.value)) missing++;
  return missing;
}

namespace {

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
