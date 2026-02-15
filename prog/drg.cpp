// Copyright 2025 Global Phasing Ltd.

#include <cstdio>
#include <iostream>
#include "gemmi/read_cif.hpp"     // for read_cif_gz
#include "gemmi/chemcomp.hpp"     // for ChemComp, make_chemcomp_from_block
#include "gemmi/mmcif.hpp"        // for ChemCompModel, make_residue_from_chemcomp_block
#include "gemmi/acedrg_tables.hpp"  // for AcedrgTables
#include "gemmi/ace_cc.hpp"        // for prepare_chemcomp
#include "gemmi/to_cif.hpp"       // for write_cif_to_stream
#include "gemmi/to_chemcomp.hpp"  // for add_chemcomp_to_block
#include "gemmi/fstream.hpp"      // for Ofstream

#define GEMMI_PROG drg
#include "options.h"
#include "timer.h"

using namespace gemmi;

namespace {

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle, OutputDir, OnlyBonds, CoordModel
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n " EXE_NAME " [options] -o DIR INPUT1.cif INPUT2.cif ..."
    "\n\nFill missing restraint values (bonds, angles) in a monomer CIF file"
    "\nusing COD/CSD statistical data from AceDRG tables."
    "\nIf OUTPUT.cif is -, the output is printed to stdout."
    "\nWith -o, multiple input files are processed and written to DIR."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Tables, 0, "t", "tables", Arg::Required,
    "  -t, --tables=DIR  \tDirectory with AceDRG tables (default: $ACEDRG_TABLES"
    "\n\t\tor $CCP4/share/acedrg/tables)." },
  { OutputDir, 0, "o", "output-dir", Arg::Required,
    "  -o, --output-dir=DIR  \tOutput directory for batch processing." },
  { Sigma, 0, "", "sigma", Arg::Float,
    "  --sigma=NUM  \tMaximum sigma for bond restraints (default: 0.02)." },
  { Timing, 0, "", "timing", Arg::None,
    "  --timing  \tPrint timing information." },
  { CifStyle, 0, "", "style", Arg::CifStyle,
    "  --style=STYLE  \tOutput style: default, pdbx, aligned." },
  { OnlyBonds, 0, "", "only-bonds", Arg::None,
    "  --only-bonds  \tOnly calculate bond restraints (skip angles, torsions, etc.)." },
  { CoordModel, 0, "", "coord-model", Arg::Required,
    "  --coord-model=KIND  \tCoordinates to use from _chem_comp_atom:"
    "\n\t\txyz | example | ideal | first | auto (default: auto)." },
  { 0, 0, 0, 0, 0, 0 }
};

std::string get_filename(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  return pos == std::string::npos ? path : path.substr(pos + 1);
}

enum class CoordModelChoice { Auto, Xyz, Example, Ideal, First };

const char* coord_model_name(CoordModelChoice c) {
  switch (c) {
    case CoordModelChoice::Auto: return "auto";
    case CoordModelChoice::Xyz: return "xyz";
    case CoordModelChoice::Example: return "example";
    case CoordModelChoice::Ideal: return "ideal";
    case CoordModelChoice::First: return "first";
  }
  return "first";
}

bool apply_coord_model(const cif::Block& block, ChemComp& cc, ChemCompModel kind) {
  Residue res = make_residue_from_chemcomp_block(block, kind);
  if (res.atoms.empty())
    return false;

  std::map<std::string, Position> pos_by_name;
  for (const Atom& atom : res.atoms)
    pos_by_name[atom.name] = atom.pos;

  int updated = 0;
  for (ChemComp::Atom& atom : cc.atoms) {
    auto it = pos_by_name.find(atom.id);
    if (it != pos_by_name.end()) {
      atom.xyz = it->second;
      ++updated;
    }
  }
  return updated > 0;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  bool batch_mode = p.options[OutputDir];
  if (batch_mode) {
    if (p.nonOptionsCount() < 1) {
      std::fprintf(stderr, "ERROR: At least one input file required with --output-dir.\n");
      return 1;
    }
  } else {
    p.require_positional_args(2);
  }

  int verbose = p.options[Verbose].count();
  bool only_bonds = p.options[OnlyBonds];
  CoordModelChoice coord_choice = CoordModelChoice::Auto;
  if (p.options[CoordModel]) {
    std::string kind = p.options[CoordModel].arg;
    if (kind == "auto")
      coord_choice = CoordModelChoice::Auto;
    else if (kind == "xyz")
      coord_choice = CoordModelChoice::Xyz;
    else if (kind == "example")
      coord_choice = CoordModelChoice::Example;
    else if (kind == "ideal")
      coord_choice = CoordModelChoice::Ideal;
    else if (kind == "first")
      coord_choice = CoordModelChoice::First;
    else {
      std::fprintf(stderr, "ERROR: Invalid --coord-model=%s (allowed: auto|xyz|example|ideal|first)\n",
                   kind.c_str());
      return 1;
    }
  }

  // Get tables directory
  std::string tables_dir;
  if (p.options[Tables]) {
    tables_dir = p.options[Tables].arg;
  } else {
    // Try ACEDRG_TABLES environment variable
    const char* env = std::getenv("ACEDRG_TABLES");
    if (env) {
      tables_dir = env;
    } else {
      // Fallback to $CCP4/share/acedrg/tables
      const char* ccp4 = std::getenv("CCP4");
      if (ccp4)
        tables_dir = std::string(ccp4) + "/share/acedrg/tables";
    }
  }

  if (tables_dir.empty()) {
    std::fprintf(stderr, "ERROR: No tables directory specified.\n"
                         "Use --tables=DIR or set ACEDRG_TABLES or CCP4 environment variable.\n");
    return 1;
  }

  Timer timer(p.options[Timing]);

  try {
    // Load tables
    if (verbose)
      std::fprintf(stderr, "Loading tables from %s ...\n", tables_dir.c_str());
    timer.start();
    AcedrgTables tables;
    tables.verbose = verbose;
    tables.load_tables(tables_dir, only_bonds);
    timer.print("Tables loaded in");

    if (p.options[Sigma])
      tables.lower_bond_sigma = std::strtod(p.options[Sigma].arg, nullptr);

    // Build list of (input, output) pairs
    std::vector<std::pair<std::string, std::string>> files;
    if (batch_mode) {
      std::string output_dir = p.options[OutputDir].arg;
      // Ensure output_dir ends with separator
      if (!output_dir.empty() && output_dir.back() != '/' && output_dir.back() != '\\')
        output_dir += '/';
      for (int i = 0; i < p.nonOptionsCount(); ++i) {
        std::string input = p.nonOption(i);
        std::string output = output_dir + get_filename(input);
        files.emplace_back(input, output);
      }
    } else {
      files.emplace_back(p.nonOption(0), p.nonOption(1));
    }

    int total_filled = 0;
    for (const auto& file_pair : files) {
      const std::string& input = file_pair.first;
      const std::string& output = file_pair.second;

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

        std::map<std::string, std::string> atom_stereo;
        cif::Table stereo_tab = block.find("_chem_comp_atom.",
                                           {"atom_id", "pdbx_stereo_config"});
        if (stereo_tab.ok()) {
          for (auto row : stereo_tab) {
            std::string stereo = row.str(1);
            if (stereo.empty() || stereo[0] == '.' || stereo[0] == '?')
              continue;
            atom_stereo.emplace(row.str(0), std::move(stereo));
          }
        }

        ChemComp cc = make_chemcomp_from_block(block);
        if (coord_choice != CoordModelChoice::Auto) {
          bool ok = false;
          switch (coord_choice) {
            case CoordModelChoice::Xyz:
              ok = apply_coord_model(block, cc, ChemCompModel::Xyz);
              break;
            case CoordModelChoice::Example:
              ok = apply_coord_model(block, cc, ChemCompModel::Example);
              break;
            case CoordModelChoice::Ideal:
              ok = apply_coord_model(block, cc, ChemCompModel::Ideal);
              break;
            case CoordModelChoice::First:
              ok = apply_coord_model(block, cc, ChemCompModel::First);
              break;
            case CoordModelChoice::Auto:
              break;
          }
          if (!ok) {
            std::fprintf(stderr, "ERROR: --coord-model=%s not available for block %s in %s\n",
                         coord_model_name(coord_choice), block.name.c_str(), input.c_str());
            return 1;
          }
        }

        // Count missing values before processing
        int missing_bonds = 0;
        for (const auto& b : cc.rt.bonds)
          if (std::isnan(b.value)) missing_bonds++;
        int missing_angles = 0;
        if (!only_bonds)
          for (const auto& a : cc.rt.angles)
            if (std::isnan(a.value)) missing_angles++;

        timer.start();
        prepare_chemcomp(cc, tables, atom_stereo, only_bonds);
        timer.print("Restraints filled in");

        // Count filled values
        int remaining_bonds = 0;
        for (const auto& b : cc.rt.bonds)
          if (std::isnan(b.value)) remaining_bonds++;
        int remaining_angles = 0;
        if (!only_bonds)
          for (const auto& a : cc.rt.angles)
            if (std::isnan(a.value)) remaining_angles++;

        int filled_bonds = missing_bonds - remaining_bonds;
        int filled_angles = missing_angles - remaining_angles;

        if (verbose)
          std::fprintf(stderr, "  Filled %d/%d bonds, %d/%d angles.\n",
                       filled_bonds, missing_bonds, filled_angles, missing_angles);
        filled_count += filled_bonds + filled_angles;

        // Compute acedrg_types on the processed form.
        std::vector<std::string> acedrg_types = tables.compute_acedrg_types(cc);

        // Update the block with new values
        add_chemcomp_to_block(cc, block, acedrg_types, only_bonds);
      }

      if (verbose)
        std::fprintf(stderr, "Writing %s ...\n", output.c_str());

      timer.start();
      Ofstream os(output, &std::cout);
      write_cif_to_stream(os.ref(), doc, cif_write_options(p.options[CifStyle]));
      timer.print("Output written in");

      if (verbose)
        std::fprintf(stderr, "Done with %s. Filled %d restraint values.\n",
                     input.c_str(), filled_count);
      total_filled += filled_count;
    }

    if (batch_mode && verbose)
      std::fprintf(stderr, "All done. Processed %zu files, filled %d total restraint values.\n",
                   files.size(), total_filled);

  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
