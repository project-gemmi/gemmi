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

void adjust_phosphate_group(ChemComp& cc) {
  for (auto& atom : cc.atoms) {
    if (atom.el == El::O &&
        (atom.id == "OP2" || atom.id == "OP3")) {
      atom.charge = -1.0f;
    }
  }
  std::vector<std::string> phos_h;
  for (const auto& atom : cc.atoms) {
    if (atom.el == El::H &&
        (atom.id == "HOP2" || atom.id == "HOP3")) {
      phos_h.push_back(atom.id);
    }
  }
  for (const std::string& atom_id : phos_h)
    remove_atom_by_id(cc, atom_id);
}

void adjust_carboxylate_group(ChemComp& cc) {
  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }

  std::vector<std::string> hydrogens_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    std::string h_id;
    std::string c_id;
    for (const std::string& nid : nb) {
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
      if (el == El::H)
        h_id = nid;
      else if (el == El::C)
        c_id = nid;
    }
    if (h_id.empty() || c_id.empty())
      continue;
    const auto& c_nb = neighbors[c_id];
    int o_count = 0;
    for (const std::string& nid : c_nb) {
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      if (cc.atoms[it->second].el == El::O)
        o_count += 1;
    }
    if (o_count < 2)
      continue;
    atom.charge = -1.0f;
    hydrogens_to_remove.push_back(h_id);
  }

  for (const std::string& h_id : hydrogens_to_remove)
    remove_atom_by_id(cc, h_id);
}

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle, OutputDir
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n " EXE_NAME " [options] --output-dir=DIR INPUT1.cif [INPUT2.cif ...]"
    "\n\nFill missing restraint values (bonds, angles) in a monomer CIF file"
    "\nusing COD/CSD statistical data from AceDRG tables."
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
  { 0, 0, 0, 0, 0, 0 }
};

std::string get_filename(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  return pos == std::string::npos ? path : path.substr(pos + 1);
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
    tables.load_tables(tables_dir);
    timer.print("Tables loaded in");

    if (p.options[Sigma])
      tables.lower_bond_sigma = std::strtod(p.options[Sigma].arg, nullptr);

    // Set verbose level: -v=atoms, -vv=lookup, -vvv=1D/2D failures
    tables.verbose = verbose;

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

        ChemComp cc = make_chemcomp_from_block(block);
        adjust_terminal_carboxylate(cc);
        add_n_terminal_h3(cc);
        adjust_phosphate_group(cc);
        adjust_carboxylate_group(cc);

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
