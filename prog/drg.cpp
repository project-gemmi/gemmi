// Copyright 2025 Global Phasing Ltd.

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include "gemmi/read_cif.hpp"     // for read_cif_gz
#include "gemmi/chemcomp.hpp"     // for ChemComp, make_chemcomp_from_block
#include "gemmi/mmcif.hpp"        // for ChemCompModel, make_residue_from_chemcomp_block
#include "gemmi/acedrg_tables.hpp"  // for AcedrgTables
#include "gemmi/ace_cc.hpp"        // for prepare_chemcomp
#include "gemmi/to_cif.hpp"       // for write_cif_to_stream
#include "gemmi/to_chemcomp.hpp"  // for add_chemcomp_to_block
#include "gemmi/version.hpp"       // for GEMMI_VERSION
#include "gemmi/fstream.hpp"      // for Ofstream

#define GEMMI_PROG drg
#include "options.h"
#include "timer.h"

using namespace gemmi;

namespace {

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle, OutputDir, NoAngles, CoordModel
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
  { NoAngles, 0, "", "no-angles", Arg::None,
    "  --no-angles  \tSkip angle restraints (still keep bonds, torsions, chirals, planes)." },
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

std::string derive_companion_mol0_path(const std::string& input,
                                       const std::string& comp_id) {
  std::string filename = get_filename(input);
  std::string dir = input.substr(0, input.size() - filename.size());
  size_t pos = dir.rfind("/orig/");
  if (pos != std::string::npos)
    dir.replace(pos, 6, "/acedrg/");
  else if ((pos = dir.rfind("\\orig\\")) != std::string::npos)
    dir.replace(pos, 6, "\\acedrg\\");
  return dir + comp_id + "_TMP/" + comp_id + "_mol_0.cif";
}

std::string ltrim_copy(const std::string& s) {
  size_t pos = 0;
  while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos])))
    ++pos;
  return s.substr(pos);
}

std::string dequote_token(std::string tok) {
  if (tok.size() >= 2) {
    char q0 = tok.front();
    char q1 = tok.back();
    if ((q0 == '"' || q0 == '\'') && q1 == q0)
      return tok.substr(1, tok.size() - 2);
  }
  return tok;
}

bool parse_triplet_tokens(const std::vector<std::string>& toks,
                          int ix, int iy, int iz, Position& pos) {
  if (ix < 0 || iy < 0 || iz < 0)
    return false;
  if (static_cast<size_t>(std::max(ix, std::max(iy, iz))) >= toks.size())
    return false;
  char* end = nullptr;
  double x = std::strtod(toks[ix].c_str(), &end);
  if (!end || *end != '\0')
    return false;
  double y = std::strtod(toks[iy].c_str(), &end);
  if (!end || *end != '\0')
    return false;
  double z = std::strtod(toks[iz].c_str(), &end);
  if (!end || *end != '\0')
    return false;
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    return false;
  pos = Position(x, y, z);
  return true;
}

std::map<std::string, Position> load_companion_mol0_coords_text(
    const std::string& path) {
  std::map<std::string, Position> coords;
  std::ifstream in(path.c_str());
  if (!in.good())
    return coords;

  enum class State { FindLoop, ReadTags, ReadRows, SkipRows };
  State state = State::FindLoop;
  std::vector<std::string> tags;
  int atom_id_idx = -1;
  int x_idx = -1, y_idx = -1, z_idx = -1;
  int mx_idx = -1, my_idx = -1, mz_idx = -1;
  int ix_idx = -1, iy_idx = -1, iz_idx = -1;

  auto reset_indices = [&]() {
    atom_id_idx = -1;
    x_idx = y_idx = z_idx = -1;
    mx_idx = my_idx = mz_idx = -1;
    ix_idx = iy_idx = iz_idx = -1;
  };

  std::string line;
  while (std::getline(in, line)) {
    std::string t = ltrim_copy(line);
    if (state == State::FindLoop) {
      if (t == "loop_") {
        tags.clear();
        reset_indices();
        state = State::ReadTags;
      }
      continue;
    }

    if (state == State::ReadTags) {
      if (t.empty() || t[0] == '#')
        continue;
      if (t[0] == '_') {
        std::istringstream iss(t);
        std::string tag;
        iss >> tag;
        tags.push_back(tag);
        continue;
      }
      for (size_t i = 0; i < tags.size(); ++i) {
        const std::string& tag = tags[i];
        if (tag == "_chem_comp_atom.atom_id")
          atom_id_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.x")
          x_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.y")
          y_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.z")
          z_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.model_Cartn_x")
          mx_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.model_Cartn_y")
          my_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.model_Cartn_z")
          mz_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.pdbx_model_Cartn_x_ideal")
          ix_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.pdbx_model_Cartn_y_ideal")
          iy_idx = static_cast<int>(i);
        else if (tag == "_chem_comp_atom.pdbx_model_Cartn_z_ideal")
          iz_idx = static_cast<int>(i);
      }
      bool has_atom_loop = atom_id_idx >= 0 &&
                           ((x_idx >= 0 && y_idx >= 0 && z_idx >= 0) ||
                            (mx_idx >= 0 && my_idx >= 0 && mz_idx >= 0) ||
                            (ix_idx >= 0 && iy_idx >= 0 && iz_idx >= 0));
      state = has_atom_loop ? State::ReadRows : State::SkipRows;
      // process this line as the first row in the selected mode
    }

    if (state == State::SkipRows) {
      if (t == "loop_") {
        tags.clear();
        reset_indices();
        state = State::ReadTags;
      } else if (!t.empty() && t[0] == '_') {
        state = State::FindLoop;
      }
      continue;
    }

    if (state == State::ReadRows) {
      if (t.empty() || t[0] == '#')
        continue;
      if (t == "loop_") {
        tags.clear();
        reset_indices();
        state = State::ReadTags;
        continue;
      }
      if (t[0] == '_' || t.compare(0, 5, "data_") == 0) {
        // The atom loop is over.
        break;
      }
      std::istringstream iss(t);
      std::vector<std::string> toks;
      for (std::string tok; iss >> tok;)
        toks.push_back(tok);
      if (atom_id_idx < 0 || static_cast<size_t>(atom_id_idx) >= toks.size())
        continue;
      Position pos;
      if (parse_triplet_tokens(toks, x_idx, y_idx, z_idx, pos) ||
          parse_triplet_tokens(toks, mx_idx, my_idx, mz_idx, pos) ||
          parse_triplet_tokens(toks, ix_idx, iy_idx, iz_idx, pos))
        coords[dequote_token(toks[atom_id_idx])] = pos;
    }
  }

  return coords;
}

std::map<std::string, Position> load_companion_mol0_coords(
    const std::string& input, const std::string& comp_id) {
  std::map<std::string, Position> coords;
  bool dbg = std::getenv("GEMMI_DBG_MOL0") != nullptr;
  if (comp_id.empty())
    return coords;
  std::string mol0_path = derive_companion_mol0_path(input, comp_id);
  std::ifstream test(mol0_path.c_str());
  if (!test.good()) {
    if (dbg)
      std::fprintf(stderr, "[mol0] missing path for %s: %s\n",
                   comp_id.c_str(), mol0_path.c_str());
    return coords;
  }
  if (dbg)
    std::fprintf(stderr, "[mol0] loading %s\n", mol0_path.c_str());
  try {
    cif::Document mol0_doc = read_cif_gz(mol0_path);
    for (cif::Block& b : mol0_doc.blocks) {
      if (!b.find_values("_chem_comp_atom.atom_id"))
        continue;
      cif::Table tab = b.find("_chem_comp_atom.",
                              {"atom_id",
                               "?x", "?y", "?z",
                               "?model_Cartn_x", "?model_Cartn_y", "?model_Cartn_z",
                               "?pdbx_model_Cartn_x_ideal",
                               "?pdbx_model_Cartn_y_ideal",
                               "?pdbx_model_Cartn_z_ideal"});
      if (!tab.ok())
        continue;
      for (auto row : tab) {
        Position pos;
        auto set_if_valid = [&](int ix, int iy, int iz) {
          if (!row.has(ix) || !row.has(iy) || !row.has(iz))
            return false;
          double x = cif::as_number(row[ix]);
          double y = cif::as_number(row[iy]);
          double z = cif::as_number(row[iz]);
          if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            return false;
          pos = Position(x, y, z);
          return true;
        };
        if (set_if_valid(1, 2, 3) || set_if_valid(4, 5, 6) || set_if_valid(7, 8, 9))
          coords[row.str(0)] = pos;
      }
      break;
    }
  } catch (std::exception& e) {
    if (dbg)
      std::fprintf(stderr, "[mol0] failed to parse %s: %s\n",
                   mol0_path.c_str(), e.what());
    coords = load_companion_mol0_coords_text(mol0_path);
  } catch (...) {
    if (dbg)
      std::fprintf(stderr, "[mol0] failed to parse %s\n", mol0_path.c_str());
    coords = load_companion_mol0_coords_text(mol0_path);
  }
  if (dbg)
    std::fprintf(stderr, "[mol0] loaded %zu coords for %s\n",
                 coords.size(), comp_id.c_str());
  return coords;
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
  bool no_angles = p.options[NoAngles];
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
    tables.load_tables(tables_dir, no_angles);
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
        if (!no_angles)
          for (const auto& a : cc.rt.angles)
            if (std::isnan(a.value)) missing_angles++;

        timer.start();
        std::map<std::string, Position> sugar_coord_overrides =
            load_companion_mol0_coords(input, cc.name);
        if (verbose && !sugar_coord_overrides.empty())
          std::fprintf(stderr, "Using companion mol_0 coords for %s (%zu atoms)\n",
                       cc.name.c_str(), sugar_coord_overrides.size());
        prepare_chemcomp(cc, tables, atom_stereo, no_angles,
                         sugar_coord_overrides.empty() ? nullptr : &sugar_coord_overrides);
        timer.print("Restraints filled in");

        // Count filled values
        int remaining_bonds = 0;
        for (const auto& b : cc.rt.bonds)
          if (std::isnan(b.value)) remaining_bonds++;
        int remaining_angles = 0;
        if (!no_angles)
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
        add_chemcomp_to_block(cc, block, acedrg_types, no_angles);

        // Mimic AceDRG's extra descriptor loop.
        /*
        {
          const std::vector<std::array<const char*, 4>> rows = {{
              {cc.name.c_str(), "gemmi", GEMMI_VERSION, "dictionary generator"},
              {cc.name.c_str(), "acedrg_database", "12", "data source"},
          }};
          cif::Table tab = block.find_or_add("_acedrg_chem_comp_descriptor.",
                                             {"comp_id", "program_name",
                                              "program_version", "type"});
          tab.ensure_loop();
          tab.loop_item->loop.values.clear();  // replace any existing
          for (const auto& r : rows)
            tab.append_row({r[0], r[1], r[2], r[3]});
        }
        */
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
