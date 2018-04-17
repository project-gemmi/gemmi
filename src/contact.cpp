// Copyright 2018 Global Phasing Ltd.
//
// Searches for contacts -- neighbouring atoms.

#include <gemmi/subcells.hpp>
#include <gemmi/calculate.hpp>  // for are_connected
#include <gemmi/pdb.hpp>  // for split_nonpolymers
#include "input.h"
#define GEMMI_PROG contact
#include "options.h"
#include <stdio.h>
#include <cstdlib>  // for strtod
#include <algorithm>  // for min, max

using namespace gemmi;

enum OptionIndex { Verbose=3, MaxDist };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nSearches for contacts in a model (PDB or mmCIF)."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { MaxDist, 0, "d", "maxdist", Arg::Float,
    "  -d, --maxdist  \tMaximal distance in A (default 4.0)" },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_contacts(const Structure& st, float max_dist, bool verbose) {
  SubCells sc(st, 5);

  if (verbose) {
    printf(" Cell grid: %d x %d x %d\n", sc.grid.nu, sc.grid.nv, sc.grid.nw);
    size_t min_count = SIZE_MAX, max_count = 0, total_count = 0;
    for (const auto& el : sc.grid.data) {
      min_count = std::min(min_count, el.size());
      max_count = std::max(max_count, el.size());
      total_count += el.size();
    }
    printf(" Items per cell: from %zu to %zu, average: %.2g\n",
           min_count, max_count, double(total_count) / sc.grid.data.size());
  }

  const Model& model = st.models.at(0);
	for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    const Chain& chain = model.chains[n_ch];
    const gemmi::Entity* ent = st.get_entity_of(chain);
    PolymerType pt = PolymerType::Unknown;
    if (ent)
      pt = ent->polymer_type;
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Atom& atom = res.atoms[n_atom];
        sc.for_each(atom.pos, max_dist, [&](const SubCells::AtomImage& a) {
            if (a.image_idx == 0 && a.chain_idx == n_ch &&
                (a.residue_idx == n_res ||
                 are_connected(res, chain.residues[a.residue_idx], pt) ||
                 are_connected(chain.residues[a.residue_idx], res, pt)))
              return;
            const Chain& chain2 = model.chains[a.chain_idx];
            const Residue& res2 = chain2.residues[a.residue_idx];
            const Atom& atom2 = res2.atoms[a.atom_idx];
            NearbyImage im = sc.grid.unit_cell.find_nearest_image(
                              atom.pos, atom2.pos, SymmetryImage::Unspecified);
            /*
            printf(" %4s %4s %4s %4s - %4s %4s %4s %4s  %.2f\n",
                   chain.name.c_str(), res.name.c_str(),
                   res.seq_id().c_str(), atom.name.c_str(),
                   chain2.name.c_str(), res2.name.c_str(),
                   res2.seq_id().c_str(), atom2.name.c_str(),
                   im.dist());
            */
            printf("            %-4s%c%3s%2s%5s   "
                   "            %-4s%c%3s%2s%5s  %6s %6s %5.2f\n",
                   atom.padded_name().c_str(),
                   atom.altloc ? std::toupper(atom.altloc) : ' ',
                   res.name.c_str(),
                   chain.name_for_pdb().c_str(),
                   res.seq_id().c_str(),
                   atom2.padded_name().c_str(),
                   atom2.altloc ? std::toupper(atom2.altloc) : ' ',
                   res2.name.c_str(),
                   chain2.name_for_pdb().c_str(),
                   res2.seq_id().c_str(),
                   "1555", im.pdb_symbol(false).c_str(), im.dist());
        });
      }
    }
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  bool verbose = p.options[Verbose];
  float max_dist = 3.0;
  if (p.options[MaxDist])
    max_dist = std::strtod(p.options[MaxDist].arg, nullptr);
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
      if (st.input_format == gemmi::CoorFormat::Pdb)
        split_nonpolymers(st);
      print_contacts(st, max_dist, verbose);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
