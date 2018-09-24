// Copyright 2018 Global Phasing Ltd.
//
// Searches for contacts -- neighbouring atoms.

#include <gemmi/subcells.hpp>
#include <gemmi/polyheur.hpp>  // for are_connected
#include <gemmi/elem.hpp>  // for is_hydrogen
#include <gemmi/gzread.hpp>
#define GEMMI_PROG contact
#include "options.h"
#include <stdio.h>
#include <cstdlib>  // for strtod
#include <algorithm>  // for min, max

using namespace gemmi;

enum OptionIndex { Verbose=3, MaxDist, Occ, Any, NoH, Count };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nSearches for contacts in a model (PDB or mmCIF)."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { MaxDist, 0, "d", "maxdist", Arg::Float,
    "  -d, --maxdist=D  \tMaximal distance in A (default 3.0)" },
  { Occ, 0, "", "occsum", Arg::Float,
    "  --occsum=MIN  \tIgnore atom pairs with summed occupancies < MIN." },
  { Any, 0, "", "any", Arg::None,
    "  --any  \tOutput any atom pair, even from the same residue." },
  { NoH, 0, "", "noh", Arg::None,
    "  --noh  \tIgnore hydrogen (and deuterium) atoms." },
  { Count, 0, "", "count", Arg::None,
    "  --count  \tPrint only a count of atom pairs." },
  { 0, 0, 0, 0, 0, 0 }
};

struct Parameters {
  float max_dist;
  float occ_sum;
  bool any;
  bool print_count;
  bool no_hydrogens;
  int verbose;
};

static void print_contacts(const Structure& st, const Parameters& params) {
  const float special_pos_cutoff = 0.8f;
  SubCells sc(st.models.at(0), st.cell, std::max(5.0f, params.max_dist));

  if (params.verbose > 0) {
    if (params.verbose > 1) {
      if (st.cell.explicit_matrices)
        printf(" Using fractionalization matrix from the file.\n");
      printf(" Each atom has %zu extra images.\n", st.cell.images.size());
    }
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

  int counter = 0;
  const Model& model = st.models.at(0);
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    const Chain& chain = model.chains[n_ch];
    PolymerType pt = check_polymer_type(chain.get_polymer());
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Atom& atom = res.atoms[n_atom];
        if (params.no_hydrogens && is_hydrogen(atom.element))
          continue;
        double min_occ = params.occ_sum - atom.occ;
        sc.for_each(atom.pos, atom.altloc, params.max_dist,
                    [&](const SubCells::Mark& m, float dist_sq) {
            if (!params.any) {
              if (m.image_idx == 0 && m.chain_idx == n_ch &&
                  (m.residue_idx == n_res ||
                   are_connected(res, chain.residues[m.residue_idx], pt) ||
                   are_connected(chain.residues[m.residue_idx], res, pt)))
                return;
            }
            if (m.chain_idx == n_ch && m.residue_idx == n_res &&
                m.atom_idx == n_atom && (m.image_idx == 0 ||
                                         dist_sq < sq(special_pos_cutoff)))
              return;
            if (params.no_hydrogens && is_hydrogen(m.element))
              return;
            const_CRA cra = m.to_cra(model);
            if (cra.atom->occ < min_occ)
              return;
            ++counter;
            if (params.print_count)
              return;
            SymImage im = st.cell.find_nearest_pbc_image(
                                        atom.pos, cra.atom->pos, m.image_idx);
            printf("            %-4s%c%3s%2s%5s   "
                   "            %-4s%c%3s%2s%5s  %6s %6s %5.2f\n",
                   atom.padded_name().c_str(),
                   atom.altloc ? std::toupper(atom.altloc) : ' ',
                   res.name.c_str(),
                   chain.name.c_str(),
                   res.seqid.str().c_str(),
                   cra.atom->padded_name().c_str(),
                   cra.atom->altloc ? std::toupper(cra.atom->altloc) : ' ',
                   cra.residue->name.c_str(),
                   cra.chain->name.c_str(),
                   cra.residue->seqid.str().c_str(),
                   "1555", im.pdb_symbol(false).c_str(), im.dist());
        });
      }
    }
  }
  if (params.print_count)
    printf("%s:%g\n", st.name.c_str(), 0.5 * counter);
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  Parameters params;
  params.verbose = p.options[Verbose].count();
  params.max_dist = 3.0;
  if (p.options[MaxDist])
    params.max_dist = std::strtof(p.options[MaxDist].arg, nullptr);
  params.occ_sum = 0;
  if (p.options[Occ])
    params.occ_sum = std::strtof(p.options[Occ].arg, nullptr);
  params.any = p.options[Any];
  params.print_count = p.options[Count];
  params.no_hydrogens = p.options[NoH];
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (params.verbose > 0 ||
          (p.nonOptionsCount() > 1 && !params.print_count))
        std::printf("%sFile: %s\n", (i > 0 ? "\n" : ""), input.c_str());
      Structure st = read_structure_gz(input);
      print_contacts(st, params);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
