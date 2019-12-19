// Copyright 2018 Global Phasing Ltd.
//
// Searches for contacts -- neighbouring atoms.

#include <cstdio>
#include <cstdlib>    // for strtod
#include <algorithm>  // for min, max
#include <gemmi/subcells.hpp>
#include <gemmi/polyheur.hpp>  // for are_connected
#include <gemmi/elem.hpp>      // for is_hydrogen
#include <gemmi/gzread.hpp>
#include <gemmi/to_pdb.hpp>    // for padded_atom_name
#define GEMMI_PROG contact
#include "options.h"

using namespace gemmi;
using std::printf;

enum OptionIndex { Cov=4, CovMult, MaxDist, Occ, Any, NoH, NoSym, Count };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nSearches for contacts in a model (PDB or mmCIF)."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { MaxDist, 0, "d", "maxdist", Arg::Float,
    "  -d, --maxdist=D  Maximal distance in A (default 3.0)" },
  { Cov, 0, "", "cov", Arg::Float,
    "  --cov=TOL  \tUse max distance = covalent radii sum + TOL [A]." },
  { CovMult, 0, "", "covmult", Arg::Float,
    "  --covmult=M  \tUse max distance = M * covalent radii sum + TOL [A]." },
  { Occ, 0, "", "occsum", Arg::Float,
    "  --occsum=MIN  \tIgnore atom pairs with summed occupancies < MIN." },
  { Any, 0, "", "any", Arg::None,
    "  --any  \tOutput any atom pair, even from the same residue." },
  { NoH, 0, "", "noh", Arg::None,
    "  --noh  \tIgnore hydrogen (and deuterium) atoms." },
  { NoSym, 0, "", "nosym", Arg::None,
    "  --nosym  \tIgnore contacts with symmetry mates." },
  { Count, 0, "", "count", Arg::None,
    "  --count  \tPrint only a count of atom pairs." },
  { 0, 0, 0, 0, 0, 0 }
};

struct Parameters {
  bool use_cov_radius;
  bool any;
  bool print_count;
  bool no_hydrogens;
  bool no_symmetry;
  float cov_tol = 0.0f;
  float cov_mult = 1.0f;
  float max_dist = 3.0f;
  float occ_sum = 0.0f;
  int verbose;
};

static void print_contacts(const Structure& st, const Parameters& params) {
  const float special_pos_cutoff = 0.8f;
  float max_r = params.use_cov_radius ? params.max_dist : 4.f + params.cov_tol;
  SubCells sc(st.models.at(0), st.cell, std::max(5.0f, max_r));
  sc.populate(/*include_h=*/!params.no_hydrogens);

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
        float d_part = 0.f;
        if (params.use_cov_radius) {
          d_part = params.cov_mult * atom.element.covalent_r() + params.cov_tol;
          max_r = d_part + 2.4f;
        }
        sc.for_each(atom.pos, atom.altloc, max_r,
                    [&](const SubCells::Mark& m, float dist_sq) {
            if (!params.any) {
              // do not consider connections inside a residue
              // or between this and previous/next residue
              if (m.image_idx == 0 && m.chain_idx == n_ch &&
                  (m.residue_idx == n_res ||
                   are_connected(res, chain.residues[m.residue_idx], pt) ||
                   are_connected(chain.residues[m.residue_idx], res, pt)))
                return;
            }
            // atom can be linked with its image, but if the image
            // is too close the atom is likely on special position.
            if (m.chain_idx == n_ch && m.residue_idx == n_res &&
                m.atom_idx == n_atom && (m.image_idx == 0 ||
                                         dist_sq < sq(special_pos_cutoff)))
              return;
            const_CRA cra = m.to_cra(model);
            if (cra.atom->occ < min_occ)
              return;
            if (params.use_cov_radius) {
              float limit = d_part + params.cov_mult *
                                     cra.atom->element.covalent_r();
              if (limit < 0 || dist_sq > sq(limit))
                return;
            }
            ++counter;
            if (params.print_count)
              return;
            const char* sym1;
            std::string sym2;
            double dist;
            if (params.no_symmetry) {
              sym1 = "";
              dist = std::sqrt(dist_sq);
            } else {
              SymImage im = st.cell.find_nearest_pbc_image(
                                          atom.pos, cra.atom->pos, m.image_idx);
              sym1 = "1555";
              sym2 = im.pdb_symbol(false);
              dist = im.dist();
            }
            printf("            %-4s%c%3s%2s%5s   "
                   "            %-4s%c%3s%2s%5s  %6s %6s %5.2f\n",
                   padded_atom_name(atom).c_str(),
                   atom.altloc ? std::toupper(atom.altloc) : ' ',
                   res.name.c_str(),
                   chain.name.c_str(),
                   res.seqid.str().c_str(),
                   padded_atom_name(*cra.atom).c_str(),
                   cra.atom->altloc ? std::toupper(cra.atom->altloc) : ' ',
                   cra.residue->name.c_str(),
                   cra.chain->name.c_str(),
                   cra.residue->seqid.str().c_str(),
                   sym1, sym2.c_str(), dist);
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
  params.use_cov_radius = (p.options[Cov] || p.options[CovMult]);
  if (p.options[Cov])
    params.cov_tol = std::strtof(p.options[Cov].arg, nullptr);
  if (p.options[CovMult])
    params.cov_mult = std::strtof(p.options[CovMult].arg, nullptr);
  if (p.options[MaxDist])
    params.max_dist = std::strtof(p.options[MaxDist].arg, nullptr);
  if (p.options[Occ])
    params.occ_sum = std::strtof(p.options[Occ].arg, nullptr);
  params.any = p.options[Any];
  params.print_count = p.options[Count];
  params.no_hydrogens = p.options[NoH];
  params.no_symmetry = p.options[NoSym];
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (params.verbose > 0 ||
          (p.nonOptionsCount() > 1 && !params.print_count))
        std::printf("%sFile: %s\n", (i > 0 ? "\n" : ""), input.c_str());
      Structure st = read_structure_gz(input);
      if (params.no_symmetry)
        st.cell = UnitCell();
      print_contacts(st, params);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
