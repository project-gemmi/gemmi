// Copyright 2018 Global Phasing Ltd.
//
// Searches for contacts -- neighbouring atoms.

#include <cstdio>
#include <cstdlib>    // for strtod
#include <algorithm>  // for min, max
#include <gemmi/contact.hpp>
#include <gemmi/neighbor.hpp>
#include "gemmi/assembly.hpp"  // for transform_to_assembly
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/sprintf.hpp>   // for snprintf_z
#define GEMMI_PROG contact
#include "options.h"

namespace {

using namespace gemmi;
using std::printf;

enum OptionIndex { Cov=4, CovMult, MaxDist, Occ, Ignore, NoSym, AsAssembly,
                   NoH, NoWater, NoLigand, Count, Twice, Sort };

const option::Descriptor Usage[] = {
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
  { Occ, 0, "", "minocc", Arg::Float,
    "  --minocc=MIN  \tIgnore atoms with occupancy < MIN." },
  { Ignore, 0, "", "ignore", Arg::Int,
    "  --ignore=N  \tIgnores atom pairs from the same: 0=none, 1=residue, "
    "2=same or adjacent residue, 3=chain, 4=asu." },
  { NoSym, 0, "", "nosym", Arg::None,
    "  --nosym  \tIgnore contacts between symmetry mates." },
  { AsAssembly, 0, "", "assembly", Arg::Required,
    "  --assembly=ID  \tOutput bioassembly with given ID (1, 2, ...)." },
  { NoH, 0, "", "noh", Arg::None,
    "  --noh  \tIgnore hydrogen (and deuterium) atoms." },
  { NoWater, 0, "", "nowater", Arg::None,
    "  --nowater  \tIgnore water." },
  { NoLigand, 0, "", "noligand", Arg::None,
    "  --noligand  \tIgnore ligands and water." },
  { Count, 0, "", "count", Arg::None,
    "  --count  \tPrint only a count of atom pairs." },
  { Twice, 0, "", "twice", Arg::None,
    "  --twice  \tPrint each atom pair A-B twice (A-B and B-A)." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tSort output by distance." },
  { 0, 0, 0, 0, 0, 0 }
};

struct ContactParameters {
  bool use_cov_radius;
  ContactSearch::Ignore ignore = ContactSearch::Ignore::AdjacentResidues;
  bool print_count;
  bool no_hydrogens;
  bool no_symmetry;
  bool twice;
  bool sort;
  float cov_tol = 0.0f;
  float cov_mult = 1.0f;
  float max_dist = 3.0f;
  float min_occ = 0.0f;
  int verbose;
};

void print_contacts(Structure& st, const ContactParameters& params) {
  float max_r = params.use_cov_radius ? 4.f + params.cov_tol : params.max_dist;
  NeighborSearch ns(st.first_model(), st.cell, std::max(5.0f, max_r));
  ns.populate(/*include_h=*/!params.no_hydrogens);

  if (params.verbose > 0) {
    if (params.verbose > 1) {
      if (st.cell.explicit_matrices)
        printf(" Using fractionalization matrix from the file.\n");
      printf(" Each atom has %zu extra images.\n", st.cell.images.size());
      // cf. Structure::setup_cell_images()
      if (const SpaceGroup* sg = st.find_spacegroup()) {
        gemmi::GroupOps group_ops = sg->operations();
        int n = 0;
        for (gemmi::Op op : group_ops)
          printf("  %2d %s\n", ++n, op.triplet().c_str());
        for (const NcsOp& ncs_op : st.ncs)
          if (!ncs_op.given) {
            for (gemmi::Op op : group_ops)
              printf("  %2d NCS %s and %s\n", ++n, ncs_op.id.c_str(), op.triplet().c_str());
          }
      }
    }
    printf(" Cell grid: %d x %d x %d\n", ns.grid.nu, ns.grid.nv, ns.grid.nw);
    size_t min_count = SIZE_MAX, max_count = 0, total_count = 0;
    for (const auto& el : ns.grid.data) {
      min_count = std::min(min_count, el.size());
      max_count = std::max(max_count, el.size());
      total_count += el.size();
    }
    printf(" Items per cell: from %zu to %zu, average: %.2g\n",
           min_count, max_count, double(total_count) / ns.grid.data.size());
  }

  // the code here is similar to LinkHunt::find_possible_links()
  int counter = 0;
  ContactSearch contacts(max_r);
  contacts.twice = params.twice;
  contacts.ignore = params.ignore;
  if (params.use_cov_radius)
    contacts.setup_atomic_radii(params.cov_mult, params.cov_tol);
  std::multimap<double, std::string> lines;
  char buf[256];
  contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                    int image_idx, double dist_sq) {
      ++counter;
      if (params.print_count)
        return;
      std::string sym1, sym2;
      if (!params.no_symmetry) {
        NearestImage im = st.cell.find_nearest_pbc_image(cra1.atom->pos,
                                                         cra2.atom->pos, image_idx);
        sym1 = "1555";
        sym2 = im.symmetry_code(false);
      }
      std::string conn_info;
      if (Connection* conn = st.find_connection_by_cra(cra1, cra2))
        conn_info = conn->name.empty() ? "(link)" : conn->name;
      snprintf_z(buf, 255, "%-11s %-4s%c%3s%2s%4s%c         "
                           "      %-4s%c%3s%2s%4s%c  %6s %6s %5.2f\n",
             conn_info.c_str(),
             cra1.atom->padded_name().c_str(),
             cra1.atom->altloc ? std::toupper(cra1.atom->altloc) : ' ',
             cra1.residue->name.c_str(),
             cra1.chain->name.c_str(),
             cra1.residue->seqid.num.str().c_str(), cra1.residue->seqid.icode,
             cra2.atom->padded_name().c_str(),
             cra2.atom->altloc ? std::toupper(cra2.atom->altloc) : ' ',
             cra2.residue->name.c_str(),
             cra2.chain->name.c_str(),
             cra2.residue->seqid.num.str().c_str(), cra2.residue->seqid.icode,
             sym1.c_str(), sym2.c_str(), std::sqrt(dist_sq));
      if (params.sort)
        lines.emplace(dist_sq, buf);
      else
        printf("%s", buf);
  });
  if (params.sort)
    for (const auto& it : lines)
      printf("%s", it.second.c_str());
  if (params.print_count)
    printf("%s:%g\n", st.name.c_str(), 0.5 * counter);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  ContactParameters params;
  params.verbose = p.options[Verbose].count();
  params.use_cov_radius = (p.options[Cov] || p.options[CovMult]);
  if (p.options[Cov])
    params.cov_tol = std::strtof(p.options[Cov].arg, nullptr);
  if (p.options[CovMult])
    params.cov_mult = std::strtof(p.options[CovMult].arg, nullptr);
  if (p.options[MaxDist])
    params.max_dist = std::strtof(p.options[MaxDist].arg, nullptr);
  if (p.options[Occ])
    params.min_occ = std::strtof(p.options[Occ].arg, nullptr);
  if (p.options[Ignore]) {
    int ignore_level = std::atoi(p.options[Ignore].arg);
    if (ignore_level < 0 || ignore_level > 4) {
      std::fprintf(stderr, "Error: value of --ignore is out of range.\n");
      return 1;
    }
    params.ignore = (ContactSearch::Ignore) ignore_level;
  }
  params.print_count = p.options[Count];
  params.no_hydrogens = p.options[NoH];
  params.no_symmetry = p.options[NoSym];
  params.twice = p.options[Twice];
  params.sort = p.options[Sort];
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (params.verbose > 0 ||
          (p.nonOptionsCount() > 1 && !params.print_count))
        std::printf("%sFile: %s\n", (i > 0 ? "\n" : ""), input.c_str());
      Structure st = read_structure_gz(input);
      if (params.ignore == ContactSearch::Ignore::AdjacentResidues)
        setup_entities(st);
      if (p.options[NoWater])
        remove_waters(st);
      if (p.options[NoLigand])
        remove_ligands_and_waters(st);
      if (p.options[AsAssembly])
        transform_to_assembly(st, p.options[AsAssembly].arg,
                              HowToNameCopiedChain::Short, nullptr);
      if (params.no_symmetry || p.options[AsAssembly])
        st.cell = UnitCell();
      print_contacts(st, params);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
