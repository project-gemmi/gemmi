// This is a manually run test. It searches for disulfide (SG-SG only)
// bonds in two ways:
// - using a simple search - find_disulfide_bonds1()
// - using cell-linked lists method from neighbor.hpp - find_disulfide_bonds2()
// and compares the results.
// It was run for all PDB entries to test the cell lists implementation.

#include <gemmi/neighbor.hpp>
#include <gemmi/contact.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/dirwalk.hpp>
#include <stdexcept>  // for runtime_error
#include <chrono>

static int verbose = false;
using namespace gemmi;

struct BondInfo {
  const_CRA cra1, cra2;
  int image_idx;
  double dist_sq;

  void print(const UnitCell& cell) const {
    NearestImage im = cell.find_nearest_pbc_image(cra1.atom->pos, cra2.atom->pos, image_idx);
    assert(fabs(dist_sq - im.dist_sq) < 1e-3);
    std::printf("%s - %s  im:%s  %.3f\n",
                atom_str(cra1).c_str(), atom_str(cra2).c_str(),
                im.symmetry_code(true).c_str(), std::sqrt(dist_sq));
  }

};

// Number of SG atoms is relatively small; checking all pairs should be fast.
std::vector<BondInfo> find_disulfide_bonds1(const Model& model,
                                              const UnitCell& cell) {
  const double max_dist = 3.0;
  // Find all SG sulfur atoms.
  std::vector<const_CRA> atoms;
  const std::string sg = "SG";
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg)
          atoms.push_back(const_CRA{&chain, &res, &atom});
  // Check distances.
  std::vector<BondInfo> ret;
  for (size_t i = 0; i < atoms.size(); ++i) {
    const Atom* a1 = atoms[i].atom;
    for (size_t j = i; j < atoms.size(); ++j) {
      const Atom* a2 = atoms[j].atom;
      if (a1->same_conformer(*a2)) {
        Asu asu = (i != j ? Asu::Any : Asu::Different);
        NearestImage im = cell.find_nearest_image(a1->pos, a2->pos, asu);
        // if i == j and the image is nearby the atom is on special position
        if (im.dist_sq < max_dist * max_dist && (i != j || im.dist_sq > 1.0))
          ret.push_back({atoms[i], atoms[j], im.sym_idx, im.dist_sq});
      }
    }
  }
  return ret;
}

// use NeighborSearch (cell list method) to find disulfide SG-SG bonds
static std::vector<BondInfo> find_disulfide_bonds2(Model& model,
                                                   const UnitCell& cell) {
  const double max_dist = 3.0;
  const std::string sg = "SG";
  NeighborSearch ns(model, cell, 5.0);
#if 0
  ns.populate();
#else
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg) {
          auto indices = model.get_indices(&chain, &res, &atom);
          ns.add_atom(atom, indices[0], indices[1], indices[2]);
        }
#endif
  std::vector<BondInfo> ret;
#if 0  // faster, but requires more code
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg) {
          auto indices = model.get_indices(&chain, &res, &atom);
          ns.for_each(atom.pos, atom.altloc, max_dist,
                      [&](const NeighborSearch::Mark& m, double dist_sq) {
              if (m.element != El::S)
                return;
              if (indices[0] > m.chain_idx || (indices[0] == m.chain_idx &&
                    (indices[1] > m.residue_idx ||
                     (indices[1] == m.residue_idx && dist_sq < 1))))
                return;
              const_CRA cra1 = {&chain, &res, &atom};
              const_CRA cra2 = m.to_cra(model);
              if (cra2.atom->name == sg)
                ret.push_back({cra1, cra2, m.image_idx, dist_sq});
          });
        }
#else  // slower, but simpler
  ContactSearch contacts(max_dist);
  contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                    int image_idx, double dist_sq) {
      if (cra1.atom->element == El::S && cra1.atom->name == sg &&
          cra2.atom->element == El::S && cra2.atom->name == sg)
        ret.push_back({cra1, cra2, image_idx, dist_sq});
  });
#endif
  return ret;
}

static void check_disulf(const std::string& path) {
  if (verbose)
    printf("path: %s\n", path.c_str());
  Structure st = read_structure(MaybeGzipped(path));
  Model& model = st.first_model();
  using Clock = std::chrono::steady_clock;
  auto start = Clock::now();
  std::vector<BondInfo> c1 = find_disulfide_bonds1(model, st.cell);
  std::chrono::duration<double> elapsed1 = Clock::now() - start;
  start = Clock::now();
  std::vector<BondInfo> c2 = find_disulfide_bonds2(model, st.cell);
  std::chrono::duration<double> elapsed2 = Clock::now() - start;
  printf("%10s  %zu %zu  (%.3g %.3g s)\n", st.name.c_str(),
         c1.size(), c2.size(), elapsed1.count(), elapsed2.count());
  if (c1.size() != c2.size() || verbose) {
    for (const BondInfo& bi : c1)
      bi.print(st.cell);
    printf("---\n");
    for (const BondInfo& bi : c2)
      bi.print(st.cell);
  }
}

int main(int argc, char* argv[]) {
  int pos = 1;
  while (pos < argc && argv[pos] == std::string("-v")) {
    ++pos;
    verbose = true;
  }
  if (pos >= argc) {
    std::fprintf(stderr, "No paths given.\n");
    return 1;
  }
  try {
    for (; pos != argc; ++pos)
      for (std::string path : CoorFileWalk(argv[pos], 'M'))
        check_disulf(path);
  } catch (std::runtime_error& err) {
    std::fprintf(stderr, "Error: %s\n", err.what());
    return 1;
  }
  return 0;
}
