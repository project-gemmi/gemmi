// This is a manually run test. It searches for disulfide (SG-SG only)
// bonds in two ways:
// - using a simple search - find_disulfide_bonds()
// - using cell-linked lists method from subcells.hpp
// and compares the results.
// It was run for all PDB entries to test the cell lists implementation.

#include <gemmi/subcells.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/dirwalk.hpp>
#include <stdexcept>  // for runtime_error

static int verbose = false;
using namespace gemmi;

struct BondInfo {
  const_CRA ref;
  SubCells::Mark mark;
  float dist_sq;
};

// Number of SG atoms is relatively small; checking all pairs should be fast.
std::vector<Connection> find_disulfide_bonds(const Model& model,
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
  std::vector<Connection> ret;
  for (size_t i = 0; i < atoms.size(); ++i) {
    const Atom* a1 = atoms[i].atom;
    for (size_t j = i; j < atoms.size(); ++j) {
      const Atom* a2 = atoms[j].atom;
      if (a1->same_conformer(*a2)) {
        Asu asu = (i != j ? Asu::Any : Asu::Different);
        SymImage im = cell.find_nearest_image(a1->pos, a2->pos, asu);
        // if i == j and the image is nearby the atom is on special position
        if (im.dist_sq < max_dist * max_dist && (i != j || im.dist_sq > 1.0)) {
          Connection c;
          c.name = "disulf" + std::to_string(ret.size() + 1);
          c.type = Connection::Disulf;
          c.asu = im.same_asu() ? Asu::Same : Asu::Different;
          c.atom_addr1 = AtomAddress(*atoms[i].chain, *atoms[i].residue, *a1);
          c.atom_addr2 = AtomAddress(*atoms[j].chain, *atoms[j].residue, *a2);
          ret.push_back(c);
        }
      }
    }
  }
  return ret;
}

// use SubCells (cell list method) to find disulfide SG-SG bonds
static std::vector<BondInfo> find_disulfide_bonds_cl(const Model& model,
                                                     const UnitCell& cell) {
  const double max_dist = 3.0;
  SubCells sc(model, cell, 5.0);
  const std::string sg = "SG";
  std::vector<BondInfo> ret;
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg) {
          auto indices = model.get_indices(&chain, &res, &atom);
          sc.for_each(atom.pos, atom.altloc, max_dist,
                      [&](const SubCells::Mark& m, float dist_sq) {
              if (m.element != El::S ||
                  indices[0] > m.chain_idx ||
                  (indices[0] == m.chain_idx && indices[1] > m.residue_idx) ||
                  (indices[0] == m.chain_idx && indices[1] == m.residue_idx &&
                   dist_sq < 1))
                return;
              if (m.to_cra(model).atom->name == sg)
                ret.push_back({{&chain, &res, &atom}, m, dist_sq});
          });
        }
  return ret;
}

static void print_connection(const AtomAddress& a1, const AtomAddress& a2,
                             const SymImage& im) {
  std::printf("%s - %s  im:%s  %.3f\n",
              a1.str().c_str(), a2.str().c_str(),
              im.pdb_symbol(true).c_str(), im.dist());
}

static void check_disulf(const std::string& path) {
  Structure st = read_structure(MaybeGzipped(path));
  const Model& model = st.models.at(0);
  std::vector<Connection> c1 = find_disulfide_bonds(model, st.cell);
  std::vector<BondInfo> c2 = find_disulfide_bonds_cl(model, st.cell);
  printf("%10s  %zu %zu\n", st.name.c_str(), c1.size(), c2.size());
  if (c1.size() != c2.size() || verbose) {
    for (const Connection& con : c1) {
      const Atom* a1 = st.models[0].find_atom(con.atom_addr1);
      const Atom* a2 = st.models[0].find_atom(con.atom_addr2);
      if (!a1 || !a2)
        fail("Ooops, cannot find atom.");
      SymImage im = st.cell.find_nearest_image(a1->pos, a2->pos, con.asu);
      print_connection(con.atom_addr1, con.atom_addr2, im);
    }
    printf("---\n");
    for (const BondInfo& bi : c2) {
      AtomAddress a1(*bi.ref.chain, *bi.ref.residue, *bi.ref.atom);
      CRA cra = bi.mark.to_cra(st.models[0]);
      AtomAddress a2(*cra.chain, *cra.residue, *cra.atom);
      SymImage im = st.cell.find_nearest_pbc_image(bi.ref.atom->pos,
                                                   cra.atom->pos,
                                                   bi.mark.image_idx);
      assert(fabs(bi.dist_sq - im.dist_sq) < 1e-3);
      print_connection(a1, a2, im);
    }
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
      for (std::string path : CoorFileWalk(expand_if_pdb_code(argv[pos])))
        check_disulf(path);
  } catch (std::runtime_error& err) {
    std::fprintf(stderr, "Error: %s\n", err.what());
    return 1;
  }
  return 0;
}
