// This is a manually run test. It searches for disulfide (SG-SG only)
// bonds in two ways:
// - using a simple search from find_disulfide_bonds() from conn.hpp
// - using cell-linked lists method from subcells.hpp
// and compares the results.
// It was run for all PDB entries to test the cell lists implementation.

#include <gemmi/subcells.hpp>
#include <gemmi/model.hpp>
#include <gemmi/conn.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/dirwalk.hpp>
#include <stdexcept>  // for runtime_error

static int verbose = false;
using namespace gemmi;

// use SubCells (cell list method) to find disulfide SG-SG bonds
static std::vector<Connection> find_disulfide_bonds_cl(const Model& model,
                                                       const UnitCell& cell) {
  const double max_dist = 3.0;
  SubCells sc(model, cell, 5.0);
  const std::string sg = "SG";
  std::vector<Connection> ret;
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg) {
          auto indices = model.get_indices(&chain, &res, &atom);
          sc.for_each(atom.pos, atom.altloc, max_dist,
                      [&](const SubCells::AtomImage& a, float dist_sq) {
              if (a.element != El::S ||
                  indices[0] > a.chain_idx ||
                  (indices[0] == a.chain_idx && indices[1] > a.residue_idx) ||
                  (indices[0] == a.chain_idx && indices[1] == a.residue_idx &&
                   dist_sq < 1))
                return;
              const Chain& chain2 = model.chains[a.chain_idx];
              const Residue& res2 = chain2.residues[a.residue_idx];
              const Atom& atom2 = res2.atoms[a.atom_idx];
              if (atom2.name != sg)
                return;
              Connection c;
              //c.name = "disulf" + std::to_string(ret.size() + 1);
              c.type = Connection::Disulf;
              c.asu = a.image_idx == 0 ? SameAsu::Yes : SameAsu::No;
              c.atom[0] = AtomAddress(chain, res, atom);
              c.atom[1] = AtomAddress(chain2, res2, atom2);
              ret.push_back(c);
          });
        }
  return ret;
}

static void print_connection(const Connection& con, Structure& st) {
  const Atom* a1 = st.models[0].find_atom(con.atom[0]);
  const Atom* a2 = st.models[0].find_atom(con.atom[1]);
  if (a1 && a2) {
    SymImage im = st.cell.find_nearest_image(a1->pos, a2->pos, con.asu);
    std::printf("%s - %s  im:%s  %.3f\n",
                con.atom[0].str().c_str(), con.atom[1].str().c_str(),
                im.pdb_symbol(true).c_str(), im.dist());
  } else {
    std::printf("Ooops, cannot find atom.\n");
  }
}

static void check_disulf(const char* path) {
  Structure st = read_structure(MaybeGzipped(path));
  const Model& model = st.models.at(0);
  std::vector<Connection> c1 = find_disulfide_bonds(model, st.cell);
  std::vector<Connection> c2 = find_disulfide_bonds_cl(model, st.cell);
  printf("%10s  %zu %zu\n", st.name.c_str(), c1.size(), c2.size());
  if (c1.size() != c2.size() || verbose) {
    for (const Connection& con : c1)
      print_connection(con, st);
    printf("---\n");
    for (const Connection& con : c2)
      print_connection(con, st);
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
      for (const char* path : CoorFileWalk(expand_if_pdb_code(argv[pos])))
        check_disulf(path);
  } catch (std::runtime_error& err) {
    std::fprintf(stderr, "Error: %s\n", err.what());
    return 1;
  }
  return 0;
}
