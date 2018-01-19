// check disulfide bonds
#include <gemmi/util.hpp> // for giends_with
#include <gemmi/gz.hpp>   // for MaybeGzipped
#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp> // for as_number
#include <gemmi/mmcif.hpp>
#include <gemmi/dirwalk.hpp> // for CifWalk
#include <cstdio>
#include <map>

using namespace gemmi;

int verbose = false;

inline const Atom* find_ssbond_atom(Connection& con, int n) {
  if (!con.res[n])
    return nullptr;
  return con.res[n]->find_by_name_altloc_elem(con.atom[n], con.altloc[n],
                                              El::S);
}

void check_ssbond(gemmi::cif::Block& block) {
  cif::Table struct_conn = block.find("_struct_conn.", {"id", "conn_type_id",
                                                        "ptnr2_symmetry",
                                                        "pdbx_dist_value" });
  Structure st = read_atoms_from_block(block);
  //TODO: check that no atom is in 2 connections?
  for (Connection& con : st.models[0].connections)
    if (con.type == Connection::Disulf) {
      const Atom* atom[2];
      for (int n : {0, 1}) {
        atom[n] = find_ssbond_atom(con, n);
        if (!atom[n])
          std::printf("%s: atom not found in %s\n", block.name.c_str(),
                      con.res[n]->ident().c_str());
      }
      if (!atom[0] || !atom[1])
        continue;
      NearestImage im = st.cell.find_nearest_image(atom[0]->pos,
                                                   atom[1]->pos, true);
      double dist = std::sqrt(im.dist_sq);
      cif::Table::Row row = struct_conn.find_row(con.id);
      assert(row.str(1) == "disulf");
      std::string ref_sym = row.str(2);
      double ref_dist = cif::as_number(row[3]);
      if (verbose || std::abs(dist - ref_dist) > 0.002) {
        std::printf("%s %s  im:%s  %.3f != %.3f (%s)\n",
                    block.name.c_str(), con.id.c_str(),
                    im.pdb_symbol(true).c_str(), dist,
                    ref_dist, ref_sym.c_str());
      }
    }
}

int main(int argc, char* argv[]) {
  if (argc == 3 && argv[1] == std::string("-v"))
    verbose = true;
  else if (argc != 2)
    return 1;
  int counter = 0;
  try {
    for (const char* path : CifWalk(argv[argc-1])) {
      cif::Document doc = cif::read(gemmi::MaybeGzipped(path));
      check_ssbond(doc.sole_block());
      if (++counter % 1000 == 0) {
        std::printf("[progress: %d files]\n", counter);
        std::fflush(stdout);
      }
    }
  } catch (std::runtime_error& err) {
    std::fprintf(stderr, "Error: %s\n", err.what());
    return 1;
  }
  return 0;
}
