// This program re-calculates _struct_conn.pdbx_dist_value values
// and prints message if it differs by more than 0.002A from the value in file.

#include <gemmi/util.hpp> // for giends_with
#include <gemmi/gz.hpp>   // for MaybeGzipped
#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp> // for as_number
#include <gemmi/mmcif.hpp>
#include <gemmi/conn.hpp>
#include <gemmi/dirwalk.hpp> // for CifWalk
#include <cstdio>
#include <map>

using namespace gemmi;

int verbose = false;

void check_struct_conn(cif::Block& block) {
  cif::Table struct_conn = block.find("_struct_conn.", {"id", "conn_type_id",
                                                        "ptnr2_symmetry",
                                                        "pdbx_dist_value" });
  Structure st = make_structure_from_block(block);
  for (Connection& con : st.models[0].connections) {
    const Atom* atom[2] = {nullptr, nullptr};
    for (int n : {0, 1}) {
      atom[n] = st.models[0].find_atom(con.atom[n]);
      if (!atom[n])
        std::printf("%s: %s atom not found in res. %s\n", block.name.c_str(),
                    con.name.c_str(), con.atom[n].str().c_str());
    }
    if (!atom[0] || !atom[1])
      continue;
    NearbyImage im = st.cell.find_nearest_image(atom[0]->pos,
                                                atom[1]->pos, con.image);
    double dist = std::sqrt(im.dist_sq);
    cif::Table::Row row = struct_conn.find_row(con.name);
    if (!starts_with(con.name, row.str(1)))
      std::printf("%s: Unexpected connection name: %s for %s\n",
                  block.name.c_str(), con.name.c_str(), row.str(1).c_str());
    if (dist > 5)
      std::printf("%s: Long connection %s: %g\n",
                  block.name.c_str(), con.name.c_str(), dist);
    std::string ref_sym = row.str(2);
    double ref_dist = cif::as_number(row[3]);
    bool differs = std::abs(dist - ref_dist) > 0.002;
    if (verbose || differs) {
      std::printf("%s %s  im:%s  %.3f %c= %.3f (%s)%s\n",
                  block.name.c_str(), con.name.c_str(),
                  im.pdb_symbol(true).c_str(), dist, (differs ? '!' : '='),
                  ref_dist, ref_sym.c_str(),
                  st.cell.explicit_matrices ? "  {fract}" : "");
    }
  }
  for (cif::Table::Row row : struct_conn)
    if (st.models[0].find_connection_by_name(row.str(0)) == nullptr)
      std::printf("%s: connection not read: %s\n", block.name.c_str(),
                  row.str(0).c_str());
  auto ssbonds = find_disulfide_bonds(st.models[0], st.cell);
  if (verbose) {
    std::printf("\nSearch for disulfide bonds gives:\n");
    for (const Connection& con : ssbonds) {
      bool in_file = false;
      for (Connection& fc : st.models[0].connections) {
        if (fc.type == Connection::Disulf &&
            ((fc.atom[0] == con.atom[0] && fc.atom[1] == con.atom[1]) ||
             (fc.atom[0] == con.atom[1] && fc.atom[1] == con.atom[0]))) {
          in_file = true;
          break;
        }
      }
      if (!in_file)
        std::printf("Connection not in the file.\n");
      const Atom* a1 = st.models[0].find_atom(con.atom[0]);
      const Atom* a2 = st.models[0].find_atom(con.atom[1]);
      if (!a1 || !a2) {
        std::printf("Ooops, cannot find atom.\n");
        continue;
      }
      NearbyImage im = st.cell.find_nearest_image(a1->pos, a2->pos, con.image);
      std::printf("%s  %s  im:%s  %.3f\n",
                  con.atom[0].str().c_str(), con.atom[1].str().c_str(),
                  im.pdb_symbol(true).c_str(), im.dist());
    }
    if (ssbonds.empty())
      std::printf("nothing\n");
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
      cif::Document doc = cif::read(MaybeGzipped(path));
      check_struct_conn(doc.sole_block());
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
