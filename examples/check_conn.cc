// This program re-calculates _struct_conn.pdbx_dist_value values
// and prints message if it differs by more than 0.002A from the value in file.

#include <gemmi/util.hpp> // for giends_with
#include <gemmi/gz.hpp>   // for MaybeGzipped
#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp> // for as_number
#include <gemmi/mmcif.hpp>
#include <gemmi/conn.hpp>
#include <gemmi/dirwalk.hpp> // for CifWalk
#include <gemmi/mmread.hpp> // for read_structure
#include <cstdio>
#include <map>

using namespace gemmi;

int verbose = false;


static bool has_connection(std::vector<Connection> vec, const Connection& con) {
  for (const Connection& item : vec)
    if (item.type == con.type &&
        (item.image == con.image || item.image == SymmetryImage::Unspecified ||
                                    con.image == SymmetryImage::Unspecified) &&
        ((item.atom[0] == con.atom[0] && item.atom[1] == con.atom[1]) ||
         (item.atom[0] == con.atom[1] && item.atom[1] == con.atom[0])))
      return true;
  return false;
}

static void print_connection(const Connection& con, Structure& st) {
  const Atom* a1 = st.models[0].find_atom(con.atom[0]);
  const Atom* a2 = st.models[0].find_atom(con.atom[1]);
  if (a1 && a2) {
    NearbyImage im = st.cell.find_nearest_image(a1->pos, a2->pos, con.image);
    std::printf("%s - %s  im:%s  %.3f\n",
                con.atom[0].str().c_str(), con.atom[1].str().c_str(),
                im.pdb_symbol(true).c_str(), im.dist());
  } else {
    std::printf("Ooops, cannot find atom.\n");
  }
}

static void check_struct_conn(cif::Block& block) {
  cif::Table struct_conn = block.find("_struct_conn.", {"id", "conn_type_id",
                                                        "ptnr2_symmetry",
                                                        "pdbx_dist_value" });
  Structure st = make_structure_from_block(block);
  size_t disulf_count = 0;
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
      std::printf("%s %-9s %-14s %-14s im:%s  %.3f %c= %.3f (%s)%s\n",
                  block.name.c_str(), con.name.c_str(),
                  con.atom[0].str().c_str(), con.atom[1].str().c_str(),
                  im.pdb_symbol(true).c_str(), dist, (differs ? '!' : '='),
                  ref_dist, ref_sym.c_str(),
                  st.cell.explicit_matrices ? "  {fract}" : "");
    }
    if (con.type == Connection::Disulf)
        ++disulf_count;
  }
  for (cif::Table::Row row : struct_conn)
    if (st.models[0].find_connection_by_name(row.str(0)) == nullptr)
      std::printf("%s: connection not read: %s\n", block.name.c_str(),
                  row.str(0).c_str());
  auto ssbonds = find_disulfide_bonds(st.models[0], st.cell);
  if (disulf_count != ssbonds.size()) {
      printf("%s: S-S bonds: %zu annotated, %zu found by gemmi.\n",
             block.name.c_str(), disulf_count, ssbonds.size());
  }
  for (const Connection& con : ssbonds)
    if (!has_connection(st.models[0].connections, con)) {
      std::printf("%s: missing S-S: ", block.name.c_str());
      print_connection(con, st);
    }
  if (verbose) {
    std::printf("Search for disulfide bonds gives:\n");
    for (const Connection& con : ssbonds)
      print_connection(con, st);
    if (ssbonds.empty())
      std::printf("nothing\n");
    std::printf("\n");
  }
}

static void check_disulf(const char* start) {
  for (const char* path : CoorFileWalk(expand_if_pdb_code(start))) {
    Structure st = read_structure(MaybeGzipped(path));
    if (!st.cell.is_crystal())
      continue;
    const Model& model = st.models.at(0);
    std::vector<Connection> c1 = find_disulfide_bonds(model, st.cell);
    std::vector<Connection> c2 = find_disulfide_bonds2(model, st.cell);
    printf("%10s  %zu %zu\n", st.name.c_str(), c1.size(), c2.size());
    if (c1.size() != c2.size() || verbose) {
      for (const Connection& con : c1)
        print_connection(con, st);
      printf("---\n");
      for (const Connection& con : c2)
        print_connection(con, st);
    }
  }
}

int main(int argc, char* argv[]) {
  int pos = 1;
  if (argc >= 3 && argv[1] == std::string("-v")) {
    ++pos;
    verbose = true;
  } else if (argc < 2) {
    return 1;
  }
  int counter = 0;
  try {
    for (; pos != argc; ++pos) {
      check_disulf(argv[pos]);
      continue;
      for (const char* path : CifWalk(expand_if_pdb_code(argv[pos]))) {
        cif::Document doc = cif::read(MaybeGzipped(path));
        check_struct_conn(doc.sole_block());
        if (++counter % 1000 == 0) {
          std::printf("[progress: %d files]\n", counter);
          std::fflush(stdout);
        }
      }
    }
  } catch (std::runtime_error& err) {
    std::fprintf(stderr, "Error: %s\n", err.what());
    return 1;
  }
  return 0;
}
