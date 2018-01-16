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

void check_ssbond(gemmi::cif::Block& block) {
  std::map<std::string, double> distances;
  for (auto row : block.find("_struct_conn.", {"id", "conn_type_id",
                                               "pdbx_dist_value" }))
    if (row.str(1) == "disulf")
      distances[row.str(0)] = cif::as_number(row[2]);
  Structure st = read_atoms_from_block(block);
  for (Connection& con : st.models[0].connections)
    if (con.type == Connection::Disulf) {
      const Atom* cg1 = con.res1->find_by_name_and_elem("SG", El::S);
      if (!cg1)
        std::printf("%s: no SG atom in %s\n", block.name.c_str(),
                    con.res1->ident().c_str());
      const Atom* cg2 = con.res2->find_by_name_and_elem("SG", El::S);
      if (!cg2)
        std::printf("%s: no SG atom in %s\n", block.name.c_str(),
                    con.res2->ident().c_str());
      if (!cg1 || !cg2)
          continue;
      NearestImage near = st.cell.find_nearest_image(cg1->pos, cg2->pos, true);
      double dist = std::sqrt(near.dist_sq);
      double ref_dist = distances[con.id];
      if (std::abs(dist - ref_dist) > 1e-3)
        std::printf("%s %s  %g != %g\n", block.name.c_str(), con.id.c_str(),
                    dist, ref_dist);
    }
}

int main(int argc, char* argv[]) {
  if (argc != 2)
    return 1;
  int counter = 0;
  for (const char* path : CifWalk(argv[1])) {
      cif::Document doc = cif::read(gemmi::MaybeGzipped(path));
      check_ssbond(doc.sole_block());
      if (++counter % 1000 == 0)
        std::printf("[progress: %d files]\n", counter);
  }
  return 0;
}
