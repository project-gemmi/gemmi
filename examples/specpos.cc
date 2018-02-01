// This program re-calculates _struct_conn.pdbx_dist_value values
// and prints message if it differs by more than 0.002A from the value in file.

#include <gemmi/gz.hpp>   // for MaybeGzipped
#include <gemmi/mmread.hpp>
#include <gemmi/dirwalk.hpp> // for CifWalk
#include <cstdio>
#include <map>

using namespace gemmi;

void check(const char* path) {
  Structure st = gemmi::read_structure(MaybeGzipped(path));
  for (Chain& chain : st.models.at(0).chains)
    for (Residue& res : chain.residues)
      for (Atom& atom : res.atoms)
        if (int n = st.cell.is_special_position(atom.pos)) {
          NearbyImage im = st.cell.find_nearest_image(atom.pos, atom.pos,
                                                      SymmetryImage::Different);
          std::printf("%s: %s %4d %3s %-3s %c *%d  occ=%.2f  B=%5.1f  d=%.4f\n",
              st.name.c_str(),
              chain.name_for_pdb().c_str(), res.seq_id_for_pdb(),
              res.name.c_str(), atom.name.c_str(), (atom.altloc | 0x20),
              n+1, atom.occ, atom.b_iso, std::sqrt(im.dist_sq));
        }
}

int main(int argc, char* argv[]) {
  if (argc != 2)
    return 1;
  int counter = 0;
  try {
    for (const char* path : CoorFileWalk(argv[1])) {
      check(path);
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
// vim:sw=2:ts=2:et
