// This program reads PDB files and reports discrepancies between
// symmetry from CRYST1 and symmetry from REMARK 290.

#include <cstdio>
#include <gemmi/dirwalk.hpp> // for PdbWalk
#include <gemmi/gz.hpp>      // for MaybeGzipped
#include <gemmi/pdb.hpp>     // for read_pdb
#include <gemmi/remarks.hpp> // for read_remark_290

void check_remark290(const std::string& path) {
  using namespace gemmi;
  Structure st = read_pdb(MaybeGzipped(path));
  std::vector<Op> ops = read_remark_290(st.raw_remarks);
  if (ops.empty()) {
    if (st.cell.is_crystal())
      printf("no remark 290: %s\n", path.c_str());
    return;
  }
  const SpaceGroup* sg = find_spacegroup_by_ops(split_centering_vectors(ops));
  if (sg) {
    const SpaceGroup* cryst1_sg = st.find_spacegroup();
    if (sg != cryst1_sg) {
      printf("%s\n", path.c_str());
      printf("CRYST1: %s -> %s\n", st.spacegroup_hm.c_str(),
             cryst1_sg ? cryst1_sg->xhm().c_str() : "n/a");
      printf("REMARK 290: %zu symops -> %s\n", ops.size(), sg->xhm().c_str());
    }
  } else {
    printf("failed remark 290 ops to spacegroup: %s\n", path.c_str());
  }
}

int main(int argc, char* argv[]) {
  int counter = 0;
  try {
    for (int pos = 1; pos < argc; ++pos) {
      for (const std::string& path : gemmi::PdbWalk(argv[pos])) {
        check_remark290(path);
        if (++counter % 10000 == 0) {
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
