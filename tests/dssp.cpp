#include "doctest.h"

// Topo::ChainInfo's constructor is not exported from gemmi DLL on Windows
// (dllexport of Topo does not extend to nested classes).
#if !(defined(_WIN32) && defined(GEMMI_SHARED))

#include <gemmi/dssp.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/pdb.hpp>     // for read_pdb_string
#include <gemmi/topo.hpp>

namespace {

// backbone of the first helical turns of 5i55 (residues 1-8),
// with N-terminal hydrogens H1, H2, H3 and no atom named H
const char* HELIX = R"(
ATOM      1  N   MET A   1       3.333   3.447  27.186  1.00 17.89           N
ATOM      2  CA  MET A   1       3.736   4.598  26.363  1.00 18.57           C
ATOM      3  C   MET A   1       2.932   4.565  25.111  1.00 15.97           C
ATOM      4  O   MET A   1       2.725   3.521  24.522  1.00 15.37           O
ATOM      5  H1  MET A   1       2.590   3.640  27.860  1.00 17.89           H
ATOM      6  H2  MET A   1       4.140   3.190  27.740  1.00 17.89           H
ATOM      7  H3  MET A   1       3.080   2.650  26.590  1.00 17.89           H
ATOM      8  N   GLU A   2       2.464   5.718  24.671  1.00 14.40           N
ATOM      9  CA  GLU A   2       1.798   5.810  23.368  1.00 13.26           C
ATOM     10  C   GLU A   2       2.699   5.494  22.176  1.00 13.21           C
ATOM     11  O   GLU A   2       3.918   5.557  22.268  1.00 12.46           O
ATOM     12  N   PHE A   3       2.082   5.185  21.028  1.00 12.84           N
ATOM     13  CA  PHE A   3       2.801   4.791  19.855  1.00 12.83           C
ATOM     14  C   PHE A   3       3.938   5.753  19.462  1.00 12.67           C
ATOM     15  O   PHE A   3       5.049   5.349  19.179  1.00 12.95           O
ATOM     16  N   VAL A   4       3.609   7.022  19.391  1.00 12.01           N
ATOM     17  CA  VAL A   4       4.591   8.006  18.921  1.00 12.19           C
ATOM     18  C   VAL A   4       5.772   8.060  19.907  1.00 11.62           C
ATOM     19  O   VAL A   4       6.915   8.159  19.520  1.00 12.68           O
ATOM     20  N   ALA A   5       5.483   8.006  21.176  1.00 12.16           N
ATOM     21  CA  ALA A   5       6.519   8.006  22.229  1.00 12.33           C
ATOM     22  C   ALA A   5       7.370   6.745  22.142  1.00 11.41           C
ATOM     23  O   ALA A   5       8.602   6.812  22.191  1.00 13.20           O
ATOM     24  N   LYS A   6       6.734   5.602  21.981  1.00 11.49           N
ATOM     25  CA  LYS A   6       7.539   4.376  21.698  1.00 11.34           C
ATOM     26  C   LYS A   6       8.500   4.495  20.515  1.00 11.09           C
ATOM     27  O   LYS A   6       9.645   3.993  20.534  1.00 11.85           O
ATOM     28  N   LEU A   7       7.975   5.035  19.404  1.00 11.11           N
ATOM     29  CA  LEU A   7       8.742   5.141  18.195  1.00 12.14           C
ATOM     30  C   LEU A   7       9.990   6.021  18.401  1.00 10.85           C
ATOM     31  O   LEU A   7      11.103   5.677  17.988  1.00 12.05           O
ATOM     32  N   PHE A   8       9.782   7.143  19.064  1.00 12.36           N
ATOM     33  CA  PHE A   8      10.901   8.086  19.329  1.00 12.51           C
ATOM     34  C   PHE A   8      11.920   7.433  20.294  1.00 12.15           C
ATOM     35  O   PHE A   8      13.123   7.537  20.053  1.00 12.00           O
)";

// hydrogen bonds found by the DSSP hydrogen-bond pass:
// acceptors[0] and acceptors[1] of each residue, as indices (-1 for none)
std::vector<int> best_acceptors(gemmi::Structure& st, gemmi::HydrogenMode mode) {
  gemmi::Model& model = st.models[0];
  gemmi::NeighborSearch ns(model, st.cell, 9.0);
  gemmi::Chain& chain = model.chains[0];
  for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
    gemmi::Residue& res = chain.residues[n_res];
    if (gemmi::Atom* ca = res.find_atom("CA", '*', gemmi::El::C))
      ns.add_atom(*ca, 0, n_res, int(ca - &res.atoms[0]));
  }
  gemmi::ResidueSpan span = chain.whole();
  gemmi::Topo::ChainInfo ci(span, chain, nullptr);
  // links to the previous residue, as prepare_topology() would add
  for (size_t i = 1; i < ci.res_infos.size(); ++i) {
    ci.res_infos[i].prev.emplace_back();
    ci.res_infos[i].prev.back().res1 = ci.res_infos[i-1].res;
    ci.res_infos[i].prev.back().res2 = ci.res_infos[i].res;
  }
  gemmi::DsspOptions opt;
  opt.hydrogen_mode = mode;
  gemmi::calculate_dssp(ns, ci, opt);
  std::vector<int> result;
  for (const gemmi::Topo::ResInfo& ri : ci.res_infos) {
    result.push_back(ri.acceptors[0] ? int(ri.acceptors[0] - &ci.res_infos[0]) : -1);
    result.push_back(ri.acceptors[1] ? int(ri.acceptors[1] - &ci.res_infos[0]) : -1);
  }
  return result;
}

} // namespace

TEST_CASE("DSSP with existing hydrogens: residue without H (#404)") {
  gemmi::Structure st = gemmi::read_pdb_string(HELIX, "helix");
  // used to dereference a null pointer for residue 1 (H1, H2, H3) and
  // for every other residue (no hydrogens)
  std::vector<int> existing = best_acceptors(st, gemmi::HydrogenMode::Existing);
  std::vector<int> calculated = best_acceptors(st, gemmi::HydrogenMode::Calculate);
  // no atom named H anywhere, so the positions are calculated in both modes
  CHECK(existing == calculated);
  // the N-H of residue 6 binds to the O of residue 2 (i -> i+4 in a helix)
  CHECK_EQ(existing.at(2 * 5 + 1), 1);
}

#endif
