// prints torsion restraints that are part of a ring

#include <gemmi/cif.hpp>
#include <gemmi/chemcomp.hpp>
//#include <gemmi/to_cif.hpp>  // for write_cif_to_stream
#include <stdio.h>

namespace cif = gemmi::cif;
using AtomId = gemmi::Restraints::AtomId;

void check_torsions(const char* file_path) {
  cif::Document doc = cif::read_file(file_path);
  for (cif::Block& block : doc.blocks) {
    if (block.name.empty() && block.name == "comp_list")
      continue;
    gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
    for (gemmi::Restraints::Torsion& tor : cc.rt.torsions) {
      std::vector<AtomId> ring = cc.rt.find_shortest_path(tor.id4, tor.id1,
                                                          {tor.id2, tor.id3});
      if (ring.empty())
        continue;
      printf("[%s] torsion %s-%s-%s-%s   angle %g +/- %g period %d\n",
             cc.name.c_str(),
             tor.id1.atom.c_str(), tor.id2.atom.c_str(),
             tor.id3.atom.c_str(), tor.id4.atom.c_str(),
             tor.value, tor.esd, tor.period);
      printf("  rest of the ring: %s\n",
             join_str(ring, '-',
                      [](const AtomId& id) { return id.atom; }).c_str());

      // example modifications of torsion angle and esd
      tor.value += 3.5;
      tor.esd = 0.3;
    }
  }
  // writing
  //std::ofstream of("output.cif");
  //of << "# " << file_path << '\n';
  //of << "# modified by ...\n";
  //cif::write_cif_to_stream(of, doc);
  //of.close();
}

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i)
    check_torsions(argv[i]);
  return 0;
}
