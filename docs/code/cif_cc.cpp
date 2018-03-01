#include <iostream>
#include <gemmi/cif.hpp>

namespace cif = gemmi::cif;

int main() {
  cif::Document doc = cif::read_file("1mru.cif");
  for (cif::Block& block : doc.blocks)
    for (auto cc : block.find("_chem_comp.", {"id", "formula_weight"}))
      std::cout << cc[0] << " weights " << cc[1] << std::endl;
}

// the next example is in docs from line 15

#include <iostream>
#include <gemmi/cif.hpp>
#include <gemmi/util.hpp>  // for gemmi::join_str
//#include <boost/algorithm/string/join.hpp>  // for boost::algorithm::join

void convert_to_xyz(cif::Document doc) {
  cif::Table atoms = doc.sole_block().find("atom_site_.",
                          {"type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"});
  std::cout << atoms.length() << "\n\n";
  for (auto row : atoms)
    std::cout << gemmi::join_str(row, " ") << "\n";
    //std::cout << boost::algorithm::join(row, " ") << "\n";
}
