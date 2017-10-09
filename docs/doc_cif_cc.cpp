#include <iostream>
#include <gemmi/cif.hpp>

namespace cif = gemmi::cif;

int main() {
  cif::Document doc = cif::read_file("1mru.cif");
  for (const cif::Block& block : doc.blocks)
    for (const auto& cc : block.find("_chem_comp.", {"id", "formula_weight"}))
      std::cout << cc[0] << " weights " << cc[1] << std::endl;
}
