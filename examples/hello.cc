#include <iostream>
#include <set>
#include <gemmi/cif.hpp>

namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
  std::set<std::string> greeted;
  if (argc != 2) return 1;
  try {
    cif::Document doc = cif::read_file(argv[1]);
    cif::Block& block = doc.sole_block(); // mmCIF has exactly one block
    for (const std::string &s : block.find_loop("_atom_site.type_symbol"))
      if (greeted.insert(s).second) // insert() returns pair<iterator,bool>
        std::cout << "Hello " << s << std::endl;
  } catch (std::exception& e) {
    std::cerr << "Oops. " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
