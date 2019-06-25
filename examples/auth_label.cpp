// Compare pairs of columns from the _atom_site table.
// Compiled with: g++-6 -O2 -Iinclude auth_label.cc -lstdc++fs -lz
#include <gemmi/gz.hpp>
#include <gemmi/cif.hpp>
#include <iostream>
#include <experimental/filesystem>  // just <filesystem> in C++17

namespace fs = std::experimental::filesystem;
namespace cif = gemmi::cif;

void print_differences(cif::Block& block, const std::string& name) {
  for (auto row : block.find("_atom_site.", {"auth_" + name, "label_" + name}))
    if (row[0] != row[1])
      std::cout << block.name << ": " << name << "  "
                << row[0] << " -> " << row[1] << std::endl;
}

bool ends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl && str.compare(str.length() - sl, sl, suffix) == 0;
}

int main(int argc, char* argv[]) {
  if (argc != 2)
    return 1;
  for (auto& p : fs::recursive_directory_iterator(argv[1])) {
    std::string path = p.path().u8string();
    if (ends_with(path, ".cif") || ends_with(path, ".cif.gz")) {
      cif::Document doc = cif::read(gemmi::MaybeGzipped(path));
      // What author's atom names were changed?
      print_differences(doc.sole_block(), "atom_id");
      // What author's residue names were changed?
      print_differences(doc.sole_block(), "comp_id");
    }
  }
  return 0;
}
