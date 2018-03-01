#include <iostream>
#include <gemmi/mmread.hpp>
#include <gemmi/gz.hpp>

int main(int argc, char** argv) {
  for (int i = 1; i < argc; ++i)
    try {
      auto st = gemmi::read_structure(gemmi::MaybeGzipped(argv[i]));
      std::cout << "This file has " << st.models.size() << " models.\n";
    } catch (std::runtime_error& e) {
      std::cout << "Oops. " << e.what() << std::endl;
    }
}
