// Copyright 2017 Global Phasing Ltd.

#include "to_json.hh"
#include "write_cif.hh"

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3)
    return 1;
  const char* filename = argv[1];
  gemmi::cif::Document d;
  try {
    d.read_file(filename);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  // temporary - for testing write_cif
  if (argc > 2) {
    std::string output = argv[2];
    std::ofstream of(output);
    if (!of)
      throw std::runtime_error("Failed to open " + output);
    if (output.size() > 4 && output.substr(output.size() - 4) == ".cif")
      of << d;
    else
      gemmi::cif::JsonWriter(of).write_json(d);
    of.close();
  } else {
    gemmi::cif::JsonWriter(std::cout).write_json(d);
  }
  return 0;
}

// vim:sw=2:ts=2:et
