// Copyright 2017 Global Phasing Ltd.

#include "to_json.hh"

int main(int argc, char **argv) {
  if (argc < 2)
    return 1;
  const char* filename = argv[1];
  gemmi::cif::Document d;
  try {
    d.read_file(filename);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  gemmi::cif::JsonWriter(std::cout).write_json(d);
  return 0;
}

// vim:sw=2:ts=2:et
