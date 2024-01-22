// Copyright 2017 Global Phasing Ltd.
//
// Reading CIF-JSON (COMCIFS) and mmJSON (PDBj) formats into cif::Document.

#ifndef GEMMI_JSON_HPP_
#define GEMMI_JSON_HPP_

#include <utility>      // for forward
#include "cifdoc.hpp"   // for Document, etc
#include "fileutil.hpp" // for read_file_into_buffer

namespace gemmi {
namespace cif {

// reads mmJSON file mutating the input buffer as a side effect
GEMMI_DLL Document read_mmjson_insitu(char* buffer, std::size_t size,
                                      const std::string& name="mmJSON");

inline Document read_mmjson_file(const std::string& path) {
  CharArray buffer = read_file_into_buffer(path);
  return read_mmjson_insitu(buffer.data(), buffer.size(), path);
}

template<typename T>
Document read_mmjson(T&& input) {
  std::string name = input.is_stdin() ? "stdin" : input.path();
  CharArray buffer = read_into_buffer(std::forward<T>(input));
  return read_mmjson_insitu(buffer.data(), buffer.size(), name);
}

} // namespace cif
} // namespace gemmi
#endif
