// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CIF files.

#ifndef GEMMI_READ_CIF_HPP_
#define GEMMI_READ_CIF_HPP_

#include "cifdoc.hpp" // for Document
#include "input.hpp"  // for CharArray

namespace gemmi {

GEMMI_DLL cif::Document read_cif_gz(const std::string& path);
GEMMI_DLL bool check_cif_syntax_gz(const std::string& path, std::string* msg);
GEMMI_DLL cif::Document read_mmjson_gz(const std::string& path);
GEMMI_DLL CharArray read_into_buffer_gz(const std::string& path);
GEMMI_DLL cif::Document read_cif_from_memory(const char* data, size_t size, const char* name);
GEMMI_DLL cif::Document read_first_block_gz(const std::string& path, size_t limit);

// cif::read_string() was moved here from cif.hpp to speed up compilation
namespace cif {
inline Document read_string(const std::string& data) {
  return read_cif_from_memory(data.data(), data.size(), "string");
}
}  // namespace cif

inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

} // namespace gemmi

#endif
