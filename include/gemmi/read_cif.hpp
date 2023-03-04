// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CIF files.
// Trivial wrappers that can make compilation faster
// by having a separate implementation file src/read_cif.cpp.

#ifndef GEMMI_READ_CIF_HPP_
#define GEMMI_READ_CIF_HPP_

#include "cifdoc.hpp" // for Document
#include "input.hpp"  // for CharArray

namespace gemmi {

GEMMI_DLL cif::Document read_cif_gz(const std::string& path);
GEMMI_DLL cif::Document read_mmjson_gz(const std::string& path);
GEMMI_DLL CharArray read_into_buffer_gz(const std::string& path);
GEMMI_DLL cif::Document read_cif_from_buffer(const CharArray& buffer, const char* name);
GEMMI_DLL cif::Document read_first_block_gz(const std::string& path, size_t limit);

inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

} // namespace gemmi

#endif
