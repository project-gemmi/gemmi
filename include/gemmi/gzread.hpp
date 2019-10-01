// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file.
// To make compilation faster, the implementation is in a separate file.
// If you use this header in your program, in exactly one file you need to:
// #include <gemmi/gzread_impl.hpp>
// See src/input.cpp for an example.

#ifndef GEMMI_GZREAD_HPP_
#define GEMMI_GZREAD_HPP_

#include "cifdoc.hpp" // for Document
#include "model.hpp"  // for Structure

namespace gemmi {

cif::Document read_cif_gz(const std::string& path);
cif::Document read_mmjson_gz(const std::string& path);

inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

Structure make_structure(const cif::Document& doc);

Structure read_structure_gz(const std::string& path,
                            CoorFormat format=CoorFormat::Unknown);

Structure read_pdb_gz(const std::string& path);

CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif
