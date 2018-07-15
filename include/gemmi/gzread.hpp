// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file.
// To make compilation faster, the implementation is behind an #ifdef.
// If you use this header in your program, in exactly one file you need
// to have line:
// #define GEMMI_GZREAD_IMPLEMENTATION
// before #including this file. It may be (but doesn't need to be)
// a separate file, such as the src/mmread_gz.hpp file here.

#ifndef GEMMI_GZREAD_HPP_
#define GEMMI_GZREAD_HPP_

#include "gemmi/cifdoc.hpp" // for Document
#include "gemmi/model.hpp"  // for Structure

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

CoorFormat coordinate_format_from_extension_gz(const std::string& path);

} // namespace gemmi

#ifdef GEMMI_GZREAD_IMPLEMENTATION
#include "mmread.hpp" // for read_structure
#include "cif.hpp"    // for cif::read
#include "mmcif.hpp"  // for make_structure_from_block
#include "pdb.hpp"    // for read_pdb
#include "json.hpp"   // for read_mmjson
#include "gz.hpp"     // for MaybeGzipped
#include "json.hpp"   // for read_mmjson
#include "util.hpp"   // for giends_with

namespace gemmi {

cif::Document read_cif_gz(const std::string& path) {
  return cif::read(MaybeGzipped(path));
}

cif::Document read_mmjson_gz(const std::string& path) {
  return cif::read_mmjson(MaybeGzipped(path));
}

Structure make_structure(const cif::Document& doc) {
  return make_structure_from_block(doc.sole_block());
}

Structure read_pdb_gz(const std::string& path) {
  return read_pdb(MaybeGzipped(path));
}

Structure read_structure_gz(const std::string& path, CoorFormat format) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension_gz(path);
  return read_structure(MaybeGzipped(path), format);
}

CoorFormat coordinate_format_from_extension_gz(const std::string& path) {
  if (ends_with(path, ".gz"))
    return coordinate_format_from_extension(path.substr(0, path.size() - 3));
  return coordinate_format_from_extension(path);
}

} // namespace gemmi
#endif // GEMMI_GZREAD_IMPLEMENTATION

#endif
// vim:sw=2:ts=2:et
