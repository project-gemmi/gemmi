// Copyright 2017 Global Phasing Ltd.
//
// Implementation of functions from gzread.hpp
//
// If you use these functions include this file in exactly one compilation unit.

#include "gzread.hpp"
#include "mmread.hpp" // for read_structure
#include "cif.hpp"    // for cif::read
#include "mmcif.hpp"  // for make_structure_from_block
#include "pdb.hpp"    // for read_pdb
#include "json.hpp"   // for read_mmjson
#include "gz.hpp"     // for MaybeGzipped
#include "json.hpp"   // for read_mmjson
#include "util.hpp"   // for ends_with

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
  if ((format == CoorFormat::Unknown || format == CoorFormat::UnknownAny) &&
      ends_with(path, ".gz"))
    format = coordinate_format_from_extension_gz(path);
  return read_structure(MaybeGzipped(path), format);
}

CoorFormat coordinate_format_from_extension_gz(const std::string& path) {
  if (ends_with(path, ".gz"))
    return coordinate_format_from_extension(path.substr(0, path.size() - 3));
  return coordinate_format_from_extension(path);
}

} // namespace gemmi
