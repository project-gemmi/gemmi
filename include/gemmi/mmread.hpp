// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file.

#ifndef GEMMI_MMREAD_HPP_
#define GEMMI_MMREAD_HPP_

#include "model.hpp"
#include "cif.hpp"
#include "mmcif.hpp"
#include "pdb.hpp"
#include "json.hpp"  // mmJSON
#include "util.hpp"  // for giends_with

namespace gemmi {

inline CoorFormat coordinate_format_from_extension(const std::string& path) {
  if (giends_with(path, ".pdb") || giends_with(path, ".ent"))
    return CoorFormat::Pdb;
  if (giends_with(path, ".cif"))
    return CoorFormat::Mmcif;
  if (giends_with(path, ".json") || giends_with(path, ".js"))
    return CoorFormat::Mmjson;
  return CoorFormat::Unknown;
}

inline Structure read_structure_file(const std::string& path,
                                     CoorFormat format=CoorFormat::Unknown) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension(path);
  switch (format) {
    case CoorFormat::Pdb:  return read_pdb_file(path);
    case CoorFormat::Mmcif:  return make_structure(cif::read_file(path));
    case CoorFormat::Mmjson: return make_structure(cif::read_mmjson_file(path));
    case CoorFormat::Unknown: fail("Unknown format.");
  }
}

template<typename T>
Structure read_structure(T&& input, CoorFormat format=CoorFormat::Unknown) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension(input.path());
  switch (format) {
    case CoorFormat::Pdb:  return read_pdb(input);
    case CoorFormat::Mmcif:  return make_structure(cif::read(input));
    case CoorFormat::Mmjson: return make_structure(cif::read_mmjson(input));
    case CoorFormat::Unknown: fail("Unknown format.");
  }
  fail(""); // avoid GCC5 warning
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
