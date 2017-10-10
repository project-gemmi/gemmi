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
#include "util.hpp"  // for ends_with

namespace gemmi {

inline CoorFormat coordinate_format_from_extension(const std::string& path) {
  if (iends_with(path, ".pdb") || iends_with(path, ".ent") ||
      iends_with(path, ".pdb.gz") || iends_with(path, ".ent.gz"))
    return CoorFormat::Pdb;
  if (iends_with(path, ".cif") || iends_with(path, ".cif.gz"))
    return CoorFormat::Cif;
  if (iends_with(path, ".json") || iends_with(path, ".js") ||
      iends_with(path, ".json.gz") || iends_with(path, ".js.gz"))
    return CoorFormat::Json;
  return CoorFormat::Unknown;
}

inline Structure read_structure_file(const std::string& path,
                                     CoorFormat format=CoorFormat::Unknown) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension(path);
  switch (format) {
    case CoorFormat::Pdb:  return read_pdb_file(path);
    case CoorFormat::Cif:  return read_atoms(cif::read_file(path));
    case CoorFormat::Json: return read_atoms(cif::read_mmjson_file(path));
    case CoorFormat::Unknown: fail("Unknown format.");
  }
}

template<typename T>
Structure read_structure(T&& input, CoorFormat format=CoorFormat::Unknown) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension(input.path());
  switch (format) {
    case CoorFormat::Pdb:  return read_pdb(input);
    case CoorFormat::Cif:  return read_atoms(cif::read(input));
    case CoorFormat::Json: return read_atoms(cif::read_mmjson(input));
    case CoorFormat::Unknown: fail("Unknown format.");
  }
  fail(""); // avoid GCC5 warning
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
