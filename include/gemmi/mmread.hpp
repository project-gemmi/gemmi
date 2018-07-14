// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file.

#ifndef GEMMI_MMREAD_HPP_
#define GEMMI_MMREAD_HPP_

#include "model.hpp"     // for Structure
#include "cif.hpp"       // for cif::read
#include "mmcif.hpp"     // for make_structure_from_block
#include "pdb.hpp"       // for read_pdb
#include "json.hpp"      // for read_mmjson
#include "util.hpp"      // for iends_with
#include "input.hpp"     // for JustFile

namespace gemmi {

inline CoorFormat coordinate_format_from_extension(const std::string& path) {
  if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
    return CoorFormat::Pdb;
  if (iends_with(path, ".cif"))
    return CoorFormat::Mmcif;
  if (iends_with(path, ".json"))
    return CoorFormat::Mmjson;
  return CoorFormat::Unknown;
}

template<typename T>
Structure read_structure(T&& input, CoorFormat format=CoorFormat::Unknown) {
  if (format == CoorFormat::Unknown)
    format = coordinate_format_from_extension(input.path());
  switch (format) {
    case CoorFormat::Pdb:
      return read_pdb(input);
    case CoorFormat::Mmcif:
      return make_structure_from_block(cif::read(input).sole_block());
    case CoorFormat::Mmjson:
      return make_structure_from_block(cif::read_mmjson(input).sole_block());
    case CoorFormat::Unknown:
      fail("Unknown format.");
  }
}

inline Structure read_structure_file(const std::string& path,
                                     CoorFormat format=CoorFormat::Unknown) {
  return read_structure(JustFile(path), format);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
