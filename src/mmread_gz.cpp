// Copyright 2021 Global Phasing Ltd.

#include <gemmi/mmread_gz.hpp>
#include <gemmi/mmread.hpp> // for read_structure
#include <gemmi/pdb.hpp>    // for read_pdb
#include <gemmi/gz.hpp>     // for MaybeGzipped

namespace gemmi {

Structure read_pdb_gz(const std::string& path, PdbReadOptions options) {
  return read_pdb(MaybeGzipped(path), options);
}

Structure read_structure_gz(const std::string& path, CoorFormat format,
                            cif::Document* save_doc) {
  return read_structure(MaybeGzipped(path), format, save_doc);
}

CoorFormat coor_format_from_ext_gz(const std::string& path) {
  return coor_format_from_ext(MaybeGzipped(path).basepath());
}

} // namespace gemmi
