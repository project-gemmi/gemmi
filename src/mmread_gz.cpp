// Copyright 2021 Global Phasing Ltd.

#include <gemmi/mmread_gz.hpp>
#include <gemmi/mmread.hpp> // for read_structure
#include <gemmi/pdb.hpp>    // for read_pdb
#include <gemmi/gz.hpp>     // for MaybeGzipped
#include <gemmi/read_cif.hpp>  // for read_cif_gz

namespace gemmi {

Structure read_structure_gz(const std::string& path, CoorFormat format,
                            cif::Document* save_doc) {
  return read_structure(MaybeGzipped(path), format, save_doc);
}

Structure read_pdb_gz(const std::string& path, PdbReadOptions options) {
  return read_pdb(MaybeGzipped(path), options);
}

Structure read_structure_from_chemcomp_gz(const std::string& path,
                                          cif::Document* save_doc, int which) {
  return make_structure_from_chemcomp_doc(read_cif_gz(path), save_doc, which);
}

CoorFormat coor_format_from_ext_gz(const std::string& path) {
  return coor_format_from_ext(MaybeGzipped(path).basepath());
}

} // namespace gemmi
