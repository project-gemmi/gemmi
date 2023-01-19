// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped coordinate files.
// Trivial wrappers that can make compilation faster
// by having a separate implementation file src/mmread_gz.cpp.

#ifndef GEMMI_MMREAD_GZ_HPP_
#define GEMMI_MMREAD_GZ_HPP_

#include "model.hpp"  // for Structure

namespace gemmi {

namespace cif { struct Document; }

Structure read_structure_gz(const std::string& path, CoorFormat format=CoorFormat::Unknown,
                            cif::Document* save_doc=nullptr);

Structure read_pdb_gz(const std::string& path, PdbReadOptions options=PdbReadOptions());

CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif
