// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped coordinate files.

#ifndef GEMMI_MMREAD_GZ_HPP_
#define GEMMI_MMREAD_GZ_HPP_

#include "model.hpp"  // for Structure

namespace gemmi {

namespace cif { struct Document; }

GEMMI_DLL Structure read_structure_gz(const std::string& path,
                                      CoorFormat format=CoorFormat::Unknown,
                                      cif::Document* save_doc=nullptr);

GEMMI_DLL Structure read_pdb_gz(const std::string& path,
                                PdbReadOptions options=PdbReadOptions());

GEMMI_DLL Structure read_structure_from_chemcomp_gz(const std::string& path,
                                                    cif::Document* save_doc=nullptr,
                                                    int which=7);

GEMMI_DLL CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif
