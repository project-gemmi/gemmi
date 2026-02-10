//! @file
//! @brief Functions for reading possibly gzipped coordinate files.

// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped coordinate files.

#ifndef GEMMI_MMREAD_GZ_HPP_
#define GEMMI_MMREAD_GZ_HPP_

#include "model.hpp"  // for Structure

namespace gemmi {

namespace cif { struct Document; }

//! @brief Read structure from possibly gzipped file.
//! @param path File path (.gz automatically detected)
//! @param format File format (auto-detected if Unknown)
//! @param save_doc Optional pointer to save CIF document
//! @return Structure
GEMMI_DLL Structure read_structure_gz(const std::string& path,
                                      CoorFormat format=CoorFormat::Unknown,
                                      cif::Document* save_doc=nullptr);

//! @brief Read PDB from possibly gzipped file.
//! @param path File path (.gz automatically detected)
//! @param options PDB reading options
//! @return Structure
GEMMI_DLL Structure read_pdb_gz(const std::string& path,
                                PdbReadOptions options=PdbReadOptions());

//! @brief Read structure from chemical component file.
//! @param path File path to chemical component (.gz automatically detected)
//! @param save_doc Optional pointer to save CIF document
//! @param which Which model to read (default 7=ideal coords)
//! @return Structure
GEMMI_DLL Structure read_structure_from_chemcomp_gz(const std::string& path,
                                                    cif::Document* save_doc=nullptr,
                                                    int which=7);

//! @brief Detect coordinate format from file extension.
//! @param path File path (handles .gz extension)
//! @return Detected format
GEMMI_DLL CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif
