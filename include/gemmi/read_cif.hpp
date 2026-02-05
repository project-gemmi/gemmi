//! @file
//! @brief Functions for reading possibly gzipped CIF files.

// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CIF files.

#ifndef GEMMI_READ_CIF_HPP_
#define GEMMI_READ_CIF_HPP_

#include "cifdoc.hpp"   // for Document
#include "fileutil.hpp" // for CharArray

namespace gemmi {

//! @brief Read CIF from possibly gzipped file.
//! @param path File path (.gz automatically detected)
//! @param check_level Syntax check level (0=none, 1=basic, 2=strict)
//! @return CIF document
GEMMI_DLL cif::Document read_cif_gz(const std::string& path, int check_level=1);

//! @brief Check CIF syntax in possibly gzipped file.
//! @param path File path (.gz automatically detected)
//! @param msg Optional pointer to store error message
//! @return True if syntax is valid
GEMMI_DLL bool check_cif_syntax_gz(const std::string& path, std::string* msg);

//! @brief Read mmJSON from possibly gzipped file.
//! @param path File path (.gz automatically detected)
//! @return CIF document
GEMMI_DLL cif::Document read_mmjson_gz(const std::string& path);

//! @brief Read file into buffer.
//! @param path File path (.gz automatically detected)
//! @return Character array buffer
GEMMI_DLL CharArray read_into_buffer_gz(const std::string& path);

//! @brief Read CIF from memory buffer.
//! @param data Pointer to data
//! @param size Data size
//! @param name Document name
//! @param check_level Syntax check level
//! @return CIF document
GEMMI_DLL cif::Document read_cif_from_memory(const char* data, size_t size, const char* name,
                                             int check_level=1);

//! @brief Read first block from possibly gzipped CIF.
//! @param path File path (.gz automatically detected)
//! @param limit Size limit for reading
//! @return CIF document with first block
GEMMI_DLL cif::Document read_first_block_gz(const std::string& path, size_t limit);

//! cif::read_string() was moved here from cif.hpp to speed up compilation
namespace cif {
//! @brief Read CIF from string.
//! @param data CIF data string
//! @param check_level Syntax check level
//! @return CIF document
inline Document read_string(const std::string& data, int check_level=1) {
  return read_cif_from_memory(data.data(), data.size(), "string", check_level);
}
}  // namespace cif

//! @brief Read CIF or mmJSON based on extension.
//! @param path File path
//! @return CIF document
inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

} // namespace gemmi

#endif
