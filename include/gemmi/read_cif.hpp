// Copyright 2021 Global Phasing Ltd.
//
/// @file
/// @brief Reading possibly gzip-compressed CIF and JSON files.

#ifndef GEMMI_READ_CIF_HPP_
#define GEMMI_READ_CIF_HPP_

#include "cifdoc.hpp"   // for Document
#include "fileutil.hpp" // for CharArray

namespace gemmi {

/// Read a CIF file, optionally gzip-compressed, from disk.
///
/// @param path    Path to the CIF file (may end with .gz for gzip compression)
/// @param check_level Syntax checking level (0=none, 1=moderate, 2=strict)
/// @return Parsed CIF document
GEMMI_DLL cif::Document read_cif_gz(const std::string& path, int check_level=1);

/// Check CIF syntax without fully parsing the file.
///
/// Performs a quick syntax validation pass on a CIF file (optionally gzipped).
///
/// @param path    Path to the CIF file (may end with .gz for gzip compression)
/// @param msg     If non-null, receives an error message if validation fails
/// @return true if file syntax is valid, false otherwise
GEMMI_DLL bool check_cif_syntax_gz(const std::string& path, std::string* msg);

/// Read an mmJSON file (optionally gzip-compressed) from disk.
///
/// mmJSON is the JSON format used by PDBj for macromolecular CIF data.
///
/// @param path    Path to the mmJSON file (may end with .gz for gzip compression)
/// @return Parsed CIF document
GEMMI_DLL cif::Document read_mmjson_gz(const std::string& path);

/// Read a file into a buffer, optionally decompressing if gzip-compressed.
///
/// @param path    Path to the file (may end with .gz for gzip compression)
/// @return Buffer containing decompressed file contents
GEMMI_DLL CharArray read_into_buffer_gz(const std::string& path);

/// Parse a CIF document from a memory buffer.
///
/// @param data        Pointer to buffer containing CIF data
/// @param size        Number of bytes to read
/// @param name        Optional name for the source (used in error messages)
/// @param check_level Syntax checking level (0=none, 1=moderate, 2=strict)
/// @return Parsed CIF document
GEMMI_DLL cif::Document read_cif_from_memory(const char* data, size_t size, const char* name,
                                             int check_level=1);

/// Read only the first block from a CIF file, optionally gzip-compressed.
///
/// Useful for reading CIF files where only the first block is needed,
/// potentially saving memory and parsing time.
///
/// @param path    Path to the CIF file (may end with .gz for gzip compression)
/// @param limit   Maximum number of bytes to read from the file
/// @return CIF document containing only the first block
GEMMI_DLL cif::Document read_first_block_gz(const std::string& path, size_t limit);

/// Read CIF data from a string.
///
/// This function was moved here from cif.hpp to speed up compilation.
///
/// @param data        CIF-formatted string
/// @param check_level Syntax checking level (0=none, 1=moderate, 2=strict)
/// @return Parsed CIF document
namespace cif {
inline Document read_string(const std::string& data, int check_level=1) {
  return read_cif_from_memory(data.data(), data.size(), "string", check_level);
}
}  // namespace cif

/// Auto-detect and read either CIF or mmJSON format from a file.
///
/// Determines format by file extension (.json, .js for JSON; otherwise CIF).
/// Handles gzip-compressed files transparently.
///
/// @param path Path to the file (may end with .gz for gzip compression)
/// @return Parsed CIF document
inline cif::Document read_cif_or_mmjson_gz(const std::string& path) {
  if (giends_with(path, "json") || giends_with(path, "js"))
    return read_mmjson_gz(path);
  return read_cif_gz(path);
}

} // namespace gemmi

#endif
