//! @file
//! @brief Read the PDB file format and store it in Structure.
//!
//! Based on the format spec:
//! https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
//!
//! Extensions beyond standard PDB format:
//! - support for two-character chain IDs (columns 21 and 22)
//! - read segment ID (columns 73-76)
//! - read hybrid-36 serial numbers (http://cci.lbl.gov/hybrid_36/)
//! - hybrid-36 sequence id for sequences longer than 9999 (no such examples)
//! - allow for longer REMARK lines (up to 120 characters)

#ifndef GEMMI_PDB_HPP_
#define GEMMI_PDB_HPP_

#include "input.hpp"    // for FileStream
#include "model.hpp"    // for Structure, ...

namespace gemmi {

//! @brief Compare first 4 letters of string with record type (case-insensitive).
//! @param s String to check (must have at least 4 characters)
//! @param record Record type to compare against (uppercase, 4 letters)
//! @return True if match
//!
//! Both args must have at least 3+1 chars. ' ' and NUL are equivalent in s.
inline bool is_record_type4(const char* s, const char* record) {
  return ialpha4_id(s) == ialpha4_id(record);
}

//! @brief Compare first 3 letters of string with record type (case-insensitive).
//! @param s String to check (must have at least 3 characters)
//! @param record Record type to compare against (uppercase, 3 letters + space)
//! @return True if match
//!
//! For record "TER": "TER ", TER\n, TER\r, TER\t match, but TERE, TER1 don't.
inline bool is_record_type3(const char* s, const char* record) {
  return (ialpha4_id(s) & ~0xf) == ialpha4_id(record);
}

//! @brief Parse REMARK 290 to get symmetry operations.
//! @param raw_remarks Vector of REMARK lines
//! @return Operations corresponding to 1555, 2555, ... N555
GEMMI_DLL std::vector<Op> read_remark_290(const std::vector<std::string>& raw_remarks);

//! @brief Populate Structure from PDB stream.
//! @param line_reader Input stream providing PDB lines
//! @param source Source name/path (for error messages)
//! @param st Structure to populate
//! @param options Reading options (skip hydrogen, etc.)
//!
//! Main function for parsing PDB format into a Structure object.
GEMMI_DLL void populate_structure_from_pdb_stream(AnyStream& line_reader,
                                                  const std::string& source,
                                                  Structure& st,
                                                  PdbReadOptions options);

//! @brief Populate Structure from PDB data in memory.
//! @param data Pointer to PDB file content in memory
//! @param size Size of data in bytes
//! @param name Source name (for error messages)
//! @param st Structure to populate
//! @param options Reading options
inline void populate_structure_from_pdb_memory(const char* data, size_t size,
                                               const std::string& name,
                                               Structure& st,
                                               PdbReadOptions options={}) {
  MemoryStream stream{data, size};
  populate_structure_from_pdb_stream(stream, name, st, options);
}

//! @brief Populate Structure from PDB string.
//! @param str PDB file content as string
//! @param name Source name (for error messages)
//! @param st Structure to populate
//! @param options Reading options
inline void populate_structure_from_pdb_string(const std::string& str,
                                              const std::string& name,
                                              Structure& st,
                                              PdbReadOptions options={}) {
  populate_structure_from_pdb_memory(str.c_str(), str.length(), name, st, options);
}

//! @brief Read PDB from stream and return Structure.
//! @param line_reader Input stream providing PDB lines
//! @param source Source name/path (for error messages)
//! @param options Reading options
//! @return Parsed Structure
inline Structure read_pdb_from_stream(AnyStream& line_reader,
                                      const std::string& source,
                                      PdbReadOptions options) {
  gemmi::Structure st;
  populate_structure_from_pdb_stream(line_reader, source, st, options);
  return st;
};

//! @brief Read PDB file and return Structure.
//! @param path Path to PDB file
//! @param options Reading options
//! @return Parsed Structure
inline Structure read_pdb_file(const std::string& path,
                               PdbReadOptions options={}) {
  FileStream stream(path.c_str(), "rb");
  return read_pdb_from_stream(stream, path, options);
}

//! @brief Read PDB from memory buffer and return Structure.
//! @param data Pointer to PDB file content in memory
//! @param size Size of data in bytes
//! @param name Source name (for error messages)
//! @param options Reading options
//! @return Parsed Structure
inline Structure read_pdb_from_memory(const char* data, size_t size,
                                      const std::string& name,
                                      PdbReadOptions options={}) {
  MemoryStream stream{data, size};
  return read_pdb_from_stream(stream, name, options);
}

//! @brief Read PDB from string and return Structure.
//! @param str PDB file content as string
//! @param name Source name (for error messages)
//! @param options Reading options
//! @return Parsed Structure
inline Structure read_pdb_string(const std::string& str,
                                 const std::string& name,
                                 PdbReadOptions options={}) {
  return read_pdb_from_memory(str.c_str(), str.length(), name, options);
}

//! @brief Read PDB from input source (file, stream, etc.) and return Structure.
//! @tparam T Input type (e.g., BasicInput, MaybeGzipped)
//! @param input Input source providing file content
//! @param options Reading options
//! @return Parsed Structure
template<typename T>
inline Structure read_pdb(T&& input, PdbReadOptions options={}) {
  return read_pdb_from_stream(*input.create_stream(), input.path(), options);
}

} // namespace gemmi
#endif
