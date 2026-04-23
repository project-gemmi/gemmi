// Copyright 2018 Global Phasing Ltd.
//
// Handling PDB ID and $PDB_DIR: is_pdb_code(), expand_pdb_code_to_path(), ...

#ifndef GEMMI_PDB_ID_HPP_
#define GEMMI_PDB_ID_HPP_

#include <cctype>    // for isdigit, isalnum
#include <cstdlib>   // getenv
#include <string>
#include "fail.hpp"  // for fail
#include "util.hpp"  // for to_lower

namespace gemmi {

/// @brief Check if a string consists entirely of alphanumeric characters.
/// @param p Null-terminated string pointer
/// @return True if all characters are alphanumeric
inline bool all_alnums(const char* p) {
  for (;;++p)
    if (!std::isalnum(*p))
      return *p == '\0';
  unreachable();
}

/// @brief Check if a string is a valid PDB code.
/// @param str Code to validate
/// @return True if str is a valid 4-character or extended PDB code
inline bool is_pdb_code(const std::string& str) {
  return (str.length() == 4 && std::isdigit(str[0]) && all_alnums(&str[1]))
      || (str.length() == 12 && str.compare(0, 4, "pdb_") == 0
                             && std::isdigit(str[4]) && all_alnums(&str[5]));
}

/// @brief Build a PDB_DIR-style path component for a given code and type.
/// @param code PDB code (4 characters only; extended codes not yet supported)
/// @param type File type: 'P' for PDB, 'M' for mmCIF, 'S' for structure factors
/// @return Path component relative to $PDB_DIR (e.g., "/structures/divided/pdb/...")
inline std::string path_in_pdb_dir(const std::string& code, char type) {
  if (code.size() == 12)
    fail("extended PDB codes are not supported yet: " + code);
  std::string lc = to_lower(code);
  int n = 0;
  if (type == 'M')
    n = 1;
  else if (type == 'S')
    n = 2;
  std::string path = "/structures/divided/";
  const char* dir[] = {"pdb/", "mmCIF/", "structure_factors/"};
  path += dir[n];
  path += lc.substr(1, 2);
  const char* prefix[] = {"/pdb", "/", "/r"};
  path += prefix[n];
  path += lc;
  const char* suffix[] = {".ent.gz", ".cif.gz", "sf.ent.gz"};
  path += suffix[n];
  return path;
}

/// @brief Expand a PDB code to a full file path using $PDB_DIR.
/// @details Call this after checking the code with gemmi::is_pdb_code(code).
/// The convention for $PDB_DIR is the same as in BioJava.
/// @param code Valid PDB code (4 or 12 characters)
/// @param type File type: 'P' for PDB, 'M' for mmCIF, 'S' for structure factors
/// @param throw_if_unset If true, fail() if $PDB_DIR is not set; if false, return empty string
/// @return Full file path (empty if $PDB_DIR is not set and throw_if_unset is false)
inline std::string expand_pdb_code_to_path(const std::string& code, char type,
                                           bool throw_if_unset=false) {
  std::string path;
  if (const char* pdb_dir = std::getenv("PDB_DIR")) {
    path = pdb_dir + path_in_pdb_dir(code, type);
  } else if (throw_if_unset) {
    fail(code + " is a PDB code, but $PDB_DIR is not set.");
  }
  return path;
}

/// @brief Expand a PDB code to a path, or return the input unchanged if not a code.
/// @param input Either a PDB code or a file path
/// @param type File type: 'P' for PDB, 'M' for mmCIF, 'S' for structure factors
/// @return Expanded path if input is a PDB code; input itself otherwise
inline std::string expand_if_pdb_code(const std::string& input, char type='M') {
  if (is_pdb_code(input))
    return expand_pdb_code_to_path(input, type, true);
  return input;
}

} // namespace gemmi
#endif
