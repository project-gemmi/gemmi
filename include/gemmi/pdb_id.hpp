// Copyright 2018 Global Phasing Ltd.
//
// handling PDB ID and $PDB_DIR: is_pdb_code(), expand_pdb_code_to_path()

#ifndef GEMMI_PDB_ID_HPP_
#define GEMMI_PDB_ID_HPP_

#include <cctype>    // for isdigit, isalnum
#include <cstdlib>   // getenv
#include <string>
#include "fail.hpp"  // for fail
#include "util.hpp"  // for to_lower

namespace gemmi {

inline bool is_pdb_code(const std::string& str) {
  return str.length() == 4 && std::isdigit(str[0]) && std::isalnum(str[1]) &&
                              std::isalnum(str[2]) && std::isalnum(str[3]);
}

/// Call it after checking the code with gemmi::is_pdb_code(code).
/// The convention for $PDB_DIR is the same as in BioJava, see the docs.
/// \par type is the requested file type: 'M' for mmCIF or 'P' for PDB, 'S' for SF-mmCIF.
inline std::string expand_pdb_code_to_path(const std::string& code, char type) {
  std::string path;
  if (const char* pdb_dir = std::getenv("PDB_DIR")) {
    int n = 0;
    if (type == 'M')
      n = 1;
    else if (type == 'S')
      n = 2;
    std::string lc = to_lower(code);
    path = pdb_dir;
    path += "/structures/divided/";
    const char* dir[] = {"pdb/", "mmCIF/", "structure_factors/"};
    path += dir[n];
    path += lc.substr(1, 2);
    const char* prefix[] = {"/pdb", "/", "/r"};
    path += prefix[n];
    path += lc;
    const char* suffix[] = {".ent.gz", ".cif.gz", "sf.ent.gz"};
    path += suffix[n];
  }
  return path;
}

/// \par type is: 'M' for mmCIF or 'P' for PDB, 'S' for SF-mmCIF.
inline std::string expand_if_pdb_code(const std::string& input, char type='M') {
  std::string path;
  if (is_pdb_code(input)) {
    path = expand_pdb_code_to_path(input, type);
    if (path.empty())
      fail(input + " is a PDB code, but $PDB_DIR is not set.");
  } else {
    path = input;
  }
  return path;
}

} // namespace gemmi
#endif
