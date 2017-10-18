// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.
#pragma once
#include "gemmi/cifdoc.hpp" // for Document
#include "gemmi/model.hpp"  // for Structure
#include "gemmi/util.hpp"   // for to_lower
#include "gemmi/gz.hpp"     // for expand_pdb_code_to_path
#include <cstdlib>          // for exit
#include <cstdio>           // for fprintf

gemmi::cif::Document cif_read_any(const std::string& path);

gemmi::Structure mmcif_read_atoms(const gemmi::cif::Document& doc);

gemmi::Structure read_structure(const std::string& path,
                    gemmi::CoorFormat format=gemmi::CoorFormat::Unknown);

gemmi::CoorFormat coordinate_format_from_extension(const std::string& path);

inline std::string expand_pdb_code_to_path_or_fail(const std::string& code) {
  std::string path = gemmi::expand_pdb_code_to_path(code);
  if (path.empty()) {
    std::fprintf(stderr,
        "The argument %s is a PDB code, but $PDB_DIR is not set.\n"
        "(To use a file or directory with such a name use: ./%s)\n",
        path.c_str(), path.c_str());
    std::exit(2);
  }
  return path;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
