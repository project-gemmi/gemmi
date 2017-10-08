// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.
#pragma once
#include "gemmi/cifdoc.hpp" // for Document
#include "gemmi/model.hpp"  // for Structure
#include "gemmi/util.hpp"   // for to_lower

gemmi::cif::Document cif_read_any(const std::string& path);

gemmi::Structure mmcif_read_atoms(const gemmi::cif::Document& doc);

gemmi::Structure pdb_read_any(const std::string& path);

inline bool is_pdb_code(const std::string& str) {
  return str.length() == 4 && std::isdigit(str[0]) && std::isalnum(str[1]) &&
                              std::isalnum(str[2]) && std::isalnum(str[3]);
}

inline std::string mmcif_subpath(const std::string& code) {
  std::string lc = gemmi::to_lower(code);
  return "/structures/divided/mmCIF/" + lc.substr(1, 2) + "/" + lc + ".cif.gz";
}


// vim:sw=2:ts=2:et:path^=../include,../third_party
