// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.

#include "gemmi/cifgz.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/pdbgz.hpp"

gemmi::cif::Document cif_read_any(const std::string& path) {
  return gemmi::cif::read_any(path);
}

gemmi::mol::Structure mmcif_read_atoms(const gemmi::cif::Document& doc) {
  return gemmi::mol::read_atoms(doc);
}

gemmi::mol::Structure pdb_read_any(const std::string& path) {
  return gemmi::mol::read_pdb_any(path);
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
