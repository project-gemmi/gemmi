// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.

#include "gemmi/gz.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/pdbgz.hpp"
#include "gemmi/json.hpp"
#include "gemmi/util.hpp"

gemmi::cif::Document cif_read_any(const std::string& path) {
  if (gemmi::iends_with(path, "json") || gemmi::iends_with(path, "js") ||
      gemmi::iends_with(path, "json.gz") || gemmi::iends_with(path, "js.gz"))
    return gemmi::cif::read_mmjson(gemmi::MaybeGzipped(path));
  return gemmi::cif::read(gemmi::MaybeGzipped(path));
}

gemmi::mol::Structure mmcif_read_atoms(const gemmi::cif::Document& doc) {
  return gemmi::mol::read_atoms(doc);
}

gemmi::mol::Structure pdb_read_any(const std::string& path) {
  return gemmi::mol::read_pdb_any(path);
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
