// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.

#include "input.h"
#include "gemmi/gz.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/json.hpp"
#include "gemmi/mmread.hpp"
#include "gemmi/util.hpp"    // for giends_with

gemmi::cif::Document cif_read_any(const std::string& path) {
  if (gemmi::giends_with(path, "json") || gemmi::giends_with(path, "js"))
    return gemmi::cif::read_mmjson(gemmi::MaybeGzipped(path));
  return gemmi::cif::read(gemmi::MaybeGzipped(path));
}

gemmi::Structure mmcif_read_atoms(const gemmi::cif::Document& doc) {
  return gemmi::read_atoms(doc);
}

gemmi::Structure read_structure(const std::string& path,
                                gemmi::CoorFormat format) {
  return gemmi::read_structure(gemmi::MaybeGzipped(path), format);
}

gemmi::CoorFormat coordinate_format_from_extension(const std::string& path) {
  return gemmi::coordinate_format_from_extension(path);
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
