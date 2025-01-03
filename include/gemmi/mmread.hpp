// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file. Usually, mmread_gz.hpp is preferred.

#ifndef GEMMI_MMREAD_HPP_
#define GEMMI_MMREAD_HPP_

#include "cif.hpp"       // for cif::read
#include "fail.hpp"      // for fail
#include "input.hpp"     // for BasicInput
#include "json.hpp"      // for read_mmjson
#include "mmcif.hpp"     // for make_structure_from_block, ...
#include "model.hpp"     // for Structure
#include "pdb.hpp"       // for read_pdb
#include "util.hpp"      // for iends_with

namespace gemmi {

inline CoorFormat coor_format_from_ext(const std::string& path) {
  if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
    return CoorFormat::Pdb;
  if (iends_with(path, ".cif") || iends_with(path, ".mmcif"))
    return CoorFormat::Mmcif;
  if (iends_with(path, ".json"))
    return CoorFormat::Mmjson;
  return CoorFormat::Unknown;
}

// If it's neither CIF nor JSON nor almost empty - we assume PDB.
inline CoorFormat coor_format_from_content(const char* buf, const char* end) {
  while (buf < end - 8) {
    if (std::isspace(*buf)) {
      ++buf;
    } else if (*buf == '#') {
      while (buf < end - 8 && *buf != '\n')
        ++buf;
    } else if (*buf == '{') {
      return CoorFormat::Mmjson;
    } else if (ialpha4_id(buf) == ialpha4_id("data") && buf[4] == '_') {
      return CoorFormat::Mmcif;
    } else {
      return CoorFormat::Pdb;
    }
  }
  return CoorFormat::Unknown;
}

inline Structure make_structure_from_doc(cif::Document&& doc, bool possible_chemcomp,
                                         cif::Document* save_doc=nullptr) {
  if (possible_chemcomp) {
    // check for special case - refmac dictionary or CCD file
    int n = check_chemcomp_block_number(doc);
    if (n != -1)
      return make_structure_from_chemcomp_block(doc.blocks[n]);
  }
  return make_structure(std::move(doc), save_doc);
}

// when reading JSON, the input buffer is changed (as an optimization)
inline Structure read_structure_from_memory(char* data, size_t size,
                                            const std::string& path,
                                            CoorFormat format=CoorFormat::Unknown,
                                            cif::Document* save_doc=nullptr) {
  if (save_doc)
    save_doc->clear();
  if (format == CoorFormat::Unknown || format == CoorFormat::Detect)
    format = coor_format_from_content(data, data + size);
  if (format == CoorFormat::Pdb)
    return read_pdb_from_memory(data, size, path);
  if (format == CoorFormat::Mmcif)
    return make_structure_from_doc(cif::read_memory(data, size, path.c_str()),
                                   true, save_doc);
  if (format == CoorFormat::Mmjson)
    return make_structure(cif::read_mmjson_insitu(data, size, path), save_doc);
  fail("wrong format of coordinate file " + path);
}

// deprecated
inline Structure read_structure_from_char_array(char* data, size_t size,
                                                const std::string& path,
                                                cif::Document* save_doc=nullptr) {
  return read_structure_from_memory(data, size, path, CoorFormat::Unknown, save_doc);
}

template<typename T>
Structure read_structure(T&& input, CoorFormat format=CoorFormat::Unknown,
                         cif::Document* save_doc=nullptr) {
  if (format == CoorFormat::Detect) {
    CharArray mem = read_into_buffer(input);
    return read_structure_from_memory(mem.data(), mem.size(), input.path(), format, save_doc);
  }
  if (save_doc)
    save_doc->clear();
  if (format == CoorFormat::Unknown)
    format = coor_format_from_ext(input.basepath());
  switch (format) {
    case CoorFormat::Pdb:
      return read_pdb(input);
    case CoorFormat::Mmcif:
      return make_structure(cif::read(input), save_doc);
    case CoorFormat::Mmjson: {
      Structure st = make_structure(cif::read_mmjson(input), save_doc);
      st.input_format = CoorFormat::Mmjson;
      return st;
    }
    case CoorFormat::ChemComp:
      return make_structure_from_chemcomp_doc(cif::read(input), save_doc);
    case CoorFormat::Unknown:
    case CoorFormat::Detect:
      fail("Unknown format of " +
           (input.path().empty() ? "coordinate file" : input.path()) + ".");
  }
  unreachable();
}

inline Structure read_structure_file(const std::string& path,
                                     CoorFormat format=CoorFormat::Unknown) {
  return read_structure(BasicInput(path), format);
}

} // namespace gemmi
#endif
