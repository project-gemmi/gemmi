// Copyright 2021 Global Phasing Ltd.

#include <gemmi/read_cif.hpp>
#include <gemmi/cif.hpp>    // for cif::read
#include <gemmi/json.hpp>   // for cif::read_mmjson
#include <gemmi/gz.hpp>     // for MaybeGzipped

namespace gemmi {

cif::Document read_cif_gz(const std::string& path) {
  return cif::read(MaybeGzipped(path));
}

cif::Document read_mmjson_gz(const std::string& path) {
  return cif::read_mmjson(MaybeGzipped(path));
}

CharArray read_into_buffer_gz(const std::string& path) {
  return read_into_buffer(MaybeGzipped(path));
}

cif::Document read_cif_from_buffer(const CharArray& buffer, const char* name) {
  return cif::read_memory(buffer.data(), buffer.size(), name);
}

cif::Document read_first_block_gz(const std::string& path, size_t limit) {
  cif::Document doc;
  doc.source = path;
  cif::read_one_block(doc, MaybeGzipped(path), limit);
  return doc;
}

} // namespace gemmi
