// Copyright 2017 Global Phasing Ltd.

// Writing cif::Document or its parts as JSON (mmJSON, CIF-JSON, etc).

#ifndef GEMMI_TO_JSON_HPP_
#define GEMMI_TO_JSON_HPP_
#include <ostream>   // for ostream
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

struct JsonWriteOptions {
  bool as_comcifs = false;  // conform to the COMCIFS CIF-JSON draft
  bool group_ddl2_categories = false;  // for mmJSON
  bool with_data_keyword = false;  // for mmJSON
  bool bare_tags = false;  // "tag" instead of "_tag"
  bool values_as_arrays = false;  // "_tag": ["value"]
  bool lowercase_names = true; // write case-insensitive names as lower case
  int quote_numbers = 1;  // 0=never (no s.u.), 1=mix, 2=always
  std::string cif_dot = "null";  // how to convert '.' from CIF

  static JsonWriteOptions comcifs() {
    JsonWriteOptions opt;
    opt.as_comcifs = true;
    opt.values_as_arrays = true;
    opt.quote_numbers = 2;
    opt.cif_dot = "false";
    return opt;
  }

  static JsonWriteOptions mmjson() {
    JsonWriteOptions opt;
    opt.group_ddl2_categories = true;
    opt.with_data_keyword = true;
    opt.bare_tags = true;
    opt.values_as_arrays = true;
    opt.lowercase_names = false;
    opt.quote_numbers = 0;
    return opt;
  }
};

GEMMI_DLL void write_json_to_stream(std::ostream& os, const Document& doc,
                                    const JsonWriteOptions& options);

inline void write_mmjson_to_stream(std::ostream& os, const Document& doc) {
  write_json_to_stream(os, doc, JsonWriteOptions::mmjson());
}

} // namespace cif
} // namespace gemmi
#endif
