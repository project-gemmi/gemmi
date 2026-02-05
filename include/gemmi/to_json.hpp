//! @file
//! @brief Write CIF documents as JSON (mmJSON, CIF-JSON).

// Copyright 2017 Global Phasing Ltd.

// Writing cif::Document or its parts as JSON (mmJSON, CIF-JSON, etc).

#ifndef GEMMI_TO_JSON_HPP_
#define GEMMI_TO_JSON_HPP_
#include <ostream>   // for ostream
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

//! @brief Options for JSON output formatting.
//!
//! Supports multiple JSON variants: COMCIFS CIF-JSON, mmJSON, etc.
struct JsonWriteOptions {
  bool as_comcifs = false;              //!< Conform to COMCIFS CIF-JSON draft
  bool group_ddl2_categories = false;   //!< Group DDL2 categories (mmJSON)
  bool with_data_keyword = false;       //!< Include "data" keyword (mmJSON)
  bool bare_tags = false;               //!< "tag" instead of "_tag"
  bool values_as_arrays = false;        //!< "_tag": ["value"] format
  bool lowercase_names = true;          //!< Write case-insensitive names as lowercase
  int quote_numbers = 1;                //!< 0=never, 1=mixed, 2=always quote numbers
  std::string cif_dot = "null";         //!< How to convert CIF '.' values

  //! @brief Create options for COMCIFS CIF-JSON format.
  //! @return Options configured for COMCIFS CIF-JSON specification
  static JsonWriteOptions comcifs() {
    JsonWriteOptions opt;
    opt.as_comcifs = true;
    opt.values_as_arrays = true;
    opt.quote_numbers = 2;
    opt.cif_dot = "false";
    return opt;
  }

  //! @brief Create options for mmJSON format.
  //! @return Options configured for mmJSON specification
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

//! @brief Write CIF document to stream as JSON.
//! @param os Output stream
//! @param doc CIF document to write
//! @param options JSON formatting options
GEMMI_DLL void write_json_to_stream(std::ostream& os, const Document& doc,
                                    const JsonWriteOptions& options);

//! @brief Write CIF document to stream as mmJSON.
//! @param os Output stream
//! @param doc CIF document to write
inline void write_mmjson_to_stream(std::ostream& os, const Document& doc) {
  write_json_to_stream(os, doc, JsonWriteOptions::mmjson());
}

} // namespace cif
} // namespace gemmi
#endif
