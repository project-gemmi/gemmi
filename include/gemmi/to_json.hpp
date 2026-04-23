// Copyright 2017 Global Phasing Ltd.

/// @file
/// @brief Writing CIF documents as JSON (mmJSON and CIF-JSON formats).

#ifndef GEMMI_TO_JSON_HPP_
#define GEMMI_TO_JSON_HPP_
#include <ostream>   // for ostream
#include "cifdoc.hpp"

namespace gemmi {
namespace cif {

/// Options for writing CIF data as JSON.
///
/// Supports multiple JSON-based serialization formats for CIF data:
/// - **CIF-JSON (COMCIFS)**: Standard JSON representation of CIF documents,
///   supporting numbered values with uncertainties.
/// - **mmJSON (PDBj)**: Specialized JSON format optimized for macromolecular
///   CIF (mmCIF) data, with DDL2 category grouping and bare tags.
///
/// Choose between preset configurations (comcifs() or mmjson()) or
/// configure individual options for custom output.
struct JsonWriteOptions {
  /// Conform to the COMCIFS CIF-JSON draft specification.
  /// If true, enables values_as_arrays, sets quote_numbers=2, and cif_dot="false".
  bool as_comcifs = false;

  /// Group items by DDL2 categories (for mmJSON compatibility).
  /// Relevant mainly for mmJSON format.
  bool group_ddl2_categories = false;

  /// Include the mmJSON "data_" keyword wrapper.
  /// Used in mmJSON format output.
  bool with_data_keyword = false;

  /// Use bare tag names (e.g., "tag" instead of "_tag").
  /// Used in mmJSON and other compact formats.
  bool bare_tags = false;

  /// Represent all values as JSON arrays (e.g., "_tag": ["value"]).
  /// Used in COMCIFS CIF-JSON and mmJSON; disabled if false.
  bool values_as_arrays = false;

  /// Write case-insensitive tag names in lowercase.
  /// CIF tag names are case-insensitive; this normalizes them.
  bool lowercase_names = true;

  /// Control quoting of numeric values with uncertainty (s.u.).
  /// - 0: Never quote numbers; s.u. information is lost (used for mmJSON)
  /// - 1: Quote numbers only when they include s.u. (default, mixed mode)
  /// - 2: Always quote numbers as strings (used for COMCIFS)
  int quote_numbers = 1;

  /// How to represent the CIF '.' (not-applicable) value in JSON.
  /// Common choices: "null" (JSON null), "false" (boolean false, used in COMCIFS).
  std::string cif_dot = "null";

  /// Preset options for COMCIFS CIF-JSON format.
  ///
  /// @return JsonWriteOptions configured for standard CIF-JSON output
  static JsonWriteOptions comcifs() {
    JsonWriteOptions opt;
    opt.as_comcifs = true;
    opt.values_as_arrays = true;
    opt.quote_numbers = 2;
    opt.cif_dot = "false";
    return opt;
  }

  /// Preset options for mmJSON format (PDBj macromolecular JSON).
  ///
  /// mmJSON is used by PDBj for macromolecular structures.
  /// It groups data by DDL2 categories and uses bare tag names.
  ///
  /// @return JsonWriteOptions configured for mmJSON output
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

/// Write a CIF document as JSON to an output stream.
///
/// Serializes a CIF document in JSON format according to the specified options.
/// See JsonWriteOptions for details on supported formats and customization.
///
/// @param os      Output stream to write to
/// @param doc     The CIF document to write
/// @param options Formatting and format selection options
GEMMI_DLL void write_json_to_stream(std::ostream& os, const Document& doc,
                                    const JsonWriteOptions& options);

/// Write a CIF document as mmJSON (PDBj macromolecular JSON) to an output stream.
///
/// Convenience function equivalent to:
/// @code
/// write_json_to_stream(os, doc, JsonWriteOptions::mmjson());
/// @endcode
///
/// @param os  Output stream to write to
/// @param doc The CIF document to write
inline void write_mmjson_to_stream(std::ostream& os, const Document& doc) {
  write_json_to_stream(os, doc, JsonWriteOptions::mmjson());
}

} // namespace cif
} // namespace gemmi
#endif
