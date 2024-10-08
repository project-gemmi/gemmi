// Copyright Global Phasing Ltd.

#include <gemmi/json.hpp>
#include <utility>  // for move

#define SAJSON_UNSORTED_OBJECT_KEYS
#define SAJSON_NUMBERS_AS_STRINGS
#include "../third_party/sajson.h"

namespace gemmi {
namespace cif {
using std::size_t;

static std::string json_type_as_string(sajson::type t) {
  switch (t) {
    case sajson::TYPE_INTEGER: return "<integer>";
    case sajson::TYPE_DOUBLE:  return "<double>";
    case sajson::TYPE_NULL:    return "<null>";
    case sajson::TYPE_FALSE:   return "<false>";
    case sajson::TYPE_TRUE:    return "<true>";
    case sajson::TYPE_STRING:  return "<string>";
    case sajson::TYPE_ARRAY:   return "<array>";
    case sajson::TYPE_OBJECT:  return "<object>";
    default:           return "<unknown type>";
  }
}

static std::string as_cif_value(const sajson::value& val) {
  switch (val.get_type()) {
    case sajson::TYPE_DOUBLE:
      return val.as_string();
    case sajson::TYPE_NULL:
      return "?";
    // mmJSON files from PDBj (this format has no spec) have special support
    // for boolean YES|NO, which is used only in category _em_specimen.
    // IMO it's a bad idea, but we must handle it if we want to read mmJSON.
    case sajson::TYPE_FALSE:
      return "NO";  // "." in CIF-JSON
    case sajson::TYPE_TRUE:
      return "YES";
    case sajson::TYPE_STRING:
      return quote(val.as_string());
    // Another undocumented feature of mmJSON: arrays as values.
    // It seems that obscure types int-range and float-range are converted to
    // 2-element arrays. But not only. link_entity_pdbjplus.db_accession has
    // arrays with strings.
    case sajson::TYPE_ARRAY: {
      std::string s;
      for (size_t i = 0; i < val.get_length(); ++i) {
        if (i != 0)
          s += ' ';
        s += val.get_array_element(0).as_string();
      }
      return quote(s);
    }
    default:
      fail("Unexpected ", json_type_as_string(val.get_type()), " as value in JSON.");
      return "";
  }
}

static void fill_document_from_sajson(Document& d, const sajson::document& s) {
  // assuming mmJSON here, we'll add handling of CIF-JSON later on
  sajson::value root = s.get_root();
  if (root.get_type() != sajson::TYPE_OBJECT)
    fail("not mmJSON - the root is not of type object");
  for (size_t block_index = 0; block_index < root.get_length(); ++block_index) {
    std::string block_name = root.get_object_key(block_index).as_string();
    if (!starts_with(block_name, "data_"))
      fail("not mmJSON - top level key should start with data_\n"
           "(if you use gemmi-cif2json to write JSON, use -m for mmJSON)");
    d.blocks.emplace_back(block_name.substr(5));
    std::vector<Item>& items = d.blocks[block_index].items;
    sajson::value top = root.get_object_value(block_index);
    if (top.get_type() != sajson::TYPE_OBJECT)
      fail("");
    for (size_t i = 0; i != top.get_length(); ++i) {
      std::string category_name = "_" + top.get_object_key(i).as_string() + ".";
      sajson::value category = top.get_object_value(i);
      if (category.get_type() != sajson::TYPE_OBJECT ||
          category.get_length() == 0 ||
          category.get_object_value(0).get_type() != sajson::TYPE_ARRAY)
        fail("");
      size_t cif_cols = category.get_length();
      size_t cif_rows = category.get_object_value(0).get_length();
      if (cif_rows > 1) {
        items.emplace_back(LoopArg{});
        Loop& loop = items.back().loop;
        loop.tags.reserve(cif_cols);
        loop.values.resize(cif_cols * cif_rows);
      }
      for (size_t j = 0; j != cif_cols; ++j) {
        std::string tag = category_name + category.get_object_key(j).as_string();
        sajson::value arr = category.get_object_value(j);
        if (arr.get_type() != sajson::TYPE_ARRAY)
          fail("Expected array, got " + json_type_as_string(arr.get_type()));
        if (arr.get_length() != cif_rows)
          fail("Expected array of length ", std::to_string(cif_rows), " not ",
               std::to_string(arr.get_length()));
        if (cif_rows == 1) {
          items.emplace_back(tag, as_cif_value(arr.get_array_element(0)));
        } else if (cif_rows != 0) {
          Loop& loop = items.back().loop;
          loop.tags.emplace_back(std::move(tag));
          for (size_t k = 0; k != cif_rows; ++k)
            loop.values[j + k*cif_cols] = as_cif_value(arr.get_array_element(k));
        }
      }
    }
  }
}

Document read_mmjson_insitu(char* buffer, size_t size, const std::string& name) {
  Document doc;
  sajson::document json = sajson::parse(sajson::dynamic_allocation(),
                                    sajson::mutable_string_view(size, buffer));
  if (!json.is_valid())
    fail(name + ":", std::to_string(json.get_error_line()), " error: ",
         json.get_error_message_as_string());
  fill_document_from_sajson(doc, json);
  doc.source = name;
  return doc;
}



} // namespace cif
} // namespace gemmi
