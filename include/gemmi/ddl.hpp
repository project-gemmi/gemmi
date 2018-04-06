// Copyright 2017 Global Phasing Ltd.
//
// Class DDL that represents DDL1 or DDL2 dictionary (ontology).
// Used to validate CIF files.

#ifndef GEMMI_DDL_HPP_
#define GEMMI_DDL_HPP_

#include <algorithm>  // for find
#include <cmath>      // for INFINITY
#include <string>
#include <unordered_map>
#include <utility>    // for pair
#include <vector>
#include "cifdoc.hpp"
#include "numb.hpp"
#include "cif.hpp"    // for read_file

namespace gemmi {
namespace cif {

class DDL {
public:
  void open_file(const std::string& filename) {
    ddl_ = read_file(filename);
    if (ddl_.blocks.size() > 1) {
      version_ = 1;
      sep_ = "_";
      read_ddl1();
    } else {
      version_ = 2;
      sep_ = ".";
      read_ddl2();
    }
  }
  // does the dictionary name/version correspond to _audit_conform_dict_*
  bool check_audit_conform(const Document& doc, std::string* msg) const;
  template <class Output>
  bool validate(Document& doc, Output& out, bool quiet);

private:
  Block* find_rules(const std::string& name) {
    auto iter = name_index_.find(name);
    return iter != name_index_.end() ? iter->second : nullptr;
  }

  void read_ddl1() {
    for (Block& b : ddl_.blocks) {
      for (std::string& name : b.find_values("_name"))
        name_index_.emplace(as_string(name), &b);
      if (b.name == "on_this_dictionary") {
        const std::string* dic_name = b.find_value("_dictionary_name");
        if (dic_name)
          dict_name_ = as_string(*dic_name);
        const std::string* dic_ver = b.find_value("_dictionary_version");
        if (dic_ver)
          dict_version_ = as_string(*dic_ver);
      }
    }
  }

  void read_ddl2() {
    for (Block& block : ddl_.blocks) // a single block is expected
      for (Item& item : block.items) {
        if (item.type == ItemType::Frame) {
          for (const std::string& name : item.frame.find_values("_item.name"))
            name_index_.emplace(as_string(name), &item.frame);
        } else if (item.type == ItemType::Pair) {
          if (item.pair[0] == "_dictionary.title")
            dict_name_ = item.pair[1];
          else if (item.pair[0] == "_dictionary.version")
            dict_version_ = item.pair[1];
        }
      }
  }

  template <class Output, class TypeCheckDDL>
  bool do_validate(Document& doc, Output& out, bool quiet);

  int version_;
  Document ddl_;
  std::unordered_map<std::string, Block*> name_index_;
  std::string dict_name_;
  std::string dict_version_;
  // "_" or ".", used to unify handling of DDL1 and DDL2, for example when
  // reading _audit_conform_dict_version and _audit_conform.dict_version.
  std::string sep_;
};



inline
bool DDL::check_audit_conform(const Document& doc, std::string* msg) const {
  std::string audit_conform = "_audit_conform" + sep_;
  for (const Block& b : doc.blocks) {
    const std::string* dict_name = b.find_value(audit_conform + "dict_name");
    if (!dict_name)
      continue;
    std::string name = as_string(*dict_name);
    if (name != dict_name_) {
      if (msg)
          *msg = "Dictionary name mismatch: " + name + " vs " + dict_name_;
      return false;
    }
    const std::string* dict_ver = b.find_value(audit_conform + "dict_version");
    if (dict_ver) {
      std::string version = as_string(*dict_ver);
      if (version != dict_version_) {
        if (msg)
          *msg = "CIF conforms to " + name + " ver. " + version
                 + " while DDL has ver. " + dict_version_;
        return false;
      }
    }
  }
  if (msg)
    *msg = "The cif file is missing " + audit_conform + "dict_(name|version)";
  return true;
}

enum class Trinary : char { Unset, Yes, No };
enum class ValType : char { Unset, Numb, Any };

class TypeCheckCommon {
protected:
  ValType type_ = ValType::Unset;
  bool has_su_ = false; // _type_conditions esd|su
  bool range_inclusive_ = false;
  std::vector<std::pair<double, double>> range_;
  std::vector<std::string> enumeration_;

public:
  bool validate_value(const std::string& value, std::string* msg) const {
    if (is_null(value))
      return true;
    if (type_ == ValType::Numb && !is_numb(value)) {
      if (msg)
        *msg = "expected number, got: " + value;
      return false;
    }
    // ignoring has_su_ - not sure if we should check it
    if (!range_.empty() && !validate_range(value, msg))
      return false;
    if (!enumeration_.empty() && !validate_enumeration(value, msg))
      return false;
    return true;
  }

private:
  bool validate_range(const std::string& value, std::string *msg) const {
    const double x = as_number(value);
    for (const auto& r : range_)
      if (r.first == r.second ? x == r.first
                              : r.first < x && x < r.second)
        return true;
    if (msg)
      *msg = "value out of expected range: " + value;
    return false;
  }

  bool validate_enumeration(const std::string& val, std::string *msg) const {
    if (std::find(enumeration_.begin(), enumeration_.end(), as_string(val))
          != enumeration_.end())
      return true;
    // TODO: case-insensitive search when appropriate
    if (msg) {
      *msg = "'" + val + "' is not one of:";
      for (const std::string& e : enumeration_)
        *msg += " " + e + ",";
      (*msg)[msg->size() - 1] = '.';
    }
    return false;
  }
};

class TypeCheckDDL1 : public TypeCheckCommon {
public:
  void from_block(Block& b) {
    if (const std::string* list = b.find_value("_list")) {
      if (*list == "yes")
        is_list_ = Trinary::Yes;
      else if (*list == "no")
        is_list_ = Trinary::No;
    }
    const std::string* type = b.find_value("_type");
    if (type)
      type_ = (*type == "numb" ? ValType::Numb : ValType::Any);
    // Hypotetically _type_conditions could be a list, but it never is.
    const std::string* conditions = b.find_value("_type_conditions");
    if (conditions)
      has_su_ = (*conditions == "esd" || *conditions == "su");
    const std::string* range = b.find_value("_enumeration_range");
    if (range) {
      size_t colon_pos = range->find(':');
      if (colon_pos != std::string::npos) {
        std::string low = range->substr(0, colon_pos);
        std::string high = range->substr(colon_pos+1);
        range_inclusive_ = true;
        range_.emplace_back(low.empty() ? -INFINITY : as_number(low),
                            high.empty() ? INFINITY : as_number(high));
      }
    }
    for (const std::string& e : b.find_loop("_enumeration"))
      enumeration_.emplace_back(as_string(e));
  }

  Trinary is_list() const { return is_list_; }

private:
  Trinary is_list_ = Trinary::Unset; // _list yes
  // type_construct regex - it is rarely used, ignore for now
  // type_conditions seq - seems to be never used, ignore it
  // For now we don't check at all relational attributes, i.e.
  // _category, _list_* and _related_*
};


class TypeCheckDDL2 : public TypeCheckCommon {
public:
  void from_block(Block& b) {
    if (const std::string* code = b.find_value("_item_type.code")) {
      type_code_ = as_string(*code);
      if (type_code_ == "float" || type_code_ == "int")
        type_ = ValType::Numb;
    }
    for (auto row : b.find("_item_range.", {"minimum", "maximum"}))
      range_.emplace_back(as_number(row[0], -INFINITY),
                          as_number(row[1], +INFINITY));
    for (const std::string& e : b.find_loop("_item_enumeration.value"))
      enumeration_.emplace_back(as_string(e));
  }

  Trinary is_list() const { return Trinary::Unset; }

private:
  std::string type_code_;
};

template <class Output>
bool DDL::validate(Document& doc, Output& out, bool quiet) {
  return version_ == 1 ? do_validate<Output, TypeCheckDDL1>(doc, out, quiet)
                       : do_validate<Output, TypeCheckDDL2>(doc, out, quiet);
}

template <class Output, class TypeCheckDDL>
bool DDL::do_validate(Document& doc, Output& out, bool quiet) {
  std::string msg;
  bool ok = true;
  auto err = [&](const Block& b, const Item& item, const std::string& s) {
    ok = false;
    out << doc.source << ":" << item.line_number
        << " in data_" << b.name << ": " << s << "\n";
  };
  for (Block& b : doc.blocks) {
    for (const Item& item : b.items) {
      if (item.type == ItemType::Pair) {
        Block* dict_block = find_rules(item.pair[0]);
        if (!dict_block) {
          if (!quiet)
            out << "Note: unknown tag: " << item.pair[0] << "\n";
          continue;
        }
        TypeCheckDDL tc;
        tc.from_block(*dict_block);
        if (tc.is_list() == Trinary::Yes)
          err(b, item, item.pair[0] + " must be a list");
        if (!tc.validate_value(item.pair[1], &msg))
          err(b, item, msg);
      } else if (item.type == ItemType::Loop) {
        const size_t ncol = item.loop.tags.size();
        for (size_t i = 0; i != ncol; i++) {
          const std::string& tag = item.loop.tags[i];
          Block* dict_block = find_rules(tag);
          if (!dict_block) {
            if (!quiet)
              out << "Note: unknown tag: " << tag << "\n";
            continue;
          }
          TypeCheckDDL tc;
          tc.from_block(*dict_block);
          if (tc.is_list() == Trinary::No)
            err(b, item, tag + " in list");
          for (size_t j = i; j < item.loop.values.size(); j += ncol)
            if (!tc.validate_value(item.loop.values[j], &msg)) {
              err(b, item, tag + ": " + msg);
              break; // stop after first error to avoid clutter
            }
        }
      }
    }
  }
  return ok;
}


} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
