// Copyright 2017 Global Phasing Ltd.
#ifndef GEMMI_DDL_HH_
#define GEMMI_DDL_HH_

#include <cassert>
#include <cfloat>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <iostream> // temporary
#include "cif.hh"

namespace ddl {

class DDL {
public:
  void open_file(const std::string& filename) {
    ddl_.parse_file(filename);
    version_ = (ddl_.blocks.size() > 1 ? 1 : 2);
    sep_ = version_ == 1 ? "_" : ".";
    if (version_ == 1)
      read_ddl1();
    else
      read_ddl2();
  }
  // does the dictionary name/version correspond to _audit_conform_dict_*
  bool check_audit_conform(const cif::Document& c, std::string* msg) const;
  void validate(const cif::Document& c,
                std::vector<std::string>* unknown_tags=nullptr) const;

private:
  const cif::Block* find(const std::string& name) const {
    auto iter = name_index_.find(name);
    return iter != name_index_.end() ? iter->second : nullptr;
  }

  void add_to_index(const cif::Block& b, const std::string& name_tag) {
    const std::string* name = b.by_tag(name_tag);
    if (name)
      name_index_.emplace(cif::as_string(*name), &b);
    else
      for (const std::string& lname : b.looped_values(name_tag))
        name_index_.emplace(cif::as_string(lname), &b);
  }

  void read_ddl1() {
    for (const cif::Block& b : ddl_.blocks) {
      add_to_index(b, "_name");
      if (b.name == "on_this_dictionary") {
        const std::string* dic_name = b.by_tag("_dictionary_name");
        if (dic_name)
          dict_name_ = cif::as_string(*dic_name);
        const std::string* dic_ver = b.by_tag("_dictionary_version");
        if (dic_ver)
          dict_version_ = cif::as_string(*dic_ver);
      }
    }
  }

  void read_ddl2() {
    for (const cif::Block& block : ddl_.blocks) // a single block is expected
      for (const cif::Item& item : block.items) {
        if (item.type == cif::ItemType::Frame) {
          add_to_index(item.frame, "_item.name");
        } else if (item.type == cif::ItemType::Value) {
          if (item.tv.tag == "_dictionary.title")
            dict_name_ = item.tv.value;
          else if (item.tv.tag == "_dictionary.version")
            dict_version_ = item.tv.value;
        }
      }
  }

  int version_;
  cif::Document ddl_;
  std::unordered_map<std::string, const cif::Block*> name_index_;
  std::string dict_name_;
  std::string dict_version_;
  // "_" or ".", used to unify handling of DDL1 and DDL2, for example when
  // reading _audit_conform_dict_version and _audit_conform.dict_version.
  std::string sep_;
};



inline
bool DDL::check_audit_conform(const cif::Document& c, std::string* msg) const {
  std::string audit_conform = "_audit_conform" + sep_;
  for (const cif::Block& b : c.blocks) {
    const std::string* dict_name = b.by_tag(audit_conform + "dict_name");
    if (!dict_name)
      continue;
    std::string name = cif::as_string(*dict_name);
    if (name != dict_name_) {
      if (msg)
          *msg = "Dictionary name mismatch: " + name + " vs " + dict_name_;
      return false;
    }
    const std::string* dict_version = b.by_tag(audit_conform + "dict_version");
    if (dict_version) {
      std::string version = cif::as_string(*dict_version);
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

inline bool validate_enumeration(const std::string& val,
                                 const std::vector<std::string>& en,
                                 std::string *msg) {
  if (en.empty() || val == "." || val == "?" ||
      std::find(en.begin(), en.end(), cif::as_string(val)) != en.end())
    return true;
  // TODO: case-insensitive search when appropriate
  if (msg) {
    *msg = "'" + val + "' is not one of:";
    for (const std::string& e : en)
      *msg += " " + e + ",";
    (*msg)[msg->size() - 1] = '.';
  }
  return false;
}

class TypeCheckDDL1 {
public:
  void from_block(const cif::Block& b) {
    const std::string* list = b.by_tag("_list");
    if (list) {
      if (*list == "yes")
        is_list_ = Trinary::Yes;
      else if (*list == "no")
        is_list_ = Trinary::No;
    }
    const std::string* type = b.by_tag("_type");
    if (type)
      is_numb_ = (*type == "numb" ? Trinary::Yes : Trinary::No);
    // Hypotetically _type_conditions could be a list, but it never is.
    const std::string* conditions = b.by_tag("_type_conditions");
    if (conditions)
      has_su_ = (*conditions == "esd" || *conditions == "su");
    const std::string* range = b.by_tag("_enumeration_range");
    if (range) {
      has_range_ = true;
      size_t colon_pos = range->find(':');
      if (colon_pos == std::string::npos)
        return;
      std::string low = range->substr(0, colon_pos);
      std::string high = range->substr(colon_pos+1);
      range_low_ = low.empty() ? -INFINITY : cif::as_number(low);
      range_high_ = high.empty() ? INFINITY : cif::as_number(high);
    }
    for (const std::string& e : b.looped_values("_enumeration"))
      enumeration_.emplace_back(cif::as_string(e));
  }

  bool validate_value(const std::string& value, std::string* msg) const {
    auto fail = [msg](std::string&& t) { if (msg) *msg = t; return false; };
    if (is_numb_ == Trinary::Yes) {
      if (value != "." && value != "?" && !cif::is_numb(value))
        return fail("expected number");
      if (has_range_) {
        float x = cif::as_number(value);
        if (x < range_low_ || x > range_high_)
          return fail("value out of expected range: " + value);
      }
      // ignoring has_su_ - not sure if we should check it
    }
    if (!validate_enumeration(value, enumeration_, msg))
      return false;
    return true;
  }

  Trinary is_list() const { return is_list_; }

private:
  Trinary is_list_ = Trinary::Unset; // _list yes
  Trinary is_numb_ = Trinary::Unset; // _type numb
  bool has_su_ = false; // _type_conditions esd|su
  bool has_range_ = false; // _enumeration_range
  float range_low_;
  float range_high_;
  std::vector<std::string> enumeration_; // loop_ _enumeration
  // type_construct regex - it is rarely used, ignore for now
  // type_conditions seq - seems to be never used, ignore it
  // For now we don't check at all relational attributes, i.e.
  // _category, _list_* and _related_*
};


class TypeCheckDDL2 {
public:
  void from_block(const cif::Block& b) {
    for (const std::string& e : b.looped_values("_item_enumeration.value"))
      enumeration_.emplace_back(cif::as_string(e));
  }

  bool validate_value(const std::string& value, std::string* msg) const {
    //auto fail = [msg](std::string&& t) { if (msg) *msg = t; return false; };
    if (!validate_enumeration(value, enumeration_, msg))
      return false;
    return true;
  }

private:
  Trinary is_numb_ = Trinary::Unset; // _type numb
  bool has_range_ = false; // _enumeration_range
  float range_low_;
  float range_high_;
  std::vector<std::string> enumeration_; // loop_ _enumeration
};

inline void DDL::validate(const cif::Document& c,
                          std::vector<std::string>* unknown_tags) const {
  std::string msg;
  for (const cif::Block& b : c.blocks) {
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Value) {
        const cif::Block* dict_block = find(item.tv.tag);
        if (!dict_block) {
          if (unknown_tags)
            unknown_tags->emplace_back(item.tv.tag);
          continue;
        }
        if (version_ == 1) {
          TypeCheckDDL1 tc;
          tc.from_block(*dict_block);
          if (tc.is_list() == Trinary::Yes)
            throw_validation_err(c, b, item, item.tv.tag + " must be a list");
          if (!tc.validate_value(item.tv.value, &msg))
            throw_validation_err(c, b, item, msg);
        } else { // version_ == 2
          TypeCheckDDL2 tc;
          tc.from_block(*dict_block);
          if (!tc.validate_value(item.tv.value, &msg))
            throw_validation_err(c, b, item, msg);
        }
      } else if (item.type == cif::ItemType::Loop) {
        const int ncol = item.loop.tags.size();
        for (int i = 0; i != ncol; i++) {
          const std::string& tag = item.loop.tags[i].tag;
          const cif::Block* dict_block = find(tag);
          if (!dict_block) {
            if (unknown_tags)
              unknown_tags->emplace_back(tag);
            continue;
          }
          if (version_ == 1) {
            TypeCheckDDL1 tc;
            tc.from_block(*dict_block);
            if (tc.is_list() == Trinary::No)
              throw_validation_err(c, b, item, tag + " in list");
            for (size_t j = i; j < item.loop.values.size(); j += ncol)
              if (!tc.validate_value(item.loop.values[j], &msg))
                throw_validation_err(c, b, item, msg);
          } else { // version_ == 2
            TypeCheckDDL2 tc;
            tc.from_block(*dict_block);
            for (size_t j = i; j < item.loop.values.size(); j += ncol)
              if (!tc.validate_value(item.loop.values[j], &msg))
                throw_validation_err(c, b, item, msg);
          }
        }
      }
    }
  }
}


} // namespace ddl
#endif
// vim:sw=2:ts=2:et
