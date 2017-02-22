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

class DDL1 {
public:
  void open_file(const std::string& filename);
  // does the dictionary name/version correspond to _audit_conform_dict_*
  bool check_audit_conform(const cif::Document& c, std::string* msg) const;
  void validate(const cif::Document& c,
                std::vector<std::string>* unknown_tags=NULL) const;

private:
  const cif::Block* find(const std::string& name) const {
    auto iter = name_index_.find(name);
    return iter != name_index_.end() ? iter->second : nullptr;
  }

  cif::Document ddl_;
  std::unordered_map<std::string, const cif::Block*> name_index_;
  std::string dict_name_;
  std::string dict_version_;
};

inline void DDL1::open_file(const std::string& filename) {
  ddl_.parse_file(filename);
  for (const cif::Block& b : ddl_.blocks) {
    const std::string* name = b.by_tag("_name");
    if (name) {
      name_index_.emplace(cif::as_string(*name), &b);
    } else {
      int idx;
      const cif::Loop* loop = b.by_loop_tag("_name", &idx);
      if (loop)
        for (size_t i = idx; i < loop->values.size(); i += loop->tags.size())
          name_index_.emplace(cif::as_string(loop->values[i]), &b);
    }
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

inline
bool DDL1::check_audit_conform(const cif::Document& c, std::string* msg) const {
  assert(msg != nullptr);
  for (const cif::Block& b : c.blocks) {
    const std::string* dict_name = b.by_tag("_audit_conform_dict_name");
    if (!dict_name)
      continue;
    std::string name = cif::as_string(*dict_name);
    if (name != dict_name_) {
      *msg = "Dictionary name mismatch: " + name + " vs " + dict_name_;
      return false;
    }
    const std::string* dict_version = b.by_tag("_audit_conform_dict_version");
    if (dict_version) {
      std::string version = cif::as_string(*dict_version);
      if (version != dict_version_) {
        *msg = "CIF conforms to " + name + " ver. " + version
               + " while DDL has ver. " + dict_version_;
        return false;
      }
    }
  }
  *msg = "The cif file is missing _audit_conform_dict_(name|version)";
  return true;
}

enum class Trinary : char { Unset, Yes, No };

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
    const std::string* range = b.by_tag("_enumeration_range");
    if (range) {
      has_range_ = true;
      range_low_ = -FLT_MAX; //TODO
      range_high_ = FLT_MAX;
    }
    int idx;
    const cif::Loop* enumeration = b.by_loop_tag("_enumeration", &idx);
    if (enumeration) {
      // TODO
      //enumeration_.emplace_back();
    }
  }

  bool validate_value(const std::string& value, std::string* msg) const {
    auto fail = [=](std::string&& t) { *msg = t; return false; };
    if (is_numb_ == Trinary::Yes) {
      if (!cif::is_numb(value))
        return fail("expected number");
      if (has_range_) {
        float x = 0; //TODO cif::as_numb(value);
        if (x < range_low_ || x > range_high_)
          return fail("value out of expected range: " + value);
      }
    }
    if (!enumeration_.empty()) {
        if (std::find(enumeration_.begin(), enumeration_.end(),
                      cif::as_string(value)) == enumeration_.end())
          return fail("not one of enumeration values: " + value);
    }
    return true;
  }

  Trinary is_list() const { return is_list_; }

private:
  Trinary is_list_ = Trinary::Unset; // _list yes
  Trinary is_numb_ = Trinary::Unset; // _type numb
  bool has_range_ = false; // _enumeration_range
  float range_low_;
  float range_high_;
  std::vector<std::string> enumeration_; // loop_ _enumeration
  // type_construct regex
  // type_conditions none|esd|su|seq
};

inline void DDL1::validate(const cif::Document& c,
                           std::vector<std::string>* unknown_tags) const {
  std::string msg;
  for (const cif::Block& b : c.blocks) {
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Value) {
        //std::cout << item.tv.tag << "\n";
        const cif::Block* dict_block = find(item.tv.tag);
        if (dict_block) {
          TypeCheckDDL1 tc;
          tc.from_block(*dict_block);
          if (tc.is_list() == Trinary::Yes)
            throw_validation_err(c, b, item, item.tv.value + " not in list");
          if (!tc.validate_value(item.tv.value, &msg))
            throw_validation_err(c, b, item, msg);
        } else {
          if (unknown_tags)
            unknown_tags->emplace_back(item.tv.tag);
        }
      } else if (item.type == cif::ItemType::Loop) {
        const int ncol = item.loop.tags.size();
        for (int i = 0; i != ncol; i++) {
          const std::string& tag = item.loop.tags[i].tag;
          const cif::Block* dict_block = find(tag);
          if (dict_block) {
            TypeCheckDDL1 tc;
            tc.from_block(*dict_block);
            if (tc.is_list() == Trinary::No)
              throw_validation_err(c, b, item, tag + " in list");
            for (size_t j = i; j < item.loop.values.size(); j += ncol)
              if (!tc.validate_value(item.loop.values[j], &msg))
                throw_validation_err(c, b, item, msg);
          } else {
            if (unknown_tags)
              unknown_tags->emplace_back(item.loop.tags[i].tag);
          }
        }
      }
    }
  }
}


class DDL2 {
public:
  //void open_file(const std::string& filename);
  // does the dictionary name/version correspond to _audit_conform.dict_*
  //bool check_audit_conform(const cif::Document& c, std::string* msg) const;
  //void validate(const cif::Document& c) const;

private:
  //const cif::Block* search_frame(const std::string& name) const {
  //  auto iter = frame_index_.find(name);
  //  return iter != frame_index_.end() ? iter->second : nullptr;
  //}

  cif::Document ddl_;
  //std::unordered_map<std::string, const cif::Frame*> frame_index_;
  std::string dict_name_;
  std::string dict_version_;
};


} // namespace ddl
#endif
// vim:sw=2:ts=2:et
