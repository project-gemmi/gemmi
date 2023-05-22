// Copyright 2017-2020 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/numb.hpp"
#include <cstdio>
#include <cmath>      // for INFINITY
#include <algorithm>  // for find
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>    // for pair
#include <vector>
#include <regex>

#ifdef ANALYZE_RULES
# include <tao/pegtl/analyze.hpp>
#endif

#define GEMMI_PROG validate
#include "options.h"

namespace cif = gemmi::cif;

// defined in validate_mon.cpp
void check_monomer_doc(const cif::Document& doc);

namespace {

enum OptionIndex {
  Quiet=4, Fast, Stat, Context, Ddl, NoRegex, NoMandatory, NoUniqueKeys, Parents, Monomer
};
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None, "Usage: " EXE_NAME " [options] FILE [...]"
                                "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Ddl, 0, "d", "ddl", Arg::Required, "  -d, --ddl=PATH  \tDDL for validation." },
  { Context, 0, "c", "context", Arg::None,
    "  -c, --context  \tCheck _pdbx_{category|item}_context.type." },
  { NoRegex, 0, "", "no-regex", Arg::None,
    "  --no-regex  \tSkip regex checking (when using DDL2)" },
  { NoMandatory, 0, "", "no-mandatory", Arg::None,
    "  --no-mandatory  \tSkip checking if mandatory tags are present." },
  { NoUniqueKeys, 0, "", "no-unique", Arg::None,
    "  --no-unique  \tSkip checking if category keys (DDL2) are unique." },
  { Parents, 0, "p", "", Arg::None,
    "  -p  \tCheck if parent items (DDL2) are present." },
  { Monomer, 0, "m", "monomer", Arg::None,
    "  -m, --monomer  \tExtra checks for Refmac dictionary files." },
  { 0, 0, 0, 0, 0, 0 }
};

// basic types, used for token statistics only
enum class ValueType : unsigned char {
  NotSet,
  Char,
  Numb,
  Dot,
  QuestionMark,
};

inline std::string value_type_to_str(ValueType v) {
  switch (v) {
    case ValueType::NotSet: return "n/a";
    case ValueType::Char: return "char";
    case ValueType::Numb: return "numb";
    case ValueType::Dot: return "'.'";
    case ValueType::QuestionMark: return "'?'";
  }
  return "";
}

// For now the infer_* functions are used only here, not sure where they belong
inline ValueType infer_value_type(const std::string& val) {
  assert(!val.empty());
  if (val == ".")
    return ValueType::Dot;
  if (val == "?")
    return ValueType::QuestionMark;
  if (cif::is_numb(val))
    return ValueType::Numb;
  return ValueType::Char;
}

std::string format_7zd(size_t k) {
  char buf[64];
  snprintf(buf, 63, "%7zu", k);
  return buf;
}

std::string br(const std::string& block_name) {
  return "[" + block_name + "] ";
}

std::string token_stats(const cif::Document& d) {
  size_t nframes = 0, nvals = 0, nloops = 0, nlooptags = 0, nloopvals = 0;
  size_t vals_by_type[5] = {0};
  size_t looptags_by_type[5] = {0};
  for (const cif::Block& block : d.blocks) {
    for (const cif::Item& item : block.items) {
      if (item.type == cif::ItemType::Pair) {
        nvals++;
        ValueType vt = infer_value_type(item.pair[1]);
        vals_by_type[static_cast<int>(vt)]++;
      } else if (item.type == cif::ItemType::Frame) {
        nframes++;
      } else if (item.type == cif::ItemType::Loop) {
        nloops++;
        size_t width = item.loop.width();
        nlooptags += width;
        nloopvals += item.loop.values.size();
        for (size_t i = 0; i != width; ++i) {
          ValueType vt = ValueType::NotSet;
          // TODO: ConstColumn(const::Item*, ...)
          const cif::Column col(const_cast<cif::Item*>(&item), i);
          for (const std::string& v : col) {
            ValueType this_vt = infer_value_type(v);
            if (this_vt != vt) {
              // if we are here: vt != ValueType::Char
              if (vt == ValueType::NotSet || this_vt == ValueType::Numb) {
                vt = this_vt;
              } else if (this_vt == ValueType::Char) {
                vt = this_vt;
                break;
              }
            }
          }
          looptags_by_type[static_cast<int>(vt)]++;
        }
      }
    }
  }
  std::string info;
  gemmi::cat_to(info, format_7zd(d.blocks.size()), " block(s)\n");
  gemmi::cat_to(info, format_7zd(nframes), " frames\n");
  gemmi::cat_to(info, format_7zd(nvals), " non-loop items:");
  for (int i = 1; i != 5; ++i)
    gemmi::cat_to(info, "  ", value_type_to_str(static_cast<ValueType>(i)),
                  ':', vals_by_type[i]);
  gemmi::cat_to(info, '\n', format_7zd(nloops), " loops w/"
                "\n        ", format_7zd(nlooptags), " tags:");
  for (int i = 1; i != 5; ++i)
    gemmi::cat_to(info, "  ", value_type_to_str(static_cast<ValueType>(i)),
                  ':', looptags_by_type[i]);
  gemmi::cat_to(info, "\n        ", format_7zd(nloopvals), " values\n");
  return info;
}

// Empty loop is not a valid CIF syntax, but we parse it to accommodate
// some broken CIF files. Only validation shows an error.
void check_empty_loops(const cif::Block& block) {
  for (const cif::Item& item : block.items) {
    if (item.type == cif::ItemType::Loop) {
      if (item.loop.values.empty() && !item.loop.tags.empty())
        throw std::runtime_error("Empty loop in block " + block.name +
                                 ": " + item.loop.tags[0]);
    } else if (item.type == cif::ItemType::Frame) {
      check_empty_loops(item.frame);
    }
  }
}

enum class Trinary : char { Unset, Yes, No };

bool is_integer(const std::string& s) {
  auto b = s.begin() + (s[0] == '+' || s[0] == '-' ? 1 : 0);
  return b != s.end() && std::all_of(b, s.end(), gemmi::is_digit);
}

// Class DDL that represents DDL1 or DDL2 dictionary (ontology).
class DDL {
public:
  DDL(bool enable_regex)
    : regex_enabled_(enable_regex), use_context_(false) {}

  void set_use_context(bool use_context) { use_context_ = use_context; }

  void open_file(const std::string& filename) {
    ddl_ = cif::read_file(filename);
    if (ddl_.blocks.size() > 1) {
      version_ = 1;
      sep_ = '_';
      read_ddl1();
    } else {
      version_ = 2;
      sep_ = '.';
      read_ddl2();
    }
  }

  // use _pdbx_item_linked_group_list
  void read_parent_links() {
    cif::Table tab = ddl_.blocks.at(0).find("_pdbx_item_linked_group_list.",
                                            {"child_category_id", "link_group_id",
                                             "child_name", "parent_name"});
    std::string prev_group;
    ParentLink* it = nullptr;
    for (cif::Table::Row row : tab) {
      std::string group = row.str(0);
      group += ' ';
      group += row[1];
      // here we assume the table is ordered
      if (!it || group != prev_group) {
        parents_.emplace_back();
        it = &parents_.back();
        prev_group = group;
      }
      it->child_tags.push_back(row.str(2));
      it->parent_tags.push_back(row.str(3));
    }
  }

  // check if the dictionary name/version correspond to _audit_conform_dict_*
  void check_audit_conform(const cif::Document& doc, bool verbose) const;

  bool validate(cif::Document& doc, std::ostream& out, bool quiet);

  void check_mandatory_items(cif::Block& b, std::ostream& out) {
    if (version_ != 2)
      return;
    // make a list of items in each category in the block
    std::map<std::string, std::vector<std::string>> categories;
    auto add_category = [&](const std::string& tag) {
      size_t pos = tag.find('.');
      if (pos != std::string::npos)
        categories[tag.substr(0, pos+1)].emplace_back(tag, pos+1);
    };
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Pair)
        add_category(item.pair[0]);
      else if (item.type == cif::ItemType::Loop)
        for (const std::string& tag : item.loop.tags)
          add_category(tag);
    }
    // go over categories and check if nothing is missing
    for (const auto& cat : categories) {
      std::string cat_name = cat.first.substr(1, cat.first.size()-2);
      cif::Block* cat_block = find_rules(cat_name);
      if (!cat_block) { // should not happen
        out << br(b.name) << "category not in the dictionary: " << cat_name << std::endl;
        continue;
      }
      // check context type
      if (use_context_)
        if (const std::string* ct = cat_block->find_value("_pdbx_category_context.type"))
          out << br(b.name) << "category indicated as "
              << *ct << ": " << cat_name << std::endl;
      // check key items
      for (const std::string& v : cat_block->find_values("_category_key.name")) {
        std::string key = cif::as_string(v);
        assert(gemmi::starts_with(key, cat.first));
        if (!gemmi::in_vector(key.substr(cat.first.size()), cat.second))
          out << br(b.name) << "missing category key: " << key << std::endl;
      }
      // check mandatory items
      for (auto i = name_index_.lower_bound(cat.first);
           i != name_index_.end() && gemmi::starts_with(i->first, cat.first);
           ++i) {
        cif::Table items = i->second->find("_item.", {"name", "mandatory_code"});
        if (items.find_row(i->first).str(1)[0] == 'y')
          if (!gemmi::in_vector(i->first.substr(cat.first.size()), cat.second))
            out << br(b.name) << "missing mandatory tag: " << i->first << std::endl;
      }
    }
  }

  void check_unique_keys_in_loop(const cif::Loop& loop, std::ostream& out,
                                 const std::string& block_name) {
    const std::string& tag1 = loop.tags[0];
    size_t dot_pos = tag1.find('.');
    std::string cat_name = tag1.substr(1, dot_pos - 1);
    if (cif::Block* cat_block = find_rules(cat_name)) {
      std::vector<int> key_positions;
      for (const std::string& v : cat_block->find_values("_category_key.name")) {
        int key_pos = loop.find_tag(cif::as_string(v));
        if (key_pos < 0)  // missing keys are reported elsewhere
          return;
        key_positions.push_back(key_pos);
      }
      std::unordered_set<std::string> seen_keys;
      size_t dup_row = 0;
      int dup_counter = 0;
      for (size_t row_pos = 0; row_pos < loop.values.size(); row_pos += loop.width()) {
        std::string row_key;
        for (int k : key_positions) {
          row_key += cif::as_string(loop.values[row_pos + k]);
          row_key += '\1';
        }
        auto r = seen_keys.insert(row_key);
        if (!r.second) {
          dup_counter++;
          if (dup_row == 0)
            dup_row = row_pos;
        }
      }
      if (dup_counter != 0) {
        out << br(block_name) << "category " << cat_name << " has ";
        if (dup_counter == 1)
          out << "1 duplicated key: ";
        else
          out << dup_counter << " duplicated keys, example: ";
        bool first = true;
        for (int k : key_positions) {
          if (first)
            first = false;
          else
            out << " and ";
          out << loop.tags[k].substr(dot_pos+1) << '=' << loop.values[dup_row + k];
        }
        out << std::endl;
      }
    }
  }

  void check_unique_keys(const cif::Block& b, std::ostream& out) {
    if (version_ != 2)
      return;
    for (const cif::Item& item : b.items)
      if (item.type == cif::ItemType::Loop)
        check_unique_keys_in_loop(item.loop, out, b.name);
  }

  void check_parents(cif::Block& b, std::ostream& out) {
    if (version_ != 2)
      return;
    std::unordered_set<std::string> present;
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Pair) {
        if (!cif::is_null(item.pair[1]))
          present.insert(item.pair[0]);
      } else if (item.type == cif::ItemType::Loop) {
        for (const std::string& tag : item.loop.tags)
          present.insert(tag);
      }
    }
    for (const ParentLink& link : parents_)
      if (present.find(link.child_tags[0]) != present.end())
        check_parent_for(link, b, out);
  }

private:
  struct ParentLink {
    std::vector<std::string> child_tags;
    std::vector<std::string> parent_tags;
  };

  static bool equal_rows(cif::Table::Row r1, cif::Table::Row r2) {
    assert(r1.size() == r2.size());
    for (int i = 0; i != (int) r1.size(); ++i) {
      if (cif::is_null(r1[i]) != cif::is_null(r2[i]))
        return false;
      if (cif::as_string(r1[i]) != cif::as_string(r2[i]))
        return false;
    }
    return true;
  }

  static std::string row_as_string(cif::Table::Row row) {
    return gemmi::join_str(row, '\v', [](const std::string& v) {
        return cif::is_null(v) ? std::string(1, '\0') : cif::as_string(v);
    });
  }

  void check_parent_for(const ParentLink& link, cif::Block& b, std::ostream& out) {
    cif::Table child_tab = b.find(link.child_tags);
    if (!child_tab.ok())
      return;
    cif::Table parent_tab = b.find(link.parent_tags);
    if (!parent_tab.ok()) {
      out << br(b.name) << "missing "
          << gemmi::join_str(link.parent_tags, '+') << "\n  parent of "
          << gemmi::join_str(link.child_tags, '+') << std::endl;
      return;
    }
    std::unordered_set<std::string> parent_hashes;
    //int dup_counter = 0;
    for (cif::Table::Row row : parent_tab) {
      parent_hashes.insert(row_as_string(row));
      /* apparently parent group doesn't need to be unique
      auto ret = parent_hashes.insert(row_as_string(row));
      if (!ret.second) {
        ++dup_counter;
        if (dup_counter < 2)
          out << br(b.name) << "duplicated parent group "
              << gemmi::join_str(link.parent_tags, '+') << ":\n  "
              << gemmi::join_str(row, '+') << std::endl;
      }
      */
    }
    int miss_counter = 0;
    for (cif::Table::Row row : child_tab) {
      if (std::all_of(row.begin(), row.end(), cif::is_null))
        continue;
      if (parent_hashes.count(row_as_string(row)) == 0) {
        ++miss_counter;
        if (miss_counter < 2)
          out << br(b.name)
              << gemmi::join_str(row, '+') << " from "
              << gemmi::join_str(link.child_tags, '+') << "\n  not in "
              << gemmi::join_str(link.parent_tags, '+') << std::endl;
      }
    }
    if (miss_counter > 1)
      out << "  [total " << miss_counter << " missing parents in this group]\n";
  }

  cif::Block* find_rules(const std::string& name) {
    auto iter = name_index_.find(name);
    return iter != name_index_.end() ? iter->second : nullptr;
  }

  void read_ddl1() {
    for (cif::Block& b : ddl_.blocks) {
      for (std::string& name : b.find_values("_name"))
        name_index_.emplace(cif::as_string(name), &b);
      if (b.name == "on_this_dictionary") {
        const std::string* dic_name = b.find_value("_dictionary_name");
        if (dic_name)
          dict_name_ = cif::as_string(*dic_name);
        const std::string* dic_ver = b.find_value("_dictionary_version");
        if (dic_ver)
          dict_version_ = cif::as_string(*dic_ver);
      }
    }
  }

  void read_ddl2() {
    for (cif::Block& block : ddl_.blocks) { // a single block is expected
      for (cif::Item& item : block.items) {
        if (item.type == cif::ItemType::Frame) {
          for (const char* tag : {"_item.name", "_category.id"}) {
            if (cif::Column col = item.frame.find_values(tag)) {
              for (const std::string& name : col)
                name_index_.emplace(cif::as_string(name), &item.frame);
              break;
            }
          }
        } else if (item.type == cif::ItemType::Pair) {
          if (item.pair[0] == "_dictionary.title")
            dict_name_ = item.pair[1];
          else if (item.pair[0] == "_dictionary.version")
            dict_version_ = item.pair[1];
        }
      }
      if (regex_enabled_)
        for (auto row : block.find("_item_type_list.", {"code", "construct"})) {
          if (cif::is_text_field(row[1]))
            // text field is problematic, but it's used only for "binary"
            // which in turn is never used
            continue;
          try {
            std::string re_str = row.str(1);
            // mmcif_pdbx_v50.dic uses custom flavour of regex:
            // character classes have unescaped \, but recognize \n, \t, etc.
            // Here is a quick fix:
            std::string::size_type pos = re_str.find("/\\{}");
            if (pos != std::string::npos)
              re_str.replace(pos, 4, "/\\\\{}");
            auto flag = std::regex::awk | std::regex::optimize;
            regexes_.emplace(row.str(0), std::regex(re_str, flag));
          } catch (const std::regex_error& e) {
            std::cout << "Note: DDL has invalid regex for " << row[0] << ":\n      "
                      << row.str(1) << "\n      "
                      << e.what() << '\n';
          }
        }
    }
  }

  int version_;
  bool regex_enabled_;
  bool use_context_;
  cif::Document ddl_;
  std::map<std::string, cif::Block*> name_index_;
  std::string dict_name_;
  std::string dict_version_;
  std::map<std::string, std::regex> regexes_;
  std::vector<ParentLink> parents_;
  // "_" or ".", used to unify handling of DDL1 and DDL2, for example when
  // reading _audit_conform_dict_version and _audit_conform.dict_version.
  char sep_;
};


class Validator1 {
public:
  Validator1(cif::Block& b) {
    if (const std::string* list = b.find_value("_list")) {
      if (*list == "yes")
        is_list_ = Trinary::Yes;
      else if (*list == "no")
        is_list_ = Trinary::No;
    }
    if (const std::string* type = b.find_value("_type"))
      if (*type == "numb")
        numb_only_ = true;
    /*
    // Hypotetically _type_conditions could be a list, but it never is.
    // It is commented out b/c we don't check esd/su - should we?
    const std::string* conditions = b.find_value("_type_conditions");
    if (conditions)
      has_su_ = (*conditions == "esd" || *conditions == "su");
    */
    const std::string* range = b.find_value("_enumeration_range");
    if (range) {
      size_t colon_pos = range->find(':');
      if (colon_pos != std::string::npos) {
        std::string low = range->substr(0, colon_pos);
        std::string high = range->substr(colon_pos+1);
        range_low_ = low.empty() ? -INFINITY : cif::as_number(low);
        range_high_ = high.empty() ? INFINITY : cif::as_number(high);
        has_range_ = true;
      }
    }
    for (const std::string& e : b.find_values("_enumeration"))
      enumeration_.emplace_back(cif::as_string(e));
  }

  // takes raw value
  bool validate_value(const std::string& value, std::string* msg) const {
    if (cif::is_null(value))
      return true;

    if (numb_only_ && !cif::is_numb(value)) {
      if (msg)
        *msg = "expected number, got: " + value;
      return false;
    }

    if (has_range_) {
      const double x = cif::as_number(value);
      if (x < range_low_ || x > range_high_) {
        if (msg)
          *msg = "value out of expected range: " + value;
        return false;
      }
    }

    if (!enumeration_.empty()) {
      std::string v = cif::as_string(value);
      if (!gemmi::in_vector(v, enumeration_)) {
        if (msg) {
          *msg = value + " is not one of the allowed values:";
          for (const std::string& e : enumeration_)
            *msg += "\n\t" + e;
        }
        return false;
      }
    }
    return true;
  }

  Trinary is_list() const { return is_list_; }

private:
  Trinary is_list_ = Trinary::Unset; // _list yes
  bool has_range_ = false;
  bool numb_only_ = false;
  double range_low_;
  double range_high_;
  std::vector<std::string> enumeration_;
  // type_construct regex - it is rarely used, ignore for now
  // type_conditions seq - seems to be never used, ignore it
  // For now we don't check at all relational attributes, i.e.
  // _category, _list_* and _related_*
};


class Validator2 {
public:
  enum class Type : char { Unset, Int, Float };
  enum class ItemContext { Default, Local, Deprecated };

  Validator2(cif::Block& b, const std::map<std::string, std::regex>& regexes) {
    if (const std::string* code = b.find_value("_item_type.code")) {
      type_code_ = cif::as_string(*code);
      if (type_code_ == "float") {
        type_ = Type::Float;
      } else if (type_code_ == "int") {
        type_ = Type::Int;
      } else {  // to make it faster, we don't use regex for int and float
        auto it = regexes.find(*code);
        if (it != regexes.end())
          re_ = &it->second;
      }
    }
    for (auto row : b.find("_item_range.", {"minimum", "maximum"}))
      range_.emplace_back(cif::as_number(row[0], -INFINITY),
                          cif::as_number(row[1], +INFINITY));
    for (const std::string& e : b.find_values("_item_enumeration.value"))
      enumeration_.emplace_back(cif::as_string(e));
    icase_ = (type_code_[0] == 'u');
    /* we could check for esd without value
    for (auto row : b.find("_item_related.", {"related_name", "function_code"})) {
      if (row[1] == "associated_value")
        associated_value_ = row.str(0);
    }
    */
    if (const std::string* context = b.find_value("_pdbx_item_context.type")) {
      if (*context == "WWPDB_LOCAL")
        context_ = ItemContext::Local;
      else if (*context == "WWPDB_DEPRECATED")
        context_ = ItemContext::Deprecated;
    }
  }

  // takes raw value
  bool validate_value(const std::string& value, std::string* msg) const {
    if (cif::is_null(value))
      return true;
    // int and float have hard-coded checks to avoid more expensive regex checks
    if (type_ == Type::Float && !cif::is_numb(value)) {
      if (msg)
        *msg = "expected number, got: " + value;
      return false;
    }
    if (type_ == Type::Int && !is_integer(value)) {
      if (msg)
        *msg = "expected integer, got: " + value;
    }
    if (!range_.empty()) {
      const double x = cif::as_number(value);
      bool matches = false;
      for (const auto& r : range_)
        if (r.first == r.second ? x == r.first
                                : r.first < x && x < r.second) {
          matches = true;
          break;
        }
      if (!matches) {
        *msg = "value out of expected range: " + value;
        return false;
      }
    }
    if (!enumeration_.empty() && !validate_enumeration(value, msg))
      return false;
    if (re_ && !std::regex_match(cif::as_string(value), *re_)) {
      *msg = value + " does not match the " + type_code_ + " regex";
      return false;
    }
    return true;
  }

  bool validate_enumeration(const std::string& val, std::string *msg) const {
    std::string v = cif::as_string(val);
    if (gemmi::in_vector(v, enumeration_))
      return true;
    if (icase_) {
      v = gemmi::to_lower(v);
      if (gemmi::in_vector_f([&v](const std::string& x) { return gemmi::iequal(x, v); },
                             enumeration_))
      return true;
    }
    if (msg) {
      *msg = val + " is not one of the allowed values:";
      for (const std::string& e : enumeration_)
        *msg += "\n\t" + e;
    }
    return false;
  }

  bool check_context_type(std::string* msg) const {
    if (context_ == ItemContext::Deprecated) {
      *msg = " is deprecated";
      return false;
    }
    if (context_ == ItemContext::Local) {
      *msg = " is for pdb internal use";
      return false;
    }
    return true;
  }

private:
  Type type_ = Type::Unset;
  bool icase_ = false;
  ItemContext context_ = ItemContext::Default;
  std::vector<std::string> enumeration_;
  std::string type_code_;
  std::vector<std::pair<double, double>> range_;
  const std::regex* re_ = nullptr;
};


void DDL::check_audit_conform(const cif::Document& doc, bool verbose) const {
  std::string audit_conform = "_audit_conform";
  audit_conform += sep_;
  for (const cif::Block& b : doc.blocks) {
    const std::string* raw_name = b.find_value(audit_conform + "dict_name");
    if (!raw_name)
      continue;
    std::string name = cif::as_string(*raw_name);
    if (name == dict_name_) {
      if (verbose)
        if (const std::string* dict_ver = b.find_value(audit_conform + "dict_version")) {
          std::string version = cif::as_string(*dict_ver);
          if (version != dict_version_)
            std::cout << "Note: " << br(b.name) << "conforms to " << name
                      << " ver. " << version << " while DDL has ver. " << dict_version_
                      << '\n';
        }
    } else {
      std::cout << "Note: " << br(b.name) << "dictionary name mismatch: " << name
                << " vs " << dict_name_ << '\n';
    }
  }
}

bool DDL::validate(cif::Document& doc, std::ostream& out, bool quiet) {
  std::string msg;
  bool ok = true;
  auto err = [&](const cif::Block& b, const cif::Item& item,
                 const std::string& s) {
    ok = false;
    out << doc.source << ":" << item.line_number
        << ' ' << br(b.name) << s << "\n";
  };
  for (cif::Block& b : doc.blocks) {
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Pair) {
        cif::Block* dict_block = find_rules(item.pair[0]);
        if (!dict_block) {
          if (!quiet)
            out << "Note: " << br(b.name) << "unknown tag " << item.pair[0] << '\n';
          continue;
        }
        // validate pair
        if (version_ == 1) {
          Validator1 tc(*dict_block);
          if (tc.is_list() == Trinary::Yes)
            err(b, item, item.pair[0] + " must be a list");
          if (!tc.validate_value(item.pair[1], &msg))
            err(b, item, msg);
        } else {
          Validator2 tc(*dict_block, regexes_);
          if (use_context_ && !tc.check_context_type(&msg))
            err(b, item, item.pair[0] + msg);
          if (!tc.validate_value(item.pair[1], &msg))
            err(b, item, msg);
        }
      } else if (item.type == cif::ItemType::Loop) {
        const size_t ncol = item.loop.tags.size();
        for (size_t i = 0; i != ncol; i++) {
          const std::string& tag = item.loop.tags[i];
          cif::Block* dict_block = find_rules(tag);
          if (!dict_block) {
            if (!quiet)
              out << "Note: " << br(b.name) << "unknown tag " << tag << '\n';
            continue;
          }
          // validate column in loop
          if (version_ == 1) {
            Validator1 tc(*dict_block);
            if (tc.is_list() == Trinary::No)
              err(b, item, tag + " in list");
            for (size_t j = i; j < item.loop.values.size(); j += ncol)
              if (!tc.validate_value(item.loop.values[j], &msg)) {
                err(b, item, tag + ": " + msg);
                break; // stop after first error to avoid clutter
              }
          } else {
            Validator2 tc(*dict_block, regexes_);
            if (use_context_ && !tc.check_context_type(&msg))
              err(b, item, tag + msg);
            for (size_t j = i; j < item.loop.values.size(); j += ncol)
              if (!tc.validate_value(item.loop.values[j], &msg)) {
                err(b, item, tag + ": " + msg);
                break; // stop after first error to avoid clutter
              }
          }
        }
      }
    }
  }
  return ok;
}

} // anonymous namespace


int GEMMI_MAIN(int argc, char **argv) {
#ifdef ANALYZE_RULES // for debugging only
  tao::pegtl::analyze<cif::rules::file>();
  tao::pegtl::analyze<cif::numb_rules::numb>();
#endif
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();

  bool quiet = p.options[Quiet];
  bool total_ok = true;
#if __GNUC__+0 == 4 && __GNUC_MINOR__+0 < 9 && !defined(__clang__)
  DDL dict(false);
  if (!p.options[NoRegex])
    std::cerr << "Note: regex support disabled on GCC 4.8\n";
#else
  DDL dict(!p.options[NoRegex]);
#endif
  dict.set_use_context(p.options[Context]);
  if (p.options[Ddl]) {
    try {
      for (option::Option* ddl = p.options[Ddl]; ddl; ddl = ddl->next()) {
        dict.open_file(ddl->arg);
        if (p.options[Parents])
          dict.read_parent_links();
      }
    } catch (std::runtime_error& e) {
      std::cerr << "Error when reading dictionary: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  }
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* path = p.nonOption(i);
    std::string msg;
    bool ok = true;
    try {
      if (p.options[Fast]) {
        ok = cif::check_syntax_any(gemmi::MaybeGzipped(path), &msg);
      } else {
        cif::Document d = cif::read(gemmi::MaybeGzipped(path));
        for (const cif::Block& block : d.blocks) {
          if (block.name == " ")
            std::cout << d.source << ": missing block name (bare data_)\n";
          check_empty_loops(block);
        }
        if (p.options[Stat])
          msg = token_stats(d);
        if (p.options[Ddl]) {
          dict.check_audit_conform(d, p.options[Verbose]);
          ok = dict.validate(d, std::cout, quiet);
          if (!p.options[NoMandatory])
            for (cif::Block& block : d.blocks)
              dict.check_mandatory_items(block, std::cout);
          if (!p.options[NoUniqueKeys])
            for (cif::Block& block : d.blocks)
              dict.check_unique_keys(block, std::cout);
          if (p.options[Parents])
            for (cif::Block& block : d.blocks)
              dict.check_parents(block, std::cout);
        }
        if (p.options[Monomer])
          check_monomer_doc(d);
      }
    } catch (std::runtime_error& e) {
      ok = false;
      msg = e.what();
    }
    if (!msg.empty())
      std::cout << msg << std::endl;

    if (p.options[Verbose])
      std::cout << (ok ? "OK" : "FAILED") << std::endl;
    total_ok = total_ok && ok;
  }
  return total_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
