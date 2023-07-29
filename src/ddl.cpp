// Copyright Global Phasing Ltd.

#include "gemmi/ddl.hpp"
#include "gemmi/numb.hpp"  // for is_numb
#include <cmath>      // for INFINITY
#include <algorithm>  // for find
#include <utility>    // for pair

namespace gemmi { namespace cif {

namespace {

std::string br(const std::string& block_name) {
  return "[" + block_name + "] ";
}

enum class Trinary : char { Unset, Yes, No };

bool is_integer(const std::string& s) {
  auto b = s.begin() + (s[0] == '+' || s[0] == '-' ? 1 : 0);
  return b != s.end() && std::all_of(b, s.end(), gemmi::is_digit);
}

// Returns after-dot position if all tags have the same category, otherwise 0.
size_t common_category(const std::vector<std::string>& v) {
  if (v.empty())
    return 0;
  size_t dot = v[0].find('.');
  if (dot == std::string::npos)
    return 0;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i].compare(0, dot+1, v[0], 0, dot+1) != 0)
      return 0;
  return dot+1;
}

std::string tags_as_str(const std::vector<std::string>& v) {
  if (v.empty())
    return "";
  std::string s = v[0];
  size_t pos = common_category(v);
  for (size_t i = 1; i < v.size(); ++i) {
    s += '+';
    s += v[i].substr(pos);
  }
  return s;
}

std::string row_as_string(cif::Table::Row row) {
  return gemmi::join_str(row, '\v', [](const std::string& v) {
      return cif::is_null(v) ? std::string(1, '\0') : cif::as_string(v);
  });
}

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

} // anonymous namespace

// check if the dictionary name/version correspond to _audit_conform_dict_*
void Ddl::check_audit_conform(const cif::Document& doc, std::ostream& out) const {
  std::string audit_conform = "_audit_conform.";
  if (major_version == 1)
    audit_conform.back() = '_';  // for _audit_conform_dict_version, etc
  for (const cif::Block& b : doc.blocks) {
    const std::string* raw_name = b.find_value(audit_conform + "dict_name");
    if (!raw_name)
      continue;
    std::string name = cif::as_string(*raw_name);
    if (name == dict_name) {
      if (print_extra_diagnostics)
        if (const std::string* dict_ver = b.find_value(audit_conform + "dict_version")) {
          std::string version = cif::as_string(*dict_ver);
          if (version != dict_version)
            out << "Note: " << br(b.name) << "conforms to " << name
                << " ver. " << version << " while DDL has ver. " << dict_version << '\n';
        }
    } else {
      out << "Note: " << br(b.name) << "dictionary name mismatch: " << name
          << " vs " << dict_name << '\n';
    }
  }
}

void Ddl::check_mandatory_items(const cif::Block& b, std::ostream& out) const {
  // make a list of items in each category in the block
  std::map<std::string, std::vector<std::string>> categories;
  auto add_category = [&](const std::string& tag) {
    size_t pos = tag.find('.');
    if (pos != std::string::npos)
      categories[to_lower(tag.substr(0, pos+1))].push_back(to_lower(tag.substr(pos+1)));
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
    size_t n = cat.first.size();
    std::string cat_name = cat.first.substr(1, n-2);
    cif::Block* cat_block = find_rules(cat_name);
    if (!cat_block) { // should not happen
      out << br(b.name) << "category not in the dictionary: " << cat_name << std::endl;
      continue;
    }
    // check context type
    if (use_context)
      if (const std::string* ct = cat_block->find_value("_pdbx_category_context.type"))
        out << br(b.name) << "category indicated as "
            << *ct << ": " << cat_name << std::endl;
    // check key items
    for (const std::string& v : cat_block->find_values("_category_key.name")) {
      std::string key = cif::as_string(v);
      if (!gemmi::istarts_with(key, cat.first))  // inconsistent dictionary
        out << "DDL2: wrong _category_key for " << cat_name << std::endl;
      if (!gemmi::in_vector(to_lower(key.substr(n)), cat.second))
        out << br(b.name) << "missing category key: " << key << std::endl;
    }
    // check mandatory items
    for (auto i = name_index_.lower_bound(cat.first);
         i != name_index_.end() && gemmi::starts_with(i->first, cat.first);
         ++i) {
      for (auto row : i->second->find("_item.", {"name", "mandatory_code"}))
        if (row.str(1)[0] == 'y' && iequal(row.str(0), i->first) &&
            !gemmi::in_vector(i->first.substr(n), cat.second))
          out << br(b.name) << "missing mandatory tag: " << i->first << std::endl;
    }
  }
}

void Ddl::check_unique_keys_in_loop(const cif::Loop& loop, std::ostream& out,
                                    const std::string& block_name) const {
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

static void check_parent_for(const std::vector<std::string>& child_tags,
                             const std::vector<std::string>& parent_tags,
                             const cif::Block& b,
                             std::ostream& out) {
  cif::Table child_tab = const_cast<cif::Block&>(b).find(child_tags);
  if (!child_tab.ok())
    return;
  cif::Table parent_tab = const_cast<cif::Block&>(b).find(parent_tags);
  if (!parent_tab.ok()) {
    out << br(b.name) << "missing " << tags_as_str(parent_tags)
        << "\n  parent of " << tags_as_str(child_tags) << std::endl;
    return;
  }
  std::unordered_set<std::string> parent_hashes;
  //int dup_counter = 0;
  for (const cif::Table::Row row : parent_tab) {
    parent_hashes.insert(row_as_string(row));
    /* apparently parent group doesn't need to be unique
    auto ret = parent_hashes.insert(row_as_string(row));
    if (!ret.second) {
      ++dup_counter;
      if (dup_counter < 2)
        out << br(b.name) << "duplicated parent group "
            << tags_as_str(parent_tags) << ":\n  "
            << gemmi::join_str(row, '+') << std::endl;
    }
    */
  }
  int miss_counter = 0;
  for (const cif::Table::Row row : child_tab) {
    if (std::all_of(row.begin(), row.end(), cif::is_null))
      continue;
    if (parent_hashes.count(row_as_string(row)) == 0) {
      ++miss_counter;
      if (miss_counter < 2)
        out << br(b.name) << gemmi::join_str(row, '+')
            << " from " << tags_as_str(child_tags)
            << "\n  not in " << tags_as_str(parent_tags) << std::endl;
    }
  }
  if (miss_counter > 1)
    out << "  [total " << miss_counter << " missing parents in this group]\n";
}

void Ddl::check_parents(const cif::Block& b, std::ostream& out) const {
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
  for (const std::string& child_tag : present) {
    std::string first_missing;
    auto it = item_parents_.find(child_tag);
    if (it != item_parents_.end()) {
      const std::string& parent_tag = it->second;
      if (present.count(parent_tag) == 0) {
        out << br(b.name) << "parent tag of " << child_tag << " is absent: "
            << parent_tag << std::endl;
      } else {
        std::unordered_set<std::string> parent_hashes;
        for (const std::string& parent : const_cast<cif::Block&>(b).find_values(parent_tag))
          if (!cif::is_null(parent))
            parent_hashes.insert(cif::as_string(parent));
        int missing_counter = 0;
        for (const std::string& child : const_cast<cif::Block&>(b).find_values(child_tag))
          if (!cif::is_null(child) && parent_hashes.count(cif::as_string(child)) == 0)
            if (missing_counter++ == 0)
              first_missing = cif::as_string(child);
        if (missing_counter != 0)
          out << br(b.name) << missing_counter << " missing parent(s) of " << child_tag
              << " in " << parent_tag << ", first one: " << first_missing << std::endl;
      }
    }
  }
  for (const ParentLink& link : parents_)
    if (present.find(link.child_tags[0]) != present.end())
      check_parent_for(link.child_tags, link.parent_tags, b, out);
}

void Ddl::read_ddl1_block(cif::Block& block) {
  for (std::string& name : block.find_values("_name"))
    name_index_.emplace(to_lower(cif::as_string(name)), &block);
  if (block.name == "on_this_dictionary") {
    const std::string* dic_name = block.find_value("_dictionary_name");
    if (dic_name)
      dict_name = cif::as_string(*dic_name);
    const std::string* dic_ver = block.find_value("_dictionary_version");
    if (dic_ver)
      dict_version = cif::as_string(*dic_ver);
  }
}

void Ddl::read_ddl2_block(cif::Block& block, std::ostream& out) {
  for (cif::Item& item : block.items) {
    if (item.type == cif::ItemType::Frame) {
      for (const char* tag : {"_item.name", "_category.id"}) {
        if (cif::Column col = item.frame.find_values(tag)) {
          for (const std::string& name : col)
            name_index_.emplace(to_lower(cif::as_string(name)), &item.frame);
          break;
        }
      }
    } else if (item.type == cif::ItemType::Pair) {
      if (item.pair[0] == "_dictionary.title")
        dict_name = item.pair[1];
      else if (item.pair[0] == "_dictionary.version")
        dict_version = item.pair[1];
    }
  }

  if (use_regex)
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
        out << "Note: Ddl has invalid regex for " << row[0] << ":\n      "
            << row.str(1) << "\n      "
            << e.what() << '\n';
      }
    }

  if (use_parents) {
    // read child-parent relations from _pdbx_item_linked_group_list
    cif::Table tab = block.find("_pdbx_item_linked_group_list.",
                                {"child_category_id", "link_group_id",
                                 "child_name", "parent_name"});
    std::string prev_group;
    ParentLink* it = nullptr;
    for (cif::Table::Row row : tab) {
      std::string group = row.str(0);
      group += ' ';
      group += row[1];
      // here we assume the table is ordered
      if (!it || group != it->group) {
        parents_.emplace_back();
        it = &parents_.back();
        it->group = group;
      }
      it->child_tags.push_back(row.str(2));
      it->parent_tags.push_back(row.str(3));
    }

    // sanity check
    for (ParentLink& link : parents_) {
      bool ok = true;
      if (common_category(link.child_tags) == 0) {
        if (print_extra_diagnostics)
          out << "Bad DDL2: linked group [" << link.group
              << "] has children in different categories" << std::endl;
        ok = false;
      }
      if (common_category(link.parent_tags) == 0) {
        if (print_extra_diagnostics)
          out << "Bad DDL2: linked group [" << link.group
              << "] has parents in different categories" << std::endl;
        ok = false;
      }
      if (!ok) {
        // the simplest fix: leave only the first relation
        link.child_tags.resize(1);
        link.parent_tags.resize(1);
      }
    }
  }

  if (use_parents)
    // check child-parent relations from _item_linked
    for (cif::Item& item : block.items) {
      if (item.type == cif::ItemType::Frame)
        for (auto row : item.frame.find("_item_linked.", {"child_name", "parent_name"})) {
          std::string child_name = row.str(0);
          std::string parent_name = row.str(1);
          auto it = item_parents_.find(child_name);
          if (it == item_parents_.end())
            item_parents_.emplace(child_name, parent_name);
          else if (it->second != parent_name && print_extra_diagnostics)
            out << "Bad DDL2: different parents for " << child_name << ": "
                << it->second << " and " << parent_name << std::endl;
        }
    }
}

bool Ddl::validate_cif(const cif::Document& doc, std::ostream& out) const {
  std::string msg;
  bool ok = true;
  auto err = [&](const cif::Block& b, const cif::Item& item,
                 const std::string& s) {
    ok = false;
    out << doc.source << ":" << item.line_number
        << ' ' << br(b.name) << s << "\n";
  };

  for (const cif::Block& b : doc.blocks) {
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Pair) {
        cif::Block* dict_block = find_rules(item.pair[0]);
        if (!dict_block) {
          if (print_unknown_tags)
            out << "Note: " << br(b.name) << "unknown tag " << item.pair[0] << '\n';
          continue;
        }
        // validate pair
        if (major_version == 1) {
          Validator1 tc(*dict_block);
          if (tc.is_list() == Trinary::Yes)
            err(b, item, item.pair[0] + " must be a list");
          if (!tc.validate_value(item.pair[1], &msg))
            err(b, item, msg);
        } else {
          Validator2 tc(*dict_block, regexes_);
          if (use_context && !tc.check_context_type(&msg))
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
            if (print_unknown_tags)
              out << "Note: " << br(b.name) << "unknown tag " << tag << '\n';
            continue;
          }
          // validate column in loop
          if (major_version == 1) {
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
            if (use_context && !tc.check_context_type(&msg))
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

    if (major_version == 2) {
      if (use_mandatory)
        check_mandatory_items(b, out);
      if (use_unique_keys) {
        for (const cif::Item& item : b.items)
          if (item.type == cif::ItemType::Loop)
            check_unique_keys_in_loop(item.loop, out, b.name);
      }
      if (use_parents)
        check_parents(b, out);
    }
  }

  return ok;
}

void Ddl::read_ddl(cif::Document&& doc, std::ostream& out) {
  ddl_docs_.emplace_back(new cif::Document(std::move(doc)));
  cif::Document& ddl_doc = *ddl_docs_.back();
  // perhaps we should check the content instead
  if (major_version == 0)
    major_version = (ddl_doc.blocks.size() > 1 ? 1 : 2);

  for (cif::Block& b : ddl_doc.blocks)
    if (major_version == 1)
      read_ddl1_block(b);
    else
      read_ddl2_block(b, out);
}

}} // namespace gemmi::cif
