// Copyright Global Phasing Ltd.

#include "gemmi/ddl.hpp"
#include "gemmi/numb.hpp"  // for is_numb
#include <cmath>      // for INFINITY
#include <algorithm>  // for find
#include <utility>    // for pair

namespace gemmi { namespace cif {

namespace {

std::string br(const cif::Block& block) {
  return "[" + block.name + "] ";
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

class Ddl1Rules {
public:
  Ddl1Rules(cif::Block& b) {
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
    // Hypothetically, _type_conditions could be a list, but it never is.
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
            *msg += "\n  " + e;
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


class Ddl2Rules {
public:
  enum class Type : char { Unset, Int, Float };

  Ddl2Rules(cif::Block& b, const Ddl* ddl, const std::string& tag) {
    if (const std::string* code = b.find_value("_item_type.code")) {
      type_code_ = cif::as_string(*code);
      if (type_code_ == "float") {
        type_ = Type::Float;
      } else if (type_code_ == "int") {
        type_ = Type::Int;
      } else {  // to make it faster, we don't use regex for int and float
        auto it = ddl->regexes().find(type_code_);
        if (it != ddl->regexes().end())
          re_ = &it->second;
        else
          ddl->logger.mesg("Bad DDL2: ", tag, " has undefined type: ", type_code_);
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

private:
  Type type_ = Type::Unset;
  bool icase_ = false;
  std::vector<std::string> enumeration_;
  std::string type_code_;
  std::vector<std::pair<double, double>> range_;
  const std::regex* re_ = nullptr;
};

std::string major_ver(const std::string &s) {
  return s.substr(0, s.find('.'));
}

} // anonymous namespace

// check if the dictionary name/version correspond to _audit_conform_dict_*
void Ddl::check_audit_conform(const cif::Document& doc) const {
  std::string audit_conform = "_audit_conform.";
  if (major_version == 1)
    audit_conform.back() = '_';  // for _audit_conform_dict_version, etc
  for (const cif::Block& b : doc.blocks) {
    const std::string* raw_name = b.find_value(audit_conform + "dict_name");
    if (!raw_name)
      continue;
    std::string name = cif::as_string(*raw_name);
    if (name != dict_name) {
      logger.note(br(b), "dictionary name mismatch: ", name, " vs ", dict_name);
    } else if (const std::string* dict_ver = b.find_value(audit_conform + "dict_version")) {
      std::string version = cif::as_string(*dict_ver);
      if (version != dict_version) {
        if (logger.threshold >= 7 || major_ver(version) != major_ver(dict_version))
          logger.note(br(b), "conforms to ", name, " ver. ", version,
                      " while DDL has ver. ", dict_version);
      }
    }
  }
}

void Ddl::check_mandatory_items(const cif::Block& b) const {
  // make a list of items in each category in the block
  std::map<std::string, std::vector<std::string>> categories;
  auto add_category = [&](const std::string& tag) {
    if (!find_rules(tag))  // ignore unknown tags
      return;
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
      warn(b, "category not in the dictionary: ", cat_name);
      continue;
    }
    // check context type
    if (use_context)
      if (const std::string* ct = cat_block->find_value("_pdbx_category_context.type"))
        logger.mesg(br(b), "category indicated as ", *ct, ": ", cat_name);
    std::vector<std::string> implicit_items;
    // check mandatory items
    for (auto i = name_index_.lower_bound(cat.first);
         i != name_index_.end() && gemmi::starts_with(i->first, cat.first);
         ++i) {
      for (auto row : i->second->find("_item.", {"name", "mandatory_code"})) {
        // mandatory_code can be one of: yes, not, implicit, implicit-ordinal
        char mc0 = row.str(1)[0];
        if (mc0 == 'y' && iequal(row.str(0), i->first) &&
            !gemmi::in_vector(i->first.substr(n), cat.second))
          warn(b, "missing mandatory tag: ", i->first);
        else if (mc0 == 'i')
          implicit_items.push_back(gemmi::to_lower(row.str(0)));
      }
    }
    // check key items
    for (const std::string& v : cat_block->find_values("_category_key.name")) {
      std::string key = gemmi::to_lower(cif::as_string(v));
      if (!gemmi::starts_with(key, cat.first))
        logger.level<3>("inconsistent dictionary: wrong _category_key for ", cat_name);
      if (!gemmi::in_vector(key.substr(n), cat.second) &&
          // check if the key item is implicit (a feature is used in mmcif_ddl.dic)
          !in_vector(key, implicit_items))
        warn(b, "missing category key: ", key);
    }
  }
}

void Ddl::check_unique_keys_in_loop(const cif::Loop& loop, const Block& block) const {
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
      warn(block, "category ", cat_name, " has ", dup_counter, " duplicated key",
           dup_counter == 1 ? ":\n  " : "s, first one:\n  ",
           gemmi::join_str(key_positions, " + ", [&](int k) {
             return gemmi::cat(loop.tags[k].substr(dot_pos+1), '=', loop.values[dup_row + k]);
           }));
    }
  }
}

void Ddl::check_parent_link(const ParentLink& link, const cif::Block& b) const {
  cif::Table child_tab = const_cast<cif::Block&>(b).find(link.child_tags);
  if (!child_tab.ok())
    return;
  cif::Table parent_tab = const_cast<cif::Block&>(b).find(link.parent_tags);
  if (!parent_tab.ok()) {
    warn(b, "missing ", tags_as_str(link.parent_tags),
         "\n  parent of ", tags_as_str(link.child_tags));
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
        warn(b, "duplicated parent group ", tags_as_str(link.parent_tags), ":\n  ",
             gemmi::join_str(row, '+'));
    }
    */
  }
  int miss_counter = 0;
  std::string first_miss;
  for (const cif::Table::Row row : child_tab) {
    if (std::all_of(row.begin(), row.end(), cif::is_null))
      continue;
    if (parent_hashes.count(row_as_string(row)) == 0) {
      if (miss_counter == 0)
        first_miss = gemmi::join_str(row, '+');
      ++miss_counter;
    }
  }
  if (miss_counter != 0)
    warn(b, miss_counter, " missing parent(s) in item-linked-group ", link.group,
         ":\n  child:  ", tags_as_str(link.child_tags),
         "\n  parent: ", tags_as_str(link.parent_tags),
         "\n  1st miss: ", first_miss);
}

void Ddl::check_parents(const cif::Block& b) const {
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
        warn(b, "parent tag of ", child_tag, " is absent: ", parent_tag);
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
          warn(b, missing_counter, " missing parent(s) of ", child_tag,
               " in ", parent_tag, ", first one: ", first_missing);
      }
    }
  }
  for (const ParentLink& link : parents_)
    if (present.find(link.child_tags[0]) != present.end())
      check_parent_link(link, b);
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

void Ddl::read_ddl2_block(cif::Block& block) {
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
      try {
        std::string re_str = row.str(1);
        // mmcif_pdbx_v50.dic uses custom flavour of regex:
        // character classes have unescaped \, but recognize \n, \t, etc.
        // Here is a quick fix:
        gemmi::replace_all(re_str, "/\\{}", "/\\\\{}");
        // in binary, \<newline> is apparently meant to be ignored
        gemmi::replace_all(re_str, "\\\n", "");
        gemmi::replace_all(re_str, "\\\r\n", "");
        auto flag = std::regex::awk | std::regex::optimize;
        regexes_.emplace(row.str(0), std::regex(re_str, flag));
      } catch (const std::regex_error& e) {
        logger.mesg("Bad DDL2: can't parse regex for '", row[0], "': ", e.what());
        // add an always-matching placeholder to avoid errors later
        regexes_.emplace(row.str(0), std::regex(".*"));
      }
    }

  if (use_parents) {
    // read child-parent relations from _pdbx_item_linked_group_list
    cif::Table tab = block.find("_pdbx_item_linked_group_list.",
                                {"child_category_id", "link_group_id",
                                 "child_name", "parent_name"});
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
        logger.mesg("Bad DDL2: linked group [", link.group,
                    "] has children in different categories");
        ok = false;
      }
      if (common_category(link.parent_tags) == 0) {
        logger.mesg("Bad DDL2: linked group [", link.group,
                    "] has parents in different categories");
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
          else if (it->second != parent_name)
            logger.mesg("Bad DDL2: different parents for ", child_name, ": ",
                        it->second, " and ", parent_name);
        }
    }
}

static const char* wrong_ddl2_context(const cif::Block& dict_block) {
  const std::string* context = dict_block.find_value("_pdbx_item_context.type");
  if (context && *context == "WWPDB_LOCAL")
    return " is for pdb internal use";
  if (context && *context == "WWPDB_DEPRECATED")
    return " is deprecated";
  return nullptr;
}

bool Ddl::validate_cif(const cif::Document& doc) const {
  bool ok = true;
  for (const cif::Block& b : doc.blocks)
    if (!validate_block(b, doc.source))
      ok = false;
  return ok;
}

bool Ddl::validate_block(const cif::Block& b, const std::string& source) const {
  bool ok = true;
  std::string msg;
  auto err = [&](const cif::Item& item, const std::string& s) {
    ok = false;
    logger.level<3>(source, ':', item.line_number, " [", b.name, "] ", s);
  };
  for (const cif::Item& item : b.items) {
    if (item.type == cif::ItemType::Pair) {
      const std::string& tag = item.pair[0];
      cif::Block* dict_block = find_rules(tag);
      if (!dict_block) {
        if (print_unknown_tags)
          warn(b, "unknown tag ", tag);
        continue;
      }
      // validate pair
      if (major_version == 1) {
        Ddl1Rules rules(*dict_block);
        if (rules.is_list() == Trinary::Yes)
          err(item, tag + " must be a list");
        if (!rules.validate_value(item.pair[1], &msg))
          err(item, msg);
      } else {
        if (use_context)
          if (const char* bad_ctx = wrong_ddl2_context(*dict_block))
            err(item, tag + bad_ctx);
        Ddl2Rules rules(*dict_block, this, tag);
        if (!rules.validate_value(item.pair[1], &msg))
          err(item, msg);
      }
    } else if (item.type == cif::ItemType::Loop) {
      const size_t ncol = item.loop.tags.size();
      for (size_t i = 0; i != ncol; i++) {
        const std::string& tag = item.loop.tags[i];
        cif::Block* dict_block = find_rules(tag);
        if (!dict_block) {
          if (print_unknown_tags)
            warn(b, "unknown tag ", tag);
          continue;
        }
        // validate column in loop
        if (major_version == 1) {
          Ddl1Rules rules(*dict_block);
          if (rules.is_list() == Trinary::No)
            err(item, tag + " in list");
          for (size_t j = i; j < item.loop.values.size(); j += ncol)
            if (!rules.validate_value(item.loop.values[j], &msg)) {
              err(item, cat(tag, ": ", msg));
              break; // stop after first error to avoid clutter
            }
        } else {
          if (use_context)
            if (const char* bad_ctx = wrong_ddl2_context(*dict_block))
              err(item, tag + bad_ctx);
          Ddl2Rules rules(*dict_block, this, tag);
          for (size_t j = i; j < item.loop.values.size(); j += ncol)
            if (!rules.validate_value(item.loop.values[j], &msg)) {
              err(item, cat(tag, ": ", msg));
              break; // stop after first error to avoid clutter
            }
        }
      }
    } else if (item.type == cif::ItemType::Frame) {
      validate_block(item.frame, source);
    }
  }

  if (major_version == 2) {
    if (use_mandatory)
      check_mandatory_items(b);
    if (use_unique_keys) {
      for (const cif::Item& item : b.items)
        if (item.type == cif::ItemType::Loop)
          check_unique_keys_in_loop(item.loop, b);
    }
    if (use_parents)
      check_parents(b);
  }

  return ok;
}

void Ddl::read_ddl(cif::Document&& doc) {
  ddl_docs_.emplace_back(new cif::Document(std::move(doc)));
  cif::Document& ddl_doc = *ddl_docs_.back();
  // perhaps we should check the content instead
  if (major_version == 0)
    major_version = (ddl_doc.blocks.size() > 1 ? 1 : 2);

  for (cif::Block& b : ddl_doc.blocks)
    if (major_version == 1)
      read_ddl1_block(b);
    else
      read_ddl2_block(b);
}

}} // namespace gemmi::cif
