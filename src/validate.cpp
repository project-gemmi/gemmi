// Copyright 2017-2020 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/numb.hpp"
#include "gemmi/tostr.hpp"
#include <cstdio>
#include <cmath>      // for INFINITY
#include <algorithm>  // for find
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
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

enum OptionIndex { Quiet=4, Fast, Stat, Ddl, NoRegex, Monomer };
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None, "Usage: " EXE_NAME " [options] FILE [...]"
                                "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Ddl, 0, "d", "ddl", Arg::Required,
                                   "  -d, --ddl=PATH  \tDDL for validation." },
  { NoRegex, 0, "", "no-regex", Arg::None,
    "  --no-regex  \tSkip regex checking (when using DDL2)" },
  { Monomer, 0, "m", "monomer", Arg::None,
    "  -m, --monomer  \tExtra checks for Refmac dictionary files." },
  { 0, 0, 0, 0, 0, 0 }
};

// basic types, used for token statistics only, cf. ValType
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

std::string token_stats(const cif::Document& d) {
  std::string info;
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
  info += format_7zd(d.blocks.size()) + " block(s)\n";
  info += format_7zd(nframes) + " frames\n";
  info += format_7zd(nvals) + " non-loop items:";
  for (int i = 1; i != 5; ++i)
    info += gemmi::tostr("  ", value_type_to_str(static_cast<ValueType>(i)),
                         ':', vals_by_type[i]);
  info += "\n";
  info += format_7zd(nloops) + " loops w/\n";
  info += "        " + format_7zd(nlooptags) + " tags:";
  for (int i = 1; i != 5; ++i)
    info += gemmi::tostr("  ", value_type_to_str(static_cast<ValueType>(i)),
                         ':', looptags_by_type[i]);
  info += "\n";
  info += "        " + format_7zd(nloopvals) + " values\n";
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

enum class ValType : char {
  Unset,
  Int,  // only in DDL2, in DDL1 it is also numb
  Numb,
  Any
};

bool is_integer(const std::string& s) {
  auto b = s.begin() + (s[0] == '+' || s[0] == '-' ? 1 : 0);
  return b != s.end() && std::all_of(b, s.end(), gemmi::is_digit);
}

class TypeCheckCommon {
protected:
  ValType type_ = ValType::Unset;
  bool has_su_ = false; // _type_conditions esd|su
  bool range_inclusive_ = false;
  bool icase_ = false;
  std::vector<std::pair<double, double>> range_;
  std::vector<std::string> enumeration_;

  // takes raw value
  bool validate_value_cmn(const std::string& value, std::string* msg) const {
    if (type_ == ValType::Numb && !cif::is_numb(value)) {
      if (msg)
        *msg = "expected number, got: " + value;
      return false;
    }
    if (type_ == ValType::Int && !is_integer(value)) {
      if (msg)
        *msg = "expected integer, got: " + value;
    }
    // ignoring has_su_ - not sure if we should check it
    if (!range_.empty() && !validate_range(value, msg))
      return false;
    if (!enumeration_.empty() && !validate_enumeration(value, msg))
      return false;
    return true;
  }

public:
  bool validate_tag(std::string*) const { return true; }

private:
  bool validate_range(const std::string& value, std::string *msg) const {
    const double x = cif::as_number(value);
    for (const auto& r : range_)
      if (r.first == r.second ? x == r.first
                              : r.first < x && x < r.second)
        return true;
    if (msg)
      *msg = "value out of expected range: " + value;
    return false;
  }

  // Takes raw value.
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
};

class DDL;

class TypeCheckDDL1 : public TypeCheckCommon {
public:
  void from_block(cif::Block& b, const DDL&) {
    if (const std::string* list = b.find_value("_list")) {
      if (*list == "yes")
        is_list_ = Trinary::Yes;
      else if (*list == "no")
        is_list_ = Trinary::No;
    }
    if (const std::string* type = b.find_value("_type"))
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
        range_.emplace_back(low.empty() ? -INFINITY : cif::as_number(low),
                            high.empty() ? INFINITY : cif::as_number(high));
      }
    }
    for (const std::string& e : b.find_loop("_enumeration"))
      enumeration_.emplace_back(cif::as_string(e));
  }

  bool validate_value(const std::string& value, std::string* msg) const {
    if (cif::is_null(value))
      return true;
    return validate_value_cmn(value, msg);
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
  enum class ItemContext { Default, Local, Deprecated };
  void from_block(cif::Block& b, const DDL& ddl);

  bool validate_value(const std::string& value, std::string* msg) const {
    if (cif::is_null(value))
      return true;
    if (!validate_value_cmn(value, msg))
      return false;
    if (re_ && !std::regex_match(cif::as_string(value), *re_)) {
      *msg = value + " does not match the " + type_code_ + " regex";
      return false;
    }
    return true;
  }

  bool validate_tag(std::string* msg) const {
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

  Trinary is_list() const { return Trinary::Unset; }

private:
  std::string type_code_;
  ItemContext context_ = ItemContext::Default;
  const std::regex* re_ = nullptr;
};


// Class DDL that represents DDL1 or DDL2 dictionary (ontology).
class DDL {
public:
  DDL(bool enable_regex) : regex_enabled_(enable_regex) {}

  void open_file(const std::string& filename) {
    ddl_ = cif::read_file(filename);
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
  // check if the dictionary name/version correspond to _audit_conform_dict_*
  void check_audit_conform(const cif::Document& doc) const;

  bool validate(cif::Document& doc, std::ostream& out, bool quiet) {
    return version_ == 1 ? do_validate<TypeCheckDDL1>(doc, out, quiet)
                         : do_validate<TypeCheckDDL2>(doc, out, quiet);
  }

  const std::regex* get_regex_ptr(const std::string& code) const {
    auto it = regexes_.find(code);
    return it != regexes_.end() ? &it->second : nullptr;
  }

private:
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
          for (const std::string& name : item.frame.find_values("_item.name"))
            name_index_.emplace(cif::as_string(name), &item.frame);
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

  template<class TypeCheckDDL>
  bool do_validate(cif::Document& doc, std::ostream& out, bool quiet);

  int version_;
  cif::Document ddl_;
  std::unordered_map<std::string, cif::Block*> name_index_;
  std::string dict_name_;
  std::string dict_version_;
  bool regex_enabled_;
  std::map<std::string, std::regex> regexes_;
  // "_" or ".", used to unify handling of DDL1 and DDL2, for example when
  // reading _audit_conform_dict_version and _audit_conform.dict_version.
  std::string sep_;
};

void TypeCheckDDL2::from_block(cif::Block& b, const DDL& ddl) {
  if (const std::string* code = b.find_value("_item_type.code")) {
    type_code_ = cif::as_string(*code);
    if (type_code_ == "float")
      type_ = ValType::Numb;
    else if (type_code_ == "int")
      type_ = ValType::Int;
    else
      // to make it faster, we don't use regex for int and float
      re_ = ddl.get_regex_ptr(*code);
  }
  for (auto row : b.find("_item_range.", {"minimum", "maximum"}))
    range_.emplace_back(cif::as_number(row[0], -INFINITY),
                        cif::as_number(row[1], +INFINITY));
  for (const std::string& e : b.find_loop("_item_enumeration.value"))
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

void DDL::check_audit_conform(const cif::Document& doc) const {
  std::string audit_conform = "_audit_conform" + sep_;
  for (const cif::Block& b : doc.blocks) {
    const std::string* raw_name = b.find_value(audit_conform + "dict_name");
    if (!raw_name) {
      std::cout << "Note: the cif file (block " << b.name << ") is missing "
                << audit_conform << "dict_name\n";
      continue;
    }
    std::string name = cif::as_string(*raw_name);
    if (name == dict_name_) {
      const std::string* dict_ver = b.find_value(audit_conform + "dict_version");
      if (dict_ver) {
        std::string version = cif::as_string(*dict_ver);
        if (version != dict_version_)
          std::cout << "Note: CIF conforms to " << name << " ver. " << version
                    << " while DDL has ver. " << dict_version_ << '\n';
      }
    } else {
      std::cout << "Note: dictionary name mismatch in " << b.name << ": " << name
                << " vs " << dict_name_ << '\n';
    }
  }
}

template <class TypeCheckDDL>
bool DDL::do_validate(cif::Document& doc, std::ostream& out, bool quiet) {
  std::string msg;
  bool ok = true;
  auto err = [&](const cif::Block& b, const cif::Item& item,
                 const std::string& s) {
    ok = false;
    out << doc.source << ":" << item.line_number
        << " in data_" << b.name << ": " << s << "\n";
  };
  for (cif::Block& b : doc.blocks) {
    for (const cif::Item& item : b.items) {
      if (item.type == cif::ItemType::Pair) {
        cif::Block* dict_block = find_rules(item.pair[0]);
        if (!dict_block) {
          if (!quiet)
            out << "Note: unknown tag in " << b.name << ": " << item.pair[0] << '\n';
          continue;
        }
        TypeCheckDDL tc;
        tc.from_block(*dict_block, *this);
        if (tc.is_list() == Trinary::Yes)
          err(b, item, item.pair[0] + " must be a list");
        if (!tc.validate_tag(&msg))
          err(b, item, item.pair[0] + msg);
        if (!tc.validate_value(item.pair[1], &msg))
          err(b, item, msg);
      } else if (item.type == cif::ItemType::Loop) {
        const size_t ncol = item.loop.tags.size();
        for (size_t i = 0; i != ncol; i++) {
          const std::string& tag = item.loop.tags[i];
          cif::Block* dict_block = find_rules(tag);
          if (!dict_block) {
            if (!quiet)
              out << "Note: unknown tag in " << b.name << ": " << tag << '\n';
            continue;
          }
          TypeCheckDDL tc;
          tc.from_block(*dict_block, *this);
          if (tc.is_list() == Trinary::No)
            err(b, item, tag + " in list");
          if (!tc.validate_tag(&msg))
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
#if __GNUC__+0 == 4 && __GNUC_MINOR__+0 < 9
  DDL dict(false);
  if (!p.options[NoRegex])
    std::cerr << "Note: regex support disabled on GCC 4.8\n";
#else
  DDL dict(!p.options[NoRegex]);
#endif
  if (p.options[Ddl])
    for (option::Option* ddl = p.options[Ddl]; ddl; ddl = ddl->next())
      dict.open_file(ddl->arg);
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* path = p.nonOption(i);
    std::string msg;
    bool ok = true;
    try {
      if (p.options[Fast]) {
        ok = cif::check_syntax_any(gemmi::MaybeGzipped(path), &msg);
      } else {
        cif::Document d = cif::read(gemmi::MaybeGzipped(path));
        for (const cif::Block& block : d.blocks)
          check_empty_loops(block);
        if (p.options[Stat])
          msg = token_stats(d);
        if (p.options[Ddl]) {
          if (p.options[Verbose])
            dict.check_audit_conform(d);
          ok = dict.validate(d, std::cout, quiet);
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

// vim:sw=2:ts=2:et:path^=../include,../third_party
