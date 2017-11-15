// Copyright 2017 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/ddl.hpp"
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <iostream>
#ifdef ANALYZE_RULES
# include <tao/pegtl/analyze.hpp>
#endif

#define EXE_NAME "gemmi-validate"
#include "options.h"

namespace cif = gemmi::cif;

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
      if (item.type == cif::ItemType::Value) {
        nvals++;
        vals_by_type[static_cast<int>(item.valtype)]++;
      } else if (item.type == cif::ItemType::Frame) {
        nframes++;
      } else if (item.type == cif::ItemType::Loop) {
        nloops++;
        nlooptags += item.loop.tags.size();
        for (const cif::LoopTag& tag : item.loop.tags)
          looptags_by_type[static_cast<int>(tag.valtype)]++;
        nloopvals += item.loop.values.size();
      }
    }
  }
  info += format_7zd(d.blocks.size()) + " block(s)\n";
  info += format_7zd(nframes) + " frames\n";
  info += format_7zd(nvals) + " non-loop items:";
  if (vals_by_type[0] == nvals) {
    info += " (run with -t for type breakdown)";
  } else {
    for (int i = 1; i != 5; ++i)
      info += "  " + cif::value_type_to_str(static_cast<cif::ValueType>(i))
              + ":" + std::to_string(vals_by_type[i]);
  }
  info += "\n";
  info += format_7zd(nloops) + " loops w/\n";
  info += "        " + format_7zd(nlooptags) + " tags:";
  if (looptags_by_type[0] != nlooptags) {
    for (int i = 1; i != 5; ++i)
      info += "  " + cif::value_type_to_str(static_cast<cif::ValueType>(i))
              + ":" + std::to_string(looptags_by_type[i]);
  }
  info += "\n";
  info += "        " + format_7zd(nloopvals) + " values\n";
  return info;
}

// For now the infer_* functions are used only here, not sure where they belong
inline cif::ValueType infer_valtype_of_string(const std::string& val) {
  assert(!val.empty());
  if (val == ".")
    return cif::ValueType::Dot;
  if (val == "?")
    return cif::ValueType::QuestionMark;
  if (cif::is_numb(val))
    return cif::ValueType::Numb;
  return cif::ValueType::Char;
}

inline void infer_valtypes_in_items(std::vector<cif::Item>& items) {
  using namespace cif;
  for (Item& item : items)
      if (item.type == ItemType::Value) {
        item.valtype = infer_valtype_of_string(item.tv.value);
      } else if (item.type == ItemType::Loop) {
        for (size_t i = 0; i != item.loop.tags.size(); ++i) {
          ValueType& vt = item.loop.tags[i].valtype;
          for (const std::string& v : Column{&item, i}) {
            ValueType this_vt = infer_valtype_of_string(v);
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
        }
      } else if (item.type == ItemType::Frame) {
        infer_valtypes_in_items(item.frame.items);
      }
}

inline void infer_valtypes(cif::Document &d) {
  for (cif::Block& block : d.blocks)
    infer_valtypes_in_items(block.items);
}


enum OptionIndex { Fast=3, Stat, Types, Quiet, Ddl };
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None, "Usage: " EXE_NAME " [options] FILE [...]"
                                "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tDisplay version information and exit." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Types, 0, "t", "types", Arg::None, "  -t, --types  \t"
                                       "Break down token statistics by type." },
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Ddl, 0, "d", "ddl", Arg::Required,
                                   "  -d, --ddl=PATH  \tDDL for validation." },
  { 0, 0, 0, 0, 0, 0 }
};



int main(int argc, char **argv) {
#ifdef ANALYZE_RULES // for debugging only
  tao::pegtl::analyze<cif::rules::file>();
  tao::pegtl::analyze<cif::numb_rules::numb>();
#endif
  OptParser p;
  p.simple_parse(argc, argv, Usage);
  if (p.nonOptionsCount() == 0) {
    option::printUsage(std::cout, Usage);
    return 0;
  }

  bool quiet = p.options[Quiet];
  bool total_ok = true;
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    const char* path = p.nonOption(i);
    std::string msg;
    bool ok = true;
    try {
      if (p.options[Fast]) {
        ok = cif::check_syntax_any(gemmi::MaybeGzipped(path), &msg);
      } else {
        cif::Document d = cif::read(gemmi::MaybeGzipped(path));
        if (p.options[Types])
          infer_valtypes(d);
        if (p.options[Stat])
          msg = token_stats(d);
        if (p.options[Ddl]) {
          cif::DDL dict;
          for (option::Option* ddl = p.options[Ddl]; ddl; ddl = ddl->next())
            dict.open_file(ddl->arg);
          std::string ver_msg;
          dict.check_audit_conform(d, &ver_msg);
          if (!ver_msg.empty() && !quiet)
            std::cout << "Note: " << ver_msg << std::endl;
          ok = dict.validate(d, std::cout, quiet);
        }
      }
    } catch (std::runtime_error& e) {
      ok = false;
      msg = e.what();
    }
    if (!msg.empty())
      std::cout << msg << std::endl;

    if (!quiet)
      std::cout << (ok ? "OK" : "FAILED") << std::endl;
    total_ok = total_ok && ok;
  }
  return total_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
