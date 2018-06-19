// Copyright 2017 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/ddl.hpp"
#include <cstdio>
#include <stdexcept>
#include <string>
#include <iostream>
#ifdef ANALYZE_RULES
# include <tao/pegtl/analyze.hpp>
#endif

#define GEMMI_PROG validate
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { Fast=3, Stat, Verbose, Quiet, Ddl };
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None, "Usage: " EXE_NAME " [options] FILE [...]"
                                "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tDisplay version information and exit." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Ddl, 0, "d", "ddl", Arg::Required,
                                   "  -d, --ddl=PATH  \tDDL for validation." },
  { 0, 0, 0, 0, 0, 0 }
};

enum class ValueType : unsigned char {
  NotSet,
  Char, // Line/Text?
  Numb, // Int/Float?
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

static std::string format_7zd(size_t k) {
  char buf[64];
  snprintf(buf, 63, "%7zu", k);
  return buf;
}

static std::string token_stats(const cif::Document& d) {
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
    info += "  " + value_type_to_str(static_cast<ValueType>(i))
            + ":" + std::to_string(vals_by_type[i]);
  info += "\n";
  info += format_7zd(nloops) + " loops w/\n";
  info += "        " + format_7zd(nlooptags) + " tags:";
  for (int i = 1; i != 5; ++i)
    info += "  " + value_type_to_str(static_cast<ValueType>(i))
            + ":" + std::to_string(looptags_by_type[i]);
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

    if (p.options[Verbose])
      std::cout << (ok ? "OK" : "FAILED") << std::endl;
    total_ok = total_ok && ok;
  }
  return total_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
