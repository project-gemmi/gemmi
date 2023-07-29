// Copyright 2017-2020 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/ddl.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/numb.hpp"
#include <cstdio>
#include <iostream>
#include <stdexcept>  // for std::runtime_error

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

} // anonymous namespace


int GEMMI_MAIN(int argc, char **argv) {
#ifdef ANALYZE_RULES // for debugging only
  tao::pegtl::analyze<cif::rules::file>();
  tao::pegtl::analyze<cif::numb_rules::numb>();
#endif
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();

  bool total_ok = true;
  cif::Ddl dict;
  dict.print_unknown_tags = !p.options[Quiet];
  dict.use_regex = !p.options[NoRegex];
  dict.use_context = p.options[Context];
  dict.use_parents = p.options[Parents];
  dict.use_mandatory = !p.options[NoMandatory];
  dict.use_unique_keys = !p.options[NoUniqueKeys];
  dict.print_extra_diagnostics = p.options[Verbose];
  if (p.options[Ddl]) {
    try {
      for (option::Option* ddl = p.options[Ddl]; ddl; ddl = ddl->next())
        dict.read_ddl(cif::read_file(ddl->arg), std::cout);
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
          dict.check_audit_conform(d, std::cout);
          ok = dict.validate_cif(d, std::cout);
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
