// Copyright 2017 Global Phasing Ltd.

#include "cif.hh"
#include "cifgz.hh"
#include "ddl.hh"
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <string>
#ifdef ANALYZE_RULES
# include <tao/pegtl/analyze.hpp>
#endif
#include <optionparser.h>

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
  info += "        " + format_7zd(d.comments.size()) + " comments\n";
  return info;
}


struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg) {
    if (option.arg != nullptr)
      return option::ARG_OK;
    if (msg)
      std::cerr << "Option '" << option.name << "' requires an argument\n";
    return option::ARG_ILLEGAL;
  }
};

enum OptionIndex { Unknown, Help, Fast, Stat, Types, Quiet, Ddl };
const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None, "Usage: validate [options] FILE [...]\n\n"
                                   "Options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Fast, 0, "f", "fast", Arg::None, "  -f, --fast  \tSyntax-only check." },
  { Stat, 0, "s", "stat", Arg::None, "  -s, --stat  \tShow token statistics" },
  { Types, 0, "t", "types", Arg::None, "  -t, --types  \t"
                                       "Break down token statistics by type." },
  { Quiet, 0, "q", "quiet", Arg::None, "  -q, --quiet  \tShow only errors." },
  { Ddl, 0, "d", "ddl", Arg::Required, "  -d, --ddl  \tDDL for validation." },
  { 0, 0, 0, 0, 0, 0 }
};

int main(int argc, char **argv) {
#ifdef ANALYZE_RULES // for debugging only
  tao::pegtl::analyze<cif::rules::file>();
  tao::pegtl::analyze<cif::numb_rules::numb>();
#endif
  if (argc < 1)
    return 2;
  option::Stats stats(usage, argc-1, argv+1);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc-1, argv+1, options.data(), buffer.data());
  if (parse.error() || options[Unknown]) {
    option::printUsage(std::cerr, usage);
    return 1;
  }
  if (options[Help] || parse.nonOptionsCount() == 0) {
    option::printUsage(std::cerr, usage);
    return 0;
  }

  for (int i = 0; i < parse.nonOptionsCount(); ++i) {
    const char* path = parse.nonOption(i);
    std::string msg;
    bool ok = true;
    try {
      if (options[Fast]) {
        ok = cif::check_file_syntax(path, &msg);
      } else {
        cif::Document d = cif::read_any(path);
        if (options[Types])
          cif::infer_valtypes(d);
        if (options[Stat])
          msg = token_stats(d);
        for (option::Option* ddl = options[Ddl]; ddl; ddl = ddl->next()) {
          cif::DDL dict;
          dict.open_file(ddl->arg);
          std::string ver_msg;
          dict.check_audit_conform(d, &ver_msg);
          if (!ver_msg.empty())
            std::cout << "Note: " << ver_msg << std::endl;
          std::vector<std::string> unknown;
          dict.validate(d, &unknown);
          if (!unknown.empty())
            std::cout << "Note: " << unknown.size() << " unknown tags"
                      << " - first one: " << unknown[0] << std::endl;
        }
      }
    } catch (std::runtime_error& e) {
      ok = false;
      msg = e.what();
    }
    if (!msg.empty())
      std::cout << msg << std::endl;

    if (!options[Quiet])
      std::cout << (ok ? "OK" : "FAILED") << std::endl;
  }
  return 0;
}

// vim:sw=2:ts=2:et
