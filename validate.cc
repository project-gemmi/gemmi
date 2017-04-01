// Copyright 2017 Global Phasing Ltd.

#include <iostream>
#include "cif.hh"
#include "cifgz.hh"
#include "ddl.hh"
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <string>
//#include <tao/pegtl/analyze.hpp>
#define CLARA_CONFIG_MAIN
#include <clara.h>

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

struct Options {
  std::string process_name;
  bool help = false;
  bool fast = false;
  bool stats = false;
  bool type_breakdown = false;
  std::string ddl_path;
  std::vector<std::string> paths;
  void add_path(const std::string& p) { paths.push_back(p); }
};

int main(int argc, char **argv) {
  Clara::CommandLine<Options> cli;
  cli["-h"]["--help"].describe("display usage information")
    .bind(&Options::help);
  cli["-f"]["--fast"].describe("fast syntax-only check")
    .bind(&Options::fast);
  cli["-s"]["--stats"].describe("show token statistics")
    .bind(&Options::stats);
  cli["-t"]["--types"].describe("show type breakdown in token statistics")
    .bind(&Options::type_breakdown);
  cli["-d"]["--ddl"].describe("DDL for validation")
    .bind(&Options::ddl_path, "file.dic");
  cli[Clara::_].bind(&Options::add_path, "file");
  cli.setThrowOnUnrecognisedTokens(true);
  cli.bindProcessName(&Options::process_name);
  Options options;
  try {
    cli.parseInto(Clara::argsToVector(argc, argv), options);
  } catch (std::exception& e) {
    std::cerr << "Error: " << e.what()
              << "\nOption -h shows usage." << std::endl;
    return 1;
  }
  if (options.help) {
    cli.usage(std::cerr, options.process_name);
    return 0;
  }

  for (const std::string& path : options.paths) {
    std::string msg;
    bool ok = true;
    try {
      if (options.fast) {
        ok = cif::check_file_syntax(path, &msg);
      } else {
        //tao::pegtl::analyze<cif::rules::file>();
        //tao::pegtl::analyze<cif::numb_rules::numb>();
        cif::Document d = cif::read_any(path);
        if (options.type_breakdown)
          cif::infer_valtypes(d);
        if (options.stats)
          msg = token_stats(d);
        if (!options.ddl_path.empty()) {
          cif::DDL dict;
          dict.open_file(options.ddl_path);
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

    std::cout << (ok ? "OK" : "FAILED") << std::endl;
  }
  return 0;
}

// vim:sw=2:ts=2:et
