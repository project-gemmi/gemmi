
#include "cif.hh"
#include "ddl.hh"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
//#include <pegtl/analyze.hh>

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

int main(int argc, char **argv) {
  bool quick = false;
  bool stats = false;
  bool type_breakdown = false;
  const char* ddl_path = nullptr;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-q") == 0) {
      quick = true;
      continue;
    }
    if (std::strcmp(argv[i], "-s") == 0) {
      stats = true;
      continue;
    }
    if (std::strcmp(argv[i], "-t") == 0) {
      type_breakdown = true;
      continue;
    }
    if (std::strncmp(argv[i], "--ddl=", 6) == 0) {
      ddl_path = argv[i] + 6;
      continue;
    }
    std::string msg;
    bool ok = true;
    if (quick) {
      ok = cif::check_file_syntax(argv[i], &msg);
    } else {
      //pegtl::analyze<cif::rules::file>();
      //pegtl::analyze<cif::numb_rules::numb>();
      cif::Document d;
      try {
        d.parse_file(argv[i]);
        if (type_breakdown)
          d.infer_valtypes();
        if (stats)
          msg = token_stats(d);
        if (ddl_path) {
          ddl::DDL dict;
          dict.open_file(ddl_path);
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
      } catch (std::runtime_error& e) {
        ok = false;
        msg = e.what();
      }
    }
    if (!msg.empty())
      std::cout << msg << std::endl;

    std::cout << (ok ? "OK" : "FAILED") << std::endl;
  }
  return 0;
}

// vim:sw=2:ts=2:et
