
#include "cif.hh"
#include "ddl.hh"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
//#include <pegtl/analyze.hh>

std::string format_8zd(size_t k) {
  char buf[64];
  snprintf(buf, 63, "%8zd", k);
  return buf;
}

std::string token_stats(const cif::Document& d) {
  std::string info;
  bool show_per_block = d.blocks.size() < 5;
  size_t nframes = 0, nvals = 0, nloops = 0, nlooptags = 0, nloopvals = 0;
  for (const cif::Block& block : d.blocks) {
    for (const cif::Item& item : block.items) {
      if (item.type == cif::ItemType::Value) {
        nvals++;
      } else if (item.type == cif::ItemType::Frame) {
        nframes++;
      } else if (item.type == cif::ItemType::Loop) {
        nloops++;
        nlooptags += item.loop.tags.size();
        nloopvals += item.loop.values.size();
      }
    }
    if (show_per_block) {
      info += "[" + block.name + "]\n";
      info += format_8zd(nframes) + " frames\n";
      info += format_8zd(nvals) + " non-loop items\n";
      info += format_8zd(nloops) + " loops w/\n";
      info += "        " + format_8zd(nlooptags) + " tags\n";
      info += "        " + format_8zd(nloopvals) + " values\n";
      nframes = nvals = nloops = nlooptags = nloopvals = 0;
    }
  }
  if (!show_per_block) {
    info += format_8zd(d.blocks.size()) + " block(s)\n";
    info += format_8zd(nframes) + " frames\n";
    info += format_8zd(nvals) + " non-loop items\n";
    info += format_8zd(nloops) + " loops w/\n";
    info += "        " + format_8zd(nlooptags) + " tags\n";
    info += "        " + format_8zd(nloopvals) + " values\n";
  }
  info += std::to_string(d.comments.size()) + " comments in the file.\n";
  return info;
}

int main(int argc, char **argv) {
  bool quick = false;
  bool stats = false;
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
      cif::Document d;
      try {
        d.parse_file(argv[i]);
        if (stats)
          msg = "\n" + token_stats(d);
        if (ddl_path) {
          ddl::DDL1 dict;
          dict.open_file(ddl_path);
          std::string ver_msg;
          dict.check_audit_conform(d, &ver_msg);
          if (!ver_msg.empty())
            std::cout << "Note: " << ver_msg << std::endl;
          dict.validate(d);
        }
      } catch (std::runtime_error& e) {
        ok = false;
        msg = e.what();
      }
    }
    std::cout << (ok ? "OK" : "FAILED: ") << msg << std::endl;
  }
  return 0;
}

// vim:sw=2:ts=2:et
