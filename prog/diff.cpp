// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <algorithm>  // for find
#include "gemmi/read_cif.hpp"  // for read_cif_gz

#define GEMMI_PROG diff
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex {
  OnlyCategories=4
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n\nCompares categories and tags in CIF files." },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { OnlyCategories, 0, "q", "", Arg::None,
    "  -q  \tPrint only categories." },
  { 0, 0, 0, 0, 0, 0 }
};

struct DiffItem {
  char change;
  std::string str;
};

using Diff = std::vector<DiffItem>;

template<typename T>
Diff make_diff(const T& a, const T& b) {
  Diff diff;
  diff.reserve(a.size() + b.size());
  for (const std::string& s : b)
    diff.push_back(DiffItem{'+', s});
  size_t idx = 0;
  for (const std::string& s : a) {
    auto pos = std::find_if(diff.begin(), diff.end(),
                            [&](const DiffItem& d) { return d.str == s; });
    if (pos != diff.end()) {
      pos->change = ' ';
      idx = pos - diff.begin() + 1;
    } else {
      diff.insert(diff.begin() + idx, DiffItem{'-', s});
      idx++;
    }
  }
  return diff;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* path1 = p.nonOption(0);
  const char* path2 = p.nonOption(1);
  //int verbose = p.options[Verbose].count();
  try {
    // Starting like an unified diff (with "--- ") enables colordiff.
    printf("--- Reading %s\n", path1);
    cif::Document doc1 = gemmi::read_cif_gz(path1);
    printf("+++ Reading %s\n", path2);
    cif::Document doc2 = gemmi::read_cif_gz(path2);
    cif::Block& b1 = doc1.blocks.at(0);
    cif::Block& b2 = doc2.blocks.at(0);
    Diff category_diff = make_diff(b1.get_mmcif_category_names(),
                                   b2.get_mmcif_category_names());
    for (DiffItem& cat : category_diff) {
      cif::Table t1 = b1.find_mmcif_category(cat.str);
      cif::Table t2 = b2.find_mmcif_category(cat.str);
      size_t len1 = t1.length();
      size_t len2 = t2.length();
      printf("%c %-37s rows: %5zu", cat.change, cat.str.c_str(), len1);
      if (len2 != len1)
        printf("  -> %5zu", len2);
      putchar('\n');
      if (p.options[OnlyCategories])
        continue;
      size_t prefix_size = cat.str.size();
      Diff tag_diff = make_diff(t1.tags(), t2.tags());
      for (DiffItem& di : tag_diff)
        printf("%c       %s\n", di.change, di.str.c_str() + prefix_size);
    }
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

