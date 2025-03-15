// Copyright 2023 Global Phasing Ltd.

#include <cstdio>   // for printf, fprintf
#include <algorithm>  // for find
#include "gemmi/read_cif.hpp"  // for read_cif_gz
#include "gemmi/pdb_id.hpp"  // for expand_if_pdb_code
#include "gemmi/dirwalk.hpp"  // for glob_match

#define GEMMI_PROG cifdiff
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex {
  Tag=4, OnlyCategories, NoComparison,
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n"
    " " EXE_NAME " [options] FILE1.cif FILE2.cif\n"
    " " EXE_NAME " [options] -n FILE.cif\n\n"
    "Compares (or just prints) categories and tags in CIF files." },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Tag, 0, "t", "", Arg::Required,
    "  -t TAG \tCompare values of TAG." },
  { OnlyCategories, 0, "q", "", Arg::None,
    "  -q  \tPrint only categories." },
  { NoComparison, 0, "n", "", Arg::None,
    "  -n  \tNo comparison, just list categories and tags." },
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

void compare_tag_values(cif::Block& b1, cif::Block& b2, const std::string& tag) {
  printf("  checking %s ...\n", tag.c_str());
  auto col1 = b1.find_values(tag);
  auto col2 = b2.find_values(tag);
  if (!col1 || !col2) {
    if (!col1)
      std::printf("    NOT FOUND in file 1 block %s\n", b1.name.c_str());
    if (!col2)
      std::printf("    NOT FOUND in file 2 block %s\n", b2.name.c_str());
    return;
  }
  int n1 = col1.length();
  int n2 = col2.length();
  int n = std::min(n1, n2);
  if (n1 != n2) {
    std::printf("-   number of values: %d\n", n1);
    std::printf("+   number of values: %d\n", n2);
    printf("  comparing the first %d values...\n", n);
  }
  int diff_count = 0;
  for (int i = 0; i < n; ++i)
    if (col1[i] != col2[i] && col1.str(i) != col2.str(i) && diff_count++ < 4) {
      std::printf("-   value %d: %s\n", i, col1[i].c_str());
      std::printf("+   value %d: %s\n", i, col2[i].c_str());
    }
  if (diff_count >= 4)
    printf("    ...\n"
           "  %d of %d values differ\n", diff_count, n);
  else if (diff_count == 0)
    printf("  %d identical values\n", n);
}

template<typename Func>
void for_each_tag(const cif::Block& block, const Func& func) {
  for (const cif::Item& item : block.items) {
    if (item.type == cif::ItemType::Pair) {
      func(item.pair[0]);
    } else if (item.type == cif::ItemType::Loop) {
      for (const std::string& t : item.loop.tags)
        func(t);
    }
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.check_exclusive_pair(Tag, NoComparison);
  bool one_file = p.options[NoComparison];
  p.require_positional_args(one_file ? 1 : 2);
  const char* path1 = p.nonOption(0);
  const char* path2 = one_file ? nullptr : p.nonOption(1);
  //int verbose = p.options[Verbose].count();
  try {
    // Starting like an unified diff (with "--- ") enables colordiff.
    printf("%sReading %s\n", one_file ? "" : "--- ", path1);
    cif::Document doc1 = gemmi::read_cif_or_mmjson_gz(gemmi::expand_if_pdb_code(path1));
    cif::Document doc2;
    if (!one_file) {
      printf("+++ Reading %s\n", path2);
      doc2 = gemmi::read_cif_or_mmjson_gz(gemmi::expand_if_pdb_code(path2));
    }
    for (size_t i = 0; i < std::max(doc1.blocks.size(), doc2.blocks.size()); ++i) {
      cif::Block* b1 = i < doc1.blocks.size() ? &doc1.blocks[i] : nullptr;
      // NoComparison mode is implemented as comparing Block with itself
      // (inefficient, but simple).
      cif::Block* b2 = b1;
      if (!one_file)
        b2 = i < doc2.blocks.size() ? &doc2.blocks[i] : nullptr;
      if (b1 && b2 && b1->name == b2->name) {
        printf("  =========[  %s  ]=========\n", b1->name.c_str());
      } else {
        if (b1)
          printf("- =========[  %s  ]=========\n", b1->name.c_str());
        if (b2)
          printf("+ =========[  %s  ]=========\n", b2->name.c_str());
        if (!b1 || !b2)
          continue;
      }
      if (p.options[Tag]) {
        for (const option::Option* opt = p.options[Tag]; opt; opt = opt->next()) {
          std::string tag = opt->arg;
          if (tag.find_first_of("?*") != std::string::npos) {
            std::vector<std::string> matching_tags;
            for_each_tag(*b1, [&](const std::string& t) {
                if (gemmi::glob_match(tag, t))
                  matching_tags.push_back(t);
            });
            for_each_tag(*b2, [&](const std::string& t) {
                if (gemmi::glob_match(tag, t) && !gemmi::in_vector(t, matching_tags))
                  matching_tags.push_back(t);
            });
            for (const std::string& t : matching_tags)
              compare_tag_values(*b1, *b2, t);
          } else {
            compare_tag_values(*b1, *b2, tag);
          }
        }
        return 0;
      }
      Diff category_diff = make_diff(b1->get_mmcif_category_names(),
                                     b2->get_mmcif_category_names());
      for (DiffItem& cat : category_diff) {
        cif::Table t1 = b1->find_mmcif_category(cat.str);
        cif::Table t2 = b2->find_mmcif_category(cat.str);
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
    }
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

