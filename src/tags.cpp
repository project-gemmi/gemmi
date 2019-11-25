// Copyright 2018 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/numb.hpp"     // for as_number
#include "gemmi/dirwalk.hpp"  // for CifWalk
#include "gemmi/util.hpp"     // for replace_all
#include "gemmi/gz.hpp"       // for MaybeGzipped
#include <cassert>
#include <cmath>    // for INFINITY
#include <cstdio>   // for printf, fprintf
#include <utility>  // for pair
#include <string>
#include <map>
#include <vector>

#define GEMMI_PROG tags
#include "options.h"

namespace pegtl = tao::pegtl;
namespace cif = gemmi::cif;
namespace rules = gemmi::cif::rules;

enum OptionIndex { CountFiles=3, Full=4 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] FILE_OR_DIR[...]"
    "\nList CIF tags with counts of blocks and values."},
  CommonUsage[Help],
  CommonUsage[Version],
  { CountFiles, 0, "", "count-files", Arg::None,
    "  --count-files  \tCount files instead of blocks." },
  { Full, 0, "", "full", Arg::None,
    "  --full  \tGather data for project-gemmi.github.io/pdb-stats/tags.html" },
  { 0, 0, 0, 0, 0, 0 }
};


constexpr int ENUM_LIMIT = 20;

struct TagStats {
  int file_count = 0;
  int block_count = 0;
  size_t total_count = 0;
  int min_count = INT_MAX;
  int max_count = 1;
  bool in_this_file = false;
  std::map<std::string, int> values;
  size_t text_count = 0;
  size_t multi_word_count = 0;
  size_t single_word_count = 0;
  size_t number_count = 0;
  double min_number = +INFINITY;
  double max_number = -INFINITY;

  void add_value(const std::string& raw) {
    assert(!cif::is_null(raw));
    // compare with cif::as_string()
    if (raw[0] == ';' && raw.size() > 2 && *(raw.end() - 2) == '\n') {
      ++text_count;
    } else if (raw[0] == '"' || raw[0] == '\'') {
      if (raw.find(' ') == std::string::npos)
        ++single_word_count;
      else
        ++multi_word_count;
    } else {
      double d = cif::as_number(raw);
      if (std::isnan(d)) {
        ++single_word_count;
      } else {
        ++number_count;
        if (d < min_number)
          min_number = d;
        if (d > max_number)
          max_number = d;
      }
    }
    if (values.size() <= ENUM_LIMIT)
      values[cif::as_string(raw)]++;
  }
};

struct LoopInfo {
  std::string tag;
  int counter;
  TagStats* stats_ptr;
};

struct Context {
  std::map<std::string, TagStats> stats;
  int total_blocks = 0;
  int total_files = 0;
  std::string tag;
  std::vector<LoopInfo> loop_info;
  size_t column = 0;
  bool per_block = true;
  bool full_output = false;
};

template<typename Rule> struct Counter : pegtl::nothing<Rule> {};

template<> struct Counter<rules::datablockname> {
  template<typename Input> static void apply(const Input&, Context& ctx) {
    ctx.total_blocks++;
  }
};
template<> struct Counter<rules::str_global> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    Counter<rules::datablockname>::apply(in, ctx);
  }
};

// tag-value pairs
template<> struct Counter<rules::tag> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    ctx.tag = in.string();
  }
};
template<> struct Counter<rules::value> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    if (!cif::is_null(in.string())) {
      TagStats& st = ctx.stats[ctx.tag];
      st.block_count++;
      st.total_count++;
      st.min_count = 1;
      st.in_this_file = true;
      if (ctx.full_output)
        st.add_value(in.string());
    }
    ctx.tag.clear();
  }
};

// loops
template<> struct Counter<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    // map::emplace() does not invalidate references
    auto iter = ctx.stats.emplace(in.string(), TagStats()).first;
    ctx.loop_info.push_back(LoopInfo{iter->first, 0, &iter->second});
  }
};
template<> struct Counter<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    std::string raw_value = in.string();
    if (!cif::is_null(raw_value)) {
      LoopInfo& loop_info = ctx.loop_info[ctx.column];
      loop_info.counter++;
      if (ctx.full_output)
        loop_info.stats_ptr->add_value(raw_value);
    }
    ctx.column++;
    if (ctx.column == ctx.loop_info.size())
      ctx.column = 0;
  }
};
template<> struct Counter<rules::loop_end> {
  template<typename Input> static void apply(const Input&, Context& ctx) {
    for (auto& info : ctx.loop_info) {
      TagStats& st = ctx.stats[info.tag];
      int n = info.counter;
      if (n != 0) {
        st.block_count++;
        st.total_count += n;
        st.max_count = std::max(st.max_count, n);
        st.min_count = std::min(st.min_count, n);
        st.in_this_file = true;
      }
    }
    ctx.column = 0;
    ctx.loop_info.clear();
  }
};

static void print_data_for_html(const Context& ctx) {
  //std::printf("tag\tfiles\tnmin\tnavg\tnmax\n");
  for (auto& item : ctx.stats) {
    const TagStats& st = item.second;
    if (st.block_count == 0)
      continue;
    double navg = double(st.total_count) / st.block_count;
    double pc;
    if (ctx.per_block)
      pc = 100.0 * st.block_count / ctx.total_blocks;
    else
      pc = 100.0 * st.file_count / ctx.total_files;
    std::printf("%s\t%.3f\t%d\t%.2f\t%d",
                item.first.c_str(), pc, st.min_count, navg, st.max_count);
    if (st.values.size() <= ENUM_LIMIT && st.text_count == 0) {
      // sort by count
      using Pair = std::pair<std::string, int>;
      std::vector<Pair> pairs(st.values.begin(), st.values.end());
      std::sort(pairs.begin(), pairs.end(),
                [](Pair& a, Pair& b) { return a.second > b.second; });
      for (auto& value_count : pairs) {
        gemmi::replace_all(value_count.first, "&", "&amp;");
        gemmi::replace_all(value_count.first, "<", "&lt;");
        std::printf("\t%d %s", value_count.second, value_count.first.c_str());
      }
    } else {
      if (st.text_count != 0)
        std::printf("\t%zu {text}", st.text_count);
      if (st.multi_word_count != 0)
        std::printf("\t%zu {line}", st.multi_word_count);
      if (st.single_word_count != 0)
        std::printf("\t%zu {word}", st.single_word_count);
      if (st.number_count != 0)
        std::printf("\t%zu {%g - %g}",
                    st.number_count, st.min_number, st.max_number);
    }
    std::printf("\n");
  }
}

static void print_tag_list(const Context& ctx) {
  std::printf("tag\t%s-count\tvalue-count\n", ctx.per_block ? "block" : "file");
  for (auto& item : ctx.stats) {
    const TagStats& st = item.second;
    if (st.block_count == 0)
      continue;
    int groups = ctx.per_block ? st.block_count : st.file_count;
    std::printf("%s\t%d\t%zu\n", item.first.c_str(), groups, st.total_count);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  Context ctx;
  ctx.full_output = p.options[Full];
  ctx.per_block = !p.options[CountFiles];
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    for (const char* path : gemmi::CifWalk(p.nonOption(i))) {
      try {
        gemmi::MaybeGzipped input(path);
        if (input.is_stdin()) {
          pegtl::cstream_input<> in(stdin, 16*1024, "stdin");
          pegtl::parse<rules::file, Counter, cif::Errors>(in, ctx);
        } else if (input.is_compressed()) {
          std::unique_ptr<char[]> mem = input.memory();
          pegtl::memory_input<> in(mem.get(), input.memory_size(), path);
          pegtl::parse<rules::file, Counter, cif::Errors>(in, ctx);
        } else {
          pegtl::file_input<> in(path);
          pegtl::parse<rules::file, Counter, cif::Errors>(in, ctx);
        }
        for (auto& item : ctx.stats)
          if (item.second.in_this_file) {
            item.second.file_count++;
            item.second.in_this_file = false;
          }
        ctx.total_files++;
      } catch (std::runtime_error &e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        // Some files can be incorrect, continue despite of it.
        ctx.tag.clear();
        ctx.column = 0;
        ctx.loop_info.clear();
      }
    }
  }
  if (p.options[Full])
    print_data_for_html(ctx);
  else
    print_tag_list(ctx);
  std::fprintf(stderr, "Tag count: %zu\n", ctx.stats.size());
  std::fprintf(stderr, "Block count: %d\n", ctx.total_blocks);
  std::fprintf(stderr, "File count: %d\n", ctx.total_files);
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
