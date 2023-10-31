// Copyright 2018 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/numb.hpp"     // for as_number
#include "gemmi/dirwalk.hpp"  // for CifWalk
#include "gemmi/util.hpp"     // for replace_all
#include "gemmi/gz.hpp"       // for MaybeGzipped
#include <cassert>
#include <cmath>              // for INFINITY
#include <cstdio>
#include <algorithm>          // for sort
#include <fstream>
#include <utility>            // for pair
#include <string>
#include <map>
#include <vector>

#define GEMMI_PROG tags
#include "options.h"

namespace pegtl = tao::pegtl;
namespace cif = gemmi::cif;
namespace rules = gemmi::cif::rules;

namespace {

enum OptionIndex { CountFiles=4, Glob, Full, EntriesIdx, Sf };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] FILE_OR_DIR[...]"
    "\nList CIF tags with counts of blocks and values."},
  CommonUsage[Help],
  CommonUsage[Version],
  { CountFiles, 0, "", "count-files", Arg::None,
    "  --count-files  \tCount files instead of blocks." },
  { Glob, 0, "", "glob", Arg::Required,
    "  --glob=GLOB  \tProcess files matching glob pattern.\n" },
  { NoOp, 0, "", "", Arg::None,
    "Options for making https://project-gemmi.github.io/pdb-stats/tags.html" },
  { Full, 0, "", "full", Arg::None, "  --full  \tGather data for tags.html" },
  { EntriesIdx, 0, "", "entries-idx", Arg::Required,
    "  --entries-idx  \tUse entries.idx to read more recent entries first." },
  { Sf, 0, "", "sf", Arg::None,
    "  --sf  \t(for use with --entries-idx) Read SF mmCIF files." },
  { 0, 0, 0, 0, 0, 0 }
};


constexpr int ENUM_LIMIT = 20;

struct TagStats {
  int file_count = 0;
  int block_count = 0;
  unsigned long total_count = 0;
  int min_count = INT_MAX;
  int max_count = 1;
  bool in_this_file = false;
  struct CountAndExample {
    int count = 0;
    std::string example;
    void add(const std::string& block_name) {
      if (count == 0)
        example = block_name;
      ++count;
    }
  };
  std::map<std::string, CountAndExample> values;
  CountAndExample text;
  CountAndExample multi_word;
  CountAndExample single_word;
  CountAndExample number;
  double min_number = +INFINITY;
  double max_number = -INFINITY;

  void add_value(const std::string& raw, const std::string& block_name) {
    assert(!cif::is_null(raw));
    // compare with cif::as_string()
    if (raw[0] == ';' && raw.size() > 2 && *(raw.end() - 2) == '\n') {
      text.add(block_name);
    } else if (raw[0] == '"' || raw[0] == '\'') {
      if (raw.find(' ') == std::string::npos)
        single_word.add(block_name);
      else
        multi_word.add(block_name);
    } else {
      double d = cif::as_number(raw);
      if (std::isnan(d)) {
        single_word.add(block_name);
      } else {
        number.add(block_name);
        if (d < min_number)
          min_number = d;
        if (d > max_number)
          max_number = d;
      }
    }
    if (values.size() <= ENUM_LIMIT)
      values[cif::as_string(raw)].add(block_name);
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
  std::string block_name;
  std::string tag;
  std::vector<LoopInfo> loop_info;
  size_t column = 0;
  bool per_block = true;
  bool full_output = false;
};

template<typename Rule> struct Counter : pegtl::nothing<Rule> {};

template<> struct Counter<rules::datablockname> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    ctx.block_name = in.string();
    ctx.total_blocks++;
  }
};
template<> struct Counter<rules::str_global> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    ctx.block_name = "[global]";
    Counter<rules::datablockname>::apply(in, ctx);
  }
};

// tag-value pairs
template<> struct Counter<rules::item_tag> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    ctx.tag = in.string();
  }
};
template<> struct Counter<rules::item_value> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    if (!cif::is_null(in.string())) {
      TagStats& st = ctx.stats[ctx.tag];
      st.block_count++;
      st.total_count++;
      st.min_count = 1;
      st.in_this_file = true;
      if (ctx.full_output)
        st.add_value(in.string(), ctx.block_name);
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
        loop_info.stats_ptr->add_value(raw_value, ctx.block_name);
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

void print_value_count_for_tsv(const char* item, const TagStats::CountAndExample &val) {
  std::printf("\t%d %s %s", val.count, val.example.c_str(), item);
}

void print_data_for_html(const Context& ctx) {
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
    if (st.values.size() <= ENUM_LIMIT && st.text.count == 0) {
      // sort by count
      using Pair = std::pair<std::string, TagStats::CountAndExample>;
      std::vector<Pair> pairs(st.values.begin(), st.values.end());
      std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
          return a.second.count > b.second.count;
      });
      for (auto& pair : pairs) {
        gemmi::replace_all(pair.first, "&", "&amp;");
        gemmi::replace_all(pair.first, "<", "&lt;");
        print_value_count_for_tsv(pair.first.c_str(), pair.second);
      }
    } else {
      if (st.text.count != 0)
        print_value_count_for_tsv("{text}", st.text);
      if (st.multi_word.count != 0)
        print_value_count_for_tsv("{line}", st.multi_word);
      if (st.single_word.count != 0)
        print_value_count_for_tsv("{word}", st.single_word);
      if (st.number.count != 0) {
        using namespace std;
        char buf[64] = {0};
        snprintf(buf, 63, "{%g - %g}", st.min_number, st.max_number);
        print_value_count_for_tsv(buf, st.number);
      }
    }
    std::printf("\n");
  }
}

void print_tag_list(const Context& ctx) {
  std::printf("tag\t%s-count\tvalue-count\n", ctx.per_block ? "block" : "file");
  for (auto& item : ctx.stats) {
    const TagStats& st = item.second;
    if (st.block_count == 0)
      continue;
    int groups = ctx.per_block ? st.block_count : st.file_count;
    std::printf("%s\t%d\t%lu\n", item.first.c_str(), groups, st.total_count);
  }
}

void process(Context& ctx, const std::string& path) {
  try {
    gemmi::MaybeGzipped input(path);
    if (input.is_stdin()) {
      pegtl::cstream_input<> in(stdin, 16*1024, "stdin");
      pegtl::parse<rules::file, Counter, cif::Errors>(in, ctx);
    } else if (gemmi::CharArray mem = input.uncompress_into_buffer()) {
      pegtl::memory_input<> in(mem.data(), mem.size(), path);
      pegtl::parse<rules::file, Counter, cif::Errors>(in, ctx);
    } else {
      GEMMI_CIF_FILE_INPUT(in, path);
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

bool file_exists(const std::string& path) {
  if (FILE *file = std::fopen(path.c_str(), "rb")) {
    fclose(file);
    return true;
  }
  return false;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  Context ctx;
  ctx.full_output = p.options[Full];
  ctx.per_block = !p.options[CountFiles];
  if (p.options[EntriesIdx]) {
    if (p.nonOptionsCount() != 1) {
      std::fprintf(stderr, "Expected one argument with --entries-idx\n");
      return 1;
    }
    // The entries.idx is read first so we can sort pdb entries
    // in reverse chronological order.
    std::vector<std::string> pdb_dates;
    {
      std::ifstream entries(p.options[EntriesIdx].arg);
      std::string line;
      while (std::getline(entries, line)) {
        std::vector<std::string> tokens = gemmi::split_str(line, '\t');
        if (tokens.size() < 3 || tokens[2].size() < 8)
          continue;
        const std::string& us_date = tokens[2];
        std::string month = us_date.substr(0, 2);
        std::string day = us_date.substr(3, 2);
        std::string year = us_date.substr(6, 2);
        std::string date = (year[0] > '5' ? "19" : "20") + year + month + day;
        pdb_dates.push_back(date + tokens[0]);
      }
    }
    std::sort(pdb_dates.begin(), pdb_dates.end(),
              [](const std::string& a, const std::string& b) {return a > b; });
    std::string top_dir = p.nonOption(0);
    if (!top_dir.empty() && top_dir[top_dir.size() - 1] != '/')
      top_dir += '/';
    for (const std::string& str : pdb_dates) {
      std::string lc = gemmi::to_lower(str.substr(8));
      if (p.options[Sf]) {
        std::string path = top_dir + lc.substr(1, 2) + "/r" + lc + "sf.ent.gz";
        if (file_exists(path))
            process(ctx, path);
      } else {
        process(ctx, top_dir + lc.substr(1, 2) + "/" + lc + ".cif.gz");
      }
    }
  } else {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      try { // DirWalk can throw
        if (p.options[Glob])
          for (const std::string& path : gemmi::GlobWalk(p.nonOption(i), p.options[Glob].arg))
            process(ctx, path);
        else
          for (const std::string& path : gemmi::CifWalk(p.nonOption(i)))
            process(ctx, path);
      } catch (std::runtime_error &e) {
        std::fprintf(stderr, "Error. %s.\n", e.what());
        return 1;
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
