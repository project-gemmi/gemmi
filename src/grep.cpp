// Copyright 2017 Global Phasing Ltd.

// TODO: better handling of multi-line text values

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/dirwalk.hpp"
#include "gemmi/util.hpp"  // for is_pdb_code
#include "input.h"         // for expand_pdb_code_to_path_or_fail
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#define EXE_NAME "gemmi-grep"
#include "options.h"

using std::printf;
using std::fprintf;
namespace pegtl = tao::pegtl;
namespace cif = gemmi::cif;
namespace rules = gemmi::cif::rules;


enum OptionIndex { FromFile=3, Recurse, MaxCount, OneBlock, And, Delim,
                   WithFileName, NoBlockName, WithLineNumbers, WithTag,
                   Summarize, MatchingFiles, NonMatchingFiles, Count, Raw };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage: " EXE_NAME " [options] TAG FILE_OR_DIR_OR_PDBID[...]\n"
    "       " EXE_NAME " -f FILE [options] TAG\n"
    "Search for TAG in CIF files."
    "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None,
    "  -h, --help  \tdisplay this help and exit" },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tdisplay version information and exit" },
  { FromFile, 0, "f", "file", Arg::Required,
    "  -f, --file=FILE  \tobtain file (or PDB ID) list from FILE" },
  { MaxCount, 0, "m", "max-count", Arg::Int,
    "  -m, --max-count=NUM  \tprint max NUM values per file" },
  { OneBlock, 0, "O", "one-block", Arg::None,
    "  -O, --one-block  \toptimize assuming one block per file" },
  { And, 0, "a", "and", Arg::Required,
    "  -a, --and=tag  \tAppend delimiter (default ';') and the tag value" },
  { Delim, 0, "d", "delimiter", Arg::Required,
    "  -d, --delimiter=DELIM  \tCSV-like output with specified delimiter" },
  { WithLineNumbers, 0, "n", "line-number", Arg::None,
    "  -n, --line-number  \tprint line number with output lines" },
  { WithFileName, 0, "H", "with-filename", Arg::None,
    "  -H, --with-filename  \tprint the file name for each match" },
  { NoBlockName, 0, "b", "no-blockname", Arg::None,
    "  -b, --no-blockname  \tsuppress the block name on output" },
  { WithTag, 0, "t", "with-tag", Arg::None,
    "  -t, --with-tag  \tprint the tag name for each match" },
  { MatchingFiles, 0, "l", "files-with-tag", Arg::None,
    "  -l, --files-with-tag  \tprint only names of files with the tag" },
  { NonMatchingFiles, 0, "L", "files-without-tag", Arg::None,
    "  -L, --files-without-tag  \tprint only names of files without the tag" },
  { Count, 0, "c", "count", Arg::None,
    "  -c, --count  \tprint only a count of values per block or file" },
  { Recurse, 0, "r", "recursive", Arg::None,
    "  -r, --recursive  \tignored (directories are always recursed)" },
  { Raw, 0, "w", "raw", Arg::None,
    "  -w, --raw  \tinclude '?', '.', and string quotes" },
  { Summarize, 0, "s", "summarize", Arg::None,
    "  -s, --summarize  \tdisplay joint statistics for all files" },
  { 0, 0, 0, 0, 0, 0 }
};

struct Parameters {
  // options
  std::string search_tag;
  int max_count = 0;
  bool with_filename = false;
  bool with_blockname = true;
  bool with_line_numbers = false;
  bool with_tag = false;
  bool summarize = false;
  bool only_filenames = false;
  bool inverse = false;  // for now it refers to only_filenames only
  bool print_count = false;
  bool raw = false;
  std::string delim;
  std::vector<std::string> multi_tags;
  bool globbing = false;
  // working parameters
  const char* path = "";
  std::string block_name;
  int match_value = 0;
  int match_column = -1;
  int table_width = 0;
  int column = 0;
  std::vector<int> counters;
  size_t total_count = 0;
  bool last_block = false;
  std::vector<int> multi_match_columns;
  std::vector<std::vector<std::string>> multi_values;
};

template<typename Input>
void process_match(const Input& in, Parameters& par, int n) {
  if (cif::is_null(in.string()) && !par.raw)
    return;
  ++par.counters[0];
  if (par.only_filenames)
    throw true;
  if (par.print_count)
    return;
  const char* sep = par.delim.empty() ? ":" : par.delim.c_str();
  if (par.with_filename)
    printf("%s%s", par.path, sep);
  if (par.with_blockname)
    printf("%s%s", par.block_name.c_str(), sep);
  if (par.with_line_numbers)
    printf("%zu%s", in.iterator().line, sep);
  if (par.with_tag)
    printf("[%s] ", (n < 0 ? par.search_tag : par.multi_tags[n]).c_str());
  std::string value = par.raw ? in.string() : cif::as_string(in.string());
  printf("%s\n", value.c_str());
  if (par.counters[0] == par.max_count)
    throw true;
}

// Escape delim (which normally is a a single character) with backslash.
static std::string escape(const std::string& s, char delim) {
  std::string r;
  for (char c : s) {
    if (c == '\n') {
      r += "\\n";
      continue;
    }
    if (c == '\\' || c == delim)
      r += '\\';
    r += c;
  }
  return r;
}

// linear-time glob matching: https://research.swtch.com/glob
bool glob_match(const std::string& pattern, const std::string& str) {
  size_t pat_next = 0;
  size_t str_next = std::string::npos;
  size_t pat_pos = 0;
  size_t str_pos = 0;
  while (pat_pos < pattern.size() || str_pos < str.size()) {
    if (pat_pos < pattern.size()) {
      char c = pattern[pat_pos];
      if (c == '*') {
        pat_next = pat_pos;
        str_next = str_pos + 1;
        pat_pos++;
        continue;
      } else if (str_pos < str.size() && (c == '?' || c == str[str_pos])) {
        pat_pos++;
        str_pos++;
        continue;
      }
    }
    if (str_next > str.size())
      return false;
    pat_pos = pat_next;
    str_pos = str_next;
  }
  return true;
}

static void process_multi_match(Parameters& par) {
  std::string need_escaping = "\n\\";
  if (par.delim.size() < 2)
    need_escaping += par.delim.empty() ? ';' : par.delim[0];
  for (size_t i = 0; i != par.multi_values[0].size(); ++i) {
    if (cif::is_null(par.multi_values[0][i]) && !par.raw)
      continue;
    const char* sep = par.delim.empty() ? ":" : par.delim.c_str();
    if (par.with_filename)
      printf("%s%s", par.path, sep);
    if (par.with_blockname)
      printf("%s%s", par.block_name.c_str(), sep);
    if (par.with_tag)
      printf("[%s] ", par.multi_tags[0].c_str());
    for (size_t j = 0; j != par.multi_values.size(); ++j) {
      if (j != 0)
        std::fputs(par.delim.empty() ? ";" : par.delim.c_str(), stdout);
      const auto& v = par.multi_values[j];
      if (!v.empty()) {
        const std::string& raw_str = v[i < v.size() ? i : 0];
        std::string s = par.raw ? raw_str : cif::as_string(raw_str);
        if (s.find_first_of(need_escaping) != std::string::npos)
          s = escape(s, need_escaping[2]);
        printf("%s", s.c_str());
      }
    }
    std::putc('\n', stdout);
    if (par.counters[0] == par.max_count)
      return;
  }
}

static void print_count(const Parameters& par) {
  const char* sep = par.delim.empty() ? ":" : par.delim.c_str();
  if (par.with_filename)
    printf("%s%s", par.path, sep);
  if (par.with_blockname)
    printf("%s%s", par.block_name.c_str(), sep);
  bool first = true;
  for (int c : par.counters) {
    if (!first)
      std::fputs(par.delim.empty() ? ";" : par.delim.c_str(), stdout);
    printf("%d", c);
    first = false;
  }
  std::putc('\n', stdout);
}


template<typename Rule> struct Search : pegtl::nothing<Rule> {};

template<> struct Search<rules::datablockname> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (!p.block_name.empty() && p.print_count && p.with_blockname) {
      print_count(p);
      p.total_count += p.counters[0];
      for (int& c : p.counters)
        c = 0;
    }
    p.block_name = in.string();
  }
};
template<> struct Search<rules::str_global> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    Search<rules::datablockname>::apply(in, p);
  }
};
template<> struct Search<rules::tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (!p.globbing) {
      if (p.search_tag == in.string())
        p.match_value = 1;
    } else {
      if (glob_match(p.search_tag, in.string())) {
        p.multi_tags.resize(1);
        p.multi_tags[0] = in.string();
        p.match_value = 1;
      }
    }
  }
};

template<> struct Search<rules::value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_value) {
      p.match_value = 0;
      process_match(in, p, p.globbing ? 0 : -1);
      if (p.last_block)
        throw true;
    }
  }
};
template<> struct Search<rules::str_loop> {
  template<typename Input> static void apply(const Input&, Parameters& p) {
    p.table_width = 0;
    p.multi_tags.clear();
    p.multi_match_columns.clear();
  }
};
template<> struct Search<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (!p.globbing) {
      if (p.search_tag == in.string()) {
        p.match_column = p.table_width;
        p.column = 0;
      }
    } else {
      if (glob_match(p.search_tag, in.string())) {
        p.multi_tags.emplace_back(in.string());
        p.multi_match_columns.emplace_back(p.table_width);
        p.match_column = 0;
        p.column = 0;
      }
    }
    p.table_width++;
  }
};

template<> struct Search<rules::loop_end> {
  template<typename Input> static void apply(const Input&, Parameters& p) {
    if (p.match_column != -1) {
      p.match_column = -1;
      if (p.last_block)
        throw true;
    }
  }
};
template<> struct Search<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_column == -1)
      return;
    if (!p.globbing) {
      if (p.column == p.match_column)
        process_match(in, p, -1);
    } else {
      for (size_t i = 0; i != p.multi_match_columns.size(); ++i)
        if (p.column == p.multi_match_columns[i])
          process_match(in, p, i);
    }
    p.column++;
    if (p.column == p.table_width)
      p.column = 0;
  }
};


template<typename T> bool any_empty(const std::vector<T>& v) {
  for (const T& a : v)
    if (a.empty())
      return true;
  return false;
}

template<typename Rule> struct MultiSearch : Search<Rule> {};

template<> struct MultiSearch<rules::tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    const std::string s = in.string();
    for (int i = 0; i < static_cast<int>(p.multi_tags.size()); ++i)
      if (p.multi_tags[i] == s)
        p.match_value = i + 1;
  }
};
template<> struct MultiSearch<rules::value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_value) {
      if (p.raw || !cif::is_null(in.string()))
        ++p.counters[p.match_value - 1];
      p.multi_values[p.match_value - 1].emplace_back(in.string());
      p.match_value = 0;
      if (p.last_block && !any_empty(p.multi_values))
        throw true;
    }
  }
};
template<> struct MultiSearch<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    const std::string s = in.string();
    for (size_t i = 0; i != p.multi_tags.size(); ++i)
      if (p.multi_tags[i] == s) {
        p.multi_match_columns[i] = p.table_width;
        p.match_column = 0;
        p.column = 0;
      }
    p.table_width++;
  }
};
template<> struct MultiSearch<rules::loop_end> {
  template<typename Input> static void apply(const Input&, Parameters& p) {
    if (p.match_column == 0) {
      p.match_column = -1;
      for (int& c : p.multi_match_columns)
        c = -1;
      if (p.last_block && !any_empty(p.multi_values))
        throw true;
    }
  }
};
template<> struct MultiSearch<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_column == 0) {
      for (size_t i = 0; i != p.multi_values.size(); ++i)
        if (p.column == p.multi_match_columns[i]) {
          if (p.raw || !cif::is_null(in.string()))
            ++p.counters[i];
          // if it's not the loop with the main tag, we need only one value
          if (p.multi_match_columns[0] != -1 || p.multi_values[i].empty())
            p.multi_values[i].emplace_back(in.string());
        }
      p.column++;
      if (p.column == p.table_width)
        p.column = 0;
    }
  }
};

template<typename Input>
void run_parse(Input&& in, Parameters& par) {
  if (par.multi_values.empty())
    pegtl::parse<rules::file, Search, cif::Errors>(in, par);
  else
    pegtl::parse<rules::file, MultiSearch, cif::Errors>(in, par);
}

static
void grep_file(const std::string& path, Parameters& par, int& err_count) {
  par.path = path.c_str();
  par.block_name.clear();
  par.counters.clear();
  size_t n_multi = par.multi_tags.size();
  par.counters.resize(n_multi == 0 ? 1 : n_multi, 0);
  par.match_column = -1;
  par.match_value = 0;
  par.multi_match_columns.clear();
  par.multi_match_columns.resize(n_multi, -1);
  par.multi_values.clear();
  par.multi_values.resize(n_multi);
  try {
    gemmi::MaybeGzipped input(path);
    if (input.is_stdin()) {
      pegtl::cstream_input<> in(stdin, 16*1024, "stdin");
      run_parse(in, par);
    } else if (input.is_compressed()) {
      std::unique_ptr<char[]> mem = input.memory();
      pegtl::memory_input<> in(mem.get(), input.mem_size(), path);
      run_parse(in, par);
    } else {
      pegtl::file_input<> in(path);
      run_parse(in, par);
    }
  } catch (bool) {
    // ok, "throw true" is used as goto
  } catch (std::runtime_error& e) {
    std::fflush(stdout);
    fprintf(stderr, "Error when parsing %s:\n\t%s\n", path.c_str(), e.what());
    err_count++;
    return;
  }
  if (par.print_count) {
    print_count(par);
  } else if (par.only_filenames) {
    if (par.inverse == (par.counters[0] == 0))
      printf("%s\n", par.path);
  } else if (!par.multi_values.empty() && !par.multi_values[0].empty()) {
    process_multi_match(par);
  }
  par.total_count += par.counters[0];
  std::fflush(stdout);
}


static void replace_all(std::string &s,
                        const std::string &old, const std::string &new_) {
  std::string::size_type pos = 0;
  while ((pos = s.find(old, pos)) != std::string::npos) {
    s.replace(pos, old.size(), new_);
    pos += new_.size();
  }
}

int main(int argc, char **argv) {
  OptParser p;
  p.simple_parse(argc, argv, Usage);
  if (p.options[FromFile] ? p.nonOptionsCount() != 1
                          : p.nonOptionsCount() < 2) {
    option::printUsage(fwrite, stderr, Usage);
    return 2;
  }

  Parameters params;
  if (p.options[MaxCount])
    params.max_count = std::strtol(p.options[MaxCount].arg, nullptr, 10);
  if (p.options[OneBlock])
    params.last_block = true;
  if (p.options[WithFileName])
    params.with_filename = true;
  if (p.options[NoBlockName])
    params.with_blockname = false;
  if (p.options[WithLineNumbers]) {
    if (p.options[And]) {
      fprintf(stderr, "Options --line-number and --and do not work together\n");
      return 2;
    }
    params.with_line_numbers = true;
  }
  if (p.options[WithTag])
    params.with_tag = true;
  if (p.options[Summarize])
    params.summarize = true;
  if (p.options[MatchingFiles])
    params.only_filenames = true;
  if (p.options[NonMatchingFiles]) {
    params.only_filenames = true;
    params.inverse = true;
  }
  if (p.options[Count])
    params.print_count = true;
  if (p.options[Raw])
    params.raw = true;
  if (p.options[Delim]) {
    params.delim = p.options[Delim].arg;
    replace_all(params.delim, "\\t", "\t");
  }

  const char* tag = p.nonOption(0);
  if (tag[0] != '_') {
    fprintf(stderr, "CIF tag must start with \"_\": %s\n", tag);
    return 2;
  }
  if (p.options[And]) {
    params.multi_tags.emplace_back(tag);
    for (const option::Option* opt = p.options[And]; opt; opt = opt->next()) {
      if (opt->arg[0] != '_') {
        fprintf(stderr, "CIF tags start with _; not a tag: %s\n", opt->arg);
        return 2;
      }
      params.multi_tags.emplace_back(opt->arg);
    }
    for (const std::string& t : params.multi_tags)
      if (t.find_first_of("?*") != std::string::npos) {
        fprintf(stderr, "Glob patterns are not supported together with -a.\n");
        return 2;
      }
  } else {
    params.search_tag = tag;
    if (params.search_tag.find_first_of("?*") != std::string::npos)
      params.globbing = true;
  }

  std::vector<std::string> paths;
  if (p.options[FromFile]) {
    std::FILE *f = std::fopen(p.options[FromFile].arg, "r");
    if (!f) {
      std::perror(p.options[FromFile].arg);
      return 2;
    }
    char buf[512];
    while (std::fgets(buf, 512, f)) {
      std::string s = gemmi::trim_str(buf);
      if (s.length() >= 4 && std::strchr(" \t\r\n:,;|", s[4]) &&
          gemmi::is_pdb_code(s.substr(0, 4)))
        s.resize(4);
      if (!s.empty())
        paths.emplace_back(s);
    }
    std::fclose(f);
  } else {
    for (int i = 1; i < p.nonOptionsCount(); ++i)
      paths.emplace_back(p.nonOption(i));
  }

  size_t file_count = 0;
  int err_count = 0;
  for (const std::string& path : paths) {
    if (path == "-") {
      grep_file(path, params, err_count);
      file_count++;
    } else if (gemmi::is_pdb_code(path)) {
      std::string real_path = expand_pdb_code_to_path_or_fail(path);
      params.last_block = true;  // PDB code implies -O
      grep_file(real_path, params, err_count);
      params.last_block = p.options[OneBlock];
      file_count++;
    } else {
      try {
        for (const char* file : gemmi::CifWalk(path)) {
          grep_file(file, params, err_count);
          file_count++;
        }
      } catch (std::runtime_error &e) {
        fprintf(stderr, "Error: %s\n", e.what());
        return 2;
      }
    }
  }
  if (p.options[Summarize]) {
    printf("Total count in %zu files: %zu\n", file_count, params.total_count);
    if (err_count > 0)
      printf("Errors encountered when reading %d files.\n", err_count);
  }
  if (err_count > 0)
    return 2;
  return params.total_count != 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
