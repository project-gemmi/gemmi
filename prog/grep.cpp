// Copyright 2017 Global Phasing Ltd.

// TODO: better handling of multi-line text values

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include "gemmi/dirwalk.hpp"
#include "gemmi/pdb_id.hpp"    // for is_pdb_code, expand_if_pdb_code
#include "gemmi/util.hpp"      // for replace_all
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#define GEMMI_PROG grep
#include "options.h"

using std::printf;
using std::fprintf;
namespace pegtl = tao::pegtl;
namespace cif = gemmi::cif;
namespace rules = gemmi::cif::rules;

namespace {

enum OptionIndex {
  FromFile=4, NamePattern, PdbDirSf, Recurse, MaxCount, OneBlock, And,
  Delim, WithFileName, NoBlockName, WithLineNumbers, WithTag,
  OnlyTags, Summarize, MatchingFiles, NonMatchingFiles, Count, Raw
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage: " EXE_NAME " [options] TAG FILE_OR_DIR_OR_PDBID[...]\n"
    "       " EXE_NAME " -f FILE [options] TAG\n"
    "Search for TAG in CIF files."
    "\nBy default, recursive directory search checks only *.cif(.gz) files."
    "\nTo change it, specify --name=* or --name=*.hkl."
    "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None,
    "  -h, --help  \tdisplay this help and exit" },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tdisplay version information and exit" },
  { Verbose, 0, "v", "verbose", Arg::None,
    "  -v, --verbose  \tprint additional messages to stderr" },
  { FromFile, 0, "f", "file", Arg::Required,
    "  -f, --file=FILE  \tobtain file (or PDB ID) list from FILE" },
  { NamePattern, 0, "", "name", Arg::Required,
    "  --name=PATTERN  \tfilename glob pattern used in recursive grep;"
    " by default, *.cif and *.cif.gz files are searched" },
  { PdbDirSf, 0, "S", "pdb-sf", Arg::None,
    "  -S, --pdb-sf \tif PDB ID is given, search structure factor file" },
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
  { OnlyTags, 0, "T", "only-tags", Arg::None,
    "  -T, --only-tags  \tprint only matching tags, not values" },
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

struct GrepParams {
  // options
  std::string search_tag;
  int max_count = 0;
  bool verbose = false;
  bool with_filename = false;
  bool with_blockname = true;
  bool with_line_numbers = false;
  bool with_tag = false;
  bool only_tags = false;
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
void process_match(const Input& in, GrepParams& par, int n) {
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
  if (par.with_tag) {
    const std::string& tag = n < 0 ? par.search_tag : par.multi_tags[n];
    if (par.only_tags) {
      printf("%s\n", tag.c_str());
      if (n == -1)
        par.match_column = -1;
      else
        par.multi_match_columns[n] = -1;
      return;
    } else if (par.delim.empty()) {
      printf("[%s] ", tag.c_str());
    } else {
      printf("%s%s", tag.c_str(), sep);
    }
  }
  std::string value = par.raw ? in.string() : cif::as_string(in.string());
  printf("%s\n", value.c_str());
  if (par.counters[0] == par.max_count)
    throw true;
}

// Escape delim (which normally is a a single character) with backslash.
std::string escape(const std::string& s, char delim) {
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

void process_multi_match(GrepParams& par) {
  if (par.multi_values.empty())
    return;
  if (par.print_count || par.only_filenames) {
    for (auto& mv : par.multi_values)
      mv.clear();
    return;
  }
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
    if (par.with_tag) {
      if (par.delim.empty())
        printf("[%s] ", par.multi_tags[0].c_str());
      else
        printf("%s%s", par.multi_tags[0].c_str(), sep);
    }
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
      break;
  }
  for (auto& mv : par.multi_values)
    mv.clear();
}

void print_count(const GrepParams& par) {
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
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    process_multi_match(p);
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
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    Search<rules::datablockname>::apply(in, p);
  }
};
template<> struct Search<rules::framename> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    p.block_name += " " + in.string();
  }
};
template<> struct Search<rules::endframe> {
  template<typename Input> static void apply(const Input&, GrepParams& p) {
    process_multi_match(p);
    p.block_name.erase(p.block_name.rfind(' '));
  }
};
template<> struct Search<rules::item_tag> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    if (!p.globbing) {
      if (p.search_tag.size() == in.size() && p.search_tag == in.string())
        p.match_value = 1;
    } else {
      if (gemmi::glob_match(p.search_tag, in.string())) {
        p.multi_tags.resize(1);
        p.multi_tags[0] = in.string();
        p.match_value = 1;
      }
    }
  }
};

template<> struct Search<rules::item_value> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    if (p.match_value) {
      p.match_value = 0;
      process_match(in, p, p.globbing ? 0 : -1);
      if (p.last_block && !p.globbing)
        throw true;
    }
  }
};
template<> struct Search<rules::str_loop> {
  template<typename Input> static void apply(const Input&, GrepParams& p) {
    p.table_width = 0;
    if (p.globbing) {
      p.multi_tags.clear();
      p.multi_match_columns.clear();
    }
  }
};
template<> struct Search<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    if (!p.globbing) {
      if (p.search_tag == in.string()) {
        p.match_column = p.table_width;
        p.column = 0;
      }
    } else {
      if (gemmi::glob_match(p.search_tag, in.string())) {
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
  template<typename Input> static void apply(const Input&, GrepParams& p) {
    if (p.match_column != -1) {
      p.match_column = -1;
      if (p.last_block && !p.globbing)
        throw true;
    }
  }
};
template<> struct Search<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    if (p.match_column == -1)
      return;
    if (!p.globbing) {
      if (p.column == p.match_column)
        process_match(in, p, -1);
    } else {
      for (int i = 0; i != (int) p.multi_match_columns.size(); ++i)
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

template<> struct MultiSearch<rules::item_tag> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
    const std::string s = in.string();
    for (int i = 0; i < static_cast<int>(p.multi_tags.size()); ++i)
      if (p.multi_tags[i] == s)
        p.match_value = i + 1;
  }
};
template<> struct MultiSearch<rules::item_value> {
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
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
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
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
  template<typename Input> static void apply(const Input&, GrepParams& p) {
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
  template<typename Input> static void apply(const Input& in, GrepParams& p) {
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
void run_parse(Input&& in, GrepParams& par) {
  if (par.multi_values.empty())
    pegtl::parse<rules::file, Search, cif::Errors>(in, par);
  else
    pegtl::parse<rules::file, MultiSearch, cif::Errors>(in, par);
}

void grep_file(const std::string& path, GrepParams& par, int& err_count) {
  if (par.verbose)
    fprintf(stderr, "Reading %s ...\n", path.c_str());
  par.path = path.c_str();
  par.block_name.clear();
  par.counters.clear();
  if (par.globbing)
    par.multi_tags.clear();
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
    } else if (gemmi::CharArray mem = input.uncompress_into_buffer()) {
      pegtl::memory_input<> in(mem.data(), mem.size(), path);
      run_parse(in, par);
    } else {
      GEMMI_CIF_FILE_INPUT(in, path);
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
  } else {
    process_multi_match(par);
  }
  par.total_count += par.counters[0];
  std::fflush(stdout);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  GrepParams params;
  if (p.options[MaxCount])
    params.max_count = std::strtol(p.options[MaxCount].arg, nullptr, 10);
  if (p.options[OneBlock])
    params.last_block = true;
  if (p.options[Verbose])
    params.verbose = true;
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
  if (p.options[OnlyTags]) {
    if (p.options[And]) {
      fprintf(stderr, "Options --only-tags and --and do not work together\n");
      return 2;
    }
    params.with_tag = true;
    params.only_tags = true;
  }
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
    gemmi::replace_all(params.delim, "\\t", "\t");
  }

  auto paths = p.paths_from_args_or_file(FromFile, 1);
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
      if (strchr(opt->arg, '?') || strchr(opt->arg, '*')) {
        fprintf(stderr, "Glob patterns are not supported together with -a.\n");
        return 2;
      }
      params.multi_tags.emplace_back(opt->arg);
    }
  } else {
    params.search_tag = tag;
    if (params.search_tag.find_first_of("?*") != std::string::npos)
      params.globbing = true;
  }

  size_t file_count = 0;
  int err_count = 0;
  try {
    char expand_type = p.options[PdbDirSf] ? 'S' : 'M';
    for (const std::string& path : paths) {
      if (path == "-") {
        grep_file(path, params, err_count);
        file_count++;
      } else if (p.options[FromFile] ? starts_with_pdb_code(path)
                                     : gemmi::is_pdb_code(path)) {
        std::string real_path = gemmi::expand_if_pdb_code(path.substr(0, 4), expand_type);
        params.last_block = true;  // PDB code implies -O
        grep_file(real_path, params, err_count);
        params.last_block = p.options[OneBlock];
        file_count++;
      } else {
        if (p.options[NamePattern]) {
          std::string pattern = p.options[NamePattern].arg;
          for (const std::string& file : gemmi::GlobWalk(path, pattern)) {
            grep_file(file, params, err_count);
            file_count++;
          }
        } else if (!p.options[Recurse] && (gemmi::giends_with(path, ".cif") ||
                                           gemmi::giends_with(path, ".mmcif"))) {
          // Avoid tinydir_file_open (used by CifWalk) when not necessary.
          // It was reported to fail on a Mac with files on network drive.
          // Probably reading the parent directory failed, no idea why.
          grep_file(path, params, err_count);
          file_count++;
        } else {
          for (const std::string& file : gemmi::CifWalk(path)) {
            grep_file(file, params, err_count);
            file_count++;
          }
        }
      }
    }
  } catch (std::runtime_error &e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return 2;
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
