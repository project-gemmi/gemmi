// Copyright 2017 Global Phasing Ltd.

// TODO: better handling of multi-line text values

#include "gemmi/cif.hpp"
#include "gemmi/gz.hpp"
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <tinydir.h>
#include "input.h"  // for is_pdb_code, mmcif_subpath

#define EXE_NAME "gemmi-grep"
#include "options.h"

using std::printf;
using std::fprintf;
namespace pegtl = tao::pegtl;
namespace cif = gemmi::cif;
namespace rules = gemmi::cif::rules;


enum OptionIndex { FromFile=3, Recurse, MaxCount, OneBlock, And,
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
    "  -a, --and=tag  \tAppend separator ';' and the value of this tag" },
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
  std::vector<std::string> multi_tags;
  // working parameters
  const char* path = "";
  std::string block_name;
  int match_value = 0;
  int match_column = -1;
  int table_width = 0;
  int column = 0;
  int counter = 0;
  size_t total_count = 0;
  bool last_block = false;
  std::vector<int> multi_match_columns;
  std::vector<std::vector<std::string>> multi_values;
};

template<typename Input>
void process_match(const Input& in, Parameters& par) {
  if (cif::is_null(in.string()) && !par.raw)
    return;
  ++par.counter;
  if (par.only_filenames)
    throw true;
  if (par.print_count)
    return;
  if (par.with_filename)
    printf("%s:", par.path);
  if (par.with_blockname)
    printf("%s:", par.block_name.c_str());
  if (par.with_line_numbers)
    printf("%zu:", in.iterator().line);
  if (par.with_tag)
    printf("[%s] ", par.search_tag.c_str());
  printf("%s\n", (par.raw ? in.string() : cif::as_string(in.string())).c_str());
  if (par.counter == par.max_count)
    throw true;
}

void process_multi_match(Parameters& par) {
  for (size_t i = 0; i != par.multi_values[0].size(); ++i) {
    if (cif::is_null(par.multi_values[0][i]) && !par.raw)
      continue;
    ++par.counter;
    if (par.only_filenames)
      return;
    if (par.print_count)
      continue;
    if (par.with_filename)
      printf("%s:", par.path);
    if (par.with_blockname)
      printf("%s:", par.block_name.c_str());
    if (par.with_tag)
      printf("[%s] ", par.multi_tags[0].c_str());
    for (size_t j = 0; j != par.multi_values.size(); ++j) {
      if (j != 0)
        std::putc(';', stdout);
      const auto& v = par.multi_values[j];
      if (!v.empty()) {
        const std::string& s = v[i < v.size() ? i : 0];
        printf("%s", (par.raw ? s : cif::as_string(s)).c_str());
      }
    }
    std::putc('\n', stdout);
    if (par.counter == par.max_count)
      return;
  }
}

static void print_count(const Parameters& par) {
  if (par.with_filename)
    printf("%s:", par.path);
  if (par.with_blockname)
    printf("%s:", par.block_name.c_str());
  printf("%d\n", par.counter);
}


template<typename Rule> struct Search : pegtl::nothing<Rule> {};

template<> struct Search<rules::datablockname> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (!p.block_name.empty() && p.print_count && p.with_blockname) {
      print_count(p);
      p.total_count += p.counter;
      p.counter = 0;
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
    if (p.search_tag == in.string())
      p.match_value = 1;
  }
};
template<> struct Search<rules::value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_value) {
      p.match_value = 0;
      process_match(in, p);
      if (p.last_block)
        throw true;
    }
  }
};
template<> struct Search<rules::str_loop> {
  template<typename Input> static void apply(const Input&, Parameters& p) {
    p.table_width = 0;
  }
};
template<> struct Search<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.search_tag == in.string()) {
      p.match_column = p.table_width;
      p.column = 0;
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
    if (p.match_column != -1) {
      if (p.column == p.match_column)
        process_match(in, p);
      p.column++;
      if (p.column == p.table_width)
        p.column = 0;
    }
  }
};


template<typename T> bool any_empty(const std::vector<T>& v) {
  for (const T& a : v)
    if (a.empty())
      return false;
  return true;
}
template<typename Rule> struct MultiSearch : Search<Rule> {};

template<> struct MultiSearch<rules::tag> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    const std::string s = in.string();
    for (size_t i = 0; i != p.multi_tags.size(); ++i)
      if (p.multi_tags[i] == s)
        p.match_value = i + 1;
  }
};
template<> struct MultiSearch<rules::value> {
  template<typename Input> static void apply(const Input& in, Parameters& p) {
    if (p.match_value) {
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
        if (p.column == p.multi_match_columns[i]
            // if it's not the loop with the main tag, we need only one value
            && (p.multi_match_columns[0] != -1 || p.multi_values[i].empty()))
          p.multi_values[i].emplace_back(in.string());
      p.column++;
      if (p.column == p.table_width)
        p.column = 0;
    }
  }
};

template<typename Input>
void run_parse(Input&& in, Parameters& par) {
  if (par.multi_tags.empty())
    pegtl::parse<rules::file, Search, cif::Errors>(in, par);
  else
    pegtl::parse<rules::file, MultiSearch, cif::Errors>(in, par);
}

static
void grep_file(const std::string& path, Parameters& par) {
  par.path = path.c_str();
  par.counter = 0;
  par.match_column = -1;
  par.match_value = 0;
  par.multi_match_columns.clear();
  par.multi_match_columns.resize(par.multi_tags.size(), -1);
  par.multi_values.clear();
  par.multi_values.resize(par.multi_tags.size());
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
  } catch (bool) {}
  if (!par.multi_values.empty() && !par.multi_values[0].empty())
    process_multi_match(par);
  par.total_count += par.counter;
  if (par.print_count) {
    print_count(par);
  } else if (par.only_filenames && par.inverse == (par.counter == 0)) {
    printf("%s\n", par.path);
  }
  std::fflush(stdout);
}


class DirWalker {
public:
  explicit DirWalker(const char* path) {
    if (tinydir_file_open(&top_, path) == -1) {
      //std::perror(nullptr);
      throw std::runtime_error("Cannot open file or directory: " +
                               std::string(path));
    }
  }
  ~DirWalker() {
    for (auto& d : dirs_)
      tinydir_close(&d.second);
  }
  void push_dir(int cur_pos, const char* path) {
    dirs_.emplace_back();
    dirs_.back().first = cur_pos;
    if (tinydir_open_sorted(&dirs_.back().second, path) == -1)
      throw std::runtime_error("Cannot open directory: " + std::string(path));
  }
  int pop_dir() {
    assert(!dirs_.empty());
    int old_pos = dirs_.back().first;
    tinydir_close(&dirs_.back().second);
    dirs_.pop_back();
    return old_pos;
  }

  struct Iter {
    DirWalker& walker;
    size_t cur;

    const tinydir_dir& get_dir() const { return walker.dirs_.back().second; }

    const tinydir_file& operator*() const {
      if (walker.dirs_.empty())
        return walker.top_;
      assert(cur < get_dir().n_files);
      return get_dir()._files[cur];
    }

    bool is_special(const char* name) const {
      return strcmp(name, ".") == 0 || strcmp(name, "..") == 0;
    }

    void operator++() { // depth first
      const tinydir_file& tf = **this;
      if (tf.is_dir) {
        walker.push_dir(cur, tf.path);
        cur = 0;
      } else {
        cur++;
      }
      while (!walker.dirs_.empty()) {
        if (cur == get_dir().n_files)
          cur = walker.pop_dir() + 1;
        else if (is_special(get_dir()._files[cur].name))
          cur++;
        else
          break;
      }
    }

    bool operator!=(const Iter& o) const {
      return !(walker.dirs_.empty() && cur == o.cur);
    }
  };

  Iter begin() { return Iter{*this, 0}; }
  Iter end() { return Iter{*this, 1}; }
  bool is_file() { return !top_.is_dir; }

private:
  friend struct Iter;
  tinydir_file top_;
  std::vector<std::pair<int, tinydir_dir>> dirs_;
};

static bool is_cif_file(const tinydir_file& f) {
  return !f.is_dir && (gemmi::iends_with(f.path, ".cif") ||
                       gemmi::iends_with(f.path, ".cif.gz"));
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

  const char* tag = p.nonOption(0);
  if (tag[0] != '_') {
    fprintf(stderr, "CIF tags start with _; not a tag: %s\n", tag);
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
  } else {
    params.search_tag = tag;
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
          is_pdb_code(s.substr(0, 4)))
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
  for (const std::string& path : paths) {
    try {
      if (path == "-") {
        grep_file(path, params);
        file_count++;
      } else if (is_pdb_code(path)) {
        if (const char* pdb_dir = getenv("PDB_DIR")) {
          params.last_block = true;  // PDB code implies -O
          grep_file(pdb_dir + mmcif_subpath(path), params);
          params.last_block = p.options[OneBlock];
        } else {
          fprintf(stderr,
                  "The argument %s is a PDB code, but $PDB_DIR is not set.\n"
                  "(To use a file or directory with such a name use: ./%s)\n",
                  path.c_str(), path.c_str());
          return 2;
        }
      } else {
        DirWalker walker(path.c_str());
        for (const tinydir_file& f : walker) {
          if (walker.is_file() || is_cif_file(f)) {
            grep_file(f.path, params);
            file_count++;
          }
        }
      }
    } catch (std::runtime_error& e) {
      std::fflush(stdout);
      fprintf(stderr, "Error when parsing %s:\n\t%s\n", path.c_str(), e.what());
      return 2;
    }
  }
  if (p.options[Summarize])
    printf("Total count in %zu files: %zu\n", file_count, params.total_count);
  return params.total_count != 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
