// Copyright 2017 Global Phasing Ltd.

#include "gemmi/cif.hpp"
#include "gemmi/cifgz.hpp"
#include <cstdio>
#include <stdexcept>
#include <string>
#include <optionparser.h>

#define EXE_NAME "gemmi-grep"

struct Arg: public option::Arg {
  static option::ArgStatus Int(const option::Option& option, bool msg) {
    if (option.arg) {
      char* endptr = nullptr;
      std::strtol(option.arg, &endptr, 10);
      if (endptr != option.arg && *endptr == '\0')
        return option::ARG_OK;
    }
    if (msg)
      fprintf(stderr, "Option '%s' requires a numeric argument\n", option.name);
    return option::ARG_ILLEGAL;
  }
};

enum OptionIndex { Unknown, Help, MaxCount, WithFileName, NoBlockName,
                   Summarize, MatchingFiles, NonMatchingFiles, Count };

const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage: " EXE_NAME " [options] TAG FILE_OR_DIR[...]\n"
    "Search for TAG in CIF files."
    "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None,
    "  -h, --help  \tprint usage and exit" },
  { MaxCount, 0, "m", "max-count", Arg::Int,
    "  -m, --max-count=NUM  \tprint max NUM values per block (default: 10)" },
  { WithFileName, 0, "H", "with-filename", Arg::None,
    "  -H, --with-filename  \tprint the file name for each match" },
  { NoBlockName, 0, "b", "no-blockname", Arg::None,
    "  -b, --no-blockname  \tsuppress the block name on output" },
  { Summarize, 0, "s", "summarize", Arg::None,
    "  -s, --summarize  \tdisplay only statistics" },
  { MatchingFiles, 0, "l", "files-with-tag", Arg::None,
    "  -l, --files-with-tag  \tprint only names of files with the tag" },
  { NonMatchingFiles, 0, "L", "files-without-tag", Arg::None,
    "  -L, --files-without-tag  \tprint only names of files without the tag" },
  { Count, 0, "c", "count", Arg::None,
    "  -c, --count  \tprint only a count of matching lines per file" },
  { 0, 0, 0, 0, 0, 0 }
};

struct ParsedOptions {
  int max_count = 10;
  bool with_filename = false;
  bool with_blockname = true;
  bool summarize = false;
  bool only_filenames = false;
  bool inverse = false;
  bool count = false;
};

static
void grep_file(const std::string& tag, const char* path,
               const ParsedOptions& opt) {
  using namespace gemmi::cif;
  Document d = read_any(path);
  for (Block& block : d.blocks) {
    TableView t = block.find(tag);
    if (opt.only_filenames) {
      if (t.ok() != opt.inverse)
        printf("%s\n", path);
      continue;
    }
    if (opt.count) {
      if (opt.with_filename)
        printf("%s: ", path);
      if (opt.with_blockname)
        printf("%s: ", block.name.c_str());
      printf(" %d\n", (int) t.length());
      continue;
    }
    int cnt = 0;
    for (const TableView::Row& row : t) {
      if (opt.with_filename)
        printf("%s: ", path);
      if (opt.with_blockname)
        printf("%s: ", block.name.c_str());
      printf("%s\n", row.as_str(0).c_str());
      if (++cnt == opt.max_count)
        break;
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 1)
    return 2;
  option::Stats stats(usage, argc-1, argv+1);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc-1, argv+1, options.data(), buffer.data());
  if (parse.error() || options[Unknown] ||
      (!options[Help] && parse.nonOptionsCount() < 2)) {
    option::printUsage(fwrite, stderr, usage);
    return 1;
  }
  if (options[Help]) {
    option::printUsage(fwrite, stdout, usage);
    return 0;
  }

  ParsedOptions parsed_options;
  if (options[MaxCount])
    parsed_options.max_count = std::strtol(options[MaxCount].arg, nullptr, 10);
  if (options[WithFileName])
    parsed_options.with_filename = true;
  if (options[NoBlockName])
    parsed_options.with_blockname = false;
  if (options[Summarize])
    parsed_options.summarize = true;
  if (options[MatchingFiles])
    parsed_options.only_filenames = true;
  if (options[NonMatchingFiles]) {
    parsed_options.only_filenames = true;
    parsed_options.inverse = true;
  }
  if (options[Count])
    parsed_options.count = true;

  std::string tag = parse.nonOption(0);

  for (int i = 1; i < parse.nonOptionsCount(); ++i) {
    const char* path = parse.nonOption(i);
    try {
      grep_file(tag, path, parsed_options);
    } catch (std::runtime_error& e) {
      fprintf(stderr, "Error when parsing %s:\n\t%s\n", path, e.what());
      return 1;
    }
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=include,third_party
