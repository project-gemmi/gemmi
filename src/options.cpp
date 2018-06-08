// Copyright 2018 Global Phasing Ltd.

#define GEMMI_PROG na
#include "options.h"
#include <cstdio>   // for fprintf, fopen
#include <cstdlib>  // for strtol, strtod, exit
#include <cstring>  // for strcmp, strchr
#include <gemmi/fileutil.hpp>  // for expand_if_pdb_code
#include <gemmi/version.hpp>   // for GEMMI_VERSION
#include <gemmi/util.hpp>   // for trim_str

using std::fprintf;

std::vector<int> parse_comma_separated_ints(const char* arg) {
  std::vector<int> result;
  char* endptr = nullptr;
  do {
    result.push_back(std::strtol(endptr ? endptr + 1 : arg, &endptr, 10));
  } while (*endptr == ',');
  if (*endptr != '\0')
    result.clear();
  return result;
}

option::ArgStatus Arg::Required(const option::Option& option, bool msg) {
  if (option.arg != nullptr)
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Option '%s' requires an argument\n", option.name);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Choice(const option::Option& option, bool msg,
                              std::vector<const char*> choices) {
  if (Required(option, msg) == option::ARG_ILLEGAL)
    return option::ARG_ILLEGAL;
  for (const char* a : choices)
    if (std::strcmp(option.arg, a) == 0)
      return option::ARG_OK;
  if (msg)
    // option.name here is a string "--option=arg"
    fprintf(stderr, "Invalid argument for %.*s: %s\n",
            option.namelen, option.name, option.arg);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Int(const option::Option& option, bool msg) {
  if (option.arg) {
    char* endptr = nullptr;
    std::strtol(option.arg, &endptr, 10);
    if (endptr != option.arg && *endptr == '\0')
      return option::ARG_OK;
  }
  if (msg)
    fprintf(stderr, "Option '%s' requires an integer argument\n", option.name);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Int3(const option::Option& option, bool msg) {
  if (option.arg && parse_comma_separated_ints(option.arg).size() == 3)
      return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Option '%.*s' requires three comma-separated integers "
                    "as an argument,\n for example: %.*s=11,12,13",
                    option.namelen, option.name, option.namelen, option.name);
  return option::ARG_ILLEGAL;
}


option::ArgStatus Arg::Float(const option::Option& option, bool msg) {
  if (option.arg) {
    char* endptr = nullptr;
    std::strtod(option.arg, &endptr);
    if (endptr != option.arg && *endptr == '\0')
      return option::ARG_OK;
  }
  if (msg)
    fprintf(stderr, "Option '%s' requires a numeric argument\n", option.name);
  return option::ARG_ILLEGAL;
}

void OptParser::simple_parse(int argc, char** argv,
                             const option::Descriptor usage[]) {
  if (argc < 1)
    std::exit(2);
  option::Stats stats(/*reordering*/true, usage, argc-1, argv+1);
  options.resize(stats.options_max);
  buffer.resize(stats.buffer_max);
  parse(usage, argc-1, argv+1, options.data(), buffer.data());
  if (error())
    std::exit(2);
  if (options[Help]) {
    option::printUsage(fwrite, stdout, usage);
    std::exit(0);
  }
  if (options[Version]) {
    printf("%s %s\n", program_name, GEMMI_VERSION);
    std::exit(0);
  }
  if (options[NoOp]) {
    fprintf(stderr, "Invalid option.\n");
    option::printUsage(fwrite, stderr, usage);
    std::exit(2);
  }
  for (const auto& group : exclusive_groups) {
    int first = 0;
    for (int opt : group)
      if (options[opt]) {
        if (first == 0) {
          first = opt;
        } else {
          fprintf(stderr, "Options -%s and -%s cannot be used together.\n",
                  given_name(first), given_name(opt));
          std::exit(2);
        }
      }
  }
}

void OptParser::print_try_help_and_exit(const char* msg) {
  fprintf(stderr, "%s\nTry '%s --help' for more information.\n",
                  msg, program_name);
  std::exit(2);
}

void OptParser::require_positional_args(int n) {
  if (nonOptionsCount() != n) {
    fprintf(stderr, "%s requires %d arguments but got %d.",
                    program_name, n, nonOptionsCount());
    print_try_help_and_exit("");
  }
}

void OptParser::require_input_files_as_args(int other_args) {
  if (nonOptionsCount() <= other_args)
    print_try_help_and_exit("No input files. Nothing to do.");
}

std::string OptParser::coordinate_input_file(int n) {
  return gemmi::expand_if_pdb_code(nonOption(n));
}

std::vector<std::string>
OptParser::paths_from_args_or_file(int opt, int other, bool expand) {
  std::vector<std::string> paths;
  const option::Option& file_option = options[opt];
  if (file_option) {
    if (nonOptionsCount() > other)
      print_try_help_and_exit("Error: File arguments together with option -f.");
    std::FILE *f = std::fopen(file_option.arg, "r");
    if (!f) {
      std::perror(file_option.arg);
      std::exit(2);
    }
    char buf[512];
    while (std::fgets(buf, 512, f)) {
      std::string s = gemmi::trim_str(buf);
      if (s.length() > 4 && std::strchr(" \t\r\n:,;|", s[4]) &&
          gemmi::is_pdb_code(s.substr(0, 4)))
        s.resize(4);
      if (!s.empty())
        paths.emplace_back(s);
    }
    std::fclose(f);
  } else {
    require_input_files_as_args(other);
    for (int i = other; i < nonOptionsCount(); ++i)
      paths.emplace_back(nonOption(i));
  }
  if (expand)
    for (std::string& path : paths)
      path = gemmi::expand_if_pdb_code(path);
  return paths;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
