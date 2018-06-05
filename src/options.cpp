// Copyright 2018 Global Phasing Ltd.

#define GEMMI_PROG na
#include "options.h"
#include <cstdio>   // for fprintf
#include <cstdlib>  // for strtol, strtod, exit
#include <cstring>  // for strcmp
#include <gemmi/fileutil.hpp>  // for expand_if_pdb_code
#include <gemmi/version.hpp>   // for GEMMI_VERSION

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

void OptParser::require_positional_args(int n) {
  if (nonOptionsCount() != n) {
    fprintf(stderr, "%s requires %d arguments but got %d.\n"
                    "Try '%s --help' for more information.\n",
                    program_name, n, nonOptionsCount(), program_name);
    std::exit(2);
  }
}

void OptParser::require_input_files_as_args() {
  if (nonOptionsCount() == 0) {
    std::fprintf(stderr, "No input files. Nothing to do.\n");
    std::exit(0);
  }
}

std::string OptParser::coordinate_input_file(int n) {
  return gemmi::expand_if_pdb_code(nonOption(n));
}
// vim:sw=2:ts=2:et:path^=../include,../third_party
