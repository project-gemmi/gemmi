// Copyright 2017 Global Phasing Ltd.

// Common code for parsing options in gemmi utilities,
// on top of The Lean Mean C++ Option Parser.

#include <cstdio>   // for fprintf
#include <cstdlib>  // for strtol, strtod, exit
#include <cstring>  // for strcmp
#include <vector>
#include <optionparser.h>
#include "gemmi/version.hpp"

enum { NoOp=0, Help=1, Version=2 };

inline std::vector<int> parse_comma_separated_ints(const char* arg) {
  std::vector<int> result;
  char* endptr = nullptr;
  do {
    result.push_back(std::strtol(endptr ? endptr + 1 : arg, &endptr, 10));
  } while (*endptr == ',');
  if (*endptr != '\0')
    result.clear();
  return result;
}

struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg) {
    if (option.arg != nullptr)
      return option::ARG_OK;
    if (msg)
      fprintf(stderr, "Option '%s' requires an argument\n", option.name);
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus Choice(const option::Option& option, bool msg,
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

  static option::ArgStatus Int(const option::Option& option, bool msg) {
    if (option.arg) {
      char* endptr = nullptr;
      std::strtol(option.arg, &endptr, 10);
      if (endptr != option.arg && *endptr == '\0')
        return option::ARG_OK;
    }
    if (msg)
      fprintf(stderr,
              "Option '%s' requires an integer argument\n", option.name);
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Int3(const option::Option& option, bool msg) {
    if (option.arg && parse_comma_separated_ints(option.arg).size() == 3)
        return option::ARG_OK;
    if (msg)
      fprintf(stderr, "Option '%.*s' requires three comma-separated integers "
                      "as an argument,\n for example: %.*s=11,12,13",
                      option.namelen, option.name, option.namelen, option.name);
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus Float(const option::Option& option, bool msg) {
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
};

struct OptParser : option::Parser {
  const std::vector<option::Option>&
      simple_parse(int argc, char** argv, const option::Descriptor usage[]) {
    if (argc < 1)
      std::exit(2);
    option::Stats stats(usage, argc-1, argv+1);
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
      printf("%s %s\n", EXE_NAME, GEMMI_VERSION);
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
    return options;
  }

  const char* given_name(int opt) const {  // sans one dash
    return options[opt].namelen > 1 ? options[opt].name + 1
                                    : options[opt].desc->shortopt;
  }

  void require_positional_args(int n) {
    if (nonOptionsCount() != n) {
      fprintf(stderr, EXE_NAME " requires %d arguments but got %d.\n"
                      "Try '" EXE_NAME " --help' for more information.\n",
                      n, nonOptionsCount());
      std::exit(2);
    }
  }

  std::vector<option::Option> options;
  std::vector<option::Option> buffer;
  std::vector<std::vector<int>> exclusive_groups;
};

// vim:sw=2:ts=2:et:path^=../include,../third_party
