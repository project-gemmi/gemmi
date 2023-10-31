// Copyright 2018 Global Phasing Ltd.

#define GEMMI_PROG na
#include "options.h"
#include <cstdio>   // for fprintf, fopen
#include <cstdlib>  // for strtol, strtod, exit
#include <cstring>  // for strcmp, strchr
#include <gemmi/fileutil.hpp>  // for file_open_or
#include <gemmi/pdb_id.hpp>    // for expand_if_pdb_code
#include <gemmi/version.hpp>   // for GEMMI_VERSION
#include <gemmi/util.hpp>   // for trim_str
#include <gemmi/model.hpp>  // for gemmi::CoorFormat
#include <gemmi/to_cif.hpp>  // for gemmi::cif::WriteOptions
#include <gemmi/atox.hpp>   // for skip_blank

using std::fprintf;

const option::Descriptor CommonUsage[] = {
  { 0, 0, 0, 0, 0, 0 }, // this makes CommonUsage[Help] return Help item, etc
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None,
    "  -v, --verbose  \tVerbose output." }
};

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

std::vector<double> parse_blank_separated_numbers(const char* arg) {
  std::vector<double> results;
  for (char* endptr; *arg != '\0'; arg = endptr) {
    results.push_back(std::strtod(arg, &endptr));
    if (endptr == arg || (*endptr != ' ' && *endptr != '\0')) {
      results.clear();
      break;
    }
  }
  return results;
}

option::ArgStatus Arg::Required(const option::Option& option, bool msg) {
  if (option.arg != nullptr)
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Option '%s' requires an argument\n", option.name);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Char(const option::Option& option, bool msg) {
  if (Required(option, msg) == option::ARG_ILLEGAL)
    return option::ARG_ILLEGAL;
  if (option.arg[0] != '\0' || option.arg[1] == '\0')
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Argument of '%s' must be one character\n", option.name);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::YesNo(const option::Option& option, bool msg) {
  if (Required(option, msg) == option::ARG_ILLEGAL)
    return option::ARG_ILLEGAL;
  char first_letter = gemmi::alpha_up(option.arg[0]);
  if (first_letter == 'Y' || first_letter == 'N')
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Argument of '%s' must be Y (yes) or N (no)\n", option.name);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Choice(const option::Option& option, bool msg,
                              const std::vector<const char*>& choices) {
  if (Required(option, msg) == option::ARG_ILLEGAL)
    return option::ARG_ILLEGAL;
  for (const char* a : choices)
    if (std::strcmp(option.arg, a) == 0)
      return option::ARG_OK;
  if (msg) {
    // option.name here is a string "--option=arg"
    fprintf(stderr, "Invalid argument for %.*s: %s\nAllowed arguments:",
            option.namelen, option.name, option.arg);
    for (const char* a : choices)
      fprintf(stderr, " %s", a);
    fprintf(stderr, "\n");
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::ColonPair(const option::Option& option, bool msg) {
  if (Required(option, msg) == option::ARG_ILLEGAL)
    return option::ARG_ILLEGAL;
  const char* sep = std::strchr(option.arg, ':');
  if (sep != nullptr && std::strchr(sep+1, ':') == nullptr)
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Option '%.*s' requires two colon-separated names "
                    "as an argument,\n for example: %.*s=A:B\n",
                    option.namelen, option.name, option.namelen, option.name);
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
                    "as an argument,\n for example: %.*s=11,12,13\n",
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

option::ArgStatus Arg::Float3(const option::Option& option, bool msg) {
  if (option.arg && parse_blank_separated_numbers(option.arg).size() == 3)
    return option::ARG_OK;
  if (msg)
    fprintf(stderr, "Option '%.*s' requires three numbers as an argument,\n"
                    " for example: %.*s='1.1 2.2 3'\n",
                    option.namelen, option.name, option.namelen, option.name);
  return option::ARG_ILLEGAL;
}

// we wrap fwrite because passing it directly may cause warning
// "ignoring attributes on template argument" [-Wignored-attributes]
static
size_t write_func(const void *ptr, size_t size, size_t nmemb, FILE *stream) {
  return fwrite(ptr, size, nmemb, stream);
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
    option::printUsage(write_func, stdout, usage);
    std::exit(0);
  }
  if (options[Version]) {
    print_version(program_name, options[Verbose]);
    std::exit(0);
  }
  if (options[NoOp]) {
    fprintf(stderr, "Invalid option.\n");
    option::printUsage(write_func, stderr, usage);
    std::exit(2);
  }
}

void OptParser::check_exclusive_group(const std::vector<int>& group) {
  int first = -1;
  for (int opt : group)
    if (options[opt]) {
      if (first == -1)
        first = opt;
      else
        exit_exclusive(first, opt);
    }
}

int OptParser::integer_or(int opt, int default_) const {
  if (options[opt])
    return std::atoi(options[opt].arg);
  return default_;
}

void OptParser::print_try_help_and_exit(const char* msg) const {
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

void OptParser::exit_exclusive(int opt1, int opt2) const {
  std::fprintf(stderr, "Options -%s and -%s cannot be used together.\n",
               given_name(opt1), given_name(opt2));
  std::exit(1);
}

std::string OptParser::coordinate_input_file(int n, char pdb_code_type) {
  return gemmi::expand_if_pdb_code(nonOption(n), pdb_code_type);
}

bool starts_with_pdb_code(const std::string& s) {
  return (s.length() == 4 && gemmi::is_pdb_code(s)) ||
         (s.length() > 4 && std::strchr(" \t\r\n:,;|", s[4]) &&
          gemmi::is_pdb_code(s.substr(0, 4)));
}

std::vector<std::string>
OptParser::paths_from_args_or_file(int opt, int other) {
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
      if (!s.empty())
        paths.emplace_back(s);
    }
    std::fclose(f);
  } else {
    require_input_files_as_args(other);
    for (int i = other; i < nonOptionsCount(); ++i)
      paths.emplace_back(nonOption(i));
  }
  return paths;
}

gemmi::CoorFormat coor_format_as_enum(const option::Option& format_in) {
  auto eq = [&](const char* s) { return std::strcmp(format_in.arg, s) == 0; };
  gemmi::CoorFormat format = gemmi::CoorFormat::Unknown;
  if (format_in) {
    if (eq("cif") || eq("mmcif"))
      format = gemmi::CoorFormat::Mmcif;
    else if (eq("pdb"))
      format = gemmi::CoorFormat::Pdb;
    else if (eq("json") || eq("mmjson"))
      format = gemmi::CoorFormat::Mmjson;
    else if (eq("chemcomp"))
      format = gemmi::CoorFormat::ChemComp;
  }
  return format;
}

gemmi::cif::WriteOptions cif_write_options(const option::Option& cif_style) {
  gemmi::cif::WriteOptions options;
  options.prefer_pairs = true;
  if (cif_style)
    switch (cif_style.arg[0]) {
      // value for 'd' (default) is returned at end of this function
      case 'p'/*pdbx*/:
        options.misuse_hash = true;
        break;
      case 'a'/*aligned*/:
        options.align_pairs = 33;
        options.align_loops = 30;
        break;
    }
  return options;
}

void print_version(const char* program_name, bool verbose) {
  std::printf("%s " GEMMI_VERSION
#ifdef GEMMI_VERSION_INFO
         " (" GEMMI_XSTRINGIZE(GEMMI_VERSION_INFO) ")"
#endif
         "\n", program_name);
  if (verbose) {
#if defined(_MSC_VER)
    std::printf("Compiler: MSVC %d (C++ %ld)\n", _MSC_FULL_VER, _MSVC_LANG);
#else
#  if defined(__clang__)
    std::printf("Compiler: Clang %d.%d.%d (C++ %ld)\n",
                __clang_major__, __clang_minor__, __clang_patchlevel__, __cplusplus);
#  elif defined(__GNUC__)
    std::printf("Compiler: GCC %d.%d.%d (C++ %ld)\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, __cplusplus);
#  endif
#endif
  }
}

void read_spec_file(const char* path, std::vector<std::string>& output) {
  char buf[256];
  gemmi::fileptr_t f_spec = gemmi::file_open_or(path, "r", stdin);
  while (std::fgets(buf, sizeof(buf), f_spec.get()) != NULL) {
    const char* start = gemmi::skip_blank(buf);
    if (*start != '\0' && *start != '\r' && *start != '\n' && *start != '#')
      output.emplace_back(start);
  }
}
