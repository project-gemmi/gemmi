// Copyright 2017 Global Phasing Ltd.

// Thin, leaky wrapper around The Lean Mean C++ Option Parser.

#pragma once

#include <vector>
#include <string>
#include <optionparser.h>

#ifndef GEMMI_PROG
# error "define GEMMI_PROG before including options.h"
#endif

#define GEMMI_XSTRINGIZE(s) GEMMI_STRINGIZE(s)
#define GEMMI_STRINGIZE(s) #s
#define GEMMI_XCONCAT(a, b) GEMMI_CONCAT(a, b)
#define GEMMI_CONCAT(a, b) a##b
#ifdef GEMMI_ALL_IN_ONE  // GEMMI_MAIN=foo_main EXE_NAME="gemmi foo"
# define GEMMI_MAIN GEMMI_XCONCAT(GEMMI_PROG, _main)
# define EXE_NAME "gemmi " GEMMI_XSTRINGIZE(GEMMI_PROG)
#else                    // GEMMI_MAIN=main     EXE_NAME="gemmi-foo"
# define GEMMI_MAIN main
# define EXE_NAME "gemmi-" GEMMI_XSTRINGIZE(GEMMI_PROG)
#endif

#if defined(_MSC_VER) && _MSC_VER+0 <= 1900
// warning C4800: 'option::Option *': forcing value to bool 'true' or 'false'
#pragma warning(disable: 4800)
#endif

enum { NoOp=0, Help=1, Version=2, Verbose=3 };

extern const option::Descriptor CommonUsage[];

std::vector<int> parse_comma_separated_ints(const char* arg);
std::vector<double> parse_blank_separated_numbers(const char* arg);

struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg);
  static option::ArgStatus Char(const option::Option& option, bool msg);
  static option::ArgStatus YesNo(const option::Option& option, bool msg);
  static option::ArgStatus Choice(const option::Option& option, bool msg,
                                  const std::vector<const char*>& choices);
  static option::ArgStatus ColonPair(const option::Option& option, bool msg);
  static option::ArgStatus Int(const option::Option& option, bool msg);
  static option::ArgStatus Int3(const option::Option& option, bool msg);
  static option::ArgStatus Float(const option::Option& option, bool msg);
  static option::ArgStatus Float3(const option::Option& option, bool msg);
  static option::ArgStatus CoorFormat(const option::Option& option, bool msg) {
    return Choice(option, msg, {"cif", "mmcif", "pdb", "json", "mmjson", "chemcomp"});
  }
  static option::ArgStatus CifStyle(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"plain", "pdbx", "aligned"});
  }
  static option::ArgStatus AsuChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"ccp4", "tnt"});
  }
};

struct OptParser : option::Parser {
  const char* program_name;
  std::vector<option::Option> options;
  std::vector<option::Option> buffer;

  explicit OptParser(const char* prog) : program_name(prog) {}
  void simple_parse(int argc, char** argv, const option::Descriptor usage[]);
  void require_positional_args(int n);
  void require_input_files_as_args(int other_args=0);
  // returns the nth arg except that a PDB code gets replaced by $PDB_DIR/...
  std::string coordinate_input_file(int n, char pdb_code_type='M');
  std::vector<std::string> paths_from_args_or_file(int opt, int other);
  [[noreturn]] void print_try_help_and_exit(const char* msg) const;
  [[noreturn]] void exit_exclusive(int opt1, int opt2) const;
  void check_exclusive_pair(int opt1, int opt2) {
    if (options[opt1] && options[opt2])
      exit_exclusive(opt1, opt2);
  }
  void check_exclusive_group(const std::vector<int>& group);
  const char* given_name(int opt) const {  // sans one dash
    return options[opt].namelen > 1 ? options[opt].name + 1
                                    : options[opt].desc->shortopt;
  }
  // for Arg::YesNo
  bool is_yes(int opt, bool default_) const {
    if (options[opt])
      return (options[opt].arg[0] & ~0x20) == 'Y';
    return default_;
  }
  int integer_or(int opt, int default_) const;
};

namespace gemmi { enum class CoorFormat; }
namespace gemmi { namespace cif { struct WriteOptions; } }

// to be used with Arg::CoorFormat
gemmi::CoorFormat coor_format_as_enum(const option::Option& format_in);

// to be used with Arg::CifStyle
gemmi::cif::WriteOptions cif_write_options(const option::Option& cif_style);

// can be used with paths_from_args_or_file()
bool starts_with_pdb_code(const std::string& s);

void print_version(const char* program_name, bool verbose=false);

void read_spec_file(const char* path, std::vector<std::string>& output);
