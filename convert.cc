// Copyright 2017 Global Phasing Ltd.
#include "to_json.hh"
#define STB_SPRINTF_IMPLEMENTATION
#include "to_pdb.hh"
#include "write_cif.hh"
#include "cifgz.hh"
#include "mmcif.hh"

#include <cstring>
#include <iostream>
#include <optionparser.h>

#define EXE_NAME "gemmi-convert"

struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg) {
    if (option.arg != nullptr)
      return option::ARG_OK;
    if (msg)
      std::cerr << "Option '" << option.name << "' requires an argument\n";
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Choice(const option::Option& option, bool msg,
                                  std::vector<const char*> choices) {
    if (Required(option, msg) == option::ARG_ILLEGAL)
      return option::ARG_ILLEGAL;
    for (const char* a : choices)
      if (strcmp(option.arg, a) == 0)
        return option::ARG_OK;
    if (msg)
      std::cerr << "Invalid argument for "
                << std::string(option.name, option.namelen) << ": "
                << option.arg << "\n";
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"json", "pdb", "cif"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }
};

enum OptionIndex { Unknown, Help, FormatOut, Bare, Numb, QMark };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] file.cif file.json"
    "\n " EXE_NAME " [options] file.cif file.pdb"
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { FormatOut, 0, "", "to", Arg::FileFormat,
    "  --to=json|pdb  \tOutput format (default: the file extension)." },
  { Unknown, 0, "", "", Arg::None, "\nCIF output options:" },
  { Bare, 0, "b", "bare-tags", Arg::None,
    "  -b, --bare-tags  \tOutput tags without the first underscore." },
  { Numb, 0, "", "numb", Arg::NumbChoice,
    "  --numb=quote|nosu|mix  \tConvert the CIF numb type to one of:"
                             "\v  quote - string in quotes,"
                             "\v  nosu - number without s.u.,"
                             "\v  mix (default) - quote only numbs with s.u." },
  { QMark, 0, "", "unknown", Arg::Required,
    "  --unknown=STRING  \tJSON representation of CIF's '?' (default: null)." },
  { Unknown, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

int main(int argc, char **argv) {
  if (argc < 1)
    return 2;
  option::Stats stats(usage, argc-1, argv+1);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc-1, argv+1, options.data(), buffer.data());
  if (parse.error()) {
    option::printUsage(std::cerr, usage);
    return 1;
  }
  if (options[Help]) {
    option::printUsage(std::cout, usage);
    return 0;
  }
  if (options[Unknown]) {
    std::cerr << "Invalid option.\n";
    option::printUsage(std::cerr, usage);
    return 1;
  }
  if (parse.nonOptionsCount() != 2) {
    std::cerr << "This program requires 2 arguments (input and output), "
              << parse.nonOptionsCount() << " given.\n"
                 "Try '" EXE_NAME " --help' for more information.\n";
    return 1;
  }

  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  char output_format = 0;
  if (options[FormatOut]) {
    output_format = options[FormatOut].arg[0];
  } else {
    const char* dot = std::strrchr(output, '.');
    if (dot) {
      if (strcmp(dot+1, "pdb") == 0 || strcmp(dot+1, "PDB") == 0 ||
          strcmp(dot+1, "ent") == 0 || strcmp(dot+1, "ENT") == 0) {
        output_format = 'p';
      } else if ((dot[1] == 'j' || dot[1] == 'J') &&
                 (dot[2] == 's' || dot[2] == 'S')) { // .json, .js, .jswhatever
        output_format = 'j';
      } else if (strcmp(dot+1, "cif") == 0 || strcmp(dot+1, "CIF") == 0) {
        output_format = 'c';
      }
    }
    if (output_format == 0) {
      std::cerr << "The output format cannot be determined from output"
                   " filename. Use option --to.\n";
      return 1;
    }
  }

  gemmi::cif::Document d;
  try {
    d = gemmi::cif::read_any(input);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  std::ostream* os;
  bool new_os = false;
  if (output != std::string("-")) {
    os = new std::ofstream(output);
    if (*os) {
      new_os = true;
    } else {
      std::cerr << "Failed to open for writing: " << output;
      delete os;
      return 1;
    }
  } else {
    os = &std::cout;
  }

  if (output_format == 'j') {
    gemmi::cif::JsonWriter writer(*os);
    writer.use_bare_tags = options[Bare];
    if (options[Numb]) {
      char first_letter = options[Numb].arg[0];
      if (first_letter == 'q')
        writer.quote_numbers = 2;
      else if (first_letter == 'n')
        writer.quote_numbers = 0;
    }
    if (options[QMark])
      writer.unknown = options[QMark].arg;
    writer.write_json(d);
  } else if (output_format == 'p') {
    try {
      gemmi::mol::Structure st = gemmi::mol::read_atoms(d);
      gemmi::mol::write_pdb(st, *os);
    } catch (std::runtime_error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return 2;
    }
  } else if (output_format == 'c') {
    *os << d;
  }

  if (new_os)
    delete os;
  return 0;
}

// vim:sw=2:ts=2:et
