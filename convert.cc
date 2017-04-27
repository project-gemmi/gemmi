// Copyright 2017 Global Phasing Ltd.
#include "to_json.hh"
#include "write_cif.hh"

#include <cstring>
#include <optionparser.h>

struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg) {
    if (option.arg != nullptr)
      return option::ARG_OK;
    if (msg)
      std::cerr << "Option '" << option.name << "' requires an argument\n";
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    if (Required(option, msg) == option::ARG_ILLEGAL)
      return option::ARG_ILLEGAL;
    for (const char* a : {"quote", "nosu", "mix"})
      if (strcmp(option.arg, a) == 0)
        return option::ARG_OK;
    if (msg)
      std::cerr << "Invalid argument for " << option.name << ": "
                << option.arg << "\n";
    return option::ARG_ILLEGAL;
  }
};

enum OptionIndex { Unknown, Help, Format, Bare, Numb, QMark };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:"
    "\n gemmi-convert [options] file.cif file.json"
    "\n gemmi-convert [options] file.cif file.pdb"
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Format, 0, "f", "", Arg::None,
    "  -f json|cif  \tOutput format (default: the file extension)." },
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
    "\nWhen input file is -, read standard input."
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
                 "Try 'gemmi-convert --help' for more information.\n";
    return 1;
  }

  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  char output_format;
  if (options[Format]) {
    output_format = options[Format].arg[0];
  } else {
    const char* dot = std::strrchr(output, '.');
    if (strcmp(dot, "pdb") == 0 || strcmp(dot, "PDB") == 0 ||
        strcmp(dot, "ent") == 0 || strcmp(dot, "ENT") == 0) {
      output_format = 'p';
    } else if ((dot[1] == 'j' || dot[1] == 'J') &&
               (dot[2] == 's' || dot[2] == 'S')) {  // .json, .js, .jswhatever
      output_format = 'j';
    } else {
      std::cerr << "The output format cannot be determined from the output"
                   " filename. Use option -f.\n";
      return 1;
    }
  }
  gemmi::cif::Document d;
  try {
    d.read_file(input);
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
    //gemmi::mol::write_pdb(d, *os);
  }
  if (new_os)
    delete os;
  return 0;
}

// vim:sw=2:ts=2:et
