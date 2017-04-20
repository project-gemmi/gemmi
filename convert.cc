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

enum OptionIndex { Unknown, Help, Bare, Numb, QMark };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage: to_json [options] file.cif file.json\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Bare, 0, "b", "bare-tags", Arg::None,
    "  -b, --bare-tags  \tOutput tags without the first underscore." },
  { Numb, 0, "", "numb", Arg::NumbChoice,
    "  --numb=quote|nosu|mix  \tConvert the CIF numb type to one of:"
                             "\v  quote - string in quotes,"
                             "\v  nosu - number without s.u.,"
                             "\v  mix (default) - quote only numbs with s.u." },
  { QMark, 0, "", "unknown", Arg::Required,
    "  --unknown=STRING  \tJSON representation of CIF's '?' (default: null)." },
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
              << parse.nonOptionsCount() << " given.\n";
    return 1;
  }
  const char* filename = parse.nonOption(0);
  gemmi::cif::Document d;
  try {
    d.read_file(filename);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  const char* output = parse.nonOption(1);
  std::ostream* os = nullptr;
  if (output != std::string("-")) {
    os = new std::ofstream(output);
    if (!*os) {
      std::cerr << "Failed to open for writing: " << output;
      delete os;
      return 1;
    }
  }
  gemmi::cif::JsonWriter writer(os ? *os : std::cout);
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
  delete os;
  return 0;
}

// vim:sw=2:ts=2:et
