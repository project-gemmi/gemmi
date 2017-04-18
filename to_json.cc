// Copyright 2017 Global Phasing Ltd.
#include "to_json.hh"
#include "write_cif.hh"
#include <optionparser.h>

struct Arg: public option::Arg {
  // option checkers here
};

enum OptionIndex { Unknown, Help, Bare, QuoteNum };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage: to_json [options] file.cif file.json\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Bare, 0, "b", "--bare-tags", Arg::None,
    "  -b, --bare-tags  \tOutput tags without the first underscore." },
  { QuoteNum, 0, "q", "--quote-numbers", Arg::None,
    "  -q, --quote-numbers  \tOutput all numbers as strings." },
  { 0, 0, 0, 0, 0, 0 }
};

int main(int argc, char **argv) {
  if (argc < 1)
    return 2;
  option::Stats stats(usage, argc-1, argv+1);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc-1, argv+1, options.data(), buffer.data());
  if (parse.error() || options[Unknown] || parse.nonOptionsCount() != 2) {
    option::printUsage(std::cerr, usage);
    return 1;
  }
  if (options[Help]) {
    option::printUsage(std::cerr, usage);
    return 0;
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
  writer.quote_numbers = options[QuoteNum];
  writer.write_json(d);
  delete os;
  return 0;
}

// vim:sw=2:ts=2:et
