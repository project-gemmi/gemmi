// Copyright 2017 Global Phasing Ltd.

#define STB_SPRINTF_IMPLEMENTATION  // must be set before to_pdb.hh
#include <iostream> // temporary, for debugging
#include "gemmi/cifgz.hh"
#include "gemmi/mmcif.hh"
#include "gemmi/pdb.hh"
#include "gemmi/to_cif.hh"
#include "gemmi/to_json.hh"
#include "gemmi/to_mmcif.hh"
#include "gemmi/to_pdb.hh"

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
    // the hidden option "none" is for testing only
    return Arg::Choice(option, msg, {"json", "pdb", "cif", "none"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }
};

enum OptionIndex { Unknown, Help, FormatIn, FormatOut, Bare, Numb, QMark };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions: cif->json and cif<->pdb."
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { FormatIn, 0, "", "from", Arg::FileFormat,
    "  --from=pdb|cif  \tInput format (default: from the file extension)." },
  { FormatOut, 0, "", "to", Arg::FileFormat,
    "  --to=json|pdb  \tOutput format (default: from the file extension)." },
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


char get_format_from_extension(const std::string& path) {
  using gemmi::iends_with;
  if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
    return 'p';
  if (iends_with(path, ".js") || iends_with(path, ".json"))
    return 'j';
  if (iends_with(path, ".cif") || iends_with(path, ".cif.gz"))
    return 'c';
  if (path == "/dev/null")
    return 'n';
  return 0;
}

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

void convert(const char* input, char input_format,
             const char* output, char output_format,
             const std::vector<option::Option>& options) {
  gemmi::cif::Document cif_in;
  gemmi::mol::Structure st;
  // hidden feature for testing cif -> Structure -> cif
  bool force_structure = gemmi::ends_with(output, ".ciF");
  if (input_format == 'c') {
    cif_in = gemmi::cif::read_any(input);
    if (output_format == 'p' || output_format == 'n' || force_structure) {
      st = gemmi::mol::read_atoms(cif_in);
      if (st.models.empty())
        fail("No atoms in the input file. Is it mmCIF?");
    }
  } else if (input_format == 'p') {
    st = gemmi::mol::read_pdb(input);
  } else {
    fail("Unexpected input format.");
  }

  std::ostream* os;
  std::unique_ptr<std::ostream> os_deleter;
  if (output != std::string("-")) {
    os_deleter.reset(new std::ofstream(output));
    os = os_deleter.get();
    if (!os || !*os)
      fail("Failed to open for writing: " + std::string(output));
  } else {
    os = &std::cout;
  }

  if (output_format == 'j') {
    if (input_format != 'c')
      fail("Conversion to JSON is possible only from CIF");
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
    writer.write_json(cif_in);
  }

  else if (output_format == 'p' || output_format == 'n') {
    if (output_format == 'p')
      gemmi::mol::write_pdb(st, *os);
    else {
      *os << st.name << ": " << count_atom_sites(st) << " atom locations";
      if (st.models.size() > 1)
        *os << " (total in " << st.models.size() << " models)";
      *os << ".\n";
    }
  } else if (output_format == 'c') {
    // cif to cif round trip is for testing only
    if (input_format != 'c' || force_structure) {
      cif_in.blocks.clear();  // temporary, for testing
      cif_in.blocks.resize(1);
      gemmi::mol::update_cif_block(st, cif_in.blocks[0]);
    }
    *os << cif_in;
  }
}

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

  char input_format = options[FormatIn] ? options[FormatIn].arg[0]
                                        : get_format_from_extension(input);
  if (input_format == 0) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

  char output_format = options[FormatOut] ? options[FormatOut].arg[0]
                                          : get_format_from_extension(output);
  if (output_format == 0) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }

  try {
    convert(input, input_format, output, output_format, options);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  }
  return 0;
}

// vim:sw=2:ts=2:et
