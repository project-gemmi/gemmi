// Copyright 2017 Global Phasing Ltd.

#include <iostream> // temporary, for debugging
#include "gemmi/cifgz.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/pdbgz.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/to_pdb.hpp"
// set this before only one of stb_sprintf.h includes
#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/to_mmcif.hpp"

#include <cstring>
#include <iostream>
#include <map>
#include <optionparser.h>

#define EXE_NAME "gemmi-convert"

enum class FileType : char { Json, Pdb, Cif, Null, Unknown };

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

enum OptionIndex { Unknown, Help, Verbose, FormatIn, FormatOut,
                   Bare, Numb, QMark };
static const option::Descriptor usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions: cif->json and cif<->pdb."
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
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

FileType get_format_from_extension(const std::string& path) {
  using gemmi::iends_with;
  if (iends_with(path, ".pdb") || iends_with(path, ".ent") ||
      iends_with(path, ".pdb.gz") || iends_with(path, ".ent.gz"))
    return FileType::Pdb;
  if (iends_with(path, ".js") || iends_with(path, ".json"))
    return FileType::Json;
  if (iends_with(path, ".cif") || iends_with(path, ".cif.gz"))
    return FileType::Cif;
  if (path == "/dev/null")
    return FileType::Null;
  return FileType::Unknown;
}

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

void convert(const char* input, FileType input_type,
             const char* output, FileType output_type,
             const std::vector<option::Option>& options) {
  gemmi::cif::Document cif_in;
  gemmi::mol::Structure st;
  // hidden feature for testing cif -> Structure -> cif
  bool force_structure = gemmi::ends_with(output, ".ciF");
  if (input_type == FileType::Cif) {
    cif_in = gemmi::cif::read_any(input);
    if (output_type == FileType::Pdb || output_type == FileType::Null ||
        force_structure) {
      st = gemmi::mol::read_atoms(cif_in);
      if (st.models.empty())
        fail("No atoms in the input file. Is it mmCIF?");
    }
  } else if (input_type == FileType::Pdb) {
    st = gemmi::mol::read_pdb_any(input);
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

  if (output_type == FileType::Json) {
    if (input_type != FileType::Cif)
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

  else if (output_type == FileType::Pdb || output_type == FileType::Null) {
    if (output_type == FileType::Pdb)
      gemmi::mol::write_pdb(st, *os);
    else {
      *os << st.name << ": " << count_atom_sites(st) << " atom locations";
      if (st.models.size() > 1)
        *os << " (total in " << st.models.size() << " models)";
      *os << ".\n";
    }
  } else if (output_type == FileType::Cif) {
    // cif to cif round trip is for testing only
    if (input_type != FileType::Cif || force_structure) {
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
  std::ios_base::sync_with_stdio(false);
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

  std::map<std::string, FileType> filetypes {{"json", FileType::Json},
                                             {"pdb", FileType::Pdb},
                                             {"cif", FileType::Cif},
                                             {"none", FileType::Null}};

  FileType in_type = options[FormatIn] ? filetypes[options[FormatIn].arg]
                                       : get_format_from_extension(input);
  if (in_type == FileType::Unknown) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

  FileType out_type = options[FormatOut] ? filetypes[options[FormatOut].arg]
                                         : get_format_from_extension(output);
  if (out_type == FileType::Unknown) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }
  if (options[Verbose])
    std::cerr << "Converting " << input << " ..." << std::endl;

  try {
    convert(input, in_type, output, out_type, options);
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
