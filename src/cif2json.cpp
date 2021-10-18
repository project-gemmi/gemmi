// Copyright 2017 Global Phasing Ltd.

#include <iostream>
#include "gemmi/cifdoc.hpp"    // for Document
#include "gemmi/to_json.hpp"   // for JsonWriter
#include "gemmi/fstream.hpp"   // for Ofstream
#include "gemmi/read_cif.hpp"  // for read_cif_gz
#define GEMMI_PROG cif2json
#include "options.h"
#include "cifmod.h"  // for apply_cif_doc_modifications, ...

namespace {

namespace cif = gemmi::cif;

struct ConvArg: public Arg {
  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }
};

enum OptionIndex { Comcifs=AfterCifModOptions, Mmjson, Bare, Numb, CifDot, };
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nConvert CIF file (any CIF files, including mmCIF) to JSON."
    "\nThe output can be COMCIFS CIF-JSON (-c), mmJSON (-m),"
    "\nor a custom JSON flavor (default)."
    "\n\nGeneral options:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],

  { NoOp, 0, "", "", Arg::None, "\nJSON output options:" },
  { Comcifs, 0, "c", "comcifs", Arg::None,
    "  -c, --comcifs  \tConform to the COMCIFS CIF-JSON standard draft." },
  { Mmjson, 0, "m", "mmjson", Arg::None,
    "  -m, --mmjson   \tCompatible with mmJSON from PDBj." },
  { Bare, 0, "", "bare-tags", Arg::None,
    "  --bare-tags  \tOutput tags without the first underscore." },
  { Numb, 0, "", "numb", ConvArg::NumbChoice,
    "  --numb=quote|nosu|mix  \tConvert the CIF numb type to one of:"
                             "\v  quote - string in quotes,"
                             "\v  nosu - number without s.u.,"
                             "\v  mix (default) - quote only numbs with s.u." },
  { CifDot, 0, "", "dot", Arg::Required,
    "  --dot=STRING  \tJSON representation of CIF's '.' (default: null)." },

  { NoOp, 0, "", "", Arg::None, "\nModifications:" },
  CifModUsage[SkipCat],
  CifModUsage[SortCif],

  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};


void convert(const std::string& input, const std::string& output,
             const std::vector<option::Option>& options) {
  cif::Document doc = gemmi::read_cif_gz(input);
  apply_cif_doc_modifications(doc, options);
  gemmi::Ofstream os(output, &std::cout);
  cif::JsonWriter writer(os.ref());
  if (options[Comcifs])
    writer.set_comcifs();
  if (options[Mmjson])
    writer.set_mmjson();
  if (options[Bare])
    writer.bare_tags = true;
  if (options[Numb]) {
    char first_letter = options[Numb].arg[0];
    if (first_letter == 'q')
      writer.quote_numbers = 2;
    else if (first_letter == 'n')
      writer.quote_numbers = 0;
  }
  if (options[CifDot])
    writer.cif_dot = options[CifDot].arg;
  writer.write_json(doc);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);

  if (p.options[Verbose])
    std::cerr << "Transcribing " << input << " to json ..." << std::endl;
  try {
    convert(input, output, p.options);
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  } catch (std::invalid_argument& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 3;
  }
  if (p.options[Verbose])
    std::cerr << "Done." << std::endl;
  return 0;
}
