// Copyright 2017 Global Phasing Ltd.

#include "gemmi/to_cif.hpp"   // for JsonWriter
#include "gemmi/fstream.hpp"  // for Ofstream
#include "gemmi/read_cif.hpp" // for read_cif_gz

#include <iostream>

#define GEMMI_PROG json2cif
#include "options.h"
#include "cifmod.h"  // for apply_cif_doc_modifications, ...

namespace {

enum OptionIndex { CifStyle=AfterCifModOptions, Cif2Cif };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nConvert mmJSON to mmCIF."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],

  { CifStyle, 0, "", "style", Arg::CifStyle,
    "  --style=STYLE  \tOne of: default, pdbx (categories separated with #),"
                     " aligned (left-aligned columns)." },
  { Cif2Cif, 0, "", "cif2cif", Arg::None,
    "  --cif2cif  \tRead CIF not JSON." },
  CifModUsage[SkipCat],
  CifModUsage[SortCif],

  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  const char* input = p.nonOption(0);
  const char* output = p.nonOption(1);
  namespace cif = gemmi::cif;

  if (p.options[Verbose])
    std::cerr << "Transcribing " << input << " to cif ..." << std::endl;
  try {
    cif::Document doc = p.options[Cif2Cif] ? gemmi::read_cif_gz(input)
                                           : gemmi::read_mmjson_gz(input);
    apply_cif_doc_modifications(doc, p.options);
    gemmi::Ofstream os(output, &std::cout);
    write_cif_to_stream(os.ref(), doc, cif_write_options(p.options[CifStyle]));
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
