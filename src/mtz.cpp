// Copyright 2019 Global Phasing Ltd.
//
// MTZ info

#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#define GEMMI_PROG mtz
#include "options.h"
#include <stdio.h>

using namespace gemmi;

enum OptionIndex { Verbose=3, Headers, Dump, ToggleEndian };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] MTZ_FILE[...]"
    "\nPrint informations from an mtz file."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Headers, 0, "H", "headers", Arg::None,
    "  -H, --headers  \tPrint raw headers, until the END record." },
  { Dump, 0, "d", "dump", Arg::None,
    "  -d, --dump  \tPrint a subset of CCP4 mtzdmp informations." },
  { ToggleEndian, 0, "", "toggle-endian", Arg::None,
    "  --toggle-endian  \tToggle assumed endiannes (little <-> big)." },
  { 0, 0, 0, 0, 0, 0 }
};


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  bool verbose = p.options[Verbose];


  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      const char* path = p.nonOption(i);
      if (i != 0)
        std::printf("\n\n");
      if (verbose)
        std::fprintf(stderr, "Reading %s ...\n", path);
      gemmi::fileptr_t f = gemmi::file_open(path, "rb");
      Mtz mtz;
      try {
        mtz.read_first_bytes(f.get());
        if (p.options[ToggleEndian])
          mtz.toggle_endiannes();
      } catch (std::runtime_error& e) {
        fail(std::string(e.what()) + ": " + path);
      }
      if (p.options[Headers]) {
        char buf[81] = {0};
        mtz.seek_headers(f.get());
        while (std::fread(buf, 1, 80, f.get()) != 0) {
          std::printf("%s\n", gemmi::rtrim_str(buf).c_str());
          if (gemmi::Mtz::id3u(buf) == gemmi::Mtz::id3("END"))
            break;
        }
      }
      mtz.read_headers(f.get(), verbose ? stderr : nullptr);
      if (p.options[Dump]) {
        std::printf("Title: %s\n", mtz.title.c_str());
        std::printf("Number of Columns = %d\n", mtz.ncol);
        std::printf("Number of Reflections = %d\n", mtz.nreflections);
        if (mtz.nbatches != 0)
          std::printf("Number of Batches = %d\n", mtz.nbatches);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
