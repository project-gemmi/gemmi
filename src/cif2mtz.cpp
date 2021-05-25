// Copyright 2019 Global Phasing Ltd.
//
// convert SF-mmCIF to MTZ

#include <cstdio>             // for fprintf
#include <cstdlib>            // for exit
#include <iostream>           // for cerr
#include <memory>             // for unique_ptr
#ifndef GEMMI_ALL_IN_ONE
# define GEMMI_WRITE_IMPLEMENTATION 1
#endif
#include <gemmi/read_cif.hpp> // for read_cif_gz
#include <gemmi/cif2mtz.hpp>  // for CifToMtz

#define GEMMI_PROG cif2mtz
#include "options.h"

namespace {

using std::fprintf;

enum OptionIndex { BlockName=4, Dir, Spec, PrintSpec, Title, History,
                   Unmerged, Sort };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n  " EXE_NAME " [options] CIF_FILE MTZ_FILE"
    "\n  " EXE_NAME " [options] CIF_FILE --dir=DIRECTORY"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \tmmCIF block to convert." },
  { Dir, 0, "d", "dir", Arg::Required,
    "  -d DIR, --dir=NAME  \tOutput directory." },
  { Spec, 0, "", "spec", Arg::Required,
    "  --spec=FILE  \tConversion spec." },
  { PrintSpec, 0, "", "print-spec", Arg::None,
    "  --print-spec  \tPrint default spec and exit." },
  { Title, 0, "", "title", Arg::Required,
    "  --title  \tMTZ title." },
  { History, 0, "-H", "history", Arg::Required,
    "  -H LINE, --history=LINE  \tAdd a history line." },
  { Unmerged, 0, "u", "unmerged", Arg::None,
    "  -u, --unmerged  \tWrite unmerged MTZ file(s)." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tOrder reflections according to Miller indices." },
  { NoOp, 0, "", "", Arg::None,
    "\nFirst variant: converts the first block of CIF_FILE, or the block"
    "\nspecified with --block=NAME, to MTZ file with given name."
    "\n\nSecond variant: converts each block of CIF_FILE to one MTZ file"
    "\n(block-name.mtz) in the specified DIRECTORY."
    "\n\nIf CIF_FILE is -, the input is read from stdin."
  },
  { 0, 0, 0, 0, 0, 0 }
};

gemmi::ReflnBlock& get_block_by_name(std::vector<gemmi::ReflnBlock>& rblocks,
                                     const std::string& name) {
  for (gemmi::ReflnBlock& rb : rblocks)
    if (rb.block.name == name)
      return rb;
  gemmi::fail("block not found: " + name);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[PrintSpec]) {
    std::printf("# Each line in the spec contains four words:\n"
                "# - tag (without category) from _refln or _diffrn_refln\n"
                "# - MTZ column label\n"
                "# - MTZ column type\n"
                "# - MTZ dataset for the column (must be 0 or 1)\n"
                "# The first 3 (5 in unmerged) columns are not in the spec,\n"
                "# they are always H K L (M/ISYM BATCH).\n\n");
    bool merged = !p.options[Unmerged];
    if (merged)
      std::printf("# For MERGED data only. Use --print-spec --unmerged for unmerged.\n");
    for (const char** line = gemmi::CifToMtz::default_spec(merged); *line != nullptr; ++line)
      std::printf("%s\n", *line);
    return 0;
  }
  bool convert_all = p.options[Dir];
  p.require_positional_args(convert_all ? 1 : 2);

  gemmi::CifToMtz cif2mtz;
  cif2mtz.verbose = p.options[Verbose];
  cif2mtz.force_unmerged = p.options[Unmerged];
  if (p.options[Title])
    cif2mtz.title = p.options[Title].arg;
  for (const option::Option* opt = p.options[History]; opt; opt = opt->next())
    cif2mtz.history.push_back(opt->arg);
  try {
    if (p.options[Spec])
      read_spec_file(p.options[Spec].arg, cif2mtz.spec_lines);
    const char* cif_path = p.nonOption(0);
    if (cif2mtz.verbose)
      fprintf(stderr, "Reading %s ...\n", cif_path);
    auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(cif_path).blocks);
    if (convert_all) {
      bool ok = true;
      for (gemmi::ReflnBlock& rb : rblocks) {
        std::string path = p.options[Dir].arg;
        path += '/';
        path += rb.block.name;
        path += ".mtz";
        try {
          gemmi::Mtz mtz = cif2mtz.convert_block_to_mtz(rb, std::cerr);
          if (cif2mtz.verbose)
            fprintf(stderr, "Writing %s ...\n", path.c_str());
          mtz.write_to_file(path);
        } catch (std::runtime_error& e) {
          fprintf(stderr, "ERROR: %s\n", e.what());
          ok = false;
        }
      }
      if (!ok)
        return 1;
    } else {
      const char* mtz_path = p.nonOption(1);
      const gemmi::ReflnBlock& rb = p.options[BlockName]
        ? get_block_by_name(rblocks, p.options[BlockName].arg)
        : rblocks.at(0);
      gemmi::Mtz mtz = cif2mtz.convert_block_to_mtz(rb, std::cerr);
      if (p.options[Sort]) {
        bool reordered = mtz.sort();
        if (cif2mtz.verbose)
          fprintf(stderr, "Reflection order has %schanged.\n", reordered ? "" : "not ");
      }
      if (cif2mtz.verbose)
        fprintf(stderr, "Writing %s ...\n", mtz_path);
      mtz.write_to_file(mtz_path);
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  if (cif2mtz.verbose)
    fprintf(stderr, "Done.\n");
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
