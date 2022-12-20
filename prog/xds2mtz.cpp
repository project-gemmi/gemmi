// Copyright 2022 Global Phasing Ltd.
//
// Convert reflection data from XDS_ASCII to MTZ.

#include <cmath>              // for ceil
#include <cstdio>             // for fprintf
#include <set>
#ifndef GEMMI_ALL_IN_ONE
# define GEMMI_WRITE_IMPLEMENTATION 1
#endif
#include <gemmi/gz.hpp>        // for MaybeGzipped
#include <gemmi/xds_ascii.hpp> // for XdsAscii
#include <gemmi/mtz.hpp>       // for Mtz
#include <gemmi/version.hpp>   // for GEMMI_VERSION

#define GEMMI_PROG xds2mtz
#include "options.h"

namespace {

enum OptionIndex {
  Title=4, History
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n  " EXE_NAME " [options] XDS_FILE MTZ_FILE"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Title, 0, "", "title", Arg::Required,
    "  --title  \tMTZ title." },
  { History, 0, "-H", "history", Arg::Required,
    "  -H LINE, --history=LINE  \tAdd a history line." },
  { NoOp, 0, "", "", Arg::None,
    "\nIf XDS_FILE is -, the input is read from stdin."
  },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);

  gemmi::XdsAscii xds_ascii;
  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path);
  try {
    xds_ascii.read_input(gemmi::MaybeGzipped(input_path));
    gemmi::Mtz mtz;
    if (const option::Option* opt = p.options[Title])
      mtz.title = opt->arg;
    else
      mtz.title = "Converted from XDS_ASCII";
    if (const option::Option* opt = p.options[History])
      for (; opt; opt = opt->next())
        mtz.history.emplace_back(opt->arg);
    else
      mtz.history.emplace_back("From gemmi-xds2mtz " GEMMI_VERSION);
    mtz.cell = xds_ascii.unit_cell;
    mtz.spacegroup = gemmi::find_spacegroup_by_number(xds_ascii.spacegroup_number);
    mtz.add_base();
    mtz.add_column("M/ISYM", 'Y', 0, -1, false);
    mtz.add_column("BATCH", 'B', 0, -1, false);
    mtz.add_column("I", 'B', 0, -1, false);
    mtz.add_column("SIGI", 'Q', 0, -1, false);
    mtz.add_column("XDET", 'R', 0, -1, false);
    mtz.add_column("YDET", 'R', 0, -1, false);
    mtz.add_column("ROT", 'R', 0, -1, false);
    mtz.nreflections = (int) xds_ascii.data.size();
    mtz.data.resize(mtz.columns.size() * xds_ascii.data.size());
    gemmi::UnmergedHklMover hkl_mover(mtz.spacegroup);
    std::set<int> frames;
    size_t k = 0;
    for (const gemmi::XdsAscii::Refl& refl : xds_ascii.data) {
      auto hkl = refl.hkl;
      int isym = hkl_mover.move_to_asu(hkl);
      for (size_t j = 0; j != 3; ++j)
        mtz.data[k++] = (float) hkl[j];
      mtz.data[k++] = (float) isym;
      float frame = std::ceil((float) refl.zd);
      mtz.data[k++] = frame;
      frames.insert((int) frame);
      mtz.data[k++] = (float) refl.iobs;
      mtz.data[k++] = (float) refl.sigma;
      mtz.data[k++] = (float) refl.xd;
      mtz.data[k++] = (float) refl.yd;
      mtz.data[k++] = (float) refl.zd;
    }
    for (int frame : frames) {
      gemmi::Mtz::Batch batch;
      batch.number = frame;
      batch.set_dataset_id(0);
      batch.set_cell(mtz.cell);
      //batch.set_wavelength(iset.wavelength);
      mtz.batches.push_back(batch);
    }
    if (verbose)
      std::fprintf(stderr, "Writing %d reflections to %s ...\n",
                   mtz.nreflections, output_path);
    mtz.write_to_file(output_path);
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
