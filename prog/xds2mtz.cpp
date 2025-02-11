// Copyright 2022 Global Phasing Ltd.
//
// Convert reflection data from XDS_ASCII to MTZ.

#include <gemmi/xds2mtz.hpp>   // for xds_to_mtz
#include <cstdio>              // for fprintf
#include <gemmi/gz.hpp>        // for MaybeGzipped
#include <gemmi/xds_ascii.hpp> // for XdsAscii
#include <gemmi/mtz.hpp>       // for Mtz
#include <gemmi/version.hpp>   // for GEMMI_VERSION

#define GEMMI_PROG xds2mtz
#include "options.h"

namespace {

enum OptionIndex {
  Title=4, History, Project, Crystal, Dataset, Batchmin, Polarization, Normal, Overload
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
  { Project, 0, "", "project", Arg::Required,
    "  --project=PROJECT  \tProject in MTZ hierarchy (default: 'XDSproject')" },
  { Crystal, 0, "", "crystal", Arg::Required,
    "  --crystal=CRYSTAL  \tCrystal in MTZ hierarchy (default: 'XDScrystal')" },
  { Dataset, 0, "", "dataset", Arg::Required,
    "  --dataset=DATASET  \tDataset in MTZ hierarchy (default: 'XDSdataset')" },
  { Batchmin, 0, "", "batchmin", Arg::Int,
    "  --batchmin=BATCHMIN  \tDelete reflections with BATCH<BATCHMIN (default: 1)" },
  { NoOp, 0, "", "", Arg::None,
    "\nPolarization correction and overload elimination options for INTEGRATE.HKL files:" },
  { Polarization, 0, "", "polarization", Arg::Float,
    "  --polarization=VALUE  \tXDS parameter FRACTION_OF_POLARIZATION" },
  { Normal, 0, "", "normal", Arg::Float3,
    "  --normal='Pnx Pny Pnz'  \tXDS POLARIZATION_PLANE_NORMAL (default: '0 1 0')" },
  { Overload, 0, "", "overload", Arg::Float,
    "  --overload=OVERLOAD  \tXDS parameter OVERLOAD to eliminate reflections with MAXC>OVERLOAD" },
  { NoOp, 0, "", "", Arg::None,
    "\nIf XDS_FILE is -, the input is read from stdin." },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  if (p.options[Normal] && !p.options[Polarization]) {
    std::fprintf(stderr, "Error. Option -%s without -%s.\n",
                 p.given_name(Normal), p.given_name(Polarization));
    return 1;
  }
  bool verbose = p.options[Verbose];
  const char* input_path = p.nonOption(0);
  const char* output_path = p.nonOption(1);

  gemmi::XdsAscii xds;
  if (verbose)
    std::fprintf(stderr, "Reading %s ...\n", input_path);
  try {
    xds.read_input(gemmi::MaybeGzipped(input_path));

    // batchmin handling
    int batchmin = 1;
    if (p.options[Batchmin])
      batchmin = std::atoi(p.options[Batchmin].arg);
    size_t size_before = xds.data.size();
    xds.eliminate_batchmin(batchmin);
    size_t nbatchmin = size_before - xds.data.size();
    if (verbose || nbatchmin != 0)
      std::printf("Number of deleted reflections with BATCH < %d (i.e. ZD < %d) = %zu\n",
                  batchmin, batchmin-1, nbatchmin);

    // polarization correction
    if (p.options[Polarization]) {
      if (xds.generated_by != "INTEGRATE") {
        std::fprintf(stderr,
                     "Error: --polarization given for data from %s (not from INTEGRATE).\n",
                     xds.generated_by.c_str());
        return 1;
      }
      if (gemmi::likely_in_house_source(xds.wavelength))
        std::fprintf(stderr, "WARNING: likely in-house source (wavelength %g)\n"
                             "         polarization correction can be inappropriate.\n",
                     xds.wavelength);
      if (verbose)
        std::fprintf(stderr, "Applying polarization correction...\n");
      double fraction = std::atof(p.options[Polarization].arg);
      gemmi::Vec3 pn(0., 1., 0.);
      if (p.options[Normal]) {
        auto v = parse_blank_separated_numbers(p.options[Normal].arg);
        pn = gemmi::Vec3(v[0], v[1], v[2]);
      }
      xds.apply_polarization_correction(fraction, pn);
    }

    // overload handling
    if (p.options[Overload]) {
      if (xds.generated_by != "INTEGRATE") {
        std::fprintf(stderr,
                     "Error: --overload given for data from %s (not from INTEGRATE).\n",
                     xds.generated_by.c_str());
        return 1;
      }
      if (verbose)
        std::fprintf(stderr, "Eliminating overloads...\n");
      double overload = std::atof(p.options[Overload].arg);
      size_before = xds.data.size();
      xds.eliminate_overloads(overload);
      size_t nover = size_before - xds.data.size();
      std::printf("Number of eliminated reflections with MAXC > %g = %zu\n", overload, nover);
    }

    bool merged = xds.is_merged();
    if (verbose && merged)
      std::fprintf(stderr, "Preparing merged MTZ file...\n");

    gemmi::Mtz mtz = gemmi::xds_to_mtz(xds);

    if (const option::Option* opt = p.options[Title])
      mtz.title = opt->arg;
    else
      mtz.title = "Converted from " + gemmi::path_basename(input_path, {});
    std::string xds_info = gemmi::cat("From ", xds.generated_by, ' ', xds.version_str);
    if (const option::Option* opt = p.options[History]) {
      for (; opt; opt = opt->next()) {
        mtz.history.emplace_back(opt->arg);
        if (mtz.history.back() == "/xds/")  // special option for -H
          mtz.history.back().swap(xds_info);
      }
    } else {
      mtz.history.emplace_back("From gemmi-xds2mtz " GEMMI_VERSION);
      mtz.history.push_back(std::move(xds_info));
    }

    if (const option::Option* opt = p.options[Project])
      for (size_t i = 1; i < mtz.datasets.size(); ++i)
        mtz.datasets[i].project_name = opt->arg;
    if (const option::Option* opt = p.options[Crystal])
      for (size_t i = 1; i < mtz.datasets.size(); ++i)
        mtz.datasets[i].crystal_name = opt->arg;
    if (const option::Option* opt = p.options[Dataset])
      for (size_t i = 1; i < mtz.datasets.size(); ++i)
        mtz.datasets[i].dataset_name = opt->arg;

    if (verbose)
      std::fprintf(stderr, "Writing %d reflections to (%smerged) %s ...\n",
                   mtz.nreflections, (merged ? "" : "un"), output_path);
    mtz.write_to_file(output_path);
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
