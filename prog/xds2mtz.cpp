// Copyright 2022 Global Phasing Ltd.
//
// Convert reflection data from XDS_ASCII to MTZ.

#include <cstdio>             // for fprintf
#include <map>
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

    gemmi::Mtz mtz;
    if (const option::Option* opt = p.options[Title])
      mtz.title = opt->arg;
    else
      mtz.title = "Converted from " + gemmi::path_basename(input_path, {});
    if (const option::Option* opt = p.options[History]) {
      for (; opt; opt = opt->next())
        mtz.history.emplace_back(opt->arg);
    } else {
      mtz.history.emplace_back("From gemmi-xds2mtz " GEMMI_VERSION);
      mtz.history.push_back(gemmi::cat("From ", xds.generated_by, ' ', xds.version_str));
    }
    mtz.cell = xds.unit_cell;
    mtz.spacegroup = gemmi::find_spacegroup_by_number(xds.spacegroup_number);
    mtz.add_base();
    const char* pxd[3] = {"XDSproject", "XDScrystal", "XDSdataset"};
    if (const option::Option* opt = p.options[Project])
      pxd[0] = opt->arg;
    if (const option::Option* opt = p.options[Crystal])
      pxd[1] = opt->arg;
    if (const option::Option* opt = p.options[Dataset])
      pxd[2] = opt->arg;
    mtz.datasets.push_back({1, pxd[0], pxd[1], pxd[2], mtz.cell, xds.wavelength});
    mtz.add_column("M/ISYM", 'Y', 0, -1, false);
    mtz.add_column("BATCH", 'B', 0, -1, false);
    mtz.add_column("I", 'J', 0, -1, false);
    mtz.add_column("SIGI", 'Q', 0, -1, false);
    mtz.add_column("XDET", 'R', 0, -1, false);
    mtz.add_column("YDET", 'R', 0, -1, false);
    mtz.add_column("ROT", 'R', 0, -1, false);
    if (xds.read_columns >= 11) {
      mtz.add_column("FRACTIONCALC", 'R', 0, -1, false);
      mtz.add_column("LP", 'R', 0, -1, false);
      mtz.add_column("CORR", 'R', 0, -1, false);
      if (xds.read_columns > 11)
        mtz.add_column("MAXC", 'I', 0, -1, false);
    }
    mtz.add_column("FLAG", 'I', 0, -1, false);
    mtz.nreflections = (int) xds.data.size();
    mtz.data.resize(mtz.columns.size() * xds.data.size());
    gemmi::UnmergedHklMover hkl_mover(mtz.spacegroup);
    int max_frame = 0;
    for (const gemmi::XdsAscii::Refl& refl : xds.data)
      max_frame = std::max(max_frame, refl.frame());
    int iset_offset = (max_frame + 11000) / 10000 * 10000;
    // iset,frame -> batch
    std::map<std::pair<int,int>, int> frames;
    size_t k = 0;
    for (const gemmi::XdsAscii::Refl& refl : xds.data) {
      auto hkl = refl.hkl;
      int isym = hkl_mover.move_to_asu(hkl);
      for (size_t j = 0; j != 3; ++j)
        mtz.data[k++] = (float) hkl[j];
      mtz.data[k++] = (float) isym;
      int frame = refl.frame();
      int batch = frame + iset_offset * std::max(refl.iset - 1, 0);
      frames.emplace(std::make_pair(refl.iset, frame), batch);
      mtz.data[k++] = (float) batch;
      mtz.data[k++] = (float) refl.iobs;  // I
      mtz.data[k++] = (float) std::fabs(refl.sigma);  // SIGI
      mtz.data[k++] = (float) refl.xd;
      mtz.data[k++] = (float) refl.yd;
      mtz.data[k++] = (float) xds.rot_angle(refl);  // ROT
      if (xds.read_columns >= 11) {
        mtz.data[k++] = float(0.01 * refl.peak);  // FRACTIONCALC
        mtz.data[k++] = (float) refl.rlp;
        mtz.data[k++] = float(0.01 * refl.corr);
        if (xds.read_columns > 11)
          mtz.data[k++] = (float) refl.maxc;
      }
      mtz.data[k++] = refl.sigma < 0 ? 64.f : 0.f;  // FLAG
    }
    // Prepare a similar batch header as Pointless.
    gemmi::Mtz::Batch batch;
    batch.set_dataset_id(1);

    // We don't set lbcell (refinement flags for unit cell),
    // because it's probably not used by any program anyway.
    // batch.ints[4] to [9] = left unset

    // We also skip jumpax, which is defined as:
    // reciprocal axis closest to principle goniostat axis E1
    // batch.ints[11] = left unset

    // We assume one crystal, one goniostat axis, one detector.
    batch.ints[12] = 1;  // ncryst
    batch.ints[14] = 2;  // ldtype 3D
    batch.ints[15] = 1;  // jsaxs - goniostat scan axis number
    batch.ints[17] = 1;  // ngonax - number of goniostat axes
    batch.ints[19] = 1;  // ndet

    batch.set_cell(mtz.cell);  // batch.floats[0] to [5]
    gemmi::Vec3 s0(-1, 0, 0); // will be re-set if we have geometry info
    try {
      gemmi::Mat33 Q = xds.calculate_conversion_from_cambridge().inverse();
      s0 = -Q.multiply(xds.get_s0_direction());
      // Orientation matrix U. It is calculated differently in Pointless,
      // so the results are slightly different (due to limited precision
      // of numbers in XDS file).
      gemmi::Mat33 U = Q.multiply(xds.get_orientation());
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          batch.floats[6 + 3*i + j] = (float) U[j][i];
    } catch (std::runtime_error& e) {
      // if some of the headers are absent, U is not set
      // and s0 is set to "idealised" s0
      std::fprintf(stderr, "Note: %s, orientation matrix U not set.\n", e.what());
    }
    batch.floats[21] = float(xds.reflecting_range_esd);  // crydat(0)
    // In the so-called "Cambridge" frame (as used by Mosflm),
    // the principal rotation axis is along z
    // and the incident beam S0 is along x.
    // Therefore, we set the rotation axis phi (scanax) to:
    batch.floats[38+2] = 1.f;  // scanax = [0, 0, 1]
    batch.floats[47] = float(xds.oscillation_range);  // phi range
    // E1,E2,E3 vectors are the goniostat axes in Cambridge laboratory frame.
    // E1 is set to rotation axis. E2 and E3 are not set for ngonax==1.
    batch.floats[59+2] = 1.f;  // e1 = scanax

    // Idealised source vector is -x ([-1 0 0]), antiparallel to beam.
    batch.floats[80+0] = -1.f;  // source[0]

    // s0 source vector (including tilts), CMtz::MTZBAT::so in libccp4
    batch.floats[83] = (float) s0.x;
    batch.floats[84] = (float) s0.y;
    batch.floats[85] = (float) s0.z;

    batch.set_wavelength((float)xds.wavelength);  // batch.floats[86]
    // or batch.set_wavelength(iset.wavelength);
    // Detector geometry.
    batch.floats[111] = (float) xds.detector_distance;  // dx[0]
    batch.floats[113] = 1.f;  // detlm[0][0][0]
    batch.floats[114] = (float) xds.nx;
    batch.floats[115] = 1.f;
    batch.floats[116] = (float) xds.ny;
    batch.axes.push_back("PHI");  // gonlab[0]

    for (auto& t : frames) {
      batch.set_dataset_id(t.first.first);
      batch.number = t.second;
      int frame = t.first.second;
      double phistt = xds.starting_angle +
                      xds.oscillation_range * (frame - xds.starting_frame);
      batch.floats[36] = float(phistt);
      batch.floats[37] = float(phistt + xds.oscillation_range);  // phiend
      mtz.batches.push_back(batch);
    }
    mtz.sort(5);
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
