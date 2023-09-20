// Copyright 2019 Global Phasing Ltd.
//
// convert MTZ to SF-mmCIF

#include <cstdio>
#include <cstdlib>            // for strtod
#include <iostream>           // for cout, cerr
#include <gemmi/mtz2cif.hpp>
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/fstream.hpp>  // for Ofstream
#include <gemmi/util.hpp>     // for giends_with
#include <gemmi/merge.hpp>    // for Intensities
#include <gemmi/read_cif.hpp> // for read_cif_gz
#include <gemmi/xds_ascii.hpp>
#define GEMMI_PROG mtz2cif
#include "options.h"

namespace {

enum OptionIndex {
  Spec=4, PrintSpec, BlockName, EntryId, SkipEmpty, SkipNegativeSigI,
  NoComments, NoHistory, NoStaranisoTensor, RunFrom, Wavelength, Validate,
  LessAnomalous, Separate, Deposition, NoIntensityCheck, Nfree, Trim
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] MTZ_FILE [MTZ_FILE] CIF_FILE"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Spec, 0, "", "spec", Arg::Required,
    "  --spec=FILE  \tColumn and format specification." },
  { PrintSpec, 0, "", "print-spec", Arg::None,
    "  --print-spec  \tPrint default spec and exit." },
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \tmmCIF block name: data_NAME (default: mtz)." },
  { EntryId, 0, "", "id", Arg::Required,
    "  --id=ID  \tvalue for _entry.id (default: xxxx)." },
  { SkipEmpty, 0, "", "skip-empty", Arg::Optional,
    "  --skip-empty[=COLS]  \tSkip reflections with no values. If COLS are "
    "given, eg. 'I(+),I(-)', only values in those columns are checked." },
  { SkipNegativeSigI, 0, "", "skip-negative-sigi", Arg::Optional,
    "  --skip-negative-sigi  \tSkip reflections with sigma(I)<0 in unmerged data." },
  { NoComments, 0, "", "no-comments", Arg::None,
    "  --no-comments  \tDo not write comments in the mmCIF file." },
  { NoHistory, 0, "", "no-history", Arg::None,
    "  --no-history  \tDo not write MTZ history in the mmCIF file." },
  { NoStaranisoTensor, 0, "", "no-staraniso-tensor", Arg::None,
    "  --no-staraniso-tensor  \tDo not write _reflns.pdbx_aniso_B_tensor_*." },
  { RunFrom, 0, "", "run-from", Arg::Required,
    "  --run-from=WRAPPER  \tAdd note in _software.description." },
  { Wavelength, 0, "", "wavelength", Arg::Float,
    "  --wavelength=LAMBDA  \tSet wavelength (default: from input file)." },
  { Validate, 0, "", "validate", Arg::None,
    "  --validate  \tFor two MTZ files: validate the intensities match." },
  { LessAnomalous, 0, "", "less-ano", Arg::None,
    "  --less-ano  \tSkip anomalous columns (even if they are in the spec)."
    " Used once, skips I(+)/I(-) if <I> and F(+)/F(-) are present."
    " Used twice, skips all anomalous columns."
  },
  { Separate, 0, "", "separate", Arg::None,
    "  --separate  \tWrite merged and unmerged data in separate blocks." },
  { Deposition, 0, "", "depo", Arg::None,
    "  --depo  \tPrepare merged+unmerged mmCIF file for deposition." },
  { Nfree, 0, "", "nfree", Arg::Int,
    "  --nfree=N  \tFlag value used for the free set (default: auto)" },
  { Trim, 0, "", "trim", Arg::Int,
    "  --trim=N  \t(for testing) output only reflections -N <= h,k,l <=N." },
  { NoIntensityCheck, 0, "", "no-intensity-check", Arg::None, 0 },
  { NoOp, 0, "", "", Arg::None,
    "\nOne or two MTZ files are taken as the input. If two files are given,"
    "\none must be merged and the other unmerged."
    "\nAlternatively, the two input files can be reflection mmCIF file and"
    "\nunmerged MTZ - this adds unmerged data block to the existing mmCIF file."
    "\nXDS_ASCII.HKL file can be used instead of unmerged MTZ file (also when"
    "\nconverting a single file), but then the spec file is ignored."
    "\nIf CIF_FILE is -, the output is printed to stdout."
    "\nIf spec is -, it is read from stdin."
    "\n"
    "\nLines in the spec file have format:"
    "\n  [FLAG] COLUMN TYPE TAG [FORMAT]"
    "\n or"
    "\n  $SPECIAL TAG"
    "\nfor example:"
    "\n  SIGF_native * SIGF_meas_au 12.5e"
    "\n  FREE I pdbx_r_free_flag 3.0f"
    "\n  $counter id"
    "\nFLAG (optional) is either ? or &:"
    "\n  ? = ignored if no column in the MTZ file has this name."
    "\n  & = ignored if the previous line was ignored."
    "\n  Example:"
    "\n      ? I    J intensity_meas"
    "\n      & SIGI Q intensity_sigma"
    "\nCOLUMN is MTZ column label. Columns H K L are added if not specified."
    "\n  Alternative labels can be separated with | (e.g. FREE|FreeR_flag)."
    "\nTYPE is used for checking the column type, unless it is '*'."
    "\nTAG does not include category name, it is only the part after _refln."
    "\nFORMAT (optional) is printf-like floating-point format:"
    "\n - one of e, f, g with optional flag, width and precision"
    "\n - flag is one of + - # _; '_' stands for ' ', for example '_.4f'"
    "\n - since all numbers in MTZ are stored as float, the integer columns use"
    "\n   the same format as float. The format of _refln.status is ignored."
    "\n$SPECIAL is $counter, $dataset, $image, $. or $?, as appropriate for"
    "\n  _diffrn_refln.id, _diffrn_refln.diffrn_id, image number or null (./?)"
  },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[PrintSpec]) {
    std::printf("                                for merged mtz\n");
    const char** lines = gemmi::MtzToCif::default_spec(/*for_merged=*/true);
    for (; *lines != nullptr; ++lines)
      std::printf("%s\n", *lines);
    std::printf("\n                             for unmerged mtz\n");
    lines = gemmi::MtzToCif::default_spec(/*for_merged=*/false);
    for (; *lines != nullptr; ++lines)
      std::printf("%s\n", *lines);
    return 0;
  }

  int nargs = p.nonOptionsCount();
  if (nargs != 2 && nargs != 3) {
    fprintf(stderr, "%s requires 2 or 3 arguments, got %d.", p.program_name, nargs);
    p.print_try_help_and_exit("");
  }
  for (auto opt : {Deposition, Validate})
    if (p.options[opt] && nargs != 3) {
      fprintf(stderr, "Option -%s works only with 2 input files.\n", p.given_name(opt));
      return 1;
    }
  bool verbose = p.options[Verbose];
  const char* mtz_paths[2];
  const char* xds_path = nullptr;
  const char* cif_input = nullptr;
  mtz_paths[0] = p.nonOption(0);
  mtz_paths[1] = (nargs == 3 ? p.nonOption(1) : nullptr);
  const char* cif_output = p.nonOption(nargs == 3 ? 2 : 1);
  std::unique_ptr<gemmi::Mtz> mtz[2];
  std::unique_ptr<gemmi::XdsAscii> xds_ascii;
  if (gemmi::giends_with(mtz_paths[0], ".cif") ||
      gemmi::giends_with(mtz_paths[0], ".mmcif") ||
      gemmi::giends_with(mtz_paths[0], ".ent")) {
    if (!mtz_paths[1]) {
      std::fprintf(stderr, "Error: no MTZ file was given\n");
      return 1;
    }
    std::swap(cif_input, mtz_paths[0]);
  }
  if (mtz_paths[nargs-2] && gemmi::giends_with(mtz_paths[nargs-2], ".hkl"))
    std::swap(xds_path, mtz_paths[nargs-2]);
  if (gemmi::giends_with(cif_output, ".mtz")) {
    std::fprintf(stderr, "This must be a mistake, requested output cif file"
                         " has mtz extension: %s\n", cif_output);
    return 1;
  }
  for (int i = 0; i < 2; ++i)
    if (mtz_paths[i]) {
      mtz[i].reset(new gemmi::Mtz);
      if (verbose) {
        std::fprintf(stderr, "Reading %s ...\n", mtz_paths[i]);
        mtz[i]->warnings = &std::cerr;
      }
      try {
        mtz[i]->read_input(gemmi::MaybeGzipped(mtz_paths[i]), true);
      } catch (std::runtime_error& e) {
        std::fprintf(stderr, "ERROR reading %s: %s\n", mtz_paths[i], e.what());
        return 1;
      }
    }
  if (xds_path) {
    try {
      xds_ascii.reset(new gemmi::XdsAscii);
      xds_ascii->read_input(gemmi::MaybeGzipped(xds_path));
      xds_ascii->gather_iset_statistics();
    } catch (std::runtime_error& e) {
      std::fprintf(stderr, "ERROR reading %s: %s\n", xds_path, e.what());
      return 1;
    }
  }
  if (cif_input && mtz[1] && mtz[1]->is_merged()) {
    std::fprintf(stderr, "Error: CIF file and merged MTZ files given\n");
    return 1;
  }
  if (xds_ascii && mtz[0] && !mtz[0]->is_merged()) {
    std::fprintf(stderr, "Error: Two unmerged files (MTZ and XDS_ASCII) given\n");
    return 1;
  }
  if (mtz[0] && mtz[1]) {
    if (mtz[0]->is_merged() && mtz[1]->is_merged()) {
      std::fprintf(stderr, "Error: two merged MTZ files given\n");
      return 1;
    }
    if (!mtz[0]->is_merged() && !mtz[1]->is_merged()) {
      std::fprintf(stderr, "Error: two unmerged MTZ files given\n");
      return 1;
    }
    if (mtz[1]->is_merged())
      mtz[0].swap(mtz[1]);
  }

  bool separate_blocks = p.options[Separate];

  gemmi::MtzToCif mtz_to_cif;

  if (p.options[Spec]) {
    try {
      read_spec_file(p.options[Spec].arg, mtz_to_cif.spec_lines);
    } catch (std::runtime_error& e) {
      std::fprintf(stderr, "Problem in translation spec: %s\n", e.what());
      return 2;
    }
  }

  mtz_to_cif.with_comments = !p.options[NoComments];
  mtz_to_cif.with_history = !p.options[NoHistory];
  mtz_to_cif.write_staraniso_tensor = !p.options[NoStaranisoTensor];
  if (p.options[RunFrom])
    mtz_to_cif.gemmi_run_from = p.options[RunFrom].arg;
  mtz_to_cif.less_anomalous = p.options[LessAnomalous].count();
  if (p.options[Nfree])
    mtz_to_cif.free_flag_value = std::atoi(p.options[Nfree].arg);
  bool validate = p.options[Validate];
  bool check_merged_columns = false;
  if (p.options[Deposition]) {
    mtz_to_cif.write_special_marker_for_pdb = true;
    mtz_to_cif.skip_negative_sigi = true;
    mtz_to_cif.with_history = false;
    mtz_to_cif.less_anomalous = std::max(mtz_to_cif.less_anomalous, 1);
    separate_blocks = true;
    validate = true;
    check_merged_columns = true;
  }
  if (p.options[SkipEmpty]) {
    mtz_to_cif.skip_empty = true;
    if (p.options[SkipEmpty].arg)
      mtz_to_cif.skip_empty_cols = p.options[SkipEmpty].arg;
  }
  if (p.options[SkipNegativeSigI])
    mtz_to_cif.skip_negative_sigi = true;
  if (p.options[BlockName])
    mtz_to_cif.block_name = p.options[BlockName].arg;
  if (p.options[EntryId])
    mtz_to_cif.entry_id = p.options[EntryId].arg;
  if (p.options[Wavelength])
    mtz_to_cif.wavelength = std::strtod(p.options[Wavelength].arg, nullptr);

  gemmi::CharArray cif_buf;
  if (cif_input) {
    try {
      cif_buf = gemmi::read_into_buffer_gz(cif_input);
    } catch (std::runtime_error& e) {
      fprintf(stderr, "ERROR: %s\n", e.what());
      return 1;
    }
  }

  bool ok = true;
  if (mtz[0] && mtz_to_cif.spec_lines.empty())
    remove_appendix_from_column_names(*mtz[0], std::cerr);
  if (check_merged_columns && mtz[0])
    ok = gemmi::validate_merged_mtz_deposition_columns(*mtz[0], std::cerr);
  gemmi::Intensities mi;
  if (mtz[0])
    mtz_to_cif.staraniso_version = mi.take_staraniso_b_from_mtz(*mtz[0]);
  if (validate && nargs == 3) {
    try {
      if (mtz[0]) {
        if (!mtz_to_cif.staraniso_version.empty()) {
          std::cerr << "Merged MTZ went through STARANISO " << mtz_to_cif.staraniso_version;
          if (mi.staraniso_b.ok())
            std::cerr << ". Taking into account anisotropic scaling.\n";
          else
            std::cerr << ". B tensor is unknown. Intensities won't be checked.\n";
        }
        mi.read_merged_intensities_from_mtz(*mtz[0]);
      } else {
        gemmi::ReflnBlock rblock = gemmi::get_refln_block(
            gemmi::read_cif_from_buffer(cif_buf, cif_input).blocks, {});
        if (!rblock.entry_id.empty() && rblock.entry_id != mtz_to_cif.entry_id) {
          fprintf(stderr, "Note: using _entry.id from mmCIF (%s) instead of %s.\n",
                  rblock.entry_id.c_str(), mtz_to_cif.entry_id.c_str());
          mtz_to_cif.entry_id = rblock.entry_id;
        }
        mi.take_staraniso_b_from_mmcif(rblock.block);
        mi.read_merged_intensities_from_mmcif(rblock);
      }
      gemmi::Intensities ui;
      if (mtz[1])
        ui.read_unmerged_intensities_from_mtz(*mtz[1]);
      else if (xds_ascii)
        ui.read_unmerged_intensities_from_xds(*xds_ascii);

      bool relaxed_check =
        p.options[NoIntensityCheck] ||
        // If an old StarAniso version was used that doesn't store B tensor,
        // allow intensities to differ.
        (!mtz_to_cif.staraniso_version.empty() && !mi.staraniso_b.ok());

      if (!gemmi::validate_merged_intensities(mi, ui, relaxed_check, std::cerr))
        ok = false;
    } catch (std::exception& e) {
      fprintf(stderr, "Error. Intensities could not be validated.\n%s.\n", e.what());
      ok = false;
    }
  }
  if (!ok) {
    fprintf(stderr, "\nIf you think the files are correct, contact us:\n"
            "see https://project-gemmi.github.io/depo.html\n");
    return 4;
  }

  if (p.options[Trim])
    mtz_to_cif.trim = std::atoi(p.options[Trim].arg);
  if (verbose)
    std::fprintf(stderr, "Opening output file %s ...\n", cif_output);
  try {
    gemmi::Ofstream os(cif_output, &std::cout);
    for (int i = 0; i < 2; ++i)
      if (mtz[i] && !mtz[i]->is_merged())
        mtz[i]->switch_to_original_hkl();
    gemmi::SMat33<double>* staraniso_b = mi.staraniso_b.ok() ? &mi.staraniso_b.b : nullptr;
    if (mtz[0] && mtz[1] && !separate_blocks) {
      if (verbose)
        std::fprintf(stderr, "Writing merged and unmerged MTZ data as mmCIF block...\n");
      mtz_to_cif.write_cif(*mtz[0], mtz[1].get(), staraniso_b, os.ref());
    } else {
      if (cif_input) {
        if (verbose)
          std::fprintf(stderr, "Copying input mmCIF file to output...\n");
        os.ref().write(cif_buf.data(), cif_buf.size());
      } else if (mtz[0]) {
        if (verbose)
          std::fprintf(stderr, "Writing merged MTZ data as mmCIF block...\n");
        mtz_to_cif.write_cif(*mtz[0], nullptr, staraniso_b, os.ref());
      }
      if ((cif_input || mtz[0]) && (mtz[1] || xds_ascii))
        os.ref() << "\n\n";
      mtz_to_cif.block_name = "unmerged";
      if (mtz[1]) {
        if (verbose)
          std::fprintf(stderr, "Writing unmerged MTZ data as mmCIF block...\n");
        mtz_to_cif.write_cif(*mtz[1], nullptr, nullptr, os.ref());
      } else if (xds_ascii) {
        if (verbose)
          std::fprintf(stderr, "Writing unmerged XDS data as mmCIF block...\n");
        mtz_to_cif.write_cif_from_xds(*xds_ascii, os.ref());
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR writing %s: %s\n", cif_output, e.what());
    return 3;
  }
  return 0;
}
