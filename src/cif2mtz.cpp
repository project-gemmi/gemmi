// Copyright 2019 Global Phasing Ltd.
//
// convert SF-mmCIF to MTZ

#include <algorithm>
#include <cstdlib>            // for exit
#include <stdio.h>
#ifndef GEMMI_ALL_IN_ONE
# define GEMMI_WRITE_IMPLEMENTATION 1
#endif
#include <gemmi/atox.hpp>     // for read_word
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/gzread.hpp>   // for read_cif_gz
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/refln.hpp>    // for ReflnBlock
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#define GEMMI_PROG cif2mtz
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { Verbose=3, BlockName, Dir, Title, History };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n  " EXE_NAME " [options] CIF_FILE MTZ_FILE"
    "\n  " EXE_NAME " [options] CIF_FILE --dir=DIRECTORY"
    "\nOptions:"},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \tmmCIF block to convert." },
  { Dir, 0, "d", "dir", Arg::Required,
    "  -d DIR, --dir=NAME  \tOutput directory." },
  { Title, 0, "", "title", Arg::Required,
    "  --title  \tMTZ title." },
  { History, 0, "-H", "history", Arg::Required,
    "  -H LINE, --history=LINE  \tAdd a history line." },
  { NoOp, 0, "", "", Arg::None,
    "\nFirst variant: converts the first block of CIF_FILE, or the block"
    "\nspecified with --block=NAME, to MTZ file with given name."
    "\n\nSecond variant: converts each block of CIF_FILE to one MTZ file"
    "\n(block-name.mtz) in the specified DIRECTORY."
    "\n\nIf CIF_FILE is -, the input is read from stdin."
  },
  { 0, 0, 0, 0, 0, 0 }
};

struct Entry {
  const char* refln_tag;
  const char* refln_tag2;
  const char* col_label;
  char col_type;
  unsigned char dataset_id;
};

static Entry translation_table[] = {
  {"index_h", nullptr,               "H",          'H', 0},
  {"index_k", nullptr,               "K",          'H', 0},
  {"index_l", nullptr,               "L",          'H', 0},
  {"pdbx_r_free_flag", "status",     "FreeR_flag", 'I', 0},
  {"intensity_meas", nullptr,        "I",          'J', 1},
  {"intensity_sigma", nullptr,       "SIGI",       'Q', 1},
  {"pdbx_I_plus", nullptr,           "I(+)",       'K', 1},
  {"pdbx_I_plus_sigma", nullptr,     "SIGI(+)",    'M', 1},
  {"pdbx_I_minus", nullptr,          "I(-)",       'K', 1},
  {"pdbx_I_minus_sigma", nullptr,    "SIGI(-)",    'M', 1},
  {"F_meas_au", nullptr,             "FP",         'F', 1},
  {"F_meas_sigma_au", nullptr,       "SIGFP",      'Q', 1},
  {"pdbx_F_plus", nullptr,           "F(+)",       'G', 1},
  {"pdbx_F_plus_sigma", nullptr,     "SIGF(+)",    'L', 1},
  {"pdbx_F_minus", nullptr,          "F(-)",       'G', 1},
  {"pdbx_F_minus_sigma", nullptr,    "SIGF(-)",    'L', 1},
  {"pdbx_anom_difference", nullptr,  "DP",         'D', 1},
  {"pdbx_anom_difference_sigma", nullptr, "SIGDP", 'Q', 1},
  {"F_calc", nullptr,                "FC",         'F', 1},
  {"phase_calc", nullptr,            "PHIC",       'P', 1},
  {"fom", "weight",                  "FOM",        'W', 1},
  {"pdbx_HL_A_iso", nullptr,         "HLA",        'A', 1},
  {"pdbx_HL_B_iso", nullptr,         "HLB",        'A', 1},
  {"pdbx_HL_C_iso", nullptr,         "HLC",        'A', 1},
  {"pdbx_HL_D_iso", nullptr,         "HLD",        'A', 1},
  {"pdbx_FWT", nullptr,              "FWT",        'F', 1},
  {"pdbx_PHWT", nullptr,             "PHWT",       'P', 1},
  {"pdbx_DELFWT", nullptr,           "DELFWT",     'F', 1},
  {"pdbx_DELPHWT", nullptr,          "DELPHWT",    'P', 1},
};

inline float status_to_freeflag(const std::string& str) {
  char c = str[0];
  if (c == '\'' || c == '"')
    c = str[1];
  if (c == 'o')
    return 1.f;
  if (c == 'f')
    return 0.f;
  return NAN;
}

static
void convert_cif_block_to_mtz(const gemmi::ReflnBlock& rb,
                              const std::string& mtz_path,
                              const std::vector<option::Option>& options) {
  gemmi::Mtz mtz;
  if (options[Title])
    mtz.title = options[Title].arg;
  for (const option::Option* opt = options[History]; opt; opt = opt->next())
    mtz.history.push_back(opt->arg);
  mtz.cell = rb.cell;
  mtz.spacegroup = rb.spacegroup;
  mtz.datasets.push_back({0, "HKL_base", "HKL_base", "HKL_base", mtz.cell, 0.});
  mtz.datasets.push_back({1, "unknown", "unknown", "unknown", mtz.cell,
                          rb.wavelength});
  if (!rb.refln_loop)
    gemmi::fail("_refln category not found in mmCIF block: " + rb.block.name);
  const cif::Loop& loop = *rb.refln_loop;
  std::vector<int> indices;
  std::string tag = "_refln.";
  for (const Entry& tr : translation_table) {
    tag.replace(7, std::string::npos, tr.refln_tag);
    int index = loop.find_tag(tag);
    if (index == -1 && tr.refln_tag2) {
      tag.replace(7, std::string::npos, tr.refln_tag2);
      index = loop.find_tag(tag);
    }
    if (index != -1) {
      indices.push_back(index);
      mtz.columns.emplace_back();
      gemmi::Mtz::Column& col = mtz.columns.back();
      col.dataset_id = tr.dataset_id;
      col.type = tr.col_type;
      col.label = tr.col_label;
      col.parent = &mtz;
      col.idx = mtz.columns.size() - 1;
    } else if (tr.col_type == 'H') {
      gemmi::fail("Miller index tag not found: " + tag);
    }
  }
  mtz.ncol = mtz.columns.size();
  mtz.nreflections = loop.length();
  mtz.data.resize(mtz.ncol * mtz.nreflections);
  bool uses_status = indices.size() > 3 &&
                     gemmi::iequal(loop.tags[indices[3]], "_refln.status");
  int k = 0;
  for (size_t i = 0; i < loop.values.size(); i += loop.tags.size()) {
    size_t j = 0;
    for (; j != 3; ++j)
      mtz.data[k++] = (float) cif::as_int(loop.values[i + indices[j]]);
    if (uses_status)
      mtz.data[k++] = status_to_freeflag(loop.values[i + indices[j++]]);
    for (; j != indices.size(); ++j)
      mtz.data[k++] = (float) cif::as_number(loop.values[i + indices[j]]);
  }
  if (options[Verbose])
    fprintf(stderr, "Writing %s ...\n", mtz_path.c_str());
  try {
    mtz.write_to_file(mtz_path);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR writing %s: %s\n", mtz_path.c_str(), e.what());
    std::exit(3);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  bool convert_all = p.options[Dir];
  p.require_positional_args(convert_all ? 1 : 2);
  bool verbose = p.options[Verbose];
  const char* cif_path = p.nonOption(0);
  if (verbose)
    fprintf(stderr, "Reading %s ...\n", cif_path);
  cif::Document doc = gemmi::read_cif_gz(cif_path);
  if (convert_all) {
    for (cif::Block& block : doc.blocks) {
      std::string path = p.options[Dir].arg;
      path += '/';
      path += block.name;
      gemmi::ReflnBlock rb(std::move(block));
      convert_cif_block_to_mtz(rb, path, p.options);
    }
  } else {
    const char* mtz_path = p.nonOption(1);
    cif::Block* block = &doc.blocks.at(0);
    if (p.options[BlockName]) {
      block = doc.find_block(p.options[BlockName].arg);
      if (!block) {
        fprintf(stderr, "Block not found: %s\n", p.options[BlockName].arg);
        std::exit(1);
      }
    }
    gemmi::ReflnBlock rb(std::move(*block));
    convert_cif_block_to_mtz(rb, mtz_path, p.options);
  }
  if (verbose)
    fprintf(stderr, "Done.\n");
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
