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
  const char* col_label;
  char col_type;
  unsigned char dataset_id;
};

// When we have a few alternative mmCIF tags for the same MTZ label,
// they are in consecutive rows and all but the last one have null col_label.
static Entry conv_table[] = {
  {"index_h",                    "H",          'H', 0},
  {"index_k",                    "K",          'H', 0},
  {"index_l",                    "L",          'H', 0},
  {"pdbx_r_free_flag",           nullptr,      'I', 0},
  {"status",                     "FreeR_flag", 's', 0}, // s is a special flag
  {"intensity_meas",             nullptr,      'J', 1},
  {"intensity_net",              "I",          'J', 1},
  {"intensity_sigma",            "SIGI",       'Q', 1},
  {"pdbx_I_plus",                "I(+)",       'K', 1},
  {"pdbx_I_plus_sigma",          "SIGI(+)",    'M', 1},
  {"pdbx_I_minus",               "I(-)",       'K', 1},
  {"pdbx_I_minus_sigma",         "SIGI(-)",    'M', 1},
  {"F_meas_au",                  "FP",         'F', 1},
  {"F_meas_sigma_au",            "SIGFP",      'Q', 1},
  {"pdbx_F_plus",                "F(+)",       'G', 1},
  {"pdbx_F_plus_sigma",          "SIGF(+)",    'L', 1},
  {"pdbx_F_minus",               "F(-)",       'G', 1},
  {"pdbx_F_minus_sigma",         "SIGF(-)",    'L', 1},
  {"pdbx_anom_difference",       "DP",         'D', 1},
  {"pdbx_anom_difference_sigma", "SIGDP",      'Q', 1},
  {"F_calc",                     "FC",         'F', 1},
  {"phase_calc",                 "PHIC",       'P', 1},
  {"fom",                        nullptr,      'W', 1},
  {"weight",                     "FOM",        'W', 1},
  {"pdbx_HL_A_iso",              "HLA",        'A', 1},
  {"pdbx_HL_B_iso",              "HLB",        'A', 1},
  {"pdbx_HL_C_iso",              "HLC",        'A', 1},
  {"pdbx_HL_D_iso",              "HLD",        'A', 1},
  {"pdbx_FWT",                   "FWT",        'F', 1},
  {"pdbx_PHWT",                  "PHWT",       'P', 1},
  {"pdbx_DELFWT",                "DELFWT",     'F', 1},
  {"pdbx_DELPHWT",               "DELPHWT",    'P', 1},
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
gemmi::ReflnBlock& get_block_by_name(std::vector<gemmi::ReflnBlock>& rblocks,
                                     const std::string& name) {
  for (gemmi::ReflnBlock& rb : rblocks)
    if (rb.block.name == name)
      return rb;
  gemmi::fail("block not found: " + name);
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
  const cif::Loop* loop = rb.refln_loop ? rb.refln_loop : rb.diffrn_refln_loop;
  if (!loop)
    gemmi::fail("_refln category not found in mmCIF block: " + rb.block.name);
  bool uses_status = false;
  std::vector<int> indices;
  std::string tag = loop->tags[0].substr(0, loop->tags[0].find('.') + 1);
  const size_t len = tag.length();
  for (auto c = std::begin(conv_table); c != std::end(conv_table); ++c) {
    tag.replace(len, std::string::npos, c->refln_tag);
    int index = loop->find_tag(tag);
    if (index != -1) {
      indices.push_back(index);
      mtz.columns.emplace_back();
      gemmi::Mtz::Column& col = mtz.columns.back();
      col.dataset_id = c->dataset_id;
      col.type = c->col_type;
      if (col.type == 's') {
        col.type = 'I';
        uses_status = true;
      }
      while (!c->col_label)
        ++c;
      col.label = c->col_label;
      col.parent = &mtz;
      col.idx = mtz.columns.size() - 1;
    } else if (c->col_type == 'H') {
      gemmi::fail("Miller index tag not found: " + tag);
    }
  }
  mtz.ncol = mtz.columns.size();
  mtz.nreflections = loop->length();
  mtz.data.resize(mtz.ncol * mtz.nreflections);
  int k = 0;
  for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
    size_t j = 0;
    for (; j != 3; ++j)
      mtz.data[k++] = (float) cif::as_int(loop->values[i + indices[j]]);
    if (uses_status)
      mtz.data[k++] = status_to_freeflag(loop->values[i + indices[j++]]);
    for (; j != indices.size(); ++j)
      mtz.data[k++] = (float) cif::as_number(loop->values[i + indices[j]]);
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
  auto rblocks = gemmi::as_refln_blocks(gemmi::read_cif_gz(cif_path).blocks);
  const gemmi::SpaceGroup* first_sg = nullptr;
  if (convert_all) {
    bool ok = true;
    for (gemmi::ReflnBlock& rb : rblocks) {
      std::string path = p.options[Dir].arg;
      path += '/';
      path += rb.block.name;
      path += ".mtz";
      if (!first_sg)
        first_sg = rb.spacegroup;
      else if (!rb.spacegroup)
        rb.spacegroup = first_sg;
      try {
        convert_cif_block_to_mtz(rb, path, p.options);
      } catch (std::runtime_error& e) {
        fprintf(stderr, "ERROR: %s\n", e.what());
        ok = false;
      }
    }
    if (!ok)
      return 1;
  } else {
    const char* mtz_path = p.nonOption(1);
    try {
      const gemmi::ReflnBlock& rb = p.options[BlockName]
        ? get_block_by_name(rblocks, p.options[BlockName].arg)
        : rblocks.at(0);
      convert_cif_block_to_mtz(rb, mtz_path, p.options);
    } catch (std::runtime_error& e) {
      fprintf(stderr, "ERROR: %s\n", e.what());
      return 1;
    }
  }
  if (verbose)
    fprintf(stderr, "Done.\n");
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
