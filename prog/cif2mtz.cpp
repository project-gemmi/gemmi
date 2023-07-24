// Copyright 2019 Global Phasing Ltd.
//
// convert SF-mmCIF to MTZ

#include <cstdio>             // for fprintf
#include <cstdlib>            // for exit
#include <iostream>           // for cerr
#include <memory>             // for unique_ptr
#include <gemmi/read_cif.hpp> // for read_cif_gz, read_first_block_gz
#include <gemmi/cif2mtz.hpp>  // for CifToMtz
#include <gemmi/pdb_id.hpp>   // for is_pdb_code, expand_pdb_code_to_path

#define GEMMI_PROG cif2mtz
#include "options.h"

namespace {

using std::fprintf;

struct Cif2MtzArg: public Arg {
  static option::ArgStatus ModeChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"normal", "friedel", "unmerged", "auto"});
  }
};

enum OptionIndex {
  BlockName=4, BlockNumber, Add, List, Dir, Spec, PrintSpec, Title,
  History, Wavelength, Unmerged, ReflnTo,
  Sort, Asu, SkipNegativeSigma, ZeroToMnf, Local
};

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
    "  -b NAME, --block=NAME  \tmmCIF block to convert, by name." },
  { BlockNumber, 0, "B", "", Arg::Int,
    "  -B INDEX  \tmmCIF block to convert, by index (default: 1)." },
  { Add, 0, "", "add", Arg::Required,
    "  --add CIF_FILE  \tUse additional input mmCIF file, first block." },
  { List, 0, "l", "list", Arg::None,
    "  -l, --list  \tDry run and list blocks in mmCIF file." },
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
  { Wavelength, 0, "", "wavelength", Arg::Float,
    "  --wavelength=LAMBDA  \tSet wavelength (default: from input file)." },
  { Unmerged, 0, "u", "unmerged", Arg::None,
    "  -u, --unmerged  \tWrite unmerged MTZ file(s)." },
  { ReflnTo, 0, "", "refln-to", Cif2MtzArg::ModeChoice,
    "  --refln-to=MODE  \tRead refln category as: "
    "normal, friedel (converts Fs to F+/F-), unmerged or auto (default)." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tOrder reflections according to Miller indices." },
  { Asu, 0, "", "asu", Arg::AsuChoice,
    "  --asu=ccp4|tnt  \tMove merged reflections into CCP4 or TNT ASU." },
  { SkipNegativeSigma, 0, "", "skip-negative-sigma", Arg::None,
    "  --skip-negative-sigma  \tSkip reflections with sigma<0 (in any MTZ Q column)." },
  { ZeroToMnf, 0, "", "zero-to-mnf", Arg::None,
    "  --zero-to-mnf  \tIf value and sigma are 0, set both to MNF." },
  { Local, 0, "", "local", Arg::None,
    "  --local  \tTake file from local copy of the PDB archive in "
    "$PDB_DIR/structures/divided/structure_factors/" },
  { NoOp, 0, "", "", Arg::None,
    "\nFirst variant: converts the first block of CIF_FILE, or the block"
    "\nspecified with --block=NAME, to MTZ file with given name."
    "\n\nSecond variant: converts each block of CIF_FILE to one MTZ file"
    "\n(block-name.mtz) in the specified DIRECTORY."
    "\n\nIf CIF_FILE is -, the input is read from stdin."
    "\n\nTo convert data from multiple CIF files into one MTZ file use:"
    "\n  " EXE_NAME " [options] CIF1 --add CIF2 --add CIF3 MTZ_FILE"
    "\nIt's for special cases, e.g. when map coefficients are in separate files."
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

gemmi::ReflnBlock& get_block_by_number(std::vector<gemmi::ReflnBlock>& rblocks, int n) {
  if (size_t(n - 1) >= rblocks.size())
    gemmi::fail("block index must be from 1 to ", std::to_string(rblocks.size()));
  return rblocks.at(n - 1);
}

void zero_to_mnf(gemmi::Mtz& mtz) {
  std::vector<size_t> cols;
  for (size_t i = 3; i < mtz.columns.size() - 1; ++i) {
    // find consecutive value and sigma columns
    char t1 = mtz.columns[i].type;
    char t2 = mtz.columns[i+1].type;
    if ((t1 == 'J' || t1 == 'F' || t1 == 'D' || t1 == 'G' || t1 == 'K') &&
        (t2 == 'Q' || t2 == 'L' || t2 == 'M'))
      cols.push_back(i);
  }
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size()) {
    for (size_t idx : cols)
      if (mtz.data[n+idx] == 0 && mtz.data[n+idx+1] == 0)
        mtz.data[n+idx] = mtz.data[n+idx+1] = NAN;
  }
}

void print_block_info(gemmi::ReflnBlock& rb, const gemmi::Mtz& mtz) {
  std::printf("%.*s %zu x %zu ->",
              (int)rb.tag_offset() - 1, rb.default_loop->tags.at(0).c_str(),
              rb.default_loop->width(), rb.default_loop->length());
  for (const gemmi::Mtz::Column& col : mtz.columns)
    std::printf(" %s", col.label.c_str());
  std::putchar('\n');
  for (const std::string& d : rb.block.find_values("_diffrn.details"))
    if (!gemmi::cif::is_null(d))
      std::printf("details: %s\n", gemmi::cif::as_string(d).c_str());
  std::putchar('\n');
}

bool is_column_data_identical(const gemmi::Mtz& mtz, size_t i, size_t j) {
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
    if (!gemmi::impl::is_same(mtz.data[n + i], mtz.data[n + j]))
      return false;
  return true;
}

void change_label_to_unique(gemmi::Mtz::Column &col) {
  size_t label_size = col.label.size();
  const gemmi::Mtz::Dataset& ds = col.parent->dataset(col.dataset_id);
  for (int appendix = 2; ; ++appendix) {
    col.label += std::to_string(appendix);
    if (col.parent->column_with_label(col.label, &ds) == &col)
      break;
    col.label.resize(label_size);
    assert(appendix < 1000);
  }
}

void add_columns_from_other_mtz(gemmi::Mtz& mtz, const gemmi::Mtz& mtz2) {
  size_t ncol = mtz.columns.size();
  size_t ncol2 = mtz2.columns.size();
  if (ncol2 <= 3)  // no data apart from hkl
    return;
  // first copy all columns apart from H K L;
  // reflections that are only in mtz2 are ignored
  mtz.copy_column(-1, mtz2.columns[3], std::vector<std::string>(ncol2 - 4));
  // avoid the same names: remove or rename
  for (size_t i = ncol; i < mtz.columns.size(); ++i) {
    gemmi::Mtz::Column& col = mtz.columns[i];
    for (size_t j = 3; j < ncol; ++j)
      if (col.label == mtz.columns[j].label &&
          col.dataset_id == mtz.columns[j].dataset_id) {
        if (is_column_data_identical(mtz, i, j))
          mtz.remove_column(i--);
        else
          change_label_to_unique(col);
        break;
      }
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.check_exclusive_pair(BlockName, BlockNumber);
  p.check_exclusive_pair(Add, Dir);
  p.check_exclusive_pair(Add, List);
  if (p.options[PrintSpec]) {
    std::printf("# Each line in the spec contains 4-5 words:\n"
                "# - tag (without category) from _refln or _diffrn_refln\n"
                "# - MTZ column label\n"
                "# - MTZ column type\n"
                "# - MTZ dataset for the column (must be 0 or 1)\n"
                "# - (optional) how to map mmCIF symbols to MTZ numbers\n");
    bool merged = !p.options[Unmerged];
    if (merged)
      std::printf("# For MERGED data only. Use --print-spec --unmerged for unmerged.\n"
                  "# Conversion of the first MTZ columns (H K L) is hardcoded -"
                  " not in the spec.\n");
    else
      std::printf("# For UNMERGED data only. Conversion of the first MTZ columns\n"
                  "# (H K L M/ISYM BATCH) is hardcoded - not in the spec.\n");
    for (const char** line = gemmi::CifToMtz::default_spec(merged); *line != nullptr; ++line)
      std::printf("%s\n", *line);
    return 0;
  }

  bool convert_all = p.options[Dir] || p.options[List];
  p.require_positional_args(convert_all ? 1 : 2);

  gemmi::CifToMtz cif2mtz;
  cif2mtz.verbose = p.options[Verbose];
  cif2mtz.force_unmerged = p.options[Unmerged];
  if (p.options[ReflnTo] && p.options[ReflnTo].arg[0] == 'u')
    cif2mtz.force_unmerged = true;
  if (p.options[Title])
    cif2mtz.title = p.options[Title].arg;
  for (const option::Option* opt = p.options[History]; opt; opt = opt->next())
    cif2mtz.history.push_back(opt->arg);
  if (p.options[Wavelength])
    cif2mtz.wavelength = std::strtod(p.options[Wavelength].arg, nullptr);
  try {
    if (p.options[Spec])
      read_spec_file(p.options[Spec].arg, cif2mtz.spec_lines);
    std::string cif_path = p.nonOption(0);
    if (p.options[Local]) {
      if (!gemmi::is_pdb_code(cif_path))
        gemmi::fail("Not a PDB ID: " + cif_path);
      cif_path = gemmi::expand_pdb_code_to_path(cif_path, 'S');
      if (cif_path.empty())
        gemmi::fail("To use option --local set $PDB_DIR.");
    }
    if (cif2mtz.verbose)
      fprintf(stderr, "Reading %s ...\n", cif_path.c_str());
    // If the file is gzipped and huge, reading it will result in error,
    // but if we are interested only in the first block we can recover.
    gemmi::cif::Document doc;
    // this limit is used only for when reading first block from gzipped file
    constexpr int block_limit = 1073741824;  // 1GB, used only for gz files
    if (convert_all || p.options[BlockName] || p.options[BlockNumber])
      doc = gemmi::read_cif_gz(cif_path);
    else  // optimization: ignore other blocks
      doc = gemmi::read_first_block_gz(cif_path, block_limit);
    auto rblocks = gemmi::as_refln_blocks(std::move(doc.blocks));
    char mode = p.options[ReflnTo] ? p.options[ReflnTo].arg[0] : 'a';
    if (convert_all) {
      bool ok = true;
      for (gemmi::ReflnBlock& rb : rblocks) {
        try {
          if (p.options[List])
            std::printf("--block=%s\n", rb.block.name.c_str());
          gemmi::Mtz mtz = cif2mtz.auto_convert_block_to_mtz(rb, std::cerr, mode);
          if (p.options[List]) {
            print_block_info(rb, mtz);
            continue;
          }
          std::string path = p.options[Dir].arg;
          path += '/';
          path += rb.block.name;
          path += ".mtz";
          if (cif2mtz.verbose)
            fprintf(stderr, "Writing %s ...\n", path.c_str());
          mtz.write_to_file(path);
        } catch (std::exception& e) {
          fprintf(stderr, "ERROR: %s\n", e.what());
          ok = false;
        }
      }
      if (!ok)
        return 1;
    } else {
      const char* mtz_path = p.nonOption(1);
      gemmi::ReflnBlock* rb;
      if (p.options[BlockName]) {
        rb = &get_block_by_name(rblocks, p.options[BlockName].arg);
      } else if (p.options[BlockNumber]) {
        int idx = std::atoi(p.options[BlockNumber].arg);
        rb = &get_block_by_number(rblocks, idx);
      } else {
        rb = &rblocks.at(0);
      }

      gemmi::Mtz mtz = cif2mtz.auto_convert_block_to_mtz(*rb, std::cerr, mode);
      // add data from additional cif block
      if (const option::Option* opt = p.options[BlockName])
        while ((opt = opt->next()) != nullptr) {
          gemmi::ReflnBlock& rblock = get_block_by_name(rblocks, opt->arg);
          gemmi::Mtz mtz2 = cif2mtz.auto_convert_block_to_mtz(rblock, std::cerr, mode);
          add_columns_from_other_mtz(mtz, mtz2);
        }
      if (const option::Option* opt = p.options[BlockNumber])
        while ((opt = opt->next()) != nullptr) {
          gemmi::ReflnBlock& rblock = get_block_by_number(rblocks, std::atoi(opt->arg));
          gemmi::Mtz mtz2 = cif2mtz.auto_convert_block_to_mtz(rblock, std::cerr, mode);
          add_columns_from_other_mtz(mtz, mtz2);
        }
      for (const option::Option* opt = p.options[Add]; opt; opt = opt->next()) {
        if (cif2mtz.verbose)
          fprintf(stderr, "Reading %s ...\n", opt->arg);
        auto rblocks2 = gemmi::as_refln_blocks(
                          gemmi::read_first_block_gz(opt->arg, block_limit).blocks);
        gemmi::Mtz mtz2 = cif2mtz.auto_convert_block_to_mtz(rblocks2.at(0), std::cerr, mode);
        add_columns_from_other_mtz(mtz, mtz2);
      }
      if (p.options[SkipNegativeSigma]) {
        for (const gemmi::Mtz::Column& col : mtz.columns)
          if (col.type == 'Q')  // typically, we'll find one Q column here
            mtz.remove_rows_if([&](const float* row) { return row[col.idx] < 0; });
      }
      if (p.options[ZeroToMnf])
        zero_to_mnf(mtz);
      if (p.options[Asu] && mtz.is_merged())
        mtz.ensure_asu(p.options[Asu].arg[0] == 't');
      if (p.options[Sort]) {
        bool reordered = mtz.sort();
        if (cif2mtz.verbose)
          fprintf(stderr, "Reflection order has %schanged.\n", reordered ? "" : "not ");
      }
      if (cif2mtz.verbose)
        fprintf(stderr, "Writing %s ...\n", mtz_path);
      mtz.write_to_file(mtz_path);
    }
  } catch (std::exception& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  if (cif2mtz.verbose)
    fprintf(stderr, "Done.\n");
  return 0;
}
