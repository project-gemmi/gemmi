// Copyright 2019 Global Phasing Ltd.
//
// convert SF-mmCIF to MTZ

#include <cstdio>             // for fprintf
#include <cstdlib>            // for exit
#include <iostream>           // for cerr
#include <memory>             // for unique_ptr
#include <gemmi/read_cif.hpp> // for read_cif_gz, read_first_block_gz
#include <gemmi/cif2mtz.hpp>  // for CifToMtz
#include <gemmi/to_cif.hpp>  // for CifToMtz

#define GEMMI_PROG cif2mtz
#include "options.h"

namespace {

using std::fprintf;

enum OptionIndex {
  BlockName=4, BlockNumber, Add, List, Dir, Spec, PrintSpec, Title,
  History, Wavelength, Unmerged, OldAnomalous,
  Sort, SkipNegativeSigma, ZeroToMnf, Local
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
  { OldAnomalous, 0, "", "anomalous", Arg::None,
    "  --anomalous  \tConvert _refln.F_meas_au as anomalous data." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tOrder reflections according to Miller indices." },
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
  std::printf("--block=%s - %.*s %zu x %zu ->",
              rb.block.name.c_str(), (int)rb.tag_offset() - 1,
              rb.default_loop->tags.at(0).c_str(),
              rb.default_loop->width(), rb.default_loop->length());
  for (const gemmi::Mtz::Column& col : mtz.columns)
    std::printf(" %s", col.label.c_str());
  std::putchar('\n');
  for (const std::string& d : rb.block.find_values("_diffrn.details"))
    std::printf("  details: %s\n", d.c_str());
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

// Before _refln.pdbx_F_plus/minus was introduced, anomalous data was
// stored as two F_meas_au reflections, say (1,1,3) then (-1,-1,-3),
// with F(-) directly after F(+). This function transcribes it to
// how the anomalous data is stored PDBx/mmCIF nowadays:
//  _refln.F_meas_au -> pdbx_F_plus / pdbx_F_minus,
//  _refln.F_meas_sigma_au -> pdbx_F_plus_sigma / pdbx_F_minus_sigma.
void transcript_old_anomalous_to_standard(gemmi::cif::Loop& loop) {
  using gemmi::fail;
  size_t old_width = loop.width();
  size_t length = loop.length();
  int hkl_idx = loop.find_tag_lc("_refln.index_h");
  if (hkl_idx == -1)
    fail("while reading old anomalous: _refln.index_h not found");
  if ((size_t) hkl_idx+2 >= loop.width() ||
      loop.tags[hkl_idx+1] != "_refln.index_k" ||
      loop.tags[hkl_idx+2] != "_refln.index_l")
    fail("while reading old anomalous: _refln.index_* not found");
  int f_idx = loop.find_tag("_refln.F_meas_au");
  if (f_idx == -1)
    fail("while reading old anomalous: _refln.F_meas_au not found");
  if ((size_t) f_idx+1 >= loop.width() || loop.tags[f_idx+1] != "_refln.F_meas_sigma_au")
    fail("while reading old anomalous: _refln.F_meas_sigma_au not found");
  const char* new_labels[4] = {"_refln.pdbx_F_plus", "_refln.pdbx_F_plus_sigma",
                               "_refln.pdbx_F_minus", "_refln.pdbx_F_minus_sigma"};
  for (const char* label : new_labels)
    if (loop.has_tag(label))
      fail("while reading old anomalous: _refln loop already has ", label);
  loop.tags[f_idx] = new_labels[0];
  loop.tags[f_idx+1] = new_labels[1];
  loop.tags.insert(loop.tags.begin() + f_idx + 2, {new_labels[2], new_labels[3]});
  size_t new_width = loop.width();
  std::vector<std::string> new_values;
  new_values.reserve(length * new_width); // upper bound
  auto new_f = new_values.end();
  gemmi::Miller prev_hkl = {0, 0, 0};
  for (size_t i = 0; i < loop.values.size(); i += old_width) {
    std::string* row = &loop.values[i];
    gemmi::Miller hkl;
    for (int j = 0; j < 3; ++j)
      hkl[j] = gemmi::cif::as_int(row[hkl_idx + j]);
    if (i == 0 || hkl[0] != -prev_hkl[0] || hkl[1] != -prev_hkl[1]
                                         || hkl[2] != -prev_hkl[2]) {
      new_values.insert(new_values.end(), row, row+f_idx);
      new_f = new_values.end();  // .reserve() above prevents re-allocations
      new_values.resize(new_values.size() + 4, ".");
      new_values.insert(new_values.end(), row+f_idx+2, row+old_width);
      prev_hkl = hkl;
    }
    // hkl asu for P1 is: l>0 or (l==0 and (h>0 or (h==0 and k>=0)))
    bool minus = hkl[2] < 0 || (hkl[2] == 0 && (hkl[0] < 0 || (hkl[0] == 0 && hkl[1] < 0)));
    *(new_f + 2 * minus) = row[f_idx];
    *(new_f + 2 * minus + 1) = row[f_idx+1];
  }
  loop.values = std::move(new_values);
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
                "# - (optional) how to map mmCIF symbols to MTZ numbers\n"
                "# The first 3 (5 in unmerged) columns are not in the spec,\n"
                "# they are always H K L (M/ISYM BATCH).\n\n");
    bool merged = !p.options[Unmerged];
    if (merged)
      std::printf("# For MERGED data only. Use --print-spec --unmerged for unmerged.\n");
    for (const char** line = gemmi::CifToMtz::default_spec(merged); *line != nullptr; ++line)
      std::printf("%s\n", *line);
    return 0;
  }

  bool convert_all = p.options[Dir] || p.options[List];
  p.require_positional_args(convert_all ? 1 : 2);

  gemmi::CifToMtz cif2mtz;
  cif2mtz.verbose = p.options[Verbose];
  cif2mtz.force_unmerged = p.options[Unmerged];
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
    if (convert_all) {
      bool ok = true;
      for (gemmi::ReflnBlock& rb : rblocks) {
        try {
          gemmi::Mtz mtz = cif2mtz.convert_block_to_mtz(rb, std::cerr);
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
      const gemmi::ReflnBlock* rb;
      if (p.options[BlockName]) {
        rb = &get_block_by_name(rblocks, p.options[BlockName].arg);
      } else if (p.options[BlockNumber]) {
        int idx = std::atoi(p.options[BlockNumber].arg);
        if (idx <= 0 || (size_t)idx > rblocks.size())
          gemmi::fail("block index must be from 1 to ", std::to_string(rblocks.size()));
        rb = &rblocks.at(idx - 1);
      } else {
        rb = &rblocks.at(0);
      }

      if (p.options[OldAnomalous] && rb->refln_loop != nullptr &&
          !rb->refln_loop->has_tag("_refln.pdbx_F_plus"))
        transcript_old_anomalous_to_standard(*rb->refln_loop);
      gemmi::Mtz mtz = cif2mtz.convert_block_to_mtz(*rb, std::cerr);
      for (const option::Option* opt = p.options[Add]; opt; opt = opt->next()) {
        if (cif2mtz.verbose)
          fprintf(stderr, "Reading %s ...\n", opt->arg);
        auto rblocks2 = gemmi::as_refln_blocks(
                          gemmi::read_first_block_gz(opt->arg, block_limit).blocks);
        gemmi::Mtz mtz2 = cif2mtz.convert_block_to_mtz(rblocks2.at(0), std::cerr);
        size_t ncol = mtz.columns.size();
        size_t ncol2 = mtz2.columns.size();
        if (ncol2 < 4)
          continue;
        mtz.copy_column(-1, mtz2.columns[3], std::vector<std::string>(ncol2-4));
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
      if (p.options[SkipNegativeSigma]) {
        for (const gemmi::Mtz::Column& col : mtz.columns)
          if (col.type == 'Q')  // typically, we'll find one Q column here
            mtz.remove_rows_if([&](const float* row) { return row[col.idx] < 0; });
      }
      if (p.options[ZeroToMnf])
        zero_to_mnf(mtz);
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
