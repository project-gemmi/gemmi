// Copyright 2020 Global Phasing Ltd.
//
// Reindex merged or unmerged MTZ file.
// TODO: handle operations that results in non-integral indices.

#include <cstdio>
#include <cstring>  // for strpbrk
#include <algorithm>
#ifndef GEMMI_ALL_IN_ONE
# define GEMMI_WRITE_IMPLEMENTATION 1
#endif
#include <gemmi/mtz.hpp>
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#define GEMMI_PROG reindex
#include "options.h"

using std::fprintf;

namespace {

enum OptionIndex { Hkl=4, NoHistory, NoSort };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] INPUT_MTZ OUTPUT_MTZ"
    "\nOptions:"},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Hkl, 0, "", "hkl", Arg::Required,
    "  --hkl=OP  \tReindexing transform as triplet (e.g. k,h,-l)." },
  { NoHistory, 0, "", "no-history", Arg::None,
    "  --no-history  \tDo not add 'Reindexed with...' line to mtz HISTORY." },
  { NoSort, 0, "", "no-sort", Arg::None,
    "  --no-sort  \tDo not reorder reflections." },
  { NoOp, 0, "", "", Arg::None,
    "\nInput file can be gzipped." },
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
  if (!p.options[Hkl]) {
    fprintf(stderr, "Specify transform with option --hkl\n");
    return 1;
  }
  const char* hkl_arg = p.options[Hkl].arg;
  gemmi::Op op = gemmi::parse_triplet(hkl_arg);
  if (std::strpbrk(hkl_arg, "xyzabcXYZABC"))
    gemmi::fail("specify OP in terms of h, k and l");
  if (op.tran != gemmi::Op::Tran{{0, 0, 0}})
    gemmi::fail("reindexing operator should not have a translation");
  // transpose rotation
  gemmi::Op real_space_op;
  real_space_op.rot = op.transposed_rot();
  if (verbose)
    fprintf(stderr, "real space transformation: %s\n", real_space_op.triplet().c_str());
  gemmi::Mtz mtz;
  if (verbose) {
    fprintf(stderr, "Reading %s ...\n", input_path);
    mtz.warnings = stderr;
  }
  try {
    mtz.read_input(gemmi::MaybeGzipped(input_path), true);
    mtz.switch_to_original_hkl();
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR reading %s: %s\n", input_path, e.what());
    return 1;
  }
  // reindex
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
    mtz.set_hkl(n, real_space_op.apply_to_hkl(mtz.get_hkl(n)));
  // hand change requires data modification
  if (op.det_rot() < 0)
    for (gemmi::Mtz::Column& column : mtz.columns) {
      // negate anomalous difference
      if (column.type == 'D') {
        for (float& value : column)
          value = -value;
        if (verbose)
          fprintf(stderr, "Column %s: anomalous difference negated.\n",
                  column.label.c_str());
        continue;
      }
      // swap (+) and (-)
      size_t pos = column.label.find("(+)");
      if (pos != std::string::npos) {
        std::string minus_label = column.label;
        minus_label[pos+1] = '-';
        gemmi::Mtz::Column* minus_column =
            mtz.column_with_label(minus_label, &column.dataset());
        if (minus_column) {
          for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
            std::swap(mtz.data[n + column.idx],
                      mtz.data[n + minus_column->idx]);
          if (verbose)
            fprintf(stderr, "Swapped values from %s and %s.\n",
                    column.label.c_str(), minus_label.c_str());
        } else {
          fprintf(stderr, "Warning: matching pair not found for: %s\n",
                  column.label.c_str());
        }
      }
    }
  // change space group
  mtz.spacegroup = gemmi::find_spacegroup_by_change_of_basis(mtz.spacegroup, op.inverse());
  if (mtz.spacegroup) {
    if (verbose)
      fprintf(stderr, "Space group changed from %s to %s\n.",
              mtz.spacegroup_name.c_str(), mtz.spacegroup->xhm().c_str());
  } else {
    fprintf(stderr, "WARNING: new space group name could not be determined.\n");
  }

  if (mtz.is_merged())
    mtz.ensure_asu();
  else
    mtz.switch_to_asu_hkl();
  if (!p.options[NoSort])
    mtz.sort();
  if (!p.options[NoHistory])
    mtz.history.emplace(mtz.history.begin(),
                        "Reindexed with gemmi-reindex " GEMMI_VERSION);
  if (verbose)
    fprintf(stderr, "Writing %s ...\n", output_path);
  mtz.write_to_file(output_path);
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
