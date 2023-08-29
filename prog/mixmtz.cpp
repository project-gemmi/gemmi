// Copyright 2019 Global Phasing Ltd.
//
// A subset of CAD functionality.

#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#define GEMMI_PROG mtzmix
#include "options.h"
#include <stdio.h>

using gemmi::Mtz;

namespace {

enum OptionIndex { Force=4, Asu, ToggleEndian };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] MTZ_IN1[...] MTZ_OUT"
    "\nWrite selected columns from one or more MTZ files to a new file."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Force, 0, "", "force", Arg::None,
    "  --force  \tDo not check if merged files are compatible." },
  { Asu, 0, "", "asu", Arg::None,
    "  --asu  \tMove reflections into conventional ASU." },
  { 0, 0, 0, 0, 0, 0 }
};

struct InputSpec {
  InputSpec(Mtz&& mtz_) : mtz(std::move(mtz_)) {}
  Mtz mtz;
};

Mtz merge(const std::vector<InputSpec>& input_list) {
  assert(!input_list.empty());
  const InputSpec& input0 = input_list[0];
  const Mtz& mtz0 = input0.mtz;
  std::vector<int> indices0 = mtz0.sorted_row_indices();
  Mtz out;
  out.data.reserve(mtz0.data.size());
  out.spacegroup = mtz0.spacegroup;
  out.nreflections = mtz0.nreflections;
  out.cell = mtz0.cell;
  out.sort_order = {{1, 2, 3, 0, 0}};
  out.columns = mtz0.columns;
  out.datasets = mtz0.datasets;
  out.history = mtz0.history;
  //for (int idx : indices0)
  //  out.data.insert(out.data.end(),
  //                  mtz0.data.begin() + idx * out.ncol,
  //                  mtz0.data.begin() + (idx + 1) * out.ncol);
  return out;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.nonOptionsCount() < 2)
    p.print_try_help_and_exit("Specify input and output MTZ files.");
  bool verbose = p.options[Verbose];
  std::vector<InputSpec> input_list;
  try {
    for (int i = 0; i < p.nonOptionsCount() - 1; ++i) {
      const char* path = p.nonOption(i);
      if (verbose)
        fprintf(stderr, "Reading %s ...\n", path);
      input_list.emplace_back();
      input_list.back().read_file_gz(path);
      if (p.options[Asu]) {
        // TODO
      }
    }
    Mtz output(merge(input_list));
    const char* out_path = p.nonOption(p.nonOptionsCount() - 1);
    output.write_to_file(out_path);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
