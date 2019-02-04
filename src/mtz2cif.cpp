// Copyright 2019 Global Phasing Ltd.
//
// convert MTZ to SF-mmCIF

// TODO:
//  - cell parameters may be different in CELL and DCELL records, check for it
//  - what to do with _refln.status
//  - check that the FP column is not from Refmac
//  - should we allow for repeated column name in MTZ?

#include <algorithm>
#include <stdio.h>
#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/atox.hpp>     // for read_word
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#define GEMMI_PROG mtz2cif
#include "options.h"

enum OptionIndex { Verbose=3, Spec, PrintSpec };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] MTZ_FILE CIF_FILE"
    "\nOptions:"},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Spec, 0, "", "spec", Arg::Required,
    "  --spec=FILE  \tColumn and format specification." },
  { PrintSpec, 0, "", "print-spec", Arg::None,
    "  --print-spec  \tPrint default spec and exit." },
  { NoOp, 0, "", "", Arg::None,
    "\nIf CIF_FILE is -, the output is printed to stdout."
    "\nIf spec is -, it is read from stdin."
    "\n\nLines in the spec file have format:"
    "\n  [FLAG] COLUMN TYPE TAG [FORMAT]"
    "\nfor example:"
    "\n  SIGF_native * SIGF_meas_au 12.5e"
    "\n  FREE I pdbx_r_free_flag 3d"
    "\nFLAG (optional) is either ? or &:"
    "\n  ? = column is skipped if absent in the MTZ file."
    "\n  & = column is skipped if previous line was skipped. Example:"
    "\n      ? I    J intensity_meas"
    "\n      & SIGI Q intensity_sigma"
    "\nCOLUMN is MTZ column label. Columns H K L are added if not specified."
    "\n  Alternative labels can be separated with | (e.g. FREE|FreeR_flag)."
    "\nTYPE is used for checking the columm type, unless it is '*'."
    "\nTAG does not include category name, it is only the part after _refln."
    "\nFORMAT (optional) is printf-like floating-point format:"
    "\n - one of e, f, g with optional flag, width and precision"
    "\n - flag is one of + - # _; '_' stands for ' ', for example '_.4f'"
    "\n - since all numbers in MTZ are stored as float, the integer columns use"
    "\n   the same format as float. The format of _refln.status is ignored."
  },
  { 0, 0, 0, 0, 0, 0 }
};

struct Trans {
  int col_idx;
  bool is_status;
  std::string refln_tag;
  std::string format;
};

static const char* default_spec[] = {
  "H H index_h",
  "K H index_k",
  "L H index_l",
  "? I       J intensity_meas",
  "& SIGI    Q intensity_sigma",
  "? I(+)    K pdbx_I_plus",
  "& SIGI(+) M pdbx_I_plus_sigma",
  "? I(-)    K pdbx_I_minus",
  "& SIGI(-) M pdbx_I_minus_sigma",
  "? FP      F F_meas_au", // check also if the MTZ is not from refmac
  "& SIGFP   Q F_meas_sigma_au",
  "? F(+)    G pdbx_F_plus",
  "& SIGF(+) L pdbx_F_plus_sigma",
  "? F(-)    G pdbx_F_minus",
  "& SIGF(-) L pdbx_F_minus_sigma",
  "? FREE|RFREE|FreeR_flag I status",
  "? FWT|2FOFCWT      F pdbx_FWT",
  "& PHWT|PH2FOFCWT   P pdbx_PHWT",
  "? DELFWT|FOFCWT    F pdbx_DELFWT",
  "& DELPHWT|PHFOFCWT P pdbx_DELPHWT",
};

static int find_column_index(const std::string& column, const gemmi::Mtz& mtz) {
  int idx = -1;
  for (const std::string& label : gemmi::split_str(column, '|')) {
    for (size_t i = 0; i != mtz.columns.size(); ++i) {
      if (mtz.columns[i].label == label) {
        if (idx == -1)
          idx = (int) i;
        else
          fprintf(stderr, "Warning: duplicate column %s\n", label.c_str());
      }
    }
    if (idx != -1)
      break;
  }
  return idx;
}

void check_format(const std::string& fmt) {
  // expected format: [#_+-]?\d*(\.\d+)?[fFgGeEc]
  if (fmt.find('%') != std::string::npos)
    gemmi::fail("Specify format without %. Got: " + fmt);
  size_t pos = 0;
  if (fmt[0] == '_' || fmt[0] == '+' || fmt[0] == '-' || fmt[0] == '#')
   pos = 1;
  while (gemmi::is_digit(fmt[pos]))
    ++pos;
  if (fmt[pos] == '.' && gemmi::is_digit(fmt[pos+1])) {
    pos += 2;
    while (gemmi::is_digit(fmt[pos]))
      ++pos;
  }
  if (std::isalpha(fmt[pos]))
    ++pos;
  if (fmt[pos] != '\0')
    gemmi::fail("wrong format : " + fmt + "\nCorrect examples: g, .4f, 12.5e");
  char c = gemmi::alpha_up(fmt.back());
  if (c != 'F' && c != 'G' && c != 'E')
    gemmi::fail("expected floating-point format, got: " + fmt);
}

static std::vector<Trans> parse_spec(const gemmi::Mtz& mtz,
                                     const std::vector<std::string>& lines) {
  std::vector<Trans> spec;
  size_t prev_size = 0;
  bool discard = false;
  for (const std::string& line : lines) {
    Trans tr;
    const char* p = line.c_str();
    if (*p == '&') {
      if (discard)
        continue;
    } else {
      prev_size = spec.size();
      discard = false;
    }
    bool optional = (*p == '?' || *p == '&');
    if (optional)
      ++p;
    std::string column = gemmi::read_word(p, &p);
    std::string type = gemmi::read_word(p, &p);
    if (type.size() != 1)
      gemmi::fail("Spec error: MTZ type '" + type + "' is not one character,"
                  "\nin line: " + line);
    tr.refln_tag = gemmi::read_word(p, &p);
    if (tr.refln_tag[0] == '_' || tr.refln_tag.find('.') != std::string::npos)
      gemmi::fail("Spec error: expected tag part after _refln., got: " +
                  tr.refln_tag + "\nin line: " + line);
    tr.col_idx = find_column_index(column, mtz);
    if (tr.col_idx == -1) {
      if (!optional)
        gemmi::fail("Column not found: " + column);
      spec.resize(prev_size);
      discard = true;
      continue;
    }
    const gemmi::Mtz::Column& col = mtz.columns[tr.col_idx];
    if (type[0] != '*' && col.type != type[0])
      gemmi::fail("Column " + col.label + " has type " +
                  std::string(1, col.type) + " not " + type);
    tr.is_status = gemmi::iequal(tr.refln_tag, "status");
    std::string fmt = gemmi::read_word(p, &p);
    if (fmt.empty()) {
      tr.format = "%g";
    } else if (!tr.is_status) {
      check_format(fmt);
      tr.format = "%" + fmt;
      if (tr.format[1] == '_')
        tr.format[1] = ' ';
    }
    spec.push_back(tr);
  }
  return spec;
}

#if defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif

static void write_cif(const gemmi::Mtz& mtz, const std::vector<Trans> spec,
                      FILE* out) {
  std::string id = ".";
  fprintf(out, "# Converted from MTZ by gemmi-mtz2cif " GEMMI_VERSION "\n");
  fprintf(out, "# MTZ title: %s\n", mtz.title.c_str());
  for (size_t i = 0; i != mtz.history.size(); ++i)
    fprintf(out, "# MTZ history #%zu: %s\n", i, mtz.history[i].c_str());
  fprintf(out, "data_%s\n\n", "mtz");
  fprintf(out, "_entry.id %s\n\n", id.c_str());
  const gemmi::UnitCell& cell = mtz.get_cell();
  fprintf(out, "_cell.entry_id %s\n", id.c_str());
  fprintf(out, "_cell.length_a    %8.3f\n", cell.a);
  fprintf(out, "_cell.length_b    %8.3f\n", cell.b);
  fprintf(out, "_cell.length_c    %8.3f\n", cell.c);
  fprintf(out, "_cell.angle_alpha %8.3f\n", cell.alpha);
  fprintf(out, "_cell.angle_beta  %8.3f\n", cell.beta);
  fprintf(out, "_cell.angle_gamma %8.3f\n\n", cell.gamma);
  if (const gemmi::SpaceGroup* sg = mtz.spacegroup) {
    fprintf(out, "_symmetry.entry_id %s\n", id.c_str());
    fprintf(out, "_symmetry.space_group_name_H-M %s\n", sg->hm);
    fprintf(out, "_symmetry.Int_Tables_number %d\n\n", sg->number);
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }
  fprintf(out, "loop_\n");
  for (const Trans& tr : spec) {
    const gemmi::Mtz::Column& col = mtz.columns.at(tr.col_idx);
    const gemmi::Mtz::Dataset& ds = mtz.dataset(col.dataset_number);
    fprintf(out, "_refln.%-24s # %-14s from dataset %s\n",
            tr.refln_tag.c_str(), col.label.c_str(), ds.dataset_name.c_str());
  }
  for (int i = 0; i != mtz.nreflections; ++i) {
    const float* row = &mtz.raw_data[i*mtz.ncol];
    bool first = true;
    for (const Trans& tr : spec) {
      if (first)
        first = false;
      else
        fputc(' ', out);
      float v = row[tr.col_idx];
      if (tr.is_status)
        fputc(v ? 'f' : '<', out);
      else
        fprintf(out, tr.format.c_str(), v);
    }
    fputc('\n', out);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[PrintSpec]) {
    for (const char* line : default_spec)
      printf("%s\n", line);
    return 0;
  }
  p.require_positional_args(2);
  bool verbose = p.options[Verbose];
  const char* mtz_path = p.nonOption(0);
  const char* cif_path = p.nonOption(1);
  gemmi::Mtz mtz;
  if (verbose) {
    fprintf(stderr, "Reading %s ...\n", mtz_path);
    mtz.warnings = stderr;
  }
  try {
    gemmi::fileptr_t f_in = gemmi::file_open(mtz_path, "rb");
    mtz.read_all_headers(f_in.get());
    mtz.read_raw_data(f_in.get());
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR reading %s: %s\n", mtz_path, e.what());
    return 1;
  }
  if (verbose)
    fprintf(stderr, "Writing %s ...\n", cif_path);
  std::vector<Trans> spec;
  try {
  std::vector<std::string> lines;
    if (p.options[Spec]) {
      char buf[256];
      const char* spec_path = p.options[Spec].arg;
      gemmi::fileptr_t f_spec = gemmi::file_open_or(spec_path, "r", stdin);
      while (fgets(buf, sizeof(buf), f_spec.get()) != NULL) {
        const char* start = gemmi::skip_blank(buf);
        if (*start != '\0' && *start != '#')
          lines.emplace_back(start);
      }
    } else {
      lines.reserve(sizeof(default_spec) / sizeof(default_spec[0]));
      for (const char* line : default_spec)
        lines.emplace_back(line);
    }
    spec = parse_spec(mtz, lines);
    if (spec.empty())
      gemmi::fail("Empty translation spec");
    for (size_t i = 0; i != spec.size(); ++i)
      for (size_t j = i + 1; j != spec.size(); ++j)
        if (spec[i].refln_tag == spec[j].refln_tag)
          gemmi::fail("duplicated output tag: " + spec[i].refln_tag);
    // H, K, L must be the first columns in MTZ and are required in _refln
    for (int i = 2; i != -1; --i)
      if (!gemmi::in_vector_f([&](const Trans& t) { return t.col_idx == i; },
                              spec)) {
        spec.insert(spec.begin(), Trans{i, false, "index_h", "%g "});
        spec.front().refln_tag.back() += i;
      }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "Problem in translation spec: %s\n", e.what());
    return 2;
  }
  try {
    gemmi::fileptr_t f_out = gemmi::file_open_or(cif_path, "w", stdout);
    write_cif(mtz, spec, f_out.get());
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR writing %s: %s\n", cif_path, e.what());
    return 3;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
