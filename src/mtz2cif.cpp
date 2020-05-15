// Copyright 2019 Global Phasing Ltd.
//
// convert MTZ to SF-mmCIF

// TODO:
//  - cell parameters may be different in CELL and DCELL records, check for it
//  - what to do with _refln.status
//  - check that the FP column is not from Refmac
//  - should we allow for repeated column name in MTZ?

#include <cstdio>
#include <cstdlib>            // for strtod
#include <algorithm>
#include <gemmi/mtz.hpp>
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/atox.hpp>     // for read_word
#include <gemmi/gz.hpp>       // for MaybeGzipped
#include <gemmi/version.hpp>  // for GEMMI_VERSION
#define GEMMI_PROG mtz2cif
#include "options.h"

using std::fprintf;

enum OptionIndex { Spec=4, PrintSpec, BlockName, SkipEmpty, NoComments,
                   Wavelength };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n  " EXE_NAME " [options] MTZ_FILE CIF_FILE"
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
  { SkipEmpty, 0, "", "skip-empty", Arg::None,
    "  --skip-empty  \tSkip reflections with no values." },
  { NoComments, 0, "", "no-comments", Arg::None,
    "  --no-comments  \tDo not write comments in the mmCIF file." },
  { Wavelength, 0, "", "wavelength", Arg::Float,
    "  --wavelength=LAMBDA  \tSet wavelengths (default: from input file)." },
  { NoOp, 0, "", "", Arg::None,
    "\nIf CIF_FILE is -, the output is printed to stdout."
    "\nIf spec is -, it is read from stdin."
    "\n\nLines in the spec file have format:"
    "\n  [FLAG] COLUMN TYPE TAG [FORMAT]"
    "\nfor example:"
    "\n  SIGF_native * SIGF_meas_au 12.5e"
    "\n  FREE I pdbx_r_free_flag 3.0f"
    "\nFLAG (optional) is either ? or &:"
    "\n  ? = ignored if no column in the MTZ file has this name."
    "\n  & = ignored if the previous line was ignored."
    "\n  Example:"
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
  bool is_status = false;
  std::string tag;  // excluding category
  std::string format = "%g";
  int min_width = 0;
};

struct Options {
  std::vector<Trans> spec;
  std::vector<int> value_indices;  // used for --skip_empty
  std::vector<int> sigma_indices;  // used for status 'x'
  bool with_comments;
  const char* block_name;
  const char* mtz_path;
  double wavelength = NAN;
};

static const char* default_merged_spec[] = {
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
  "& DELPHWT|PHDELWT|PHFOFCWT P pdbx_DELPHWT",
};

static const char* default_unmerged_spec[] = {
  "H H index_h",
  "K H index_k",
  "L H index_l",
  "? I       J intensity_net",
  "& SIGI    Q intensity_sigma",
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

int check_format(const std::string& fmt) {
  // expected format: [#_+-]?\d*(\.\d+)?[fFgGeEc]
  int min_width = 0;
  if (fmt.find('%') != std::string::npos)
    gemmi::fail("Specify format without %. Got: " + fmt);
  const char* p = fmt.c_str();
  if (*p == '_' || *p == '+' || *p == '-' || *p == '#')
   ++p;
  if (gemmi::is_digit(*p)) {
    min_width = *p++ - '0';
    if (gemmi::is_digit(*p)) // two digits of width number max
      min_width = min_width * 10 + (*p++ - '0');
  }
  if (*p == '.' && gemmi::is_digit(*(p+1))) {
    p += 2;
    if (gemmi::is_digit(*p)) // two digits of precision numbers max
      ++p;
  }
  if (!std::isalpha(*p) || *(p+1) != '\0')
    gemmi::fail("wrong format : " + fmt + "\nCorrect examples: g, .4f, 12.5e");
  char c = gemmi::alpha_up(*p);
  if (c != 'F' && c != 'G' && c != 'E')
    gemmi::fail("expected floating-point format, got: " + fmt);
  return min_width;
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
    tr.tag = gemmi::read_word(p, &p);
    if (tr.tag[0] == '_' || tr.tag.find('.') != std::string::npos)
      gemmi::fail("Spec error: expected tag part after _refln., got: " +
                  tr.tag + "\nin line: " + line);
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
    tr.is_status = gemmi::iequal(tr.tag, "status");
    std::string fmt = gemmi::read_word(p, &p);
    if (!fmt.empty() && !tr.is_status) {
      tr.min_width = check_format(fmt);
      tr.format = "%" + fmt;
      if (tr.format[1] == '_')
        tr.format[1] = ' ';
    }
    spec.push_back(tr);
  }
  if (spec.empty())
    gemmi::fail("Empty translation spec");
  for (size_t i = 0; i != spec.size(); ++i)
    for (size_t j = i + 1; j != spec.size(); ++j)
      if (spec[i].tag == spec[j].tag)
        gemmi::fail("duplicated output tag: " + spec[i].tag);
  // H, K, L must be the first columns in MTZ and are required in _refln
  for (int i = 2; i != -1; --i)
    if (!gemmi::in_vector_f([&](const Trans& t) { return t.col_idx == i; },
                            spec)) {
      Trans tr;
      tr.col_idx = i;
      tr.tag = "index_";
      tr.tag += ('h' + i); // h, k or l
      spec.insert(spec.begin(), tr);
    }
  return spec;
}

#if defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif

// Get the first (non-zero) DWAVEL corresponding to a sigma column from
// the template. If not found, try value columns (intensity, amplitude
// or anomalous difference).
static
double get_wavelength(const gemmi::Mtz& mtz, const std::vector<Trans>& spec) {
  for (const Trans& tr : spec) {
    const gemmi::Mtz::Column& col = mtz.columns.at(tr.col_idx);
    if (col.type == 'Q' || col.type == 'M' || col.type == 'L') { // sigma
      double wavelength = mtz.dataset(col.dataset_id).wavelength;
      if (wavelength != 0.)
        return wavelength;
    }
  }
  for (const Trans& tr : spec) {
    const gemmi::Mtz::Column& col = mtz.columns.at(tr.col_idx);
    if (col.type == 'F' || col.type == 'J' || col.type == 'D' ||
        col.type == 'K' || col.type == 'G') { // data value
      double wavelength = mtz.dataset(col.dataset_id).wavelength;
      if (wavelength != 0.)
        return wavelength;
    }
  }
  return 0.;
}

static void write_cif(const gemmi::Mtz& mtz, const Options& opt, FILE* out) {
  std::string id = ".";
  if (opt.with_comments) {
    fprintf(out, "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n");
    fprintf(out, "# from: %s\n", opt.mtz_path);
    fprintf(out, "# MTZ title: %s\n", mtz.title.c_str());
    for (size_t i = 0; i != mtz.history.size(); ++i)
      fprintf(out, "# MTZ history #%zu: %s\n", i, mtz.history[i].c_str());
  }
  fprintf(out, "data_%s\n\n", opt.block_name);
  fprintf(out, "_entry.id %s\n\n", id.c_str());

  if (!mtz.batches.empty()) {
    fprintf(out, "_exptl_crystal.id 1\n");
    fprintf(out, "_diffrn.id 1\n");
    fprintf(out, "_diffrn.crystal_id 1\n\n");
  }

  const gemmi::UnitCell& cell = mtz.get_cell();
  fprintf(out, "_cell.entry_id %s\n", id.c_str());
  fprintf(out, "_cell.length_a    %8.3f\n", cell.a);
  fprintf(out, "_cell.length_b    %8.3f\n", cell.b);
  fprintf(out, "_cell.length_c    %8.3f\n", cell.c);
  fprintf(out, "_cell.angle_alpha %8.3f\n", cell.alpha);
  fprintf(out, "_cell.angle_beta  %8.3f\n", cell.beta);
  fprintf(out, "_cell.angle_gamma %8.3f\n\n", cell.gamma);
  double wavelength = std::isnan(opt.wavelength) ? get_wavelength(mtz, opt.spec)
                                                 : opt.wavelength;
  if (wavelength != 0.) {
    fprintf(out, "_diffrn_radiation_wavelength.id 1\n");
    fprintf(out, "_diffrn_radiation_wavelength.wavelength %g\n\n", wavelength);
  }
  if (const gemmi::SpaceGroup* sg = mtz.spacegroup) {
    fprintf(out, "_symmetry.entry_id %s\n", id.c_str());
    fprintf(out, "_symmetry.space_group_name_H-M '%s'\n", sg->hm);
    fprintf(out, "_symmetry.Int_Tables_number %d\n\n", sg->number);
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }
  fprintf(out, "loop_\n");
  if (!mtz.batches.empty()) {
    fprintf(out, "_diffrn_refln.diffrn_id\n");
    fprintf(out, "_diffrn_refln.standard_code\n");
    fprintf(out, "_diffrn_refln.scale_group_code\n");
    fprintf(out, "_diffrn_refln.id\n");
  }
  for (const Trans& tr : opt.spec) {
    const gemmi::Mtz::Column& col = mtz.columns.at(tr.col_idx);
    const gemmi::Mtz::Dataset& ds = mtz.dataset(col.dataset_id);
    if (!mtz.batches.empty())
      fprintf(out, "_diffrn");
    if (opt.with_comments)
      fprintf(out, "_refln.%-26s # %-14s from dataset %s\n",
              tr.tag.c_str(), col.label.c_str(), ds.dataset_name.c_str());
    else
      fprintf(out, "_refln.%s\n", tr.tag.c_str());
  }
  for (int i = 0; i != mtz.nreflections; ++i) {
    if (!mtz.batches.empty())
      fprintf(out, "1 1 1 %d ", i + 1);
    const float* row = &mtz.data[i * mtz.columns.size()];
    if (!opt.value_indices.empty())
      if (std::all_of(opt.value_indices.begin(), opt.value_indices.end(),
                      [&](int n) { return std::isnan(row[n]); }))
        continue;
    bool first = true;
    for (const Trans& tr : opt.spec) {
      if (first)
        first = false;
      else
        std::fputc(' ', out);
      float v = row[tr.col_idx];
      if (tr.is_status) {
        char status = 'x';
        if (opt.sigma_indices.empty() ||
            !std::all_of(opt.sigma_indices.begin(), opt.sigma_indices.end(),
                         [&](int n) { return std::isnan(row[n]); }))
          status = v == 0. ? 'f' : 'o';
        std::fputc(status, out);
      } else if (std::isnan(v)) {
        for (int j = 1; j < tr.min_width; ++j)
          std::fputc(' ', out);
        std::fputc('?', out);
      } else {
        fprintf(out, tr.format.c_str(), v);
      }
    }
    std::fputc('\n', out);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  if (p.options[PrintSpec]) {
    std::printf("                                for merged mtz\n");
    for (const char* line : default_merged_spec)
      std::printf("%s\n", line);
    std::printf("\n                             for unmerged mtz\n");
    for (const char* line : default_unmerged_spec)
      std::printf("%s\n", line);
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
    mtz.read_input(gemmi::MaybeGzipped(mtz_path), true);
    mtz.switch_to_real_hkl();
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR reading %s: %s\n", mtz_path, e.what());
    return 1;
  }
  if (verbose)
    fprintf(stderr, "Writing %s ...\n", cif_path);
  Options options;
  options.mtz_path = mtz_path;
  try {
    std::vector<std::string> lines;
    if (p.options[Spec]) {
      char buf[256];
      const char* spec_path = p.options[Spec].arg;
      gemmi::fileptr_t f_spec = gemmi::file_open_or(spec_path, "r", stdin);
      while (fgets(buf, sizeof(buf), f_spec.get()) != NULL) {
        const char* start = gemmi::skip_blank(buf);
        if (*start != '\0' && *start != '\r' && *start != '\n' && *start != '#')
          lines.emplace_back(start);
      }
    } else if (mtz.batches.empty()) {
      lines.reserve(sizeof(default_merged_spec) / sizeof(char*));
      for (const char* line : default_merged_spec)
        lines.emplace_back(line);
    } else {
      lines.reserve(sizeof(default_unmerged_spec) / sizeof(char*));
      for (const char* line : default_unmerged_spec)
        lines.emplace_back(line);
    }
    options.spec = parse_spec(mtz, lines);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "Problem in translation spec: %s\n", e.what());
    return 2;
  }
  for (const Trans& tr : options.spec) {
    const gemmi::Mtz::Column& col = mtz.columns[tr.col_idx];
    if (p.options[SkipEmpty] && col.type != 'H' && col.type != 'I')
      options.value_indices.push_back(tr.col_idx);
    if (col.type != 'Q' && col.type != 'L' && col.type != 'M')
      options.sigma_indices.push_back(tr.col_idx);
  }
  options.block_name = p.options[BlockName] ? p.options[BlockName].arg : "mtz";
  options.with_comments = !p.options[NoComments];
  if (p.options[Wavelength])
    options.wavelength = std::strtod(p.options[Wavelength].arg, nullptr);
  try {
    gemmi::fileptr_t f_out = gemmi::file_open_or(cif_path, "w", stdout);
    write_cif(mtz, options, f_out.get());
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR writing %s: %s\n", cif_path, e.what());
    return 3;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
