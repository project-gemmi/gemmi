// Copyright Global Phasing Ltd.

#include <stdio.h>    // for printf, fprintf, stderr
#include <cstring>    // for memchr, memcpy
#include <exception>
#include <functional> // for function
#include <random>
#include "gemmi/atof.hpp"      // for fast_from_chars
//#include "gemmi/calculate.hpp" // for calculate_center_of_mass
#include "gemmi/cif.hpp"       // for cif::rules
#include "gemmi/pdb.hpp"       // for is_record_type4
#include "gemmi/read_cif.hpp"  // for read_into_buffer_gz
#include <gemmi/mmread.hpp>    // for read_structure_from_memory, coor_format_from_content
#include <gemmi/select.hpp>    // for Selection
#include <gemmi/sprintf.hpp>   // for snprintf_z

#define GEMMI_PROG set
#include "options.h"

namespace {

enum OptionIndex { Bfactor=4, Occupancy, Noise, Shift, Select };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nModify atom attributes in a coordinate file."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Bfactor, 0, "B", "bfactor", Arg::NumberOrRange,
    "  -B MIN[:MAX]  \tSet isotropic B-factors to a single value MIN"
    "\n\tor clamp them to MIN:MAX." },
  { Occupancy, 0, "O", "occ", Arg::NumberOrRange,
    "  -O MIN[:MAX]  \tSet occupancies to a single value or clamp them to MIN:MAX." },
  { Noise, 0, "", "noise", Arg::Float,
    "  --noise M  \tAdd random shifts, uniform in (-M,M), to x,y,z." },
  { Shift, 0, "", "shift", Arg::Float3,
    "  --shift='DX DY DZ' \tTranslate the model coordinates (units: Angstroms)." },
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL \tApply transformations only to selected atoms (MMDB syntax)." },
  { 0, 0, 0, 0, 0, 0 }
};

namespace pegtl = tao::pegtl;
namespace rules = gemmi::cif::rules;

struct Context {
  int target_column = -1;
  int current_column = -1;
  int loop_width = -1;
  int nvalues = 0;
  size_t atom_index = 0;
  std::string target_tag;
  std::size_t copied_until_byte = 0;
  std::vector<bool>* picked_ptr = nullptr;
  const gemmi::CharArray* arr;
  std::vector<char> output;
  std::function<void(double*)> func;
  double values[3];

  void append(const char* a, const char* b) { output.insert(output.end(), a, b); }
  int get_n() const { return current_column - target_column; }
};

template<typename Rule> struct Search : pegtl::nothing<Rule> {};

template<> struct Search<rules::str_loop> {
  template<typename Input> static void apply(const Input&, Context& ctx) {
    ctx.current_column = 0;
    ctx.target_column = -1;
  }
};
template<> struct Search<rules::loop_tag> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    if (ctx.target_column >= 0 && ctx.get_n() < ctx.nvalues) {
      std::string tag = ctx.target_tag;
      tag.back() += ctx.get_n();  // 'x' -> 'y'/'z'
      if (in.string() != tag)
        gemmi::fail(tag + " not found where expected");
    } else if (in.string() == ctx.target_tag) {
      ctx.target_column = ctx.current_column;
    }
    ++ctx.current_column;
  }
};
template<> struct Search<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    if (ctx.target_column < 0)
      return;
    if (ctx.loop_width < 0) {  // first value in the loop
      ctx.loop_width = ctx.current_column;
      ctx.current_column = 0;
      ctx.atom_index = 0;
    } else {
      ++ctx.current_column;
      if (ctx.current_column == ctx.loop_width) {
        ctx.current_column = 0;
        ctx.atom_index++;
      }
    }
    int n = ctx.get_n();
    if (n >= 0 && n < ctx.nvalues) {
      if (ctx.picked_ptr && !ctx.picked_ptr->at(ctx.atom_index))
        return;
      auto pos = in.iterator().byte;
      if (n == 0) {
        ctx.append(ctx.arr->data() + ctx.copied_until_byte, ctx.arr->data() + pos);
        ctx.copied_until_byte = pos;
      }
      auto result = gemmi::fast_from_chars(in.begin(), in.end(), ctx.values[n]);
      if (result.ec != std::errc() && result.ec != std::errc::result_out_of_range)
        gemmi::fail("line ", std::to_string(in.iterator().line),
                    ": not a number:" + in.string());
      if (n + 1 == ctx.nvalues) {
        ctx.func(ctx.values);
        size_t old_width = in.size() + pos - ctx.copied_until_byte;
        size_t new_width = n;
        char tmp[16] = {};
        tmp[0] = ' ';
        for (int i = 0; i <= n; ++i) {
          auto len = gemmi::snprintf_z(tmp+1, 15, "%.3f", ctx.values[i]);
          ctx.append(tmp + (i == 0 ? 1 : 0), tmp + 1 + len);
          new_width += len;
        }
        if (old_width > new_width)
          ctx.output.insert(ctx.output.end(), old_width - new_width, ' ');
      }
      ctx.copied_until_byte = pos + in.size();
    }
  }
};

// modify occupancy or tempFactor, Real(6.2)
template<typename Func> void modify_line_6_2(char* line, int column, Func& func) {
  double d = 0.;
  char* start = line + column - 1;
  auto result = gemmi::fast_from_chars(start, start + 6, d);
  if (result.ec != std::errc() && result.ec != std::errc::result_out_of_range)
    gemmi::fail("failed to parse a number in line:\n", line);
  func(&d);
  char tmp[7] = {};
  gemmi::snprintf_z(tmp, 7, "%6.2f", d);
  std::memcpy(start, tmp, 6);
}

template<typename Func> void modify_line_xyz(char* line, int column, Func& func) {
  double xyz[3];
  char* start = line + column - 1;
  for (int i = 0; i < 3; ++i, start += 8) {
    auto result = gemmi::fast_from_chars(start, start + 8, xyz[i]);
    if (result.ec != std::errc() && result.ec != std::errc::result_out_of_range)
      gemmi::fail("failed to parse a number in line:\n", line);
  }
  func(xyz);
  start -= 3 * 8;
  char tmp[9] = {};
  for (int i = 0; i < 3; ++i, start += 8) {
    gemmi::snprintf_z(tmp, 9, "%8.3f", xyz[i]);
    std::memcpy(start, tmp, 8);
  }
}

struct EditBuffer {
  gemmi::CoorFormat format = gemmi::CoorFormat::Unknown;
  std::string path;  // for error messages
  gemmi::CharArray arr;
  std::vector<bool> picked;

  template<typename Func>
  void modify_pdb(int min_length, const Func& modify_line) {
    char* line = arr.data();
    const char* buf_end = arr.data() + arr.size();
    size_t line_counter = 0;
    for (;;) {
      char* eol = (char*) std::memchr(line, '\n', size_t(buf_end - line));
      auto len = (eol ? eol : buf_end) - line;
      if (len > 4 && (gemmi::is_record_type4(line, "ATOM") ||
                      gemmi::is_record_type4(line, "HETATM"))) {
        if (picked.size() <= line_counter || picked[line_counter++])
          if (len >= min_length)
            modify_line(line);
      }
      if (!eol)
        break;
      line = eol + 1;
    }
  }

  template<typename Func>
  void modify_mmcif(int nvalues, const std::string& tag, const Func& func) {
    Context ctx;
    ctx.target_tag = tag;
    ctx.nvalues = nvalues;
    ctx.arr = &arr;
    ctx.output.reserve(arr.size() * 11 / 10);
    ctx.func = func;
    if (!picked.empty())
      ctx.picked_ptr = &picked;
    pegtl::memory_input<> in(arr.data(), arr.size(), path);
    // parse and copy (Search<>::apply() does copying) modified input to ctx.output
    pegtl::parse<rules::file, Search, gemmi::cif::Errors>(in, ctx);
    // copy the rest of input, unchanged, to ctx.output
    ctx.append(arr.data() + ctx.copied_until_byte, arr.data() + arr.size());
    // copy back: ctx.output to arr
    arr.resize(ctx.output.size());
    std::memcpy(arr.data(), ctx.output.data(), ctx.output.size());
    ctx.output.clear();
  }

  template<typename Func> void modify_xyz(Func& func) {
    if (format == gemmi::CoorFormat::Pdb)
      modify_pdb(54, [&](char* line) { modify_line_xyz(line, 31, func); });
    else
      modify_mmcif(3, "_atom_site.Cartn_x", func);
  }
};

} // anonymous namespace

struct Clamp {
  double xmin, xmax;
  void operator()(double* x) const {
    if (xmin == xmax)
      *x = xmin;
    else
      *x = gemmi::clamp(*x, xmin, xmax);
  }
};

struct NoiseAdder {
  NoiseAdder(std::random_device::result_type seed, double m) : e1(seed), dist(-m, m) {}
  void operator()(double* xyz) {
    for (int i = 0; i < 3; ++i)
      xyz[i] += dist(e1);
  }
private:
  std::default_random_engine e1;
  std::uniform_real_distribution<double> dist;
};

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  int verbosity = p.options[Verbose].count();
  EditBuffer eb;
  eb.path = p.coordinate_input_file(0);
  const char* output = p.nonOption(1);
  try {
    eb.arr = gemmi::read_into_buffer_gz(eb.path);
    eb.format = gemmi::coor_format_from_content(eb.arr.data(), eb.arr.data() + eb.arr.size());
    if (p.options[Select]) {
      gemmi::Selection sel(p.options[Select].arg);
      gemmi::Structure st = gemmi::read_structure_from_memory(eb.arr.data(), eb.arr.size(),
                                                              eb.path, eb.format);
      // Here we assume that the order of atoms in the structure is the same
      // as the order of lines in the file. This is true unless the file has
      // atoms from the same residue in non-consecutive order, or atoms from the
      // same model are non-consecutive (not possible in the PDB format).
      for (const gemmi::Model& model : st.models) {
        bool model_matches = sel.matches(model);
        for (const gemmi::Chain& chain : model.chains) {
          bool chain_matches = model_matches && sel.matches(chain);
          for (const gemmi::Residue& res : chain.residues) {
            bool res_matches = chain_matches && sel.matches(res);
            for (const gemmi::Atom& atom : res.atoms)
              eb.picked.push_back(res_matches && sel.matches(atom));
          }
        }
      }
    }
    if (p.options[Bfactor]) {
      double b_min, b_max;
      parse_number_or_range(p.options[Bfactor].arg, &b_min, &b_max);
      Clamp clamp{b_min, b_max};
      if (eb.format == gemmi::CoorFormat::Pdb)
        eb.modify_pdb(66, [&](char* line) { modify_line_6_2(line, 61, clamp); });
      else
        eb.modify_mmcif(1, "_atom_site.B_iso_or_equiv", clamp);
    }
    if (p.options[Occupancy]) {
      double occ_min, occ_max;
      parse_number_or_range(p.options[Occupancy].arg, &occ_min, &occ_max);
      Clamp clamp{occ_min, occ_max};
      if (eb.format == gemmi::CoorFormat::Pdb)
        eb.modify_pdb(60, [&](char* line) { modify_line_6_2(line, 55, clamp); });
      else
        eb.modify_mmcif(1, "_atom_site.occupancy", clamp);
    }
    if (p.options[Noise]) {
      double m = std::atof(p.options[Noise].arg);
      std::random_device r;
      NoiseAdder func(r(), m);
      eb.modify_xyz(func);
    }
    if (p.options[Shift]) {
      auto v = parse_blank_separated_numbers(p.options[Shift].arg);
      auto func = [v](double *xyz) {
        for (int i = 0; i < 3; ++i)
          xyz[i] += v[i];
      };
      eb.modify_xyz(func);
    }
    gemmi::fileptr_t f = gemmi::file_open(output, "wb");
    if (std::fwrite(eb.arr.data(), eb.arr.size(), 1, f.get()) != 1) {
      perror("Writing file failed");
      return 1;
    }
  } catch (std::exception& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  if (verbosity > 0)
    printf("Done.\n");
  return 0;
}
