// Copyright Global Phasing Ltd.

#include <stdio.h>    // for printf, fprintf, stderr
#include <cstring>    // for memchr, memcpy
#include <algorithm>  // for min, max
#include <exception>
#include "gemmi/atof.hpp"      // for fast_from_chars
//#include "gemmi/calculate.hpp" // for calculate_center_of_mass
#include "gemmi/pdb.hpp"       // for impl::is_record_type
#include "gemmi/read_cif.hpp"  // for read_into_buffer_gz
#include <gemmi/mmread.hpp>    // for coor_format_from_content
#include <gemmi/sprintf.hpp>   // for snprintf_z

#define GEMMI_PROG set
#include "options.h"

namespace {

enum OptionIndex { Bfactor=4, Occupancy, Noise, };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nModify atom attributes in a coordinate file."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Bfactor, 0, "B", "bfactor", Arg::Required,
    "  -B MIN[:MAX]  \tSet isotropic B-factors to a single value or clamp it to MIN:MAX" },
  { Occupancy, 0, "O", "occ", Arg::Required,
    "  -O MIN[:MAX]  \tSet occupancy to a single value or clamp it to MIN:MAX" },
  { 0, 0, 0, 0, 0, 0 }
};

namespace pegtl = tao::pegtl;
namespace rules = gemmi::cif::rules;

struct Context {
  int target_column = -1;
  int current_column = -1;
  int loop_width = -1;
  std::string target_tag;
  std::size_t copy_pos = 0;
  const gemmi::CharArray* arr;
  std::vector<char> output;

  void append(const char* a, const char* b) { output.insert(output.end(), a, b); }
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
    if (in.string() == ctx.target_tag)
      ctx.target_column = ctx.current_column;
    ++ctx.current_column;
  }
};
template<> struct Search<rules::loop_value> {
  template<typename Input> static void apply(const Input& in, Context& ctx) {
    if (ctx.target_column >= 0) {
      if (ctx.loop_width < 0) {
        ctx.loop_width = ctx.current_column;
        ctx.current_column = 0;
      } else {
        ++ctx.current_column;
        if (ctx.current_column == ctx.loop_width)
          ctx.current_column = 0;
      }
      if (ctx.current_column == ctx.target_column) {
        auto pos = in.iterator().byte;
        ctx.append(ctx.arr->data() + ctx.copy_pos, ctx.arr->data() + pos);
        std::string s = "{" + in.string() + "}"; // TODO
        ctx.append(s.data(), s.data() + s.size());
        ctx.copy_pos = pos + in.size();
      }
    }
  }
};

struct EditBuffer {
  gemmi::CoorFormat format = gemmi::CoorFormat::Unknown;
  std::string path;  // for error messages
  gemmi::CharArray arr;

  template<typename Func>
  void modify_values(std::pair<int,int> pdb_columns, const std::string& tag, const Func& func) {
    if (format == gemmi::CoorFormat::Pdb)
      modify_values_in_pdb(pdb_columns, func);
    else
      modify_values_in_mmcif(tag, func);
  }

  template<typename Func>
  void modify_values_in_pdb(std::pair<int,int> cols, const Func& func) {
    char* line = arr.data();
    const char* buf_end = arr.data() + arr.size();
    for (;;) {
      char* eol = (char*) std::memchr(line, '\n', size_t(buf_end - line));
      auto len = (eol ? eol : buf_end) - line;
      if (len >= cols.second &&
          (gemmi::impl::is_record_type(line, "ATOM") ||
           gemmi::impl::is_record_type(line, "HETATM"))) {
        double d = 0.;
        char* start = line + cols.first - 1;
        auto result = gemmi::fast_from_chars(start, line + cols.second, d);
        if (result.ec != std::errc())
          gemmi::fail("failed to parse a number in line:\n", line);
        d = func(d);
        char tmp[7] = {};
        gemmi::snprintf_z(tmp, 7, "%6.2f", d);
        std::memcpy(start, tmp, 6);
      }
      if (!eol)
        break;
      line = eol + 1;
    }
  }

  template<typename Func>
  void modify_values_in_mmcif(const std::string& tag, const Func& func) {
    (void) func;  // TODO
    Context ctx;
    ctx.target_tag = tag;
    ctx.arr = &arr;
    ctx.output.reserve(arr.size() * 11 / 10);
    pegtl::memory_input<> in(arr.data(), arr.size(), path);
    pegtl::parse<rules::file, Search, gemmi::cif::Errors>(in, ctx);
    ctx.append(arr.data() + ctx.copy_pos, arr.data() + arr.size());
    arr.resize(ctx.output.size());
    std::memcpy(arr.data(), ctx.output.data(), ctx.output.size());
    ctx.output.clear();
  }
};

} // anonymous namespace

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
    if (p.options[Bfactor]) {
      double b_min, b_max;
      bool ok = parse_number_or_range(p.options[Bfactor].arg, &b_min, &b_max);
      if (!ok)
        gemmi::fail("argument for -B should be a number or number:number");
      eb.modify_values({61, 66}, "_atom_site.B_iso_or_equiv", [b_min, b_max](double x) {
        return std::min(std::max(x, b_min), b_max);
      });
    }
    if (p.options[Occupancy]) {
      double occ_min, occ_max;
      bool ok = parse_number_or_range(p.options[Occupancy].arg, &occ_min, &occ_max);
      if (!ok)
        gemmi::fail("argument for -O should be a number or number:number");
      eb.modify_values({55, 60}, "_atom_site.occupancy", [occ_min, occ_max](double x) {
        return std::min(std::max(x, occ_min), occ_max);
      });
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
