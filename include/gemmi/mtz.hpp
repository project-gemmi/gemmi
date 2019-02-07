// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cassert>
#include <cstdint>   // for int32_t
#include <cstdio>    // for FILE, fread
#include <cstring>   // for memcpy
#include <array>
#include <string>
#include <vector>
#include "atox.hpp"      // for simple_atof, read_word, string_to_int
#include "iterator.hpp"  // for StrideIter
#include "fileutil.hpp"  // for file_open, is_little_endian, ...
#include "symmetry.hpp"  // for find_spacegroup_by_name, SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for fail

#ifdef  __INTEL_COMPILER
# pragma warning push
# pragma warning disable 597  // for StrideIter, the same as in cifdoc.hpp.
#endif

namespace gemmi {

struct Mtz {
  struct Column {
    int dataset_number;
    char type;
    std::string label;
    float min_value = NAN;
    float max_value = NAN;
    std::string source;  // from COLSRC
    Mtz* parent;
    std::size_t idx;

    bool has_data() const { return parent->has_data(); }
    int size() const { return has_data() ? parent->nreflections : 0; }
    int stride() const { return parent->ncol; }
    float& operator[](int n) { return parent->data[idx + n * parent->ncol]; }
    float operator[](int n) const { return parent->data[idx + n * parent->ncol]; }
    float& at(int n) { return parent->data.at(idx + n * parent->ncol); }
    float at(int n) const { return parent->data.at(idx + n * parent->ncol); }
    using iterator = StrideIter<float>;
    iterator begin() {
      assert(parent);
      assert(&parent->columns[idx] == this);
      return iterator({parent->data.data(), idx, (unsigned) parent->ncol});
    }
    iterator end() {
      return iterator({parent->data.data() + parent->data.size(),
                       idx, (unsigned) parent->ncol});
    }
    using const_iterator = StrideIter<const float>;
    const_iterator begin() const { return const_cast<Column*>(this)->begin(); }
    const_iterator end() const { return const_cast<Column*>(this)->end(); }
  };

  struct Dataset {
    int number;
    std::string project_name;
    std::string crystal_name;
    std::string dataset_name;
    UnitCell cell;
    double wavelength = NAN;
  };

  bool same_byte_order = true;
  std::int32_t header_offset = 0;
  std::string version_stamp;
  std::string title;
  int ncol = 0;
  int nreflections = 0;
  int nbatches = 0;
  std::array<int, 5> sort_order;
  double min_1_d2 = NAN;
  double max_1_d2 = NAN;
  float valm = NAN;
  int nsymop = 0;
  UnitCell cell;
  int spacegroup_number = 0;
  std::string spacegroup_name;
  std::vector<Op> symops;
  const SpaceGroup* spacegroup = nullptr;
  std::vector<Dataset> datasets;
  std::vector<Column> columns;
  std::vector<std::string> history;
  std::vector<float> data;

  FILE* warnings = nullptr;

  Mtz() = default;
  Mtz(Mtz&& o) noexcept
    : same_byte_order(o.same_byte_order),
      header_offset(o.header_offset),
      version_stamp(std::move(o.version_stamp)),
      title(std::move(o.title)),
      ncol(o.ncol),
      nreflections(o.nreflections),
      nbatches(o.nbatches),
      sort_order(std::move(o.sort_order)),
      min_1_d2(o.min_1_d2),
      max_1_d2(o.max_1_d2),
      valm(o.valm),
      nsymop(o.nsymop),
      cell(std::move(o.cell)),
      spacegroup_number(o.spacegroup_number),
      spacegroup_name(std::move(o.spacegroup_name)),
      symops(std::move(o.symops)),
      spacegroup(o.spacegroup),
      datasets(std::move(o.datasets)),
      columns(std::move(o.columns)),
      history(std::move(o.history)),
      data(std::move(o.data)),
      warnings(o.warnings) {
    for (Mtz::Column& col : columns)
      col.parent = this;
  }
  Mtz(Mtz const&) = delete;
  Mtz& operator=(Mtz const&) = delete;

  // Functions to use after MTZ headers (and data) is read.

  double resolution_high() const { return std::sqrt(1.0 / max_1_d2); }
  double resolution_low() const  { return std::sqrt(1.0 / min_1_d2); }

  UnitCell& get_cell(int dataset=-1) {
    for (Dataset& ds : datasets)
      if (ds.number == dataset && ds.cell.is_crystal())
        return ds.cell;
    return cell;
  }

  const UnitCell& get_cell(int dataset=-1) const {
    return const_cast<Mtz*>(this)->get_cell(dataset);
  }

  Dataset& last_dataset() {
    if (datasets.empty())
      fail("MTZ dataset not found (missing DATASET header line?).");
    return datasets.back();
  }
  Dataset& dataset(int number) {
    if ((size_t)number < datasets.size() && datasets[number].number == number)
      return datasets[number];
    for (Dataset& d : datasets)
      if (d.number == number)
        return d;
    fail("MTZ file has no dataset number " + std::to_string(number));
  }
  const Dataset& dataset(int number) const {
    return const_cast<Mtz*>(this)->dataset(number);
  }
  int count(const std::string& label) const {
    int n = 0;
    for (const Column& col : columns)
      if (col.label == label)
        ++n;
    return n;
  }
  Column* column_with_label(const std::string& label) {
    for (Column& col : columns)
      if (col.label == label)
        return &col;
    return nullptr;
  }
  Column* column_with_type(char type) {
    for (Column& col : columns)
      if (col.type == type)
        return &col;
    return nullptr;
  }

  bool has_data() const {
    return data.size() == (size_t) ncol * nreflections;
  }

  // Functions for reading MTZ headers and data.

  void toggle_endiannes() {
    same_byte_order = !same_byte_order;
    swap_four_bytes(&header_offset);
  }

  void read_first_bytes(std::FILE* stream) {
    char buf[12] = {0};

    if (std::fread(buf, 1, 12, stream) != 12)
      fail("Could not read the MTZ file (it it empty?)");
    if (buf[0] != 'M' || buf[1] != 'T' || buf[2] != 'Z' || buf[3] != ' ')
      fail("Not an MTZ file - it does not start with 'MTZ '");

    // Bytes 9-12 have so-called machine stamp:
    // "The first 4 half-bytes represent the real, complex, integer and
    // character formats".
    // We don't try to handle all the combinations here, only the two most
    // common: big endian (for all types) and little endian (for all types).
    // BE is denoted by 1 and LE by 4.
    // If we get a value different than 1 and 4 we assume the native byte order.
    if ((buf[9] & 0xf0) == (is_little_endian() ? 0x10 : 0x40))
      toggle_endiannes();

    std::memcpy(&header_offset, buf + 4, 4);
    if (!same_byte_order)
      swap_four_bytes(&header_offset);
  }

  static const char* skip_word(const char* line) {
    while (*line != '\0' && !std::isspace(*line))
      ++line;
    while (std::isspace(*line))
      ++line;
    return line;
  }

  static UnitCell read_cell_parameters(const char* line) {
    double a = simple_atof(line, &line);
    double b = simple_atof(line, &line);
    double c = simple_atof(line, &line);
    double alpha = simple_atof(line, &line);
    double beta = simple_atof(line, &line);
    double gamma = simple_atof(line, &line);
    return UnitCell(a, b, c, alpha, beta, gamma);
  }

  void warn(const std::string& text) {
    if (warnings)
      std::fprintf(warnings, "%s\n", text.c_str());
  }

  void parse_main_header(const char* line) {
    const char* args = skip_word(line);
    switch (ialpha4_id(line)) {
      case ialpha4_id("VERS"):
        version_stamp = rtrim_str(args);
        break;
      case ialpha4_id("TITL"):
        title = rtrim_str(args);
        break;
      case ialpha4_id("NCOL"): {
        ncol = simple_atoi(args, &args);
        nreflections = simple_atoi(args, &args);
        nbatches = simple_atoi(args, &args);
        break;
      }
      case ialpha4_id("CELL"):
        cell = read_cell_parameters(args);
        break;
      case ialpha4_id("SORT"):
        for (int& n : sort_order)
          n = simple_atoi(args, &args);
        break;
      case ialpha4_id("SYMI"): {
        nsymop = simple_atoi(args, &args);
        symops.reserve(nsymop);
        simple_atoi(args, &args); // ignore number of primitive operations
        args = skip_word(skip_blank(args)); // ignore lattice type
        spacegroup_number = simple_atoi(args, &args);
        args = skip_blank(args);
        if (*args != '\'')
          spacegroup_name = read_word(args);
        else if (const char* end = std::strchr(++args, '\''))
          spacegroup_name.assign(args, end);
        // ignore point group which is at the end of args
        break;
      }
      case ialpha4_id("SYMM"):
        symops.push_back(parse_triplet(args));
        break;
      case ialpha4_id("RESO"):
        min_1_d2 = simple_atof(args, &args);
        max_1_d2 = simple_atof(args, &args);
        break;
      case ialpha4_id("VALM"):
        if (*args != 'N') {
          const char* endptr;
          float v = (float) simple_atof(args, &endptr);
          if (*endptr == '\0' || is_space(*endptr))
            valm = v;
          else
            warn("Unexpected VALM value: " + rtrim_str(args));
        }
        break;
      case ialpha4_id("COLU"): {
        columns.emplace_back();
        Column& col = columns.back();
        col.label = read_word(args, &args);
        col.type = read_word(args, &args)[0];
        col.min_value = (float) simple_atof(args, &args);
        col.max_value = (float) simple_atof(args, &args);
        col.dataset_number = simple_atoi(args);
        col.parent = this;
        col.idx = columns.size() - 1;
        break;
      }
      case ialpha4_id("COLS"):
        if (columns.empty())
          fail("MTZ: COLSRC before COLUMN?");
        args = skip_word(args);
        columns.back().source = read_word(args);
        break;
      case ialpha4_id("COLG"):
        // Column group - not used.
        break;
      case ialpha4_id("NDIF"):
        datasets.reserve(simple_atoi(args));
        break;
      case ialpha4_id("PROJ"):
        datasets.emplace_back();
        datasets.back().number = simple_atoi(args, &args);
        datasets.back().project_name = read_word(skip_word(args));
        break;
      case ialpha4_id("CRYS"):
        if (simple_atoi(args, &args) == last_dataset().number)
          datasets.back().crystal_name = read_word(args);
        else
          warn("MTZ CRYSTAL line: unusual numbering.");
        break;
      case ialpha4_id("DATA"):
        if (simple_atoi(args, &args) == last_dataset().number)
          datasets.back().dataset_name = read_word(args);
        else
          warn("MTZ DATASET line: unusual numbering.");
        break;
      case ialpha4_id("DCEL"):
        if (simple_atoi(args, &args) == last_dataset().number)
          datasets.back().cell = read_cell_parameters(args);
        else
          warn("MTZ DCELL line: unusual numbering.");
        break;
      // case("DRES"): not in use yet
      case ialpha4_id("DWAV"):
        if (simple_atoi(args, &args) == last_dataset().number)
          datasets.back().wavelength = simple_atof(args);
        else
          warn("MTZ DWAV line: unusual numbering.");
        break;
      case ialpha4_id("BATCH"):
        // this header is not used for anything?
        break;
      default:
        warn("Unknown header: " + rtrim_str(line));
    }
  }

  void seek_headers(std::FILE* stream) {
    if (std::fseek(stream, 4 * (header_offset - 1), SEEK_SET) != 0)
      fail("Cannot rewind to the MTZ header at byte "
           + std::to_string(header_offset));
  }

  // read headers until END
  void read_main_headers(std::FILE* stream) {
    char buf[81] = {0};
    seek_headers(stream);
    while (std::fread(buf, 1, 80, stream) == 80 &&
           ialpha3_id(buf) != ialpha3_id("END"))
      parse_main_header(buf);
  }

  // read the part between END and MTZENDOFHEADERS
  void read_history_and_later_headers(std::FILE* stream) {
    char buf[81] = {0};
    int n_headers = 0;
    while (std::fread(buf, 1, 80, stream) == 80 &&
           ialpha4_id(buf) != ialpha4_id("MTZE")) {
      if (n_headers != 0) {
        const char* start = skip_blank(buf);
        const char* end = rtrim_cstr(start, start+80);;
        history.emplace_back(start, end);
        --n_headers;
      } else if (ialpha4_id(buf) == ialpha4_id("MTZH")) {
        n_headers = simple_atoi(skip_word(buf));
        if (n_headers < 0 || n_headers > 30) {
          warn("Wrong MTZ: number of headers should be between 0 and 30");
          return;
        }
        history.reserve(n_headers);
      } else if (ialpha4_id(buf) == ialpha4_id("MTZB")) {
        for (int i = 0; i < nbatches; ++i) {
          // TODO: BH, etc
        }
      }
    }
  }

  void setup_spacegroup() {
    spacegroup = find_spacegroup_by_name(spacegroup_name);
    if (!spacegroup)
      warn("MTZ: unrecognized spacegroup name: " + spacegroup_name);
    if (spacegroup->ccp4 != spacegroup_number)
      warn("MTZ: inconsistent spacegroup name and number");
  }

  void read_raw_data(std::FILE* stream) {
    data.resize(ncol * nreflections);
    if (std::fseek(stream, 80, SEEK_SET) != 0)
      fail("Cannot rewind to the MTZ data.");
    size_t n = ncol * nreflections;
    if (std::fread(data.data(), 4, n, stream) != n)
      fail("Error when reading MTZ data");
    if (!same_byte_order)
      for (float& f : data)
        swap_four_bytes(&f);
  }

  void read_all_headers(std::FILE* stream) {
    read_first_bytes(stream);
    read_main_headers(stream);
    read_history_and_later_headers(stream);
    setup_spacegroup();
  }
};


inline Mtz read_mtz_stream(std::FILE* stream, bool with_data) {
  Mtz mtz;
  mtz.read_all_headers(stream);
  if (with_data)
    mtz.read_raw_data(stream);
  return mtz;
}

inline Mtz read_mtz_file(const std::string& path) {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "rb");
  try {
    return read_mtz_stream(f.get(), true);
  } catch (std::runtime_error& e) {
    fail(std::string(e.what()) + ": " + path);
  }
}

#ifdef  __INTEL_COMPILER
# pragma warning pop
#endif

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
