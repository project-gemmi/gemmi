// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cstdint>   // for int32_t
#include <cstdio>    // for FILE, fread
#include <cstring>   // for memcpy
#include <array>
#include <string>
#include <vector>
#include "fileutil.hpp"  // for file_open, is_little_endian, ...
#include "symmetry.hpp"  // for find_spacegroup_by_name, SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for fail
#include "atox.hpp"      // for simple_atof, read_word, string_to_int

namespace gemmi {

struct Mtz {
  struct Column {
    int dataset_number;
    char type;
    std::string label;
    float min_value = NAN;
    float max_value = NAN;
    std::string source;
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
  double inv_d2_min = NAN;
  double inv_d2_max = NAN;
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

  FILE* warnings = nullptr;

  double resolution_high() const { return std::sqrt(1.0 / inv_d2_max); }
  double resolution_low() const  { return std::sqrt(1.0 / inv_d2_min); }

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
    fail("MTZ file ha no dataset number " + std::to_string(number));
  }

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

  void seek_headers(std::FILE* stream) {
    if (std::fseek(stream, 4 * (header_offset - 1), SEEK_SET) != 0)
      fail("Cannot rewinding to the MTZ header at byte "
           + std::to_string(header_offset));
  }

  static const char* skip_word(const char* line) {
    while (*line != '\0' && !std::isspace(*line))
      ++line;
    while (std::isspace(*line))
      ++line;
    return line;
  }

  static UnitCell read_cell_parameters(const char* line) {
    UnitCell cell;
    cell.a = simple_atof(line, &line);
    cell.b = simple_atof(line, &line);
    cell.c = simple_atof(line, &line);
    cell.alpha = simple_atof(line, &line);
    cell.beta = simple_atof(line, &line);
    cell.gamma = simple_atof(line, &line);
    return cell;
  }

  void warn(const std::string& text) {
    if (warnings)
      std::fprintf(warnings, "%s\n", text.c_str());
  }

  void parse_main_header(const char* line) {
    const int first_word_id = ialpha4_id(line);
    line = skip_word(line);
    switch (first_word_id) {
      case ialpha4_id("VERS"):
        version_stamp = rtrim_str(line);
        break;
      case ialpha4_id("TITL"):
        title = rtrim_str(line);
        break;
      case ialpha4_id("NCOL"): {
        ncol = simple_atoi(line, &line);
        nreflections = simple_atoi(line, &line);
        nbatches = simple_atoi(line, &line);
        break;
      }
      case ialpha4_id("CELL"):
        cell = read_cell_parameters(line);
        break;
      case ialpha4_id("SORT"):
        for (int& n : sort_order)
          n = simple_atoi(line, &line);
        break;
      case ialpha4_id("SYMI"): {
        nsymop = simple_atoi(line, &line);
        symops.reserve(nsymop);
        simple_atoi(line, &line); // ignore number of primitive operations
        line = skip_word(skip_blank(line)); // ignore lattice type
        spacegroup_number = simple_atoi(line, &line);
        line = skip_blank(line);
        if (*line != '\'')
          spacegroup_name = read_word(line);
        else if (const char* end = std::strchr(++line, '\''))
          spacegroup_name.assign(line, end);
        // ignore point group which is at the end of line
        break;
      }
      case ialpha4_id("SYMM"):
        symops.push_back(parse_triplet(line));
        break;
      case ialpha4_id("RESO"):
        inv_d2_min = simple_atof(line, &line);
        inv_d2_max = simple_atof(line, &line);
        break;
      case ialpha4_id("VALM"):
        if (*line != 'N') {
          const char* endptr;
          float v = (float) simple_atof(line, &endptr);
          if (*endptr == '\0' || isspace_c(*endptr))
            valm = v;
          else
            warn("Unexpected VALM value: " + rtrim_str(line));
        }
        break;
      case ialpha4_id("COLU"): {
        columns.emplace_back();
        Column& col = columns.back();
        col.label = read_word(line, &line);
        col.type = read_word(line, &line)[0];
        col.min_value = (float) simple_atof(line, &line);
        col.max_value = (float) simple_atof(line, &line);
        col.dataset_number = simple_atoi(line);
        break;
      }
      case ialpha4_id("COLS"):
        if (columns.empty())
          fail("MTZ: COLSRC before COLUMN?");
        line = skip_word(line);
        columns.back().source = read_word(line);
        break;
      case ialpha4_id("COLG"):
        // Column group - not used.
        break;
      case ialpha4_id("NDIF"):
        datasets.reserve(simple_atoi(line));
        break;
      case ialpha4_id("PROJ"):
        datasets.emplace_back();
        datasets.back().number = simple_atoi(line, &line);
        datasets.back().project_name = read_word(skip_word(line));
        break;
      case ialpha4_id("CRYS"):
        if (simple_atoi(line, &line) == last_dataset().number)
          datasets.back().crystal_name = read_word(line);
        else
          warn("MTZ CRYSTAL line: unusual numbering.");
        break;
      case ialpha4_id("DATA"):
        if (simple_atoi(line, &line) == last_dataset().number)
          datasets.back().dataset_name = read_word(line);
        else
          warn("MTZ DATASET line: unusual numbering.");
        break;
      case ialpha4_id("DCEL"):
        if (simple_atoi(line, &line) == last_dataset().number)
          datasets.back().cell = read_cell_parameters(line);
        else
          warn("MTZ DCELL line: unusual numbering.");
        break;
      // case("DRES"): not in use yet
      case ialpha4_id("DWAV"):
        if (simple_atoi(line, &line) == last_dataset().number)
          datasets.back().wavelength = simple_atof(line);
        else
          warn("MTZ DWAV line: unusual numbering.");
        break;
      default:
        warn("Unknown header: " + rtrim_str(line));
    }
  }

  void read_headers(std::FILE* stream) {
    char buf[81] = {0};
    seek_headers(stream);
    // main headers -- until END
    while (std::fread(buf, 1, 80, stream) == 80 &&
           ialpha3_id(buf) != ialpha3_id("END"))
      parse_main_header(buf);
    // history and batch headers -- until MTZENDOFHEADERS
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
      }
      // TODO: MTZB, BH
    }
  }

  void setup_spacegroup() {
    spacegroup = find_spacegroup_by_name(spacegroup_name);
    if (!spacegroup)
      warn("MTZ: unrecognized spacegroup name: " + spacegroup_name);
    if (spacegroup->ccp4 != spacegroup_number)
      warn("MTZ: inconsistent spacegroup name and number");
  }

  void read_data(std::FILE* stream) {
    // TODO
  }
};


inline Mtz read_mtz_stream(std::FILE* stream, const std::string& path) {
  Mtz mtz;
  try {
    mtz.read_first_bytes(stream);
    mtz.read_headers(stream);
    mtz.setup_spacegroup();
    mtz.read_data(stream);
  } catch (std::runtime_error& e) {
    fail(std::string(e.what()) + ": " + path);
  }
  return mtz;
}

inline Mtz read_mtz_file(const std::string& path) {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "rb");
  return read_mtz_stream(f.get(), path);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
