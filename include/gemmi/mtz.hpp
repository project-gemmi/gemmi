// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cstdint>   // for int32_t
#include <cstdio>    // for FILE, fread
#include <cstdlib>   // for strtol
#include <cstring>   // for memcpy
#include <string>
#include <vector>
#include <map>
#include "fileutil.hpp"  // for file_open, is_little_endian, ...
//#include "symmetry.hpp"  // for SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for fail

namespace gemmi {

struct Mtz {
  struct Dataset {
    std::string project;
    std::string dataset_name;
    UnitCell cell;
    double wavelength;
  };

  bool same_byte_order = true;
  std::int32_t header_offset = 0;
  std::string version_stamp;
  std::string title;
  int ncol = 0;
  int nreflections = 0;
  int nbatches = 0;
  UnitCell cell;
  std::map<int, Dataset> datasets;

  void toggle_endiannes() {
    same_byte_order = !same_byte_order;
    swap_four_bytes(&header_offset);
  }

  static constexpr int id3(const char* s) {
    return (s[0] << 16) | (s[1] << 8) | s[2];
  }
  static constexpr int id4(const char* s) {
    return (s[0] << 24) | (s[1] << 16) | (s[2] << 8) | s[3];
  }
  static constexpr int id3u(const char* s) { return id3(s) & ~0x202020; }
  // TODO: use ialpha4_id
  static constexpr int id4u(const char* s) { return id4(s) & ~0x20202020; }

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
    while (*line != '\0' && !isspace(*line))
      ++line;
    while (isspace(*line))
      ++line;
    return line;
  }

  bool parse_header(const char* line) {
    switch (id4u(line)) {
      case id4("VERS"):
        version_stamp = rtrim_str(skip_word(line));
        return true;
      case id4("TITL"):
        title = rtrim_str(skip_word(line));
        return true;
      case id4("NCOL"): {
        char* endptr;
        ncol = std::strtol(skip_word(line), &endptr, 10);
        nreflections = std::strtol(endptr, &endptr, 10);
        nbatches = std::strtol(endptr, &endptr, 10);
        return true;
      }
      case id4("CELL"):
        return true;
      case id4("DCEL"):
        return true;
      case id4("DWAV"):
        return true;
      default:
        return false;
    }
  }

  void read_headers(std::FILE* stream, std::FILE* warnings=nullptr) {
    char buf[81] = {0};
    seek_headers(stream);
    while (std::fread(buf, 1, 80, stream) == 80 && id3u(buf) != id3("END")) {
      bool known = parse_header(buf);
      if (warnings && !known)
        std::fprintf(warnings, "Unknown header: %s\n", rtrim_str(buf).c_str());
    }
  }
};


inline Mtz read_mtz_stream(std::FILE* stream, const std::string& path) {
  Mtz mtz;
  try {
    mtz.read_first_bytes(stream);
    mtz.read_headers(stream);
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
