// Copyright 2020 Global Phasing Ltd.
//
// Read XDS_ASCII.HKL. For now, only unmerged files are read.

#ifndef GEMMI_XDS_ASCII_HPP_
#define GEMMI_XDS_ASCII_HPP_

#include "atof.hpp"      // for fast_from_chars
#include "atox.hpp"      // for is_space
#include "input.hpp"     // for copy_line_from_stream
#include "fileutil.hpp"  // for file_open
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for fast_from_chars

namespace gemmi {

struct XdsAscii {
  struct Refl {
    Miller hkl;
    double iobs;
    double sigma;
    double xd;
    double yd;
    double zd;
    int iset = 0;
  };
  struct Iset {
    int n;
    double wavelength;
  };
  std::string source_path;
  int spacegroup_number;
  UnitCell unit_cell;
  double wavelength;
  std::vector<Iset> isets;
  std::vector<Refl> data;
};

template<typename Input>
XdsAscii read_xds_ascii_stream(Input&& infile, const std::string& source) {
  static const char* expected_columns[8] = {
    "H=1", "K=2", "L=3", "IOBS=4", "SIGMA(IOBS)=5", "XD=6", "YD=7", "ZD=8",
  };
  XdsAscii ret;
  ret.source_path = source;
  char line[96];
  size_t len0 = copy_line_from_stream(line, 95, infile);
  int iset_col = 0;
  if (len0 == 0 || !starts_with(line, "!FORMAT=XDS_ASCII    MERGE=FALSE"))
    fail("not an unmerged XDS_ASCII file: " + source);
  while (size_t len = copy_line_from_stream(line, 95, infile)) {
    if (line[0] == '!') {
      if (starts_with(line+1, "SPACE_GROUP_NUMBER=")) {
        ret.spacegroup_number = simple_atoi(line + 20);
      } else if (starts_with(line+1, "UNIT_CELL_CONSTANTS=")) {
        double par[6];
        const char* start = line + 21;
        for (int i = 0; i < 6; ++i) {
          auto result = fast_from_chars(start, line+len, par[i]);
          if (result.ec != std::errc())
            fail("failed to parse cell constants:\n", line);
          start = result.ptr;
        }
        ret.unit_cell.set(par[0], par[1], par[2], par[3], par[4], par[5]);
      } else if (starts_with(line+1, "X-RAY_WAVELENGTH=")) {
        auto result = fast_from_chars(line+18, line+len, ret.wavelength);
        if (result.ec != std::errc())
          fail("failed to parse wavelength:\n", line);
      } else if (starts_with(line+1, " ISET= ")) {
        const char* endptr;
        int iset = simple_atoi(line + 7, &endptr);
        endptr = skip_blank(endptr);
        double w;
        if (starts_with(endptr, "X-RAY_WAVELENGTH=")) {
          auto result = fast_from_chars(endptr+17, line+len, w);
          if (result.ec != std::errc())
            fail("failed to parse iset wavelength:\n", line);
          ret.isets.push_back({iset, w});
        }
      } else if (starts_with(line+1, "NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=")) {
        int num = simple_atoi(line + 37);
        if (num < 8)
          fail("expected 8+ columns, got:\n", line);
        for (const char* col : expected_columns) {
          copy_line_from_stream(line, 40, infile);
          if (std::strncmp(line, "!ITEM_", 6) != 0 ||
              std::strncmp(line+6, col, std::strlen(col)) != 0)
            fail("column !ITEM_" + std::string(col), " not found.");
        }
      } else if (starts_with(line+1, "ITEM_ISET=")) {
        iset_col = simple_atoi(line + 13);
      } else if (starts_with(line+1, "END_OF_DATA")) {
        return ret;
      }
    } else {
      ret.data.emplace_back();
      XdsAscii::Refl& r = ret.data.back();
      const char* p = line;
      for (int i = 0; i < 3; ++i)
        r.hkl[i] = simple_atoi(p, &p);
      auto result = fast_from_chars(p, line+len, r.iobs); // 4
      result = fast_from_chars(result.ptr, line+len, r.sigma); // 5
      result = fast_from_chars(result.ptr, line+len, r.xd); // 6
      result = fast_from_chars(result.ptr, line+len, r.yd); // 7
      result = fast_from_chars(result.ptr, line+len, r.zd); // 8
      if (result.ec != std::errc())
        fail("failed to parse data line:\n", line);
      if (iset_col > 8) {
        const char* iset_ptr = result.ptr;
        for (int j = 9; j < iset_col; ++j)
          iset_ptr = skip_word(skip_blank(iset_ptr));
        r.iset = simple_atoi(iset_ptr);
      }
    }
  }
  fail("incorrect or unfinished file: " + source);
}

inline XdsAscii read_xds_ascii_file(const std::string& path) {
  auto f = file_open(path.c_str(), "r");
  return read_xds_ascii_stream(FileStream{f.get()}, path);
}

template<typename T>
inline XdsAscii read_xds_ascii(T&& input) {
  if (input.is_stdin())
    return read_xds_ascii_stream(FileStream{stdin}, "stdin");
  if (input.is_compressed())
    return read_xds_ascii_stream(input.get_uncompressing_stream(), input.path());
  return read_xds_ascii_file(input.path());
}

} // namespace gemmi
#endif
