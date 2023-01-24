// Copyright 2020 Global Phasing Ltd.
//
// Read unmerged XDS files: XDS_ASCII.HKL and INTEGRATE.HKL.

#ifndef GEMMI_XDS_ASCII_HPP_
#define GEMMI_XDS_ASCII_HPP_

#include "atof.hpp"      // for fast_from_chars
#include "atox.hpp"      // for is_space
#include "input.hpp"     // for copy_line_from_stream
#include "fileutil.hpp"  // for file_open
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for starts_with

namespace gemmi {

// from Pointless docs: likely in-house source, in which case
// the unpolarised value is left unchanged (recognised wavelengths
// are CuKalpha 1.5418 +- 0.0019, Mo 0.7107 +- 0.0002, Cr 2.29 +- 0.01)
inline bool likely_in_house_source(double wavelength) {
  return std::fabs(wavelength - 1.5418) < 0.0019 ||
         std::fabs(wavelength - 0.7107) < 0.0002 ||
         std::fabs(wavelength - 2.29) < 0.01;
}

struct XdsAscii {
  struct Refl {
    Miller hkl;
    double iobs;
    double sigma;
    double xd;
    double yd;
    double zd;
    double rlp;
    double peak;
    int iset = 0;

    // I think ZD can't be negative
    int frame() const { return int(zd + 1); }
  };
  struct Iset {
    int id;
    std::string input_file;
    double wavelength = 0.;
    double cell_constants[6] = {0., 0., 0., 0., 0., 0.};
    //statistics set by gather_iset_statistics()
    int frame_number_min = -1;
    int frame_number_max = -1;
    int frame_count = -1;
    int reflection_count = -1;

    Iset(int id_) : id(id_) {}
  };
  std::string source_path;
  int spacegroup_number;
  UnitCell unit_cell;
  double wavelength;
  double incident_beam_dir[3] = {0., 0., 0.};
  double oscillation_range = 0.;
  double rotation_axis[3] = {0., 0., 0.};
  double starting_angle = 0.;
  double reflecting_range_esd = 0.;
  int starting_frame = 1;
  int nx = 0;
  int ny = 0;
  double qx = 0.;
  double qy = 0.;
  double orgx = 0.;
  double orgy = 0.;
  double detector_distance = 0.;
  std::string generated_by;
  std::vector<Iset> isets;
  std::vector<Refl> data;

  Iset& find_or_add_iset(int id) {
    for (Iset& i : isets)
      if (i.id == id)
        return i;
    isets.emplace_back(id);
    return isets.back();
  }
  template<typename Stream>
  void read_stream(Stream&& stream, const std::string& source);

  template<typename T>
  inline void read_input(T&& input) {
    if (input.is_stdin()) {
      read_stream(FileStream{stdin}, "stdin");
    } else if (input.is_compressed()) {
      read_stream(input.get_uncompressing_stream(), input.path());
    } else {
      auto f = file_open(input.path().c_str(), "r");
      read_stream(FileStream{f.get()}, input.path());
    }
  }

  void gather_iset_statistics() {
    for (Iset& iset : isets) {
      iset.frame_number_min = INT_MAX;
      iset.frame_number_max = 0;
      for (const XdsAscii::Refl& refl : data)
        if (refl.iset == iset.id) {
          ++iset.reflection_count;
          int frame = refl.frame();
          iset.frame_number_min = std::min(iset.frame_number_min, frame);
          iset.frame_number_max = std::max(iset.frame_number_max, frame);
        }
      std::vector<uint8_t> frames(iset.frame_number_max - iset.frame_number_min + 1);
      for (const XdsAscii::Refl& refl : data)
        if (refl.iset == iset.id)
          frames[refl.frame() - iset.frame_number_min] = 1;
      iset.frame_count = 0;
      for (uint8_t f : frames)
        iset.frame_count += f;
    }
  }

  double rot_angle(const Refl& refl) const {
    double z = refl.zd - starting_frame + 1;
    return starting_angle + oscillation_range * z;
  }

  void apply_polarization_correction(double fraction, Vec3 normal);
};

template<size_t N>
bool starts_with_ptr(const char* a, const char (&b)[N], const char** endptr) {
  if (std::strncmp(a, b, N-1) != 0)
    return false;
  *endptr = a + N - 1;
  return true;
}

template<size_t N>
bool starts_with_ptr_b(const char* a, const char (&b)[N], const char** endptr) {
  return starts_with_ptr<N>(skip_blank(a), b, endptr);
}

inline const char* parse_number_into(const char* start, const char* end,
                                     double& val, const char* line) {
  auto result = fast_from_chars(start, end, val);
  if (result.ec != std::errc())
    fail("failed to parse number in:\n", line);
  return result.ptr;
}

template<int N>
void parse_numbers_into_array(const char* start, const char* end,
                              double (&arr)[N], const char* line) {
  for (int i = 0; i < N; ++i) {
    auto result = fast_from_chars(start, end, arr[i]);
    if (result.ec != std::errc())
      fail("failed to parse number #", i+1, " in:\n", line);
    start = result.ptr;
  }
}

template<typename Stream>
void XdsAscii::read_stream(Stream&& stream, const std::string& source) {
  static const char* expected_columns[10] = {
    "H=1", "K=2", "L=3", "IOBS=4", "SIGMA(IOBS)=5", "XD=6", "YD=7", "ZD=8",
    "RLP=9", "PEAK=10"
  };
  source_path = source;
  char line[256];
  size_t len0 = copy_line_from_stream(line, 255, stream);
  int iset_col = 0;
  if (len0 == 0 || !(starts_with(line, "!FORMAT=XDS_ASCII    MERGE=FALSE") ||
                    (starts_with(line, "!OUTPUT_FILE=INTEGRATE.HKL"))))
    fail("not an unmerged XDS_ASCII nor INTEGRATE.HKL file: " + source_path);
  const char* rhs;
  while (size_t len = copy_line_from_stream(line, 255, stream)) {
    if (line[0] == '!') {
      if (starts_with_ptr(line+1, "Generated by ", &rhs)) {
        generated_by = read_word(rhs);
      } else if (starts_with_ptr(line+1, "SPACE_GROUP_NUMBER=", &rhs)) {
        spacegroup_number = simple_atoi(rhs);
      } else if (starts_with_ptr(line+1, "UNIT_CELL_CONSTANTS=", &rhs)) {
        double par[6];
        parse_numbers_into_array(rhs, line+len, par, line);
        unit_cell.set(par[0], par[1], par[2], par[3], par[4], par[5]);
      } else if (starts_with_ptr(line+1, "REFLECTING_RANGE_E.S.D.=", &rhs)) {
        auto result = fast_from_chars(rhs, line+len, reflecting_range_esd);
        if (result.ec != std::errc())
          fail("failed to parse mosaicity:\n", line);
      } else if (starts_with_ptr(line+1, "X-RAY_WAVELENGTH=", &rhs)) {
        auto result = fast_from_chars(rhs, line+len, wavelength);
        if (result.ec != std::errc())
          fail("failed to parse wavelength:\n", line);
      } else if (starts_with_ptr(line+1, "INCIDENT_BEAM_DIRECTION=", &rhs)) {
        parse_numbers_into_array(rhs, line+len, incident_beam_dir, line);
      } else if (starts_with_ptr(line+1, "OSCILLATION_RANGE=", &rhs)) {
        auto result = fast_from_chars(rhs, line+len, oscillation_range);
        if (result.ec != std::errc())
          fail("failed to parse:\n", line);
      } else if (starts_with_ptr(line+1, "ROTATION_AXIS=", &rhs)) {
        parse_numbers_into_array(rhs, line+len, rotation_axis, line);
      } else if (starts_with_ptr(line+1, "STARTING_ANGLE=", &rhs)) {
        auto result = fast_from_chars(rhs, line+len, starting_angle);
        if (result.ec != std::errc())
          fail("failed to parse:\n", line);
      } else if (starts_with_ptr(line+1, "STARTING_FRAME=", &rhs)) {
        starting_frame = simple_atoi(rhs);
      } else if (starts_with_ptr(line+1, " ISET= ", &rhs)) {
        const char* endptr;
        int id = simple_atoi(rhs, &endptr);
        XdsAscii::Iset& iset = find_or_add_iset(id);
        endptr = skip_blank(endptr);
        if (starts_with_ptr(endptr, "INPUT_FILE=", &rhs)) {
          iset.input_file = read_word(rhs);
        } else if (starts_with_ptr(endptr, "X-RAY_WAVELENGTH=", &rhs)) {
          double w;
          auto result = fast_from_chars(rhs, line+len, w);
          if (result.ec != std::errc())
            fail("failed to parse iset wavelength:\n", line);
          iset.wavelength = w;
        } else if (starts_with_ptr(endptr, "UNIT_CELL_CONSTANTS=", &rhs)) {
          parse_numbers_into_array(rhs, line+len, iset.cell_constants, line);
        }
      } else if (starts_with_ptr(line+1, "NX=", &rhs)) {
        const char* endptr;
        nx = simple_atoi(rhs, &endptr);
        if (starts_with_ptr_b(endptr, "NY=", &rhs))
          ny = simple_atoi(rhs, &endptr);
        if (starts_with_ptr_b(endptr, "QX=", &rhs))
          endptr = parse_number_into(rhs, line+len, qx, line);
        if (starts_with_ptr_b(endptr, "QY=", &rhs))
          parse_number_into(rhs, line+len, qy, line);
      } else if (starts_with_ptr(line+1, "ORGX=", &rhs)) {
        const char* endptr = parse_number_into(rhs, line+len, orgx, line);
        if (starts_with_ptr_b(endptr, "ORGY=", &rhs))
          endptr = parse_number_into(rhs, line+len, orgy, line);
        if (starts_with_ptr_b(endptr, "DETECTOR_DISTANCE=", &rhs))
          parse_number_into(rhs, line+len, detector_distance, line);
      } else if (starts_with_ptr(line+1, "NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=", &rhs)) {
        int num = simple_atoi(rhs);
        if (num < 10)
          fail("expected 10+ columns, got:\n", line);
        if (generated_by == "INTEGRATE") {
          copy_line_from_stream(line, 42, stream);
          if (std::strncmp(line, "!H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK", 41) != 0)
            fail("unexpected column order in INTEGRATE.HKL");
        } else {
          for (const char* col : expected_columns) {
            copy_line_from_stream(line, 42, stream);
            if (std::strncmp(line, "!ITEM_", 6) != 0 ||
                std::strncmp(line+6, col, std::strlen(col)) != 0)
              fail("column !ITEM_" + std::string(col), " not found.");
          }
        }
      } else if (starts_with_ptr(line+1, "ITEM_ISET=", &rhs)) {
        iset_col = simple_atoi(rhs);
      } else if (starts_with(line+1, "END_OF_DATA")) {
        if (isets.empty()) {
          for (XdsAscii::Refl& refl : data) {
            if (refl.iset != 0 && refl.iset != 1)
              fail("unexpected ITEM_ISET " + std::to_string(refl.iset));
            refl.iset = 1;
          }
          isets.emplace_back(1);
          isets.back().wavelength = wavelength;
        }
        return;
      }
    } else {
      data.emplace_back();
      XdsAscii::Refl& r = data.back();
      const char* p = line;
      for (int i = 0; i < 3; ++i)
        r.hkl[i] = simple_atoi(p, &p);
      auto result = fast_from_chars(p, line+len, r.iobs); // 4
      result = fast_from_chars(result.ptr, line+len, r.sigma); // 5
      result = fast_from_chars(result.ptr, line+len, r.xd); // 6
      result = fast_from_chars(result.ptr, line+len, r.yd); // 7
      result = fast_from_chars(result.ptr, line+len, r.zd); // 8
      result = fast_from_chars(result.ptr, line+len, r.rlp); // 9
      result = fast_from_chars(result.ptr, line+len, r.peak); // 10
      if (result.ec != std::errc())
        fail("failed to parse data line:\n", line);
      if (iset_col > 10) {
        const char* iset_ptr = result.ptr;
        for (int j = 11; j < iset_col; ++j)
          iset_ptr = skip_word(skip_blank(iset_ptr));
        r.iset = simple_atoi(iset_ptr);
      }
    }
  }
  fail("incorrect or unfinished file: " + source_path);
}

inline XdsAscii read_xds_ascii_file(const std::string& path) {
  auto f = file_open(path.c_str(), "r");
  XdsAscii ret;
  ret.read_stream(FileStream{f.get()}, path);
  return ret;
}

} // namespace gemmi
#endif
