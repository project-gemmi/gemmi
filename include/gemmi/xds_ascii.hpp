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

struct GEMMI_DLL XdsAscii {
  struct Refl {
    Miller hkl;
    double iobs;
    double sigma;
    double xd;
    double yd;
    double zd;
    double rlp;
    double peak;
    double corr;  // is it always integer?
    double maxc;
    int iset = 1;

    // ZD can be negative for a few reflections
    int frame() const { return (int) std::floor(zd + 1); }
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
  int read_columns = 0;  // doesn't include ITEM_ISET from XSCALE
  int spacegroup_number = 0;
  UnitCell unit_cell;
  Mat33 cell_axes{0.};
  double wavelength = 0.;
  double incident_beam_dir[3] = {0., 0., 0.};
  double oscillation_range = 0.;
  double rotation_axis[3] = {0., 0., 0.};
  double starting_angle = 0.;
  double reflecting_range_esd = 0.;
  int starting_frame = 1;
  int nx = 0;  // detector size - number of pixels
  int ny = 0;
  double qx = 0.;  // pixel size in mm
  double qy = 0.;
  double orgx = 0.;
  double orgy = 0.;
  double detector_distance = 0.;
  std::string generated_by;
  std::string version_str;
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
      if (iset.frame_number_min > iset.frame_number_max)
        continue;
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

  static Vec3 get_normalized(const double (&arr)[3], const char* name) {
    Vec3 vec(arr[0], arr[1], arr[2]);
    double length = vec.length();
    if (length == 0)
      fail("unknown ", name);
    return vec / length;
  }

  // it's already normalized, but just in case normalize it again
  Vec3 get_rotation_axis() const {
    return get_normalized(rotation_axis, "rotation axis");
  }

  // I'm not sure if always |incident_beam_dir| == 1/wavelength
  Vec3 get_s0_direction() const {
    return get_normalized(incident_beam_dir, "incident beam direction");
  }

  bool has_cell_axes() const {
    for (int i = 0; i < 3; ++i)
      if (cell_axes[i][0] == 0 && cell_axes[i][1] == 0 && cell_axes[i][2] == 0)
        return false;
    return true;
  }

  /// Return transition matrix from "Cambridge" frame to XDS frame.
  /// x_xds = M x_cam
  Mat33 calculate_conversion_from_cambridge() const {
    // Cambridge z direction is along the principal rotation axis
    Vec3 z = get_rotation_axis();
    // Cambridge z direction is along beam
    Vec3 x = get_s0_direction();
    Vec3 y = z.cross(x).normalized();
    // beam and rotation axis may not be orthogonal
    x = y.cross(z).normalized();
    return Mat33::from_columns(x, y, z);
  }

  Mat33 get_orientation() const {
    if (!has_cell_axes())
      fail("unknown unit cell axes");
    Vec3 a = cell_axes.row_copy(0);
    Vec3 b = cell_axes.row_copy(1);
    Vec3 c = cell_axes.row_copy(2);
    Vec3 ar = b.cross(c).normalized();
    Vec3 br = c.cross(a);
    Vec3 cr = ar.cross(br).normalized();
    br = cr.cross(ar);
    return Mat33::from_columns(ar, br, cr);
  }

  /// \par p is degree of polarization from range (0,1), as used in XDS.
  void apply_polarization_correction(double p, Vec3 normal);

  /// \par overload is maximally allowed pixel value in a peak (MAXC).
  void eliminate_overloads(double overload) {
    vector_remove_if(data, [&](Refl& r) { return r.maxc > overload; });
  }

  /// \par batchmin lowest allowed batch number.
  void eliminate_batchmin(int batchmin) {
    double minz = batchmin - 1;
    vector_remove_if(data, [&](Refl& r) { return r.zd < minz; });
  }
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
  source_path = source;
  read_columns = 12;
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
        generated_by = read_word(rhs, &rhs);
        version_str = trim_str(rhs);
      } else if (starts_with_ptr(line+1, "SPACE_GROUP_NUMBER=", &rhs)) {
        spacegroup_number = simple_atoi(rhs);
      } else if (starts_with_ptr(line+1, "UNIT_CELL_", &rhs)) {
        if (starts_with_ptr(rhs, "CONSTANTS=", &rhs)) {  // UNIT_CELL_CONSTANTS=
          double par[6];
          parse_numbers_into_array(rhs, line+len, par, line);
          unit_cell.set(par[0], par[1], par[2], par[3], par[4], par[5]);
        } else if (starts_with_ptr(rhs, "A-AXIS=", &rhs)) { // UNIT_CELL_A-AXIS=
          parse_numbers_into_array(rhs, line+len, cell_axes.a[0], line);
        } else if (starts_with_ptr(rhs, "B-AXIS=", &rhs)) { // UNIT_CELL_B-AXIS=
          parse_numbers_into_array(rhs, line+len, cell_axes.a[1], line);
        } else if (starts_with_ptr(rhs, "C-AXIS=", &rhs)) { // UNIT_CELL_C-AXIS=
          parse_numbers_into_array(rhs, line+len, cell_axes.a[2], line);
        }
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
        // INTEGRATE.HKL has read_columns=12, as set above
        if (generated_by == "XSCALE")
          read_columns = 8;
        else if (generated_by == "CORRECT")
          read_columns = 11;
        // check if the columns are what they always are
        if (num < read_columns)
          fail("expected ", std::to_string(read_columns), "+ columns, got:\n", line);
        if (generated_by == "INTEGRATE") {
          copy_line_from_stream(line, 52, stream);
          if (!starts_with(line, "!H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC"))
            fail("unexpected column order in INTEGRATE.HKL");
        } else {
          const char* expected_columns[12] = {
            "H=1", "K=2", "L=3", "IOBS=4", "SIGMA(IOBS)=5",
            "XD=6", "YD=7", "ZD=8", "RLP=9", "PEAK=10", "CORR=11", "MAXC=12"
          };
          for (int i = 0; i < read_columns; ++i) {
            const char* col = expected_columns[i];
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
          isets.emplace_back(1);
          isets.back().wavelength = wavelength;
        }
        for (XdsAscii::Refl& refl : data)
          if (size_t(refl.iset - 1) >= isets.size())
            fail("unexpected ITEM_ISET " + std::to_string(refl.iset));
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
      if (read_columns >= 11) {
        result = fast_from_chars(result.ptr, line+len, r.rlp); // 9
        result = fast_from_chars(result.ptr, line+len, r.peak); // 10
        result = fast_from_chars(result.ptr, line+len, r.corr); // 11
        if (read_columns > 11) {
          result = fast_from_chars(result.ptr, line+len, r.maxc); // 12
        } else {
          r.maxc = 0;
        }
      } else {
        r.rlp = r.peak = r.corr = r.maxc = 0;
      }
      if (result.ec != std::errc())
        fail("failed to parse data line:\n", line);
      if (iset_col >= read_columns) {
        const char* iset_ptr = result.ptr;
        for (int j = read_columns+1; j < iset_col; ++j)
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
