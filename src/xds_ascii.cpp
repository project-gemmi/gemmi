// Copyright 2023 Global Phasing Ltd.

#include <gemmi/xds_ascii.hpp>
#include <gemmi/atof.hpp>      // for fast_from_chars
#include <gemmi/atox.hpp>      // for skip_blank, read_word
#include <gemmi/util.hpp>      // for trim_str
#include <gemmi/gz.hpp>
#include <gemmi/math.hpp>

namespace gemmi {

void XdsAscii::gather_iset_statistics() {
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

/// Based on Phil Evans' notes and the literature, see:
/// https://github.com/project-gemmi/gemmi/discussions/248
/// \par p is defined as in XDS (p=0.5 for unpolarized beam).
void XdsAscii::apply_polarization_correction(double p, Vec3 normal) {
  if (!has_cell_axes())
    fail("unknown unit cell axes");
  Mat33 UB = cell_axes.inverse();
  Vec3 rot_axis = get_rotation_axis();
  Vec3 s0_dir = get_s0_direction();
  normal = normal.normalized();
  // The polarization normal is expected to be approx. orthogonal to the beam.
  // dot() is the same as cos_angle() for normalized vectors.
  if (normal.dot(s0_dir) > std::cos(rad(5.0)))
    fail("polarization normal is far from orthogonal to the incident beam");
  // make normal exactly orthogonal to the beam
  normal = s0_dir.cross(normal).cross(s0_dir).normalized();
  // wavevector
  Vec3 s0 = s0_dir / wavelength;
  double s0_m2 = 1. / s0.length_sq();  // s0^-2

  for (Refl& refl : data) {
    double phi = rad(rot_angle(refl));
    Vec3 h(refl.hkl[0], refl.hkl[1], refl.hkl[2]);
    Vec3 r0 = UB.multiply(h);
    Vec3 r = rotate_about_axis(r0, rot_axis, phi);
    Vec3 s = s0 + r;
#if 0
    double two_theta = s0.angle(s);
    // 2d sin(theta) = lambda
    double bragg_angle = std::asin(wavelength / (2 * unit_cell.calculate_d(refl.hkl)));
    printf("(%d %d %d) two-theta %g %g\n",
           refl.hkl[0], refl.hkl[1], refl.hkl[2], deg(two_theta), deg(2 * bragg_angle));
#endif
    // we should have |s| == |s0|, but just in case calculate it separately
    double s_m2 = 1. / s.length_sq();
    // 1 + cos^2(2theta) = 2 * correction for unpolarized beam
    double t = 1 + sq(s.dot(s0)) * s_m2 * s0_m2;
    double polariz_factor = (1 - 2*p) * (1 - sq(normal.dot(s)) * s_m2) + p * t;
    // We assume that the XDS files has polarization correction applied,
    // but for non-polarized beam. So we multiply intensities by P0=t/2
    // and divide by a hopefully more accurate polarization factor.
    double mult = 0.5 * t / polariz_factor;
    refl.iobs *= mult;
    refl.sigma *= mult;
    refl.rlp *= mult;
  }
}

namespace {

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
    fail("failed to parse a number in:\n", line);
  return result.ptr;
}

template<size_t N>
void parse_numbers_into_array(const char* start, const char* end,
                              double (&arr)[N], const char* line) {
  for (double& val : arr)
    start = parse_number_into(start, end, val, line);
}

template<size_t N>
void parse_numbers_into_array(const char* start, const char* end,
                              std::array<double,N>& arr, const char* line) {
  for (double& val : arr)
    start = parse_number_into(start, end, val, line);
}

void parse_numbers_into_vec3(const char* start, const char* end,
                             Vec3& vec, const char* line) {
  for (double* val : {&vec.x, &vec.y, &vec.z})
    start = parse_number_into(start, end, *val, line);
}


} // anonymous namespace

void XdsAscii::read_stream(AnyStream& line_reader, const std::string& source) {
  source_path = source;
  read_columns = 12;
  char line[256];
  size_t len0 = line_reader.copy_line(line, 255);
  if (len0 == 0)
    fail("empty file");
  int iset_col = 0;
  const char xds_ascii_header[] = "!FORMAT=XDS_ASCII    MERGE=";
  char xds_ascii_type = '\0';
  if (starts_with(line, xds_ascii_header)) {
    size_t n = sizeof(xds_ascii_header)-1;
    xds_ascii_type = line[n];
    // !FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=
    if (strncmp(line + n + 5, "    FRIEDEL'S_LAW=", 18) == 0)
      friedels_law = line[50];
  }
  if (!xds_ascii_type && !starts_with(line, "!OUTPUT_FILE=INTEGRATE.HKL"))
    fail("not an XDS_ASCII nor INTEGRATE.HKL file: " + source_path);
  const char* rhs;
  while (size_t len = line_reader.copy_line(line, 255)) {
    if (line[0] == '!') {
      if (starts_with_ptr(line+1, "Generated by ", &rhs)) {
        generated_by = read_word(rhs, &rhs);
        version_str = trim_str(rhs);
      } else if (starts_with_ptr(line+1, "SPACE_GROUP_NUMBER=", &rhs)) {
        spacegroup_number = simple_atoi(rhs);
      } else if (starts_with_ptr(line+1, "UNIT_CELL_", &rhs)) {
        if (starts_with_ptr(rhs, "CONSTANTS=", &rhs)) {  // UNIT_CELL_CONSTANTS=
          parse_numbers_into_array(rhs, line+len, cell_constants, line);
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
        parse_numbers_into_vec3(rhs, line+len, incident_beam_dir, line);
      } else if (starts_with_ptr(line+1, "OSCILLATION_RANGE=", &rhs)) {
        auto result = fast_from_chars(rhs, line+len, oscillation_range);
        if (result.ec != std::errc())
          fail("failed to parse:\n", line);
      } else if (starts_with_ptr(line+1, "ROTATION_AXIS=", &rhs)) {
        parse_numbers_into_vec3(rhs, line+len, rotation_axis, line);
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
        if (xds_ascii_type == 'T')  // merged file
          read_columns = 5;
        else if (generated_by == "XSCALE")
          read_columns = 8;
        else if (generated_by == "CORRECT")
          read_columns = 11;
        // check if the columns are what they always are
        if (num < read_columns)
          fail("expected ", std::to_string(read_columns), "+ columns, got:\n", line);
        if (generated_by == "INTEGRATE") {
          line_reader.copy_line(line, 52);
          if (!starts_with(line, "!H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC"))
            fail("unexpected column order in INTEGRATE.HKL");
        } else {
          const char* expected_columns[12] = {
            "H=1", "K=2", "L=3", "IOBS=4", "SIGMA(IOBS)=5",
            "XD=6", "YD=7", "ZD=8", "RLP=9", "PEAK=10", "CORR=11", "MAXC=12"
          };
          for (int i = 0; i < read_columns; ++i) {
            const char* col = expected_columns[i];
            line_reader.copy_line(line, 42);
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
      if (read_columns >= 8) {
        result = fast_from_chars(result.ptr, line+len, r.xd); // 6
        result = fast_from_chars(result.ptr, line+len, r.yd); // 7
        result = fast_from_chars(result.ptr, line+len, r.zd); // 8
        if (read_columns >= 11) {
          result = fast_from_chars(result.ptr, line+len, r.rlp); // 9
          result = fast_from_chars(result.ptr, line+len, r.peak); // 10
          result = fast_from_chars(result.ptr, line+len, r.corr); // 11
          if (read_columns >= 12) {
            result = fast_from_chars(result.ptr, line+len, r.maxc); // 12
          } else {
            r.maxc = 0;  // 12
          }
        } else {
          r.rlp = r.peak = r.corr = r.maxc = 0;  // 9-11
        }
      } else {
        r.xd = r.yd = r.zd = 0;  // 6-8
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

XdsAscii read_xds_ascii(const std::string& path) {
  XdsAscii xds_ascii;
  xds_ascii.read_input(gemmi::MaybeGzipped(path));
  return xds_ascii;
}

}  // namespace gemmi
