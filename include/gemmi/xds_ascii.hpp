// Copyright 2020 Global Phasing Ltd.
//
// Read XDS files: XDS_ASCII.HKL and INTEGRATE.HKL.

#ifndef GEMMI_XDS_ASCII_HPP_
#define GEMMI_XDS_ASCII_HPP_

#include "input.hpp"     // for AnyStream, FileStream
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for starts_with

namespace gemmi {

/// @file
/// @brief Reader for XDS_ASCII.HKL and INTEGRATE.HKL reflection files.

/// @brief Check if wavelength likely comes from in-house X-ray source.
///
/// Based on Pointless documentation; recognizes Cu, Mo, and Cr K-alpha
/// wavelengths with their typical uncertainties.
/// @param wavelength Wavelength in Angstroms.
/// @return True if wavelength matches Cu, Mo, or Cr.
inline bool likely_in_house_source(double wavelength) {
  return std::fabs(wavelength - 1.5418) < 0.0019 ||
         std::fabs(wavelength - 0.7107) < 0.0002 ||
         std::fabs(wavelength - 2.29) < 0.01;
}

/// @brief Metadata for XDS reflection data (base class).
struct XdsAsciiMetadata {
  /// @brief Properties of one integration set (dataset/sweep).
  struct Iset {
    /// Integration set ID.
    int id;
    /// Input file name for this set.
    std::string input_file;
    /// Wavelength in Angstroms (0 if not specified).
    double wavelength = 0.;
    /// Unit cell constants [a, b, c, alpha, beta, gamma].
    std::array<double,6> cell_constants = {0., 0., 0., 0., 0., 0.};
    /// Minimum frame number (set by gather_iset_statistics).
    int frame_number_min = -1;
    /// Maximum frame number (set by gather_iset_statistics).
    int frame_number_max = -1;
    /// Total number of distinct frames in set.
    int frame_count = -1;
    /// Number of reflections in this set.
    int reflection_count = -1;

    /// @brief Construct Iset with given ID.
    Iset(int id_) : id(id_) {}
  };
  /// Source file path.
  std::string source_path;
  /// Number of columns read from DATA section (0-13, excludes ITEM_ISET from XSCALE).
  int read_columns = 0;
  /// Space group number from XDS_ASCII header.
  int spacegroup_number = 0;
  /// X-ray wavelength in Angstroms.
  double wavelength = 0.;
  /// Unit cell constants [a, b, c, alpha, beta, gamma].
  std::array<double,6> cell_constants = {0., 0., 0., 0., 0., 0.};
  /// Unit cell matrix: rows are a*, b*, c* axes (reciprocal lattice).
  Mat33 cell_axes{0.};
  /// Incident beam direction (usually normalized).
  Vec3 incident_beam_dir;
  /// Oscillation range per frame in degrees.
  double oscillation_range = 0.;
  /// Rotation axis direction.
  Vec3 rotation_axis;
  /// Starting angle for rotation in degrees.
  double starting_angle = 0.;
  /// Mosaicity/reflecting range standard deviation in degrees.
  double reflecting_range_esd = 0.;
  /// Friedel's law assumption: '\0' unknown, 'T'rue, 'F'alse.
  char friedels_law = '\0';
  /// First frame number.
  int starting_frame = 1;
  /// Detector width in pixels.
  int nx = 0;
  /// Detector height in pixels.
  int ny = 0;
  /// Pixel size in x-direction (mm).
  double qx = 0.;
  /// Pixel size in y-direction (mm).
  double qy = 0.;
  /// Detector origin in x-direction (mm).
  double orgx = 0.;
  /// Detector origin in y-direction (mm).
  double orgy = 0.;
  /// Distance from sample to detector (mm).
  double detector_distance = 0.;
  /// Program name that generated the file (e.g., "XDS").
  std::string generated_by;
  /// Version string of generating program.
  std::string version_str;
  /// Integration sets (for multi-sweep data).
  std::vector<Iset> isets;
};

/// @brief Container for XDS reflection data (XDS_ASCII.HKL or INTEGRATE.HKL).
///
/// Stores integration metadata and per-reflection measurements from XDS output.
/// Supports both merged and unmerged data depending on read_columns value.
struct GEMMI_DLL XdsAscii : XdsAsciiMetadata {
  /// @brief One reflection record from XDS data.
  struct Refl {
    /// Miller indices (h, k, l).
    Miller hkl;
    /// Integration set ID (dataset number).
    int iset = 1;
    /// Integrated intensity.
    double iobs;
    /// Standard deviation of iobs.
    double sigma;
    /// Detector position x (mm).
    double xd;
    /// Detector position y (mm).
    double yd;
    /// Frame number (zd) relative to starting frame (can be negative).
    double zd;
    /// Reciprocal lattice point (RLP) value; related to partiality.
    double rlp;
    /// Peak intensity percentage (0-10000 for 0-100%).
    double peak;
    /// Correction factor (Lorentz-polarization etc; usually 100-150).
    double corr;
    /// Maximum pixel value in reflection.
    double maxc;

    /// @brief Get frame number (rounded up from zd).
    /// @return Frame number (zd is negative-friendly).
    int frame() const { return (int) std::floor(zd + 1); }
  };
  /// All reflection records.
  std::vector<Refl> data;

  /// @brief Default constructor.
  XdsAscii() = default;
  /// @brief Construct with existing metadata.
  XdsAscii(const XdsAsciiMetadata& m) : XdsAsciiMetadata(m) {}

  /// @brief Get existing or create new integration set.
  /// @param id Integration set ID.
  /// @return Reference to Iset with given ID.
  Iset& find_or_add_iset(int id) {
    for (Iset& i : isets)
      if (i.id == id)
        return i;
    isets.emplace_back(id);
    return isets.back();
  }
  /// @brief Read XDS file from stream.
  /// @param reader Input stream handler.
  /// @param source File path (for error messages).
  void read_stream(AnyStream& reader, const std::string& source);

  /// @brief Read XDS file from input object (file or stdin).
  /// @tparam T Input object with create_stream() and path() methods.
  /// @param input Input object.
  template<typename T>
  void read_input(T&& input) {
    read_stream(*input.create_stream(), input.path());
  }

  /// @brief Check if data is merged (few columns in XDS file).
  /// @return True if read_columns < 8 (no per-reflection BATCH info).
  bool is_merged() const { return read_columns < 8; }

  /// @brief Calculate frame number statistics for each integration set.
  ///
  /// Sets frame_number_min, frame_number_max, frame_count, and
  /// reflection_count for each Iset.
  void gather_iset_statistics();

  /// @brief Calculate rotation angle for a reflection.
  /// @param refl Reflection record with zd frame number.
  /// @return Rotation angle in degrees.
  double rot_angle(const Refl& refl) const {
    double z = refl.zd - starting_frame + 1;
    return starting_angle + oscillation_range * z;
  }

  /// @brief Get normalized rotation axis.
  /// @return Normalized rotation_axis vector.
  /// @throws gemmi::fail if rotation_axis is zero.
  Vec3 get_rotation_axis() const {
    double length = rotation_axis.length();
    if (length == 0)
      fail("unknown rotation axis");
    return rotation_axis / length;
  }

  /// @brief Get normalized incident beam direction (S0).
  /// @return Normalized incident_beam_dir vector.
  /// @throws gemmi::fail if incident_beam_dir is zero.
  Vec3 get_s0_direction() const {
    double length = incident_beam_dir.length();
    if (length == 0)
      fail("unknown incident beam direction");
    return incident_beam_dir / length;
  }

  /// @brief Check if reciprocal lattice vectors (cell_axes) are set.
  /// @return True if all 3 reciprocal axes are non-zero.
  bool has_cell_axes() const {
    for (int i = 0; i < 3; ++i)
      if (cell_axes[i][0] == 0 && cell_axes[i][1] == 0 && cell_axes[i][2] == 0)
        return false;
    return true;
  }

  /// @brief Calculate transformation matrix from Cambridge frame to XDS frame.
  ///
  /// Cambridge frame: z along rotation axis, x along incident beam.
  /// @return 3x3 matrix such that x_xds = M * x_cambridge.
  /// @throws gemmi::fail if geometry data missing.
  Mat33 calculate_conversion_from_cambridge() const {
    // Cambridge z direction is along the principal rotation axis
    Vec3 z = get_rotation_axis();
    // Cambridge x direction is along beam
    Vec3 x = get_s0_direction();
    Vec3 y = z.cross(x).normalized();
    // beam and rotation axis may not be orthogonal
    x = y.cross(z).normalized();
    return Mat33::from_columns(x, y, z);
  }

  /// @brief Calculate crystal orientation matrix U.
  /// @return 3x3 orientation matrix.
  /// @throws gemmi::fail if cell_axes not set.
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

  /// @brief Apply polarization correction to all intensities and sigmas.
  ///
  /// Based on incident beam direction, rotation axis, and polarization plane.
  /// Assumes XDS file has unpolarized beam correction already applied.
  /// @param p Degree of polarization in [0, 1] (0.5 for unpolarized).
  /// @param normal Normal vector to polarization plane.
  void apply_polarization_correction(double p, Vec3 normal);

  /// @brief Remove reflections with peak pixel value exceeding threshold.
  /// @param overload Maximum allowed MAXC pixel value.
  void eliminate_overloads(double overload) {
    vector_remove_if(data, [&](Refl& r) { return r.maxc > overload; });
  }

  /// @brief Remove reflections with frame number below threshold.
  /// @param batchmin Minimum frame number to keep.
  void eliminate_batchmin(int batchmin) {
    double minz = batchmin - 1;
    vector_remove_if(data, [&](Refl& r) { return r.zd < minz; });
  }
};

/// @brief Read XDS_ASCII file from path.
/// @param path File path.
/// @return Populated XdsAscii object.
inline XdsAscii read_xds_ascii_file(const std::string& path) {
  XdsAscii ret;
  FileStream stream(path.c_str(), "rb");
  ret.read_stream(stream, path);
  return ret;
}

/// @brief Read XDS_ASCII file, handling gzip compression.
/// @param path File path (may be .gz).
/// @return Populated XdsAscii object.
GEMMI_DLL XdsAscii read_xds_ascii(const std::string& path);

} // namespace gemmi
#endif
