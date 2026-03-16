// Copyright Global Phasing Ltd.

#include "common.h"
#include "iso.h"
#include <gemmi/ccp4.hpp>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using IsoVertices = std::vector<float>;
using IsoSegments = std::vector<uint32_t>;

em::val float_view(const std::vector<float>& values) {
  return em::val(emscripten::typed_memory_view(values.size(), values.data()));
}

em::val uint32_view(const std::vector<uint32_t>& values) {
  return em::val(emscripten::typed_memory_view(values.size(), values.data()));
}

template<typename Grid>
bool extract_isosurface_from_grid(Grid& grid,
                                  double radius,
                                  double x,
                                  double y,
                                  double z,
                                  double isolevel,
                                  const std::string& method,
                                  IsoVertices& vertices,
                                  IsoSegments& segments,
                                  std::string& error) {
  error.clear();
  vertices.clear();
  segments.clear();

  const gemmi::UnitCell& cell = grid.unit_cell;
  const gemmi::Fractional fc = cell.fractionalize({x, y, z});
  const std::array<double, 3> r = {{
    radius / cell.a,
    radius / cell.b,
    radius / cell.c,
  }};
  const std::array<int, 3> grid_min = {{
    static_cast<int>(std::floor((fc.x - r[0]) * grid.nu)),
    static_cast<int>(std::floor((fc.y - r[1]) * grid.nv)),
    static_cast<int>(std::floor((fc.z - r[2]) * grid.nw)),
  }};
  const std::array<int, 3> grid_max = {{
    static_cast<int>(std::floor((fc.x + r[0]) * grid.nu)),
    static_cast<int>(std::floor((fc.y + r[1]) * grid.nv)),
    static_cast<int>(std::floor((fc.z + r[2]) * grid.nw)),
  }};
  const gemmi_wasm::Size3 dims = {{
    grid_max[0] - grid_min[0] + 1,
    grid_max[1] - grid_min[1] + 1,
    grid_max[2] - grid_min[2] + 1,
  }};
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) {
    error = "Grid dimensions are zero along at least one edge";
    return false;
  }

  const size_t point_count = static_cast<size_t>(dims[0]) * dims[1] * dims[2];
  std::vector<float> points;
  std::vector<float> values;
  points.reserve(point_count * 3u);
  values.reserve(point_count);

  for (int i = grid_min[0]; i <= grid_max[0]; ++i) {
    for (int j = grid_min[1]; j <= grid_max[1]; ++j) {
      for (int k = grid_min[2]; k <= grid_max[2]; ++k) {
        gemmi::Position pos = grid.get_position(i, j, k);
        points.push_back(static_cast<float>(pos.x));
        points.push_back(static_cast<float>(pos.y));
        points.push_back(static_cast<float>(pos.z));
        values.push_back(grid.get_value(i, j, k));
      }
    }
  }

  return gemmi_wasm::calculate_isosurface(dims, values, points, isolevel,
                                          method, vertices, segments, error);
}

int16_t read_i16(const std::string& buf, size_t index, bool little_endian) {
  const size_t offset = 2 * index;
  if (offset + 1 >= buf.size())
    throw std::runtime_error("DSN6 header is truncated");
  const unsigned char b0 = static_cast<unsigned char>(buf[offset]);
  const unsigned char b1 = static_cast<unsigned char>(buf[offset + 1]);
  const uint16_t value = little_endian ? static_cast<uint16_t>(b0 | (b1 << 8))
                                       : static_cast<uint16_t>((b0 << 8) | b1);
  return static_cast<int16_t>(value);
}

}  // namespace

class Ccp4Map {
public:
  Ccp4Map(std::string buf): buf_(std::move(buf)) {}

  bool read(bool expand_symmetry=true) {
    try {
      ccp4_.read_ccp4_from_memory(buf_.data(), buf_.size(), "");
      ccp4_.setup(0.f, expand_symmetry ? gemmi::MapSetup::Full
                                       : gemmi::MapSetup::NoSymmetry);
      apply_mode0_scaling();
      if (!has_valid_header_stats())
        ccp4_.hstats = gemmi::calculate_data_statistics(ccp4_.grid.data);
    } catch (std::runtime_error& e) {
      last_error_ = "Failed to read CCP4 map: ";
      last_error_ += e.what();
      return false;
    }
    return true;
  }

  em::val data() const { return float_view(ccp4_.grid.data); }

  bool extract_isosurface(double radius, double x, double y, double z,
                          double isolevel, const std::string& method) {
    return extract_isosurface_from_grid(ccp4_.grid, radius, x, y, z, isolevel,
                                        method, iso_vertices_, iso_segments_,
                                        last_error_);
  }

  em::val isosurface_vertices() const { return float_view(iso_vertices_); }
  em::val isosurface_segments() const { return uint32_view(iso_segments_); }

  int get_nx() const { return ccp4_.grid.nu; }
  int get_ny() const { return ccp4_.grid.nv; }
  int get_nz() const { return ccp4_.grid.nw; }
  double get_mean() const { return ccp4_.hstats.dmean; }
  double get_rms() const { return ccp4_.hstats.rms; }
  std::string get_last_error() const { return last_error_; }
  gemmi::UnitCell get_cell() const { return ccp4_.grid.unit_cell; }

private:
  void apply_mode0_scaling() {
    if (ccp4_.header_i32(4) != 0)
      return;
    if (ccp4_.header_i32(40) != -128 || ccp4_.header_i32(41) != 127)
      return;
    const double min = ccp4_.hstats.dmin;
    const double max = ccp4_.hstats.dmax;
    const float b1 = static_cast<float>((max - min) / 255.0);
    const float b0 = static_cast<float>(0.5 * (min + max + b1));
    for (float& value : ccp4_.grid.data)
      value = b1 * value + b0;
  }

  bool has_valid_header_stats() const {
    const gemmi::DataStats& stats = ccp4_.hstats;
    return std::isfinite(stats.rms) && stats.rms > 0. &&
           stats.dmean >= stats.dmin && stats.dmean <= stats.dmax;
  }

  gemmi::Ccp4<float> ccp4_;
  std::string buf_;
  std::string last_error_;
  IsoVertices iso_vertices_;
  IsoSegments iso_segments_;
};

class Dsn6Map {
public:
  Dsn6Map(std::string buf): buf_(std::move(buf)) {}

  bool read() {
    try {
      parse();
    } catch (std::runtime_error& e) {
      last_error_ = "Failed to read DSN6 map: ";
      last_error_ += e.what();
      return false;
    }
    return true;
  }

  em::val data() const { return float_view(grid_.data); }

  bool extract_isosurface(double radius, double x, double y, double z,
                          double isolevel, const std::string& method) {
    return extract_isosurface_from_grid(grid_, radius, x, y, z, isolevel,
                                        method, iso_vertices_, iso_segments_,
                                        last_error_);
  }

  em::val isosurface_vertices() const { return float_view(iso_vertices_); }
  em::val isosurface_segments() const { return uint32_view(iso_segments_); }

  int get_nx() const { return grid_.nu; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nw; }
  double get_mean() const { return stats_.dmean; }
  double get_rms() const { return stats_.rms; }
  std::string get_last_error() const { return last_error_; }
  gemmi::UnitCell get_cell() const { return grid_.unit_cell; }

private:
  void parse() {
    if (buf_.size() < 512)
      throw std::runtime_error("header is truncated");

    bool little_endian = false;
    if (read_i16(buf_, 18, false) == 100) {
      little_endian = false;
    } else if (read_i16(buf_, 18, true) == 100) {
      little_endian = true;
    } else {
      throw std::runtime_error("Endian swap failed");
    }

    auto header = [this, little_endian](size_t index) {
      return read_i16(buf_, index, little_endian);
    };

    const std::array<int, 3> origin = {{header(0), header(1), header(2)}};
    const std::array<int, 3> n_real = {{header(3), header(4), header(5)}};
    const std::array<int, 3> n_grid = {{header(6), header(7), header(8)}};
    const int cell_scale = header(17);
    const int prod_word = header(15);
    if (cell_scale == 0)
      throw std::runtime_error("invalid cell scale in header");
    if (prod_word == 0)
      throw std::runtime_error("invalid density scale in header");

    const double cell_mult = 1.0 / cell_scale;
    grid_.set_unit_cell(cell_mult * header(9),
                        cell_mult * header(10),
                        cell_mult * header(11),
                        cell_mult * header(12),
                        cell_mult * header(13),
                        cell_mult * header(14));
    grid_.set_size(n_grid[0], n_grid[1], n_grid[2]);

    const float prod = static_cast<float>(prod_word) / 100.f;
    const int plus = header(16);
    size_t offset = 512;
    const std::array<int, 3> n_blocks = {{
      static_cast<int>(std::ceil(n_real[0] / 8.0)),
      static_cast<int>(std::ceil(n_real[1] / 8.0)),
      static_cast<int>(std::ceil(n_real[2] / 8.0)),
    }};

    for (int zz = 0; zz < n_blocks[2]; ++zz) {
      for (int yy = 0; yy < n_blocks[1]; ++yy) {
        for (int xx = 0; xx < n_blocks[0]; ++xx) {
          for (int k = 0; k < 8; ++k) {
            const int z = 8 * zz + k;
            for (int j = 0; j < 8; ++j) {
              const int y = 8 * yy + j;
              for (int i = 0; i < 8; ++i) {
                const int x = 8 * xx + i;
                if (offset >= buf_.size())
                  throw std::runtime_error("data is truncated");
                if (x < n_real[0] && y < n_real[1] && z < n_real[2]) {
                  const unsigned char density_byte =
                      static_cast<unsigned char>(buf_[offset]);
                  const float density = (density_byte - plus) / prod;
                  ++offset;
                  grid_.set_value(origin[0] + x,
                                  origin[1] + y,
                                  origin[2] + z,
                                  density);
                } else {
                  offset += static_cast<size_t>(8 - i);
                  if (offset > buf_.size())
                    throw std::runtime_error("data is truncated");
                  break;
                }
              }
            }
          }
        }
      }
    }

    stats_ = gemmi::calculate_data_statistics(grid_.data);
  }

  gemmi::Grid<float> grid_;
  gemmi::DataStats stats_;
  std::string buf_;
  std::string last_error_;
  IsoVertices iso_vertices_;
  IsoSegments iso_segments_;
};

void add_map() {
  em::class_<Ccp4Map>("Ccp4Map")
    .constructor<std::string>()
    .function("read", &Ccp4Map::read)
    .function("data", &Ccp4Map::data)
    .function("extract_isosurface", &Ccp4Map::extract_isosurface)
    .function("isosurface_vertices", &Ccp4Map::isosurface_vertices)
    .function("isosurface_segments", &Ccp4Map::isosurface_segments)
    .property("nx", &Ccp4Map::get_nx)
    .property("ny", &Ccp4Map::get_ny)
    .property("nz", &Ccp4Map::get_nz)
    .property("mean", &Ccp4Map::get_mean)
    .property("rms", &Ccp4Map::get_rms)
    .property("last_error", &Ccp4Map::get_last_error)
    .property("cell", &Ccp4Map::get_cell)
    ;

  em::class_<Dsn6Map>("Dsn6Map")
    .constructor<std::string>()
    .function("read", &Dsn6Map::read)
    .function("data", &Dsn6Map::data)
    .function("extract_isosurface", &Dsn6Map::extract_isosurface)
    .function("isosurface_vertices", &Dsn6Map::isosurface_vertices)
    .function("isosurface_segments", &Dsn6Map::isosurface_segments)
    .property("nx", &Dsn6Map::get_nx)
    .property("ny", &Dsn6Map::get_ny)
    .property("nz", &Dsn6Map::get_nz)
    .property("mean", &Dsn6Map::get_mean)
    .property("rms", &Dsn6Map::get_rms)
    .property("last_error", &Dsn6Map::get_last_error)
    .property("cell", &Dsn6Map::get_cell)
    ;
}
