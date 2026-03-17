// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/ccp4.hpp>
#include <gemmi/dsn6.hpp>
#include <gemmi/isosurface.hpp>
#include <cmath>
#include <emscripten/val.h>

namespace {

em::val float_view(const std::vector<float>& v) {
  return em::val(emscripten::typed_memory_view(v.size(), v.data()));
}

em::val uint32_view(const std::vector<uint32_t>& v) {
  return em::val(emscripten::typed_memory_view(v.size(), v.data()));
}

// Shared map wrapper — holds a Grid<float>, DataStats, and isosurface results.
// Both Ccp4Map and Dsn6Map delegate to this for the common interface.
struct MapData {
  gemmi::Grid<float> grid;
  gemmi::DataStats stats = {};
  gemmi::IsoSurface iso;
  std::string last_error;

  em::val data() const { return float_view(grid.data); }

  bool extract_isosurface(double radius, double x, double y, double z,
                          double isolevel, const std::string& method) {
    try {
      iso = gemmi::extract_isosurface(grid, {x, y, z}, radius,
                                      isolevel, gemmi::iso_method_from_string(method));
    } catch (std::exception& e) {
      last_error = e.what();
      return false;
    }
    return true;
  }

  em::val isosurface_vertices() const { return float_view(iso.vertices); }
  em::val isosurface_segments() const { return uint32_view(iso.triangles); }

  int get_nx() const { return grid.nu; }
  int get_ny() const { return grid.nv; }
  int get_nz() const { return grid.nw; }
  double get_mean() const { return stats.dmean; }
  double get_rms() const { return stats.rms; }
  std::string get_last_error() const { return last_error; }
  gemmi::UnitCell get_cell() const { return grid.unit_cell; }
};

bool has_valid_header_stats(const gemmi::DataStats& stats) {
  return std::isfinite(stats.rms) && stats.rms > 0. &&
         stats.dmean >= stats.dmin && stats.dmean <= stats.dmax;
}

void apply_mode0_scaling(gemmi::Ccp4<float>& ccp4) {
  if (ccp4.header_i32(4) != 0)
    return;
  if (ccp4.header_i32(40) != -128 || ccp4.header_i32(41) != 127)
    return;
  const double min = ccp4.hstats.dmin;
  const double max = ccp4.hstats.dmax;
  const float b1 = static_cast<float>((max - min) / 255.0);
  const float b0 = static_cast<float>(0.5 * (min + max + b1));
  for (float& value : ccp4.grid.data)
    value = b1 * value + b0;
}

}  // namespace

class Ccp4Map : public MapData {
public:
  Ccp4Map(std::string buf): buf_(std::move(buf)) {}

  bool read(bool expand_symmetry=true) {
    try {
      ccp4_.read_ccp4_from_memory(buf_.data(), buf_.size(), "");
      ccp4_.setup(0.f, expand_symmetry ? gemmi::MapSetup::Full
                                       : gemmi::MapSetup::NoSymmetry);
      apply_mode0_scaling(ccp4_);
      if (!has_valid_header_stats(ccp4_.hstats))
        ccp4_.hstats = gemmi::calculate_data_statistics(ccp4_.grid.data);
      grid = std::move(ccp4_.grid);
      stats = ccp4_.hstats;
    } catch (std::runtime_error& e) {
      last_error = "Failed to read CCP4 map: ";
      last_error += e.what();
      return false;
    }
    return true;
  }

private:
  gemmi::Ccp4<float> ccp4_;
  std::string buf_;
};

class Dsn6Map : public MapData {
public:
  Dsn6Map(std::string buf): buf_(std::move(buf)) {}

  bool read() {
    try {
      stats = gemmi::read_dsn6_from_memory(buf_.data(), buf_.size(), grid);
    } catch (std::runtime_error& e) {
      last_error = "Failed to read DSN6 map: ";
      last_error += e.what();
      return false;
    }
    return true;
  }

private:
  std::string buf_;
};

void add_map() {
  em::class_<MapData>("MapData")
    .function("data", &MapData::data)
    .function("extract_isosurface", &MapData::extract_isosurface)
    .function("isosurface_vertices", &MapData::isosurface_vertices)
    .function("isosurface_segments", &MapData::isosurface_segments)
    .property("nx", &MapData::get_nx)
    .property("ny", &MapData::get_ny)
    .property("nz", &MapData::get_nz)
    .property("mean", &MapData::get_mean)
    .property("rms", &MapData::get_rms)
    .property("last_error", &MapData::get_last_error)
    .property("cell", &MapData::get_cell)
    ;

  em::class_<Ccp4Map, em::base<MapData>>("Ccp4Map")
    .constructor<std::string>()
    .function("read", &Ccp4Map::read)
    ;

  em::class_<Dsn6Map, em::base<MapData>>("Dsn6Map")
    .constructor<std::string>()
    .function("read", &Dsn6Map::read)
    ;
}
