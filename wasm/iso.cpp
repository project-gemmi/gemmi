// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/isosurface.hpp>
#include <emscripten/val.h>

class Isosurface {
public:
  void resize_input(int point_count) {
    points_.resize(static_cast<size_t>(point_count) * 3u);
    values_.resize(static_cast<size_t>(point_count));
  }

  void set_size(int size_x, int size_y, int size_z) {
    dims_ = {size_x, size_y, size_z};
  }

  em::val input_points() {
    return em::val(emscripten::typed_memory_view(points_.size(), points_.data()));
  }

  em::val input_values() {
    return em::val(emscripten::typed_memory_view(values_.size(), values_.data()));
  }

  bool calculate(double isolevel, const std::string& method) {
    try {
      iso_ = gemmi::calculate_isosurface(dims_, values_, points_, isolevel,
                                         gemmi::iso_method_from_string(method));
    } catch (std::exception& e) {
      last_error_ = e.what();
      return false;
    }
    return true;
  }

  em::val vertices() const {
    return em::val(emscripten::typed_memory_view(iso_.vertices.size(), iso_.vertices.data()));
  }

  em::val triangles() const {
    return em::val(emscripten::typed_memory_view(iso_.triangles.size(), iso_.triangles.data()));
  }

  std::string get_last_error() const { return last_error_; }

private:
  std::array<int, 3> dims_ = {0, 0, 0};
  std::vector<float> points_;
  std::vector<float> values_;
  gemmi::IsoSurface iso_;
  std::string last_error_;
};

void add_iso() {
  em::class_<Isosurface>("Isosurface")
    .constructor<>()
    .function("resize_input", &Isosurface::resize_input)
    .function("set_size", &Isosurface::set_size)
    .function("input_points", &Isosurface::input_points)
    .function("input_values", &Isosurface::input_values)
    .function("calculate", &Isosurface::calculate)
    .function("vertices", &Isosurface::vertices)
    .function("triangles", &Isosurface::triangles)
    .property("last_error", &Isosurface::get_last_error)
    ;
}
