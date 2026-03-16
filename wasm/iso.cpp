// Copyright Global Phasing Ltd.

#include "common.h"
#include "isosurface_tables.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace {

using Size3 = std::array<int, 3>;

constexpr std::array<Size3, 8> kCubeVerts = {{
  {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
}};

}  // namespace

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
    last_error_.clear();
    vertices_.clear();
    segments_.clear();
    if (dims_[0] <= 0 || dims_[1] <= 0 || dims_[2] <= 0) {
      last_error_ = "Grid dimensions are zero along at least one edge";
      return false;
    }
    const size_t point_count =
        static_cast<size_t>(dims_[0]) * dims_[1] * dims_[2];
    if (values_.size() != point_count || points_.size() != point_count * 3u) {
      last_error_ = "isosurface: array size mismatch";
      return false;
    }

    const bool snap = method == "snapped MC";
    const auto& seg_offsets = method == "squarish" ?
      gemmi_wasm::kSegTable2Offsets : gemmi_wasm::kSegTableOffsets;
    const int* seg_data = method == "squarish" ?
      gemmi_wasm::kSegTable2Data.data() : gemmi_wasm::kSegTableData.data();

    std::array<int, 8> vert_offsets;
    for (int i = 0; i < 8; ++i) {
      const auto& v = kCubeVerts[i];
      vert_offsets[i] = v[0] + dims_[2] * (v[1] + dims_[1] * v[2]);
    }

    std::array<float, 8> vertex_values;
    std::array<size_t, 8> point_offsets;
    std::array<uint32_t, 12> vlist = {};
    uint32_t vertex_count = 0;

    for (int x = 0; x < dims_[0] - 1; ++x) {
      for (int y = 0; y < dims_[1] - 1; ++y) {
        for (int z = 0; z < dims_[2] - 1; ++z) {
          const int offset0 = z + dims_[2] * (y + dims_[1] * x);
          int cubeindex = 0;
          for (int i = 0; i < 8; ++i) {
            const int point_index = offset0 + vert_offsets[i];
            cubeindex |= (values_[point_index] < isolevel) ? 1 << i : 0;
          }
          if (cubeindex == 0 || cubeindex == 255)
            continue;

          for (int i = 0; i < 8; ++i) {
            const size_t point_index = static_cast<size_t>(offset0 + vert_offsets[i]);
            vertex_values[i] = values_[point_index];
            point_offsets[i] = point_index * 3u;
          }

          const int edge_mask = gemmi_wasm::kEdgeTable[cubeindex];
          for (int i = 0; i < 12; ++i) {
            if ((edge_mask & (1 << i)) == 0)
              continue;
            const auto& edge = gemmi_wasm::kEdgeIndex[i];
            double mu = (isolevel - vertex_values[edge[0]]) /
                        (vertex_values[edge[1]] - vertex_values[edge[0]]);
            if (snap) {
              if (mu > 0.85)
                mu = 1.0;
              else if (mu < 0.15)
                mu = 0.0;
            }
            const size_t p1 = point_offsets[edge[0]];
            const size_t p2 = point_offsets[edge[1]];
            vertices_.push_back(points_[p1] + (points_[p2] - points_[p1]) * mu);
            vertices_.push_back(points_[p1 + 1] + (points_[p2 + 1] - points_[p1 + 1]) * mu);
            vertices_.push_back(points_[p1 + 2] + (points_[p2 + 2] - points_[p1 + 2]) * mu);
            vlist[i] = vertex_count++;
          }

          for (int i = seg_offsets[cubeindex]; i < seg_offsets[cubeindex + 1]; ++i)
            segments_.push_back(vlist[seg_data[i]]);
        }
      }
    }
    return true;
  }

  em::val vertices() const {
    return em::val(emscripten::typed_memory_view(vertices_.size(), vertices_.data()));
  }

  em::val segments() const {
    return em::val(emscripten::typed_memory_view(segments_.size(), segments_.data()));
  }

  std::string get_last_error() const { return last_error_; }

private:
  Size3 dims_ = {0, 0, 0};
  std::vector<float> points_;
  std::vector<float> values_;
  std::vector<float> vertices_;
  std::vector<uint32_t> segments_;
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
    .function("segments", &Isosurface::segments)
    .property("last_error", &Isosurface::get_last_error)
    ;
}
