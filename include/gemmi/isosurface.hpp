// Copyright Global Phasing Ltd.
//
// Marching-cubes isosurface extraction from 3D grids.

#ifndef GEMMI_ISOSURFACE_HPP_
#define GEMMI_ISOSURFACE_HPP_

#include <array>
#include <cmath>      // for floor
#include <cstdint>    // for uint32_t
#include <string>
#include <vector>
#include "grid.hpp"   // for Grid, UnitCell, Position, Fractional

namespace gemmi {

/// Result of an isosurface extraction.
struct IsoSurface {
  std::vector<float> vertices;      // x, y, z triples
  std::vector<uint32_t> triangles;  // vertex-index triples
};

/// @brief Isosurface extraction algorithm selection.
enum class IsoMethod {
  /// @brief Standard marching cubes algorithm.
  MarchingCubes,
  /// @brief Marching cubes with snapped edge midpoints (0, 0.5, 1.0 only).
  SnappedMC
};

/// @brief Parse isosurface method name from string.
/// @param s String name: "snapped MC" returns SnappedMC, otherwise MarchingCubes.
/// @return The selected isosurface method.
inline IsoMethod iso_method_from_string(const std::string& s) {
  if (s == "snapped MC")
    return IsoMethod::SnappedMC;
  return IsoMethod::MarchingCubes;
}

namespace impl {

// Marching-cubes lookup tables.
// The tables below are defined in include/gemmi/mc_tables.hpp.
}  // namespace impl
}  // namespace gemmi

// The tables are large, keep them in a separate file.
#include "mc_tables.hpp"

namespace gemmi {

/// @brief Extract an isosurface from flat 3D arrays using marching cubes.
/// @details
/// Low-level function that operates directly on flat arrays of values and positions.
/// For convenience, use extract_isosurface() with a Grid object instead.
/// @param dims Number of grid points in each dimension [x, y, z].
/// @param values Scalar field values, length dims[0]*dims[1]*dims[2].
/// @param points Cartesian coordinates as flattened triples: x0,y0,z0,x1,y1,z1,...
///               Total length: 3 * dims[0]*dims[1]*dims[2].
/// @param isolevel The isovalue threshold for surface extraction.
/// @param method Marching cubes variant (standard or snapped).
/// @return IsoSurface containing vertices (x,y,z triples) and triangles (vertex-index triples).
inline IsoSurface calculate_isosurface(const std::array<int, 3>& dims,
                                       const std::vector<float>& values,
                                       const std::vector<float>& points,
                                       double isolevel,
                                       IsoMethod method = IsoMethod::MarchingCubes) {
  IsoSurface result;
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0)
    return result;
  const size_t point_count = static_cast<size_t>(dims[0]) * dims[1] * dims[2];
  if (values.size() != point_count || points.size() != point_count * 3u)
    fail("isosurface: array size mismatch");

  constexpr std::array<std::array<int, 3>, 8> kCubeVerts = {{
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
  }};

  const bool snap = (method == IsoMethod::SnappedMC);
  const auto& tri_offsets = impl::kTriTableOffsets;
  const int* tri_data = impl::kTriTableData.data();

  std::array<int, 8> vert_offsets;
  for (int i = 0; i < 8; ++i) {
    const auto& v = kCubeVerts[i];
    vert_offsets[i] = v[0] + dims[2] * (v[1] + dims[1] * v[2]);
  }

  std::array<float, 8> vertex_values;
  std::array<size_t, 8> point_offsets;
  std::array<uint32_t, 12> vlist = {};
  uint32_t vertex_count = 0;

  for (int x = 0; x < dims[0] - 1; ++x) {
    for (int y = 0; y < dims[1] - 1; ++y) {
      for (int z = 0; z < dims[2] - 1; ++z) {
        const int offset0 = z + dims[2] * (y + dims[1] * x);
        int cubeindex = 0;
        for (int i = 0; i < 8; ++i) {
          const int point_index = offset0 + vert_offsets[i];
          cubeindex |= (values[point_index] < isolevel) ? 1 << i : 0;
        }
        if (cubeindex == 0 || cubeindex == 255)
          continue;

        for (int i = 0; i < 8; ++i) {
          const size_t pi = static_cast<size_t>(offset0 + vert_offsets[i]);
          vertex_values[i] = values[pi];
          point_offsets[i] = pi * 3u;
        }

        const int edge_mask = impl::kEdgeTable[cubeindex];
        for (int i = 0; i < 12; ++i) {
          if ((edge_mask & (1 << i)) == 0)
            continue;
          const auto& edge = impl::kEdgeIndex[i];
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
          result.vertices.push_back(points[p1]     + (points[p2]     - points[p1])     * mu);
          result.vertices.push_back(points[p1 + 1] + (points[p2 + 1] - points[p1 + 1]) * mu);
          result.vertices.push_back(points[p1 + 2] + (points[p2 + 2] - points[p1 + 2]) * mu);
          vlist[i] = vertex_count++;
        }

        for (int i = tri_offsets[cubeindex]; i < tri_offsets[cubeindex + 1]; ++i)
          result.triangles.push_back(vlist[tri_data[i]]);
      }
    }
  }
  return result;
}

/// @brief Extract an isosurface from a Grid within a spherical region.
/// @details
/// Automatically extracts a region of the grid around a center point,
/// converts to fractional/grid coordinates, and runs marching cubes
/// on the subregion. This is the recommended high-level interface.
/// @tparam T The grid value type (float, double, etc.).
/// @param grid The 3D density map to process.
/// @param center Center point in Cartesian coordinates.
/// @param radius Radius of the sphere in Angstroms.
/// @param isolevel The isovalue threshold for surface extraction.
/// @param method Marching cubes variant (standard or snapped).
/// @return IsoSurface containing extracted vertices and triangle indices.
template<typename T>
IsoSurface extract_isosurface(const Grid<T>& grid,
                              const Position& center,
                              double radius,
                              double isolevel,
                              IsoMethod method = IsoMethod::MarchingCubes) {
  const UnitCell& cell = grid.unit_cell;
  const Fractional fc = cell.fractionalize(center);
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
  const std::array<int, 3> dims = {{
    grid_max[0] - grid_min[0] + 1,
    grid_max[1] - grid_min[1] + 1,
    grid_max[2] - grid_min[2] + 1,
  }};
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0)
    return {};

  const size_t point_count = static_cast<size_t>(dims[0]) * dims[1] * dims[2];
  std::vector<float> points;
  std::vector<float> values;
  points.reserve(point_count * 3u);
  values.reserve(point_count);

  for (int i = grid_min[0]; i <= grid_max[0]; ++i) {
    for (int j = grid_min[1]; j <= grid_max[1]; ++j) {
      for (int k = grid_min[2]; k <= grid_max[2]; ++k) {
        Position pos = grid.get_position(i, j, k);
        points.push_back(static_cast<float>(pos.x));
        points.push_back(static_cast<float>(pos.y));
        points.push_back(static_cast<float>(pos.z));
        values.push_back(grid.get_value(i, j, k));
      }
    }
  }

  return calculate_isosurface(dims, values, points, isolevel, method);
}

}  // namespace gemmi
#endif
