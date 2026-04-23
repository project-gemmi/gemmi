/// @file floodfill.hpp
/// @brief Flood fill algorithm for identifying connected regions in 3D grids.
///
/// Implements a scanline-based flood fill for marking connected regions in a grid
/// with periodic boundary conditions and 6-way (face) connectivity.
/// Useful for identifying solvent envelopes, protein/ligand masks, or other
/// binary segmentation tasks on crystallographic maps.

// Copyright 2020 Global Phasing Ltd.
//
// The flood fill (scanline fill) algorithm for Grid.
// Assumes periodic boundary conditions in the grid and 6-way connectivity.

#ifndef GEMMI_FLOODFILL_HPP_
#define GEMMI_FLOODFILL_HPP_

#include <cstdint>     // for int8_t
#include "grid.hpp"    // for Grid

namespace gemmi {

/// @brief Flood fill algorithm for finding connected regions in a 3D grid.
/// @tparam T Grid data type (typically int8_t for binary masks)
/// @tparam Land Value marking the region to fill (0 or 1); fills connected regions with this value
///
/// Implements a scanline-based flood fill that respects periodic boundary conditions
/// and uses 6-way (face) connectivity. The algorithm identifies all grid points
/// with value @p Land that are connected to a seed point.
template<typename T, int Land>
struct FloodFill {
  /// @brief Reference to the grid being processed
  Grid<T>& mask;

  /// @brief Horizontal line segment within a flood fill region
  struct Line {
    int u;           ///< Starting u (x) coordinate
    int v;           ///< v (y) coordinate
    int w;           ///< w (z) coordinate
    int ulen;        ///< Length of the line in u direction
    T* ptr;          ///< Pointer to first data element of the line
  };

  /// @brief Result of a flood fill operation: all connected line segments
  struct Result {
    std::vector<Line> lines; ///< All horizontal lines composing the filled region

    /// @brief Total number of points in the filled region
    size_t point_count() const {
      size_t n = 0;
      for (const Line& line : lines)
        n += line.ulen;
      return n;
    }
  };

  /// @brief Internal marker value (combines Land with a flag bit)
  static constexpr T this_island() { return T(Land|2); }

  /// @brief Mark all points in a line with a given value, handling periodic wraparound.
  /// @param line Line segment to mark
  /// @param value Value to write
  ///
  /// Handles periodic boundary conditions in the u-direction for lines that wrap around.
  void set_line_values(Line& line, T value) const {
    for (int i = 0; i < std::min(line.ulen, mask.nu - line.u); ++i) {
      assert(line.ptr[i] != value);
      line.ptr[i] = value;
    }
    for (int i = -line.u; i < line.ulen - mask.nu; ++i) {
      assert(line.ptr[i] != value);
      line.ptr[i] = value;
    }
  }

  /// @brief Mark all lines in a filled region with a given value.
  /// @param r Result containing all lines of the region
  /// @param value Value to write to each line
  void set_volume_values(Result& r, T value) const {
    for (Line& line : r.lines)
      set_line_values(line, value);
  }

  /// @brief Find all grid points with value Land connected to a seed point via face neighbors.
  /// @param u Starting u (x) coordinate
  /// @param v Starting v (y) coordinate
  /// @param w Starting w (z) coordinate
  /// @return Result containing all connected line segments; points are temporarily marked with this_island()
  ///
  /// Uses a scanline algorithm that expands from the seed point to find all connected regions.
  /// Respects periodic boundary conditions and 6-way (face) connectivity.
  Result find_all_connected_points(int u, int v, int w) {
    Result r;
    T* ptr = &mask.data[mask.index_q(u, v, w)];
    r.lines.push_back(line_from_point(u, v, w, ptr));
    set_line_values(r.lines.back(), this_island());
    for (size_t i = 0; i < r.lines.size()/*increasing!*/; ++i) {
      int u_ = r.lines[i].u;
      int v_ = r.lines[i].v;
      int w_ = r.lines[i].w;
      int ustart = u_ != 0 ? u_ - 1 : mask.nu - 1;
      int ulen = std::min(mask.nu, r.lines[i].ulen + 2);
      // add adjacent lines
      int vl = v_ != 0 ? v_ - 1 : mask.nv - 1;
      int vh = v_ + 1 != mask.nv ? v_ + 1 : 0;
      int wl = w_ != 0 ? w_ - 1 : mask.nw - 1;
      int wh = w_ + 1 != mask.nw ? w_ + 1 : 0;
      add_lines(ustart, vl, wl, ulen, r);
      add_lines(ustart, vl, w_, ulen, r);
      add_lines(ustart, vl, wh, ulen, r);
      add_lines(ustart, v_, wl, ulen, r);
      add_lines(ustart, v_, wh, ulen, r);
      add_lines(ustart, vh, wl, ulen, r);
      add_lines(ustart, vh, w_, ulen, r);
      add_lines(ustart, vh, wh, ulen, r);
    }
    return r;
  }

  /// @brief Iterate over all connected regions (islands) in the grid, calling a function on each.
  /// @tparam Func Callable taking a Result (one island)
  /// @param func Function invoked once per connected region
  ///
  /// Scans the entire grid, identifying connected regions with value Land.
  /// After processing all regions, resets marked points back to Land.
  template<typename Func>
  void for_each_islands(Func func) {
    size_t idx = 0;
    for (int w = 0; w != mask.nw; ++w)
      for (int v = 0; v != mask.nv; ++v)
        for (int u = 0; u != mask.nu; ++u, ++idx) {
          assert(idx == mask.index_q(u, v, w));
          if (mask.data[idx] == Land) {
            // it temporarily marks current island as this_island()
            Result r = find_all_connected_points(u, v, w);
            func(r);
          }
        }
    // set big islands (continents) back to Land
    for (T& p : mask.data)
      p = T((int)p & 1);
  }

private:
  /// @brief Internal: expand from a line to find connected lines in neighboring v and w planes.
  /// @param u Starting u coordinate
  /// @param v Starting v coordinate
  /// @param w Starting w coordinate
  /// @param ulen Length to check
  /// @param r Result accumulating all connected lines
  void add_lines(int u, int v, int w, int ulen, Result& r) {
    T* ptr = &mask.data[mask.index_q(u, v, w)];
    for (int i = 0; i < std::min(ulen, mask.nu - u); ++i)
      if (ptr[i] == Land) {
        r.lines.push_back(line_from_point(u + i, v, w, ptr + i));
        set_line_values(r.lines.back(), this_island());
      }
    for (int i = -u; i < ulen - mask.nu; ++i)
      if (ptr[i] == Land) {
        r.lines.push_back(line_from_point(u + i, v, w, ptr + i));
        set_line_values(r.lines.back(), this_island());
      }
  }

  /// @brief Internal: extract the full horizontal line containing a point.
  /// @param u Starting u coordinate
  /// @param v Starting v coordinate
  /// @param w Starting w coordinate
  /// @param ptr Pointer to grid data at (u, v, w)
  /// @return Line segment with full extent in u direction, handling periodic wraparound
  ///
  /// Finds the maximal horizontal line of Land values containing the point,
  /// handling wraparound in the periodic u-direction.
  Line line_from_point(int u, int v, int w, T* ptr) const {
    int len = 1;
    while (u + len < mask.nu && ptr[len] == Land)
      ++len;
    if (u + len == mask.nu)
      while (len < mask.nu && ptr[len - mask.nu] == Land)
        ++len;
    for (int i = 0; i > -u; --i)
      if (ptr[i-1] != Land)
        return {u + i, v, w, len - i, ptr + i};
    if (ptr[mask.nu-1-u] != Land)
      return {0, v, w, len + u, ptr - u};
    for (int i = mask.nu - 1 - u; i > 1; --i)
      if (ptr[i-1] != Land)
        return {u + i, v, w, len + mask.nu - 2 - i, ptr + i};
    return {u, v, w, mask.nu, ptr};
  }
};

/// @brief Create a binary mask of grid points above (or below) a density threshold.
/// @param mask Output binary mask grid (1 = above threshold, 0 = below)
/// @param grid Input electron density or other continuous-valued grid
/// @param threshold Density threshold value
/// @param negate If true, invert the comparison (below threshold -> 1)
///
/// Copies grid metadata to mask and sets mask data to 1 where grid > threshold (or < threshold if negate).
inline void mask_nodes_above_threshold(Grid<std::int8_t>& mask, const Grid<float>& grid,
                                       double threshold, bool negate=false) {
  mask.copy_metadata_from(grid);
  mask.data.resize(grid.data.size());
  size_t n = 0;
  for (float d : grid.data)
    mask.data[n++] = std::int8_t((negate ? -d : d) > threshold);
}

/// @brief Identify connected regions above (or below) a threshold using flood fill starting from seeds.
/// @param grid Input electron density or other continuous-valued grid
/// @param seeds Vector of seed positions (Angstrom coordinates) to start the flood fill
/// @param threshold Density threshold; regions above this value are filled
/// @param negate If true, find regions below threshold instead
/// @return Binary mask grid (1 = in connected region from seeds, 0 = background)
///
/// Creates a mask where 1 indicates points in regions connected to the seed points
/// and satisfying the threshold criterion. Useful for identifying solvent envelopes
/// (seed from solvent, find connected regions above threshold) or protein masks
/// (seed from protein region, find connected regions).
inline Grid<std::int8_t> flood_fill_above(const Grid<float>& grid,
                                          const std::vector<Position>& seeds,
                                          double threshold,
                                          bool negate=false) {
  grid.check_not_empty();
  Grid<std::int8_t> mask;
  mask_nodes_above_threshold(mask, grid, threshold, negate);
  FloodFill<std::int8_t, 1> flood_fill{mask};
  for (const Position& pos : seeds) {
    auto point = mask.get_nearest_point(pos);
    if (*point.value == 1)
      flood_fill.find_all_connected_points(point.u, point.v, point.w);
  }
  for (std::int8_t& d : mask.data)
    d = std::int8_t(d == flood_fill.this_island());
  return mask;
}

} // namespace gemmi
#endif
