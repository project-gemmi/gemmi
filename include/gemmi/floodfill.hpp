// Copyright 2020 Global Phasing Ltd.
//
// The flood fill (scanline fill) algorithm for Grid.
// Assumes periodic boundary conditions in the grid and 6-way connectivity.

#ifndef GEMMI_FLOODFILL_HPP_
#define GEMMI_FLOODFILL_HPP_

#include "grid.hpp"    // for Grid

namespace gemmi {

// Land is either 0 (when 1=sea, 0=land) or 1 (when 0=sea, 1=land)
template<typename T, int Land>
struct FloodFill {
  Grid<T>& mask;

  struct Line {
    int u, v, w, ulen;
    T* ptr;
  };

  struct Result {
    std::vector<Line> lines;

    size_t point_count() const {
      size_t n = 0;
      for (const Line& line : lines)
        n += line.ulen;
      return n;
    }
  };

  static constexpr T this_island() { return T(Land|2); }

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

  void set_volume_values(Result& r, T value) const {
    for (Line& line : r.lines)
      set_line_values(line, value);
  }

  // Find all connected points with value Land. Change them to this_island().
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

inline void mask_nodes_above_threshold(Grid<std::int8_t>& mask, const Grid<float>& grid,
                                       double threshold, bool negate=false) {
  mask.copy_metadata_from(grid);
  mask.data.resize(grid.data.size());
  size_t n = 0;
  for (float d : grid.data)
    mask.data[n++] = std::int8_t((negate ? -d : d) > threshold);
}

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
