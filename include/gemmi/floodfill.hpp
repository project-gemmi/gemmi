// Copyright 2020 Global Phasing Ltd.
//
// The flood fill (scanline fill) algorithm for Grid.
// Assumes periodic boundary conditions in the grid and 6-way connectivity.

#ifndef GEMMI_FLOODFILL_HPP_
#define GEMMI_FLOODFILL_HPP_

#include "grid.hpp"    // for Grid

namespace gemmi {

template<typename T>
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

  void set_line_values(Line& line, int value) const {
    for (int i = 0; i < std::min(line.ulen, mask.nu - line.u); ++i) {
      assert(line.ptr[i] != value);
      line.ptr[i] = value;
    }
    for (int i = -line.u; i < line.ulen - mask.nu; ++i) {
      assert(line.ptr[i] != value);
      line.ptr[i] = value;
    }
  }

  void set_volume_values(Result& r, int value) const {
    for (Line& line : r.lines)
      set_line_values(line, value);
  }

  // Find all connected points with value 0. Change them to 2.
  Result find_volume(int u, int v, int w) {
    Result r;
    T* ptr = &mask.data[mask.index_q(u, v, w)];
    r.lines.push_back(line_from_point(u, v, w, ptr));
    set_line_values(r.lines.back(), 2);
    for (size_t i = 0; i < r.lines.size()/*increasing!*/; ++i) {
      int u_ = r.lines[i].u;
      int v_ = r.lines[i].v;
      int w_ = r.lines[i].w;
      int ulen = r.lines[i].ulen;
      // add adjacent lines
      add_lines(u_, v_ != 0 ? v_ - 1 : mask.nv - 1, w_, ulen, r);
      add_lines(u_, v_ + 1 != mask.nv ? v_ + 1 : 0, w_, ulen, r);
      add_lines(u_, v_, w_ != 0 ? w_ - 1 : mask.nw - 1, ulen, r);
      add_lines(u_, v_, w_ + 1 != mask.nw ? w_ + 1 : 0, ulen, r);
    }
    return r;
  }

private:
  void add_lines(int u, int v, int w, int ulen, Result& r) {
    T* ptr = &mask.data[mask.index_q(u, v, w)];
    for (int i = 0; i < std::min(ulen, mask.nu - u); ++i)
      if (ptr[i] == 0) {
        r.lines.push_back(line_from_point(u + i, v, w, ptr + i));
        set_line_values(r.lines.back(), 2);
      }
    for (int i = -u; i < ulen - mask.nu; ++i)
      if (ptr[i] == 0) {
        r.lines.push_back(line_from_point(u + i, v, w, ptr + i));
        set_line_values(r.lines.back(), 2);
      }
  }

  Line line_from_point(int u, int v, int w, T* ptr) const {
    int len = 1;
    while (u + len < mask.nu && ptr[len] == 0)
      ++len;
    if (u + len == mask.nu)
      while (len < mask.nu && ptr[len - mask.nu] == 0)
        ++len;
    for (int i = 0; i > -u; --i)
      if (ptr[i-1] != 0)
        return {u + i, v, w, len - i, ptr + i};
    if (ptr[mask.nu-1-u] != 0)
      return {0, v, w, len + u, ptr - u};
    for (int i = mask.nu - 1 - u; i > 1; --i)
      if (ptr[i-1] != 0)
        return {u + i, v, w, len + mask.nu - 2 - i, ptr + i};
    return {u, v, w, mask.nu, ptr};
  }
};


// Removes islands of 1's in the oceans of 0's. Uses flood fill.
// cf. find_blobs_by_flood_fill()
template<typename T>
int remove_islands_in_mask(Grid<T>& mask, size_t limit) {
  using Result = typename FloodFill<T>::Result;
  int counter = 0;
  size_t idx = 0;
  FloodFill<T> flood_fill{mask};
  for (int w = 0; w != mask.nw; ++w)
    for (int v = 0; v != mask.nv; ++v)
      for (int u = 0; u != mask.nu; ++u, ++idx) {
        assert(idx == mask.index_q(u, v, w));
        if (mask.data[idx] == 0) {
          // current island is temporarily marked as 2
          Result r = flood_fill.find_volume(u, v, w);
          //printf("island %d: %zu in %zu (limit: %zu)\n",
          //       counter, r.point_count(), lines.size(), limit);
          if (r.point_count() <= limit) {
            ++counter;
            flood_fill.set_volume_values(r, 1);
          }
        }
      }
  for (T& p : mask.data)
    p &= 1;
  return counter;
}

} // namespace gemmi
#endif
