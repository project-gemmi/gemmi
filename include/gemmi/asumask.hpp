// Copyright 2022 Global Phasing Ltd.
//
// MaskedGrid used primarily as direct-space asu mask.

#ifndef GEMMI_ASUMASK_HPP_
#define GEMMI_ASUMASK_HPP_

#include "grid.hpp"

namespace gemmi {

template<typename T, typename V=std::int8_t> struct MaskedGrid {
  std::vector<V> mask;
  Grid<T>* grid;

  struct iterator {
    MaskedGrid& parent;
    size_t index;
    int u = 0, v = 0, w = 0;
    iterator(MaskedGrid& parent_, size_t index_)
      : parent(parent_), index(index_) {}
    iterator& operator++() {
      do {
        ++index;
        if (++u == parent.grid->nu) {
          u = 0;
          if (++v == parent.grid->nv) {
            v = 0;
            ++w;
          }
        }
      } while (index != parent.mask.size() && parent.mask[index] != 0);
      return *this;
    }
    typename GridBase<T>::Point operator*() {
      return {u, v, w, &parent.grid->data[index]};
    }
    bool operator==(const iterator &o) const { return index == o.index; }
    bool operator!=(const iterator &o) const { return index != o.index; }
  };
  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, mask.size()}; }

  std::array<int, 3> asu_max() const {
    std::array<int, 3> m = {{0, 0, 0}};
    for (auto pt : *const_cast<MaskedGrid*>(this))
      if (*pt.value == 0) {
        m[0] = std::max(pt.u, m[0]);
        m[1] = std::max(pt.v, m[1]);
        m[2] = std::max(pt.w, m[2]);
      }
    return m;
  }
};

template<typename T, typename V=std::int8_t>
std::vector<V> get_asu_mask(const Grid<T>& grid) {
  std::vector<V> mask(grid.data.size(), 0);
  std::vector<GridOp> ops = grid.get_scaled_ops_except_id();
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx)
        if (mask[idx] == 0)
          for (const GridOp& op : ops) {
            std::array<int, 3> t = op.apply(u, v, w);
            size_t mate_idx = grid.index_n(t[0], t[1], t[2]);
            // grid point can be on special position
            if (mate_idx != idx)
              mask[mate_idx] = 1;
          }
  return mask;
}

template<typename T>
MaskedGrid<T> masked_asu(Grid<T>& grid) {
  return {get_asu_mask(grid), &grid};
}


struct AsuBrick {
  // The brick is 0<=x<=size[0]/24, 0<=y<= size[1]/24, 0<=z<=size[2]/24
  static constexpr int denom = 24;
  std::array<int, 3> size;
  std::array<bool, 3> incl;
  int volume;

  AsuBrick(int a, int b, int c)
    // For now we don't check which boundaries are included in the asymmetric unit,
    // we assume inequalities are not strict (e.g. 0<=x<=1/2) except '<1'.
    : size({a,b,c}), incl({a < denom, b < denom, c < denom}), volume(a*b*c) {}

  bool is_in(const std::array<int,3>& p) const {
    return p[0] <= size[0] && p[1] <= size[1] && p[2] <= size[2];
  }

  std::string str() const {
    std::string s;
    for (int i = 0; i < 3; ++i) {
      if (i != 0)
        s += "; ";
      s += "0<=";
      s += "xyz"[i];
      s += incl[i] ? "<=" : "<";
      impl::append_fraction(s, impl::get_op_fraction(size[i]));
    }
    return s;
  }
};

// Returns asu brick upper bound. Lower bound is always (0,0,0).
// Currently we do not check if the boundaries are includedFor now bounds are assumed
// Brute force method that considers 8^3 sizes.
inline AsuBrick find_asu_brick(const SpaceGroup* sg) {
  if (sg == nullptr)
    fail("Missing space group");

  using Point = std::array<int, 3>;
  const int allowed_sizes[] = {3, 4, 6, 8, 12, 16, 18, 24};
  const GroupOps gops = sg->operations();
  const int n_ops = gops.order();
  constexpr int N = 24;

  Grid<std::int8_t> grid;
  grid.spacegroup = sg;
  grid.set_size(N, N, N);
  const std::vector<GridOp> ops = grid.get_scaled_ops_except_id();

  auto is_asu_brick = [&](const AsuBrick& brick) -> bool {
    // fast screening - corners and the middle
    // TODO: check which points are actually used
    static const Point checkpoints[] = {
      {N/2, N/2, N/2},
      {1, 1, N-2},
      {1, N-2, 1},
      {1, N-2, N-2},
      {N-2, 1, 1},
      {N-2, 1, N-2},
      {N-2, N-2, 1},
      {N-2, N-2, N-2},
    };
    //printf("[info] check brick %d x %d x %d\n",
    //       brick.size[0], brick.size[1], brick.size[2]);
    for (const Point& point : checkpoints)
      if (!brick.is_in(point)) {
        bool ok = false;
        for (const GridOp& op : ops) {
          Point t = op.apply(point[0], point[1], point[2]);
          grid.index_n_ref(t[0], t[1], t[2]);
          for (int i = 0; i < 3; ++i)
            assert(0 <= t[i] && t[i] < N);
          if (brick.is_in(t)) {
            ok = true;
            break;
          }
        }
        if (!ok) {
          //printf("[info] point %d %d %d\n", point[0], point[1], point[2]);
          return false;
        }
      }

    // full check
    //printf("[info] full check\n");
    grid.fill(0);
    for (int w = 0; w <= std::min(brick.size[2], 23); ++w)
      for (int v = 0; v <= std::min(brick.size[1], 23); ++v)
        for (int u = 0; u <= std::min(brick.size[0], 23); ++u) {
          int idx = grid.index_q(u, v, w);
          if (grid.data[idx] == 0) {
            grid.data[idx] = 1;
            for (const GridOp& op : ops) {
              Point t = op.apply(u, v, w);
              size_t mate_idx = grid.index_n(t[0], t[1], t[2]);
              grid.data[mate_idx] = 1;
            }
          }
        }
    return std::find(grid.data.begin(), grid.data.end(), 0) == grid.data.end();
  };

  std::vector<AsuBrick> possible_bricks;
  for (int a : allowed_sizes)
    for (int b : allowed_sizes)
      for (int c : allowed_sizes) {
        AsuBrick brick(a, b, c);
        if (brick.volume * n_ops >= N*N*N)
          possible_bricks.push_back(brick);
      }
  // the last item is the full unit cell, no need to check it
  possible_bricks.pop_back();
  // if two bricks have the same size, prefer the more cube-shaped one
  std::sort(possible_bricks.begin(), possible_bricks.end(),
            [](const AsuBrick& a, const AsuBrick& b) {
              return a.volume < b.volume || (a.volume == b.volume &&
                                             a.size[0] + a.size[1] + a.size[2] <
                                             b.size[0] + b.size[1] + b.size[2]);
  });
  for (const AsuBrick& brick : possible_bricks)
    if (is_asu_brick(brick))
      return brick;
  return AsuBrick(24, 24, 24);
}

} // namespace gemmi
#endif
