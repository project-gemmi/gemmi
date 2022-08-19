// Copyright 2022 Global Phasing Ltd.
//
// AsuBrick and MaskedGrid that is used primarily as direct-space asu mask.

#ifndef GEMMI_ASUMASK_HPP_
#define GEMMI_ASUMASK_HPP_

#include "grid.hpp"

namespace gemmi {

struct AsuBrick {
  // The brick is 0<=x<=size[0]/24, 0<=y<=size[1]/24, 0<=z<=size[2]/24
  static constexpr int denom = 24;
  std::array<int, 3> size;
  std::array<bool, 3> incl;
  int volume;

  AsuBrick(int a, int b, int c)
    // For now we don't check which boundaries are included in the asymmetric unit,
    // we assume inequalities are not strict (e.g. 0<=x<=1/2) except '<1'.
    : size({a,b,c}), incl({a < denom, b < denom, c < denom}), volume(a*b*c) {}

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

  Fractional get_upper_limit() const {
    double inv_denom = 1.0 / denom;
    return Fractional(inv_denom * size[0] + (incl[0] ? 1e-9 : -1e-9),
                      inv_denom * size[1] + (incl[1] ? 1e-9 : -1e-9),
                      inv_denom * size[2] + (incl[2] ? 1e-9 : -1e-9));
  }

  // cf. Ccp4Base::get_extent()
  Box<Fractional> get_extent() const {
    Box<Fractional> box;
    box.minimum = Fractional(-1e-9, -1e-9, -1e-9);
    box.maximum = get_upper_limit();
    return box;
  }

  std::array<int,3> uvw_end(const GridMeta& meta) const {
    if (meta.axis_order != AxisOrder::XYZ)
      fail("grid is not fully setup");
    Fractional f = get_upper_limit();
    // upper limit is positive and never exact integer
    auto iceil = [](double x) { return int(x) + 1; };
    return {iceil(f.x * meta.nu), iceil(f.y * meta.nv), iceil(f.z * meta.nw)};
  }
};

// Returns asu brick upper bound. Lower bound is always (0,0,0).
// Currently we do not check if the boundaries are includedFor now bounds are assumed
// Brute force method that considers 8^3 sizes.
inline AsuBrick find_asu_brick(const SpaceGroup* sg) {
  if (sg == nullptr)
    fail("Missing space group");

  using Point = std::array<int, 3>;
  static_assert(AsuBrick::denom == 24, "");
  const int allowed_sizes[] = {3, 4, 6, 8, 12, 16, 18, 24};
  const GroupOps gops = sg->operations();
  const int n_ops = gops.order();

  Grid<std::int8_t> grid;
  grid.spacegroup = sg;
  // Testing with grid size 24 can't distiguish x<=1/8 and x<1/6,
  // but it happens to give the same results as grid size 48 for all
  // space groups tabulated in gemmi, so it's fine.
  // M=1 -> grid size 24; M=2 -> grid size 48
  constexpr int M = 1;
  grid.set_size(M*AsuBrick::denom, M*AsuBrick::denom, M*AsuBrick::denom);
  const std::vector<GridOp> ops = grid.get_scaled_ops_except_id();

  auto is_asu_brick = [&](const AsuBrick& brick, bool check_size) -> bool {
    // The most effective screening points for grid size 24.
    // These points were determined by doing calculations for all space groups
    // from gemmi::spacegroup_tables.
    static const Point size_checkpoints[] = {
      {7, 17, 7},    // eliminates 9866 out of 14726 wrong candidates
      {11, 1, 23},   // eliminates 2208 out of the rest
      {11, 10, 11},  // eliminates 1108
      {19, 1, 1},    // eliminates 665
      {3, 7, 19},    // eliminates 305
      {13, 9, 3},
      {23, 23, 23},
      {21, 10, 7},
      {11, 22, 1},
      {9, 15, 3},
      {5, 17, 23},
      {1, 5, 23},
      {5, 7, 17},
      {7, 5, 15},
      {20, 4, 5},
      {9, 23, 23},
      {9, 13, 13},  // eliminates the last wrong candidates
    };
    static const Point boundary_checkpoints[] = {
      {6, 18, 18},
      {12, 12, 12},
      {8, 16, 0},
      {0, 12, 3},
      {12, 3, 2},
      {0, 0, 12},
      {3, 12, 6},
      {16, 8, 9},
      {9, 21, 21},
      {8, 16, 6},
      {12, 12, 7},
      {20, 15, 0},
      {12, 0, 1},
      {16, 0, 0},
      {0, 6, 18},
      // ...
    };
    if (M == 1) {
      auto is_in = [&](const Point& p) {
        return p[0] < brick.size[0] + int(brick.incl[0])
            && p[1] < brick.size[1] + int(brick.incl[1])
            && p[2] < brick.size[2] + int(brick.incl[2]);
      };
      auto check_point = [&](const Point& point) {
        if (is_in(point))
          return true;
        for (const GridOp& op : ops) {
          Point t = op.apply(point[0], point[1], point[2]);
          grid.index_n_ref(t[0], t[1], t[2]);
          if (is_in(t))
            return true;
        }
        return false;
      };

      auto it = check_size ? std::begin(size_checkpoints) : std::begin(boundary_checkpoints);
      auto end = check_size ? std::end(size_checkpoints) : std::end(boundary_checkpoints);
      for (; it != end; ++it)
        if (!check_point(*it))
          return false;
    }

    // full check (it could be skipped for M==1 and check_size)
    grid.fill(0);
    int w_lim = M * brick.size[2] + int(brick.incl[2]);
    int v_lim = M * brick.size[1] + int(brick.incl[1]);
    int u_lim = M * brick.size[0] + int(brick.incl[0]);
    for (int w = 0; w < w_lim; ++w)
      for (int v = 0; v < v_lim; ++v)
        for (int u = 0; u < u_lim; ++u) {
          size_t idx = grid.index_q(u, v, w);
          if (grid.data[idx] == 0) {
            grid.data[idx] = 1;
            for (const GridOp& op : ops) {
              Point t = op.apply(u, v, w);
              size_t mate_idx = grid.index_n(t[0], t[1], t[2]);
              grid.data[mate_idx] = 1;
            }
          }
        }
#if 0
    // this code was used for determining checkpoints
    bool found = false;
    for (size_t n = 0; n != grid.data.size(); ++n)
      if (grid.data[n] == 0) {
        auto p = grid.index_to_point(n);
        printf("[debug1] %d %d %d is missing\n", p.u, p.v, p.w);
        found = true;
      }
    if (found)
      printf("[debug2] checkpoints failed\n");
#endif
    if (std::find(grid.data.begin(), grid.data.end(), 0) != grid.data.end())
      return false;

    return true;
  };

  std::vector<AsuBrick> possible_bricks;
  for (int a : allowed_sizes)
    for (int b : allowed_sizes)
      for (int c : allowed_sizes) {
        AsuBrick brick(a, b, c);
        if (brick.volume * n_ops >= brick.denom * brick.denom * brick.denom)
          possible_bricks.push_back(brick);
      }
  // the last item is the full unit cell, no need to check it
  possible_bricks.pop_back();
  // if two bricks have the same size, prefer the more cube-shaped one
  std::stable_sort(possible_bricks.begin(), possible_bricks.end(),
      [](const AsuBrick& a, const AsuBrick& b) {
        return a.volume < b.volume ||
               (a.volume == b.volume && a.size[0] + a.size[1] + a.size[2] <
                                        b.size[0] + b.size[1] + b.size[2]);
  });
  for (AsuBrick& brick : possible_bricks)
    if (is_asu_brick(brick, true)) {
      for (int i = 0; i < 3; ++i) {
        if (brick.incl[i] && brick.size[i] != 4) {
          brick.incl[i] = false;
          if (!is_asu_brick(brick, false))
            brick.incl[i] = true;
        }
      }
      return brick;
    }
  return AsuBrick(AsuBrick::denom, AsuBrick::denom, AsuBrick::denom);
}


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
};

template<typename V=std::int8_t>
std::vector<V> get_asu_mask(const GridMeta& grid) {
  std::vector<V> mask(grid.point_count(), 2);
  std::vector<GridOp> ops = grid.get_scaled_ops_except_id();
  auto end = find_asu_brick(grid.spacegroup).uvw_end(grid);
  for (int w = 0; w < end[2]; ++w)
    for (int v = 0; v < end[1]; ++v)
      for (int u = 0; u < end[0]; ++u) {
        size_t idx = grid.index_q(u, v, w);
        if (mask[idx] == 2) {
          mask[idx] = 0;
          for (const GridOp& op : ops) {
            std::array<int, 3> t = op.apply(u, v, w);
            size_t mate_idx = grid.index_n(t[0], t[1], t[2]);
            // grid point can be on special position
            if (mate_idx != idx)
              mask[mate_idx] = 1;
          }
        }
      }
  if (std::find(mask.begin(), mask.end(), 2) != mask.end())
    fail("get_asu_mask(): internal error");
  return mask;
}

template<typename T>
MaskedGrid<T> masked_asu(Grid<T>& grid) {
  return {get_asu_mask(grid), &grid};
}


// Calculating bounding box (brick) with the data (non-zero and non-NaN).

namespace impl {
// find the shortest span (possibly wrapped) that contains all true values
inline std::pair<int, int> trim_false_values(const std::vector<bool>& vec) {
  const int n = (int) vec.size();
  assert(n != 0);
  std::pair<int, int> span{n, n};  // return value for all-true vector
  int max_trim = 0;
  // first calculate the wrapped span (initial + final non-zero values)
  if (!vec[0] || !vec[n-1]) {
    // determine trailing-false length and store it in span.first
    while (span.first != 0 && !vec[span.first-1])
      --span.first;
    if (span.first == 0)  // all-false vector
      return span;  // i.e. {0,n}
    // determine leading-false length and store it in span.second
    span.second = 0;
    while (span.second != n && !vec[span.second])
      ++span.second;
    max_trim = span.second + (n - span.first);
  }
  for (int start = 0; ;) {
    for (;;) {
      if (start == n)
        return span;
      if (!vec[start])
        break;
      ++start;
    }
    int end = start + 1;
    while (end != n && !vec[end])
      ++end;
    if (end - start > max_trim) {
      max_trim = end - start;
      span.first = start;
      span.second = end;
    }
    start = end;
  }
  unreachable();
}
}  // namespace impl

// Get the smallest box with non-zero (and non-NaN) values.
template<typename T>
Box<Fractional> get_nonzero_extent(const GridBase<T>& grid) {
  grid.check_not_empty();
  std::vector<bool> nonzero[3];
  nonzero[0].resize(grid.nu, false);
  nonzero[1].resize(grid.nv, false);
  nonzero[2].resize(grid.nw, false);
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        T val = grid.data[idx];
        if (!(impl::is_nan(val) || val == 0)) {
          nonzero[0][u] = true;
          nonzero[1][v] = true;
          nonzero[2][w] = true;
        }
      }
  Box<Fractional> box;
  for (int i = 0; i < 3; ++i) {
    auto span = impl::trim_false_values(nonzero[i]);
    double inv_n = 1.0 / nonzero[i].size();
    box.minimum.at(i) = (span.second - 0.5) * inv_n - int(span.second >= span.first);
    box.maximum.at(i) = (span.first - 0.5) * inv_n;
  }
  return box;
}

} // namespace gemmi
#endif
