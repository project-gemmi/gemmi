// Copyright 2017 Global Phasing Ltd.
//
// 3d grids used by CCP4 maps, cell-method search and hkl data.

#ifndef GEMMI_GRID_HPP_
#define GEMMI_GRID_HPP_

#include <cassert>
#include <complex>
#include <algorithm>  // for fill
#include <numeric>    // for accumulate
#include <vector>
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "fail.hpp"  // for fail

namespace gemmi {

inline int modulo(int a, int n) {
  if (a >= n)
    a %= n;
  else if (a < 0)
    a = (a + 1) % n + n - 1;
  return a;
}

inline bool has_small_factorization(int n) {
  while (n % 2 == 0)
    n /= 2;
  for (int k : {3, 5})
    while (n % k == 0)
      n /= k;
  return n == 1 || n == -1;
}

inline std::array<int, 3> good_grid_size(const std::array<double, 3>& limit,
                                         bool denser, const SpaceGroup* sg) {
  std::array<int, 3> m = {{0, 0, 0}};
  GroupOps gops = (sg ? *sg : get_spacegroup_p1()).operations();
  std::array<int, 3> sg_fac = gops.find_grid_factors();
  for (int i = 0; i != 3; ++i) {
    for (int j = 0; j < i; ++j)
      if (std::fabs(limit[i] - limit[j]) < 0.5 && sg_fac[i] == sg_fac[j]) {
        m[i] = m[j];
        break;
      }
    if (m[i] == 0) {
      // having sizes always even simplifies things
      int f = sg_fac[i] % 2 == 0 ? sg_fac[i] : 2 * sg_fac[i];
      int n;
      if (denser) {
        n = int(std::ceil(limit[i] / f));
        while (!has_small_factorization(n))
          ++n;
      } else {
        n = int(std::floor(limit[i] / f));
        if (n > 1)
          while (!has_small_factorization(n))
            --n;
        else
          n = 1;
      }
      m[i] = n * f;
    }
  }
  for (int i = 1; i != 3; ++i)
    for (int j = 0; j != i; ++j)
      if (gops.are_directions_symmetry_related(i, j) && m[i] != m[j])
        m[i] = m[j] = (denser ? std::max(m[i], m[j]) : std::min(m[i], m[j]));
  return m;
}

struct GridOp {
  Op scaled_op;

  std::array<int, 3> apply(int u, int v, int w) const {
    std::array<int, 3> t;
    const Op::Rot& rot = scaled_op.rot;
    for (int i = 0; i != 3; ++i)
      t[i] = rot[i][0] * u + rot[i][1] * v + rot[i][2] * w + scaled_op.tran[i];
    return t;
  }
};

inline void check_grid_factors(const SpaceGroup* sg, int u, int v, int w) {
  if (sg) {
    auto factors = sg->operations().find_grid_factors();
    if (u % factors[0] != 0 || v % factors[1] != 0 || w % factors[2] != 0)
      fail("Grid not compatible with the space group " + sg->xhm());
  }
}

inline double lerp_(double a, double b, double t) {
  return a + (b - a) * t;
}
template<typename T>
std::complex<T> lerp_(std::complex<T> a, std::complex<T> b, double t) {
  return a + (b - a) * (T) t;
}

template<typename T, typename V=std::int8_t> struct MaskedGrid;

// Order of grid axis. Some Grid functionality works only with the XYZ order.
// The values XYZ and XYZ are used only when the grid covers whole unit cell.
enum class AxisOrder : unsigned char {
  Unknown,
  XYZ,  // default, corresponds to CCP4 map with axis order XYZ,
        // i.e. index X (H in reciprocal space) is fast and Z (or L) is slow
  ZYX   // fast Z (or L), may not be fully supported everywhere
};

template<typename T>
struct GridBase {
  struct Point {
    int u, v, w;
    T* value;
  };

  UnitCell unit_cell;
  const SpaceGroup* spacegroup = nullptr;
  std::vector<T> data;
  int nu = 0, nv = 0, nw = 0;
  AxisOrder axis_order = AxisOrder::Unknown;


  int point_count() const { return nu * nv * nw; }

  void set_size_without_checking(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
  }

  // Quick but unsafe. assumes (for efficiency) that 0 <= u < nu, etc.
  int index_q(int u, int v, int w) const { return w * nu * nv + v * nu + u; }
  T get_value_q(int u, int v, int w) const { return data[index_q(u, v, w)]; }

  Fractional point_to_fractional(const Point& p) const {
    return {p.u * (1.0 / nu), p.v * (1.0 / nv), p.w * (1.0 / nw)};
  }
  size_t point_to_index(const Point& p) const { return p.value - data.data(); }

  void fill(T value) { std::fill(data.begin(), data.end(), value); }
  T sum() const { return std::accumulate(data.begin(), data.end(), T()); }


  struct iterator {
    GridBase& parent;
    size_t index;
    int u = 0, v = 0, w = 0;
    iterator(GridBase& parent_, size_t index_)
      : parent(parent_), index(index_) {}
    iterator& operator++() {
      ++index;
      if (++u == parent.nu) {
        u = 0;
        if (++v == parent.nv) {
          v = 0;
          ++w;
        }
      }
      return *this;
    }
    typename GridBase<T>::Point operator*() {
      return {u, v, w, &parent.data[index]};
    }
    bool operator==(const iterator &o) const { return index == o.index; }
    bool operator!=(const iterator &o) const { return index != o.index; }
  };
  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, data.size()}; }
};

// For simplicity, some operations work only if the grid covers whole unit cell
// and axes u,v,w correspond to a,b,c in the unit cell.
template<typename T=float>
struct Grid : GridBase<T> {
  using Point = typename GridBase<T>::Point;
  using GridBase<T>::nu;
  using GridBase<T>::nv;
  using GridBase<T>::nw;
  using GridBase<T>::unit_cell;
  using GridBase<T>::spacegroup;
  using GridBase<T>::data;

  double spacing[3];

  void calculate_spacing() {
    spacing[0] = 1.0 / (nu * unit_cell.ar);
    spacing[1] = 1.0 / (nv * unit_cell.br);
    spacing[2] = 1.0 / (nw * unit_cell.cr);
  }

  void set_size_without_checking(int u, int v, int w) {
    GridBase<T>::set_size_without_checking(u, v, w);
    calculate_spacing();
    this->axis_order = AxisOrder::XYZ;
  }

  void set_size(int u, int v, int w) {
    check_grid_factors(spacegroup, u, v, w);
    set_size_without_checking(u, v, w);
  }

  // The resulting spacing can be smaller (if denser=true) or greater than arg.
  void set_size_from_spacing(double approx_spacing, bool denser) {
    std::array<double, 3> limit = {{1. / (unit_cell.ar * approx_spacing),
                                    1. / (unit_cell.br * approx_spacing),
                                    1. / (unit_cell.cr * approx_spacing)}};
    auto m = good_grid_size(limit, denser, spacegroup);
    set_size_without_checking(m[0], m[1], m[2]);
  }


  void set_unit_cell(double a, double b, double c,
                     double alpha, double beta, double gamma) {
    unit_cell.set(a, b, c, alpha, beta, gamma);
    calculate_spacing();
  }

  void set_unit_cell(const UnitCell& cell) {
    unit_cell = cell;
    calculate_spacing();
  }

  // Assumes (for efficiency) that -nu <= u < 2*nu, etc.
  int index_n(int u, int v, int w) const {
    if (u >= nu) u -= nu; else if (u < 0) u += nu;
    if (v >= nv) v -= nv; else if (v < 0) v += nv;
    if (w >= nw) w -= nw; else if (w < 0) w += nw;
    return this->index_q(u, v, w);
  }

  // Safe but slower.
  int index_s(int u, int v, int w) const {
    return this->index_q(modulo(u, nu), modulo(v, nv), modulo(w, nw));
  }

  T get_value(int u, int v, int w) const {
    return data[index_s(u, v, w)];
  }

  Point get_point(int u, int v, int w) {
    return {u, v, w, &data[index_s(u, v, w)]};
  }

  Position point_to_position(const Point& p) const {
    return unit_cell.orthogonalize(this->point_to_fractional(p));
  }

  // https://en.wikipedia.org/wiki/Trilinear_interpolation
  T interpolate_value(double x, double y, double z) const {
    double tmp;
    double xd = std::modf(x, &tmp);
    int u = (int) tmp;
    double yd = std::modf(y, &tmp);
    int v = (int) tmp;
    double zd = std::modf(z, &tmp);
    int w = (int) tmp;
    assert(u >= 0 && v >= 0 && w >= 0);
    assert(u < nu && v < nv && w < nw);
    T avg[2];
    for (int i = 0; i < 2; ++i) {
      int wi = (i == 0 || w + 1 != nw ? w + i : 0);
      int idx1 = this->index_q(u, v, wi);
      int v2 = v + 1 != nv ? v + 1 : 0;
      int idx2 = this->index_q(u, v2, wi);
      int u_add = u + 1 != nu ? 1 : -u;
      avg[i] = (T) lerp_(lerp_(data[idx1], data[idx1 + u_add], xd),
                         lerp_(data[idx2], data[idx2 + u_add], xd),
                         yd);
    }
    return (T) lerp_(avg[0], avg[1], zd);
  }
  T interpolate_value(const Fractional& fctr) const {
    Fractional f = fctr.wrap_to_unit();
    return interpolate_value(f.x * nu, f.y * nv, f.z * nw);
  }
  T interpolate_value(const Position& ctr) const {
    return interpolate_value(unit_cell.fractionalize(ctr));
  }

  void set_value(int u, int v, int w, T x) { data[index_s(u, v, w)] = x; }

  template <typename Func>
  void use_points_in_box(const Fractional& fctr_, int du, int dv, int dw,
                         Func&& func, bool fail_on_too_large_radius=true) {
    if (fail_on_too_large_radius) {
      if (2 * du >= nu || 2 * dv >= nv || 2 * dw >= nw)
        fail("grid operation failed: radius bigger than half the unit cell?");
    } else {
      // If we'd use the minimum image convention the max would be (nu-1)/2.
      // The limits set here are necessary for index_n() that is used below.
      du = std::min(du, nu - 1);
      dv = std::min(dv, nv - 1);
      dw = std::min(dw, nw - 1);
    }
    const Fractional fctr = fctr_.wrap_to_unit();
    int u0 = iround(fctr.x * nu);
    int v0 = iround(fctr.y * nv);
    int w0 = iround(fctr.z * nw);
    for (int w = w0-dw; w <= w0+dw; ++w)
      for (int v = v0-dv; v <= v0+dv; ++v)
        for (int u = u0-du; u <= u0+du; ++u) {
          Fractional fdelta{fctr.x - u * (1.0 / nu),
                            fctr.y - v * (1.0 / nv),
                            fctr.z - w * (1.0 / nw)};
          Position delta = unit_cell.orthogonalize_difference(fdelta);
          func(data[index_n(u, v, w)], delta);
        }
  }

  template <typename Func>
  void use_points_around(const Fractional& fctr_, double radius, Func&& func,
                         bool fail_on_too_large_radius=true) {
    int du = (int) std::ceil(radius / spacing[0]);
    int dv = (int) std::ceil(radius / spacing[1]);
    int dw = (int) std::ceil(radius / spacing[2]);
    use_points_in_box(fctr_, du, dv, dw,
                      [&](T& point, const Position& delta) {
                        double d2 = delta.length_sq();
                        if (d2 < radius * radius)
                          func(point, d2);
                      },
                      fail_on_too_large_radius);
  }

  void set_points_around(const Position& ctr, double radius, T value) {
    Fractional fctr = unit_cell.fractionalize(ctr);
    use_points_around(fctr, radius, [&](T& point, double) { point = value; });
  }

  // used in Fortran bindings
  void mask_atom(double x, double y, double z, double radius) {
    set_points_around(Position(x, y, z), radius, 1);
  }

  void make_zeros_and_ones(double threshold) {
    for (auto& d : data)
      d = d > threshold ? 1 : 0;
  }

  // operations re-scaled for faster later calculations; identity not included
  std::vector<GridOp> get_scaled_ops_except_id() const {
    GroupOps gops = spacegroup->operations();
    std::vector<GridOp> grid_ops;
    grid_ops.reserve(gops.order());
    for (const Op& so : gops.sym_ops)
      for (const Op::Tran& co : gops.cen_ops) {
        Op op = so.add_centering(co);
        if (op != Op::identity()) {
          // Rescale. Rotations are expected to be integral.
          op.tran[0] = op.tran[0] * nu / Op::DEN;
          op.tran[1] = op.tran[1] * nv / Op::DEN;
          op.tran[2] = op.tran[2] * nw / Op::DEN;
          for (int i = 0; i != 3; ++i)
            for (int j = 0; j != 3; ++j)
              op.rot[i][j] /= Op::DEN;
          grid_ops.push_back({op});
        }
      }
    return grid_ops;
  }

  template<typename Func>
  void symmetrize_using_ops(const std::vector<GridOp>& ops, Func func) {
    std::vector<int> mates(ops.size(), 0);
    std::vector<bool> visited(data.size(), false);
    int idx = 0;
    for (int w = 0; w != nw; ++w)
      for (int v = 0; v != nv; ++v)
        for (int u = 0; u != nu; ++u, ++idx) {
          assert(idx == this->index_q(u, v, w));
          if (visited[idx])
            continue;
          for (size_t k = 0; k < ops.size(); ++k) {
            std::array<int,3> t = ops[k].apply(u, v, w);
            mates[k] = index_n(t[0], t[1], t[2]);
          }
          T value = data[idx];
          for (int k : mates) {
            assert(!visited[k]);
            value = func(value, data[k]);
          }
          data[idx] = value;
          visited[idx] = true;
          for (int k : mates) {
            data[k] = value;
            visited[k] = true;
          }
        }
    assert(idx == (int) data.size());
  }

  // Use provided function to reduce values of all symmetry mates of each
  // grid point, then assign the result to all the points.
  template<typename Func>
  void symmetrize(Func func) {
    if (spacegroup && spacegroup->number != 1 &&
        this->axis_order == AxisOrder::XYZ)
      symmetrize_using_ops(get_scaled_ops_except_id(), func);
  }

  // two most common symmetrize functions
  void symmetrize_min() {
    symmetrize([](T a, T b) { return (a < b || !(b == b)) ? a : b; });
  }
  void symmetrize_max() {
    symmetrize([](T a, T b) { return (a > b || !(b == b)) ? a : b; });
  }

  template<typename V> std::vector<V> get_asu_mask() const {
    std::vector<V> mask(data.size(), 0);
    std::vector<GridOp> ops = get_scaled_ops_except_id();
    int idx = 0;
    for (int w = 0; w != nw; ++w)
      for (int v = 0; v != nv; ++v)
        for (int u = 0; u != nu; ++u, ++idx)
          if (mask[idx] == 0)
            for (const GridOp& op : ops) {
              std::array<int, 3> t = op.apply(u, v, w);
              int mate_idx = index_n(t[0], t[1], t[2]);
              // grid point can be on special position
              if (mate_idx != idx)
                mask[mate_idx] = 1;
            }
    return mask;
  }

  MaskedGrid<T> asu();
};


template<typename T, typename V> struct MaskedGrid {
  Grid<T>* grid;
  Grid<V> mask; // should we simply store the mask as vector?

  MaskedGrid(Grid<T>& grid_, std::vector<V>&& mask_data) : grid(&grid_) {
    mask.nu = grid_.nu;
    mask.nv = grid_.nv;
    mask.nw = grid_.nw;
    mask.spacegroup = grid_.spacegroup;
    mask.data = mask_data;
  }

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
      } while (index != parent.mask.data.size() &&
               parent.mask.data[index] != 0);
      return *this;
    }
    typename GridBase<T>::Point operator*() {
      return {u, v, w, &parent.grid->data[index]};
    }
    bool operator==(const iterator &o) const { return index == o.index; }
    bool operator!=(const iterator &o) const { return index != o.index; }
  };
  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, mask.data.size()}; }
};

template<typename T>
MaskedGrid<T> Grid<T>::asu() {
  return {*this, get_asu_mask<std::int8_t>()};
}

template<typename T>
struct ReciprocalGrid : GridBase<T> {
  bool half_l = false; // hkl grid that stores only l>=0
  bool has_index(int u, int v, int w) const {
    bool half_u = (half_l && this->axis_order == AxisOrder::ZYX);
    bool half_w = (half_l && this->axis_order != AxisOrder::ZYX);
    return std::abs(half_u ? u : 2 * u) < this->nu &&
           std::abs(2 * v) < this->nv &&
           std::abs(half_w ? w : 2 * w) < this->nw;
  }
  void check_index(int u, int v, int w) const {
    if (!has_index(u, v, w))
      throw std::out_of_range("ReciprocalGrid: index out of grid.");
  }
  // Similar to Grid::index_n(), but works only for -nu <= u < nu, etc.
  int index_n(int u, int v, int w) const {
    return this->index_q(u >= 0 ? u : u + this->nu,
                         v >= 0 ? v : v + this->nv,
                         w >= 0 ? w : w + this->nw);
  }
  int index_checked(int u, int v, int w) const {
    check_index(u, v, w);
    return index_n(u, v, w);
  }
  T get_value(int u, int v, int w) const {
    return this->data[index_checked(u, v, w)];
  }
  T get_value_or_zero(int u, int v, int w) const {
    return has_index(u, v, w) ? this->data[index_n(u, v, w)] : T{};
  }
  void set_value(int u, int v, int w, T x) {
    this->data[index_checked(u, v, w)] = x;
  }
  Miller to_hkl(const typename GridBase<T>::Point& point) const {
    Miller hkl{{point.u, point.v, point.w}};
    if (2 * point.u >= this->nu &&
        !(half_l && this->axis_order == AxisOrder::ZYX))
      hkl[0] -= this->nu;
    if (2 * point.v >= this->nv)
      hkl[1] -= this->nv;
    if (2 * point.w >= this->nw &&
        !(half_l && this->axis_order != AxisOrder::ZYX))
      hkl[2] -= this->nw;
    if (this->axis_order == AxisOrder::ZYX)
      std::swap(hkl[0], hkl[2]);
    return hkl;
  }
};

template<typename T> using FPhiGrid = ReciprocalGrid<std::complex<T>>;

} // namespace gemmi
#endif
