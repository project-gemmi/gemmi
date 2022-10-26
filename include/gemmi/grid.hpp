// Copyright 2017 Global Phasing Ltd.
//
// 3d grids used by CCP4 maps, cell-method search and hkl data.

#ifndef GEMMI_GRID_HPP_
#define GEMMI_GRID_HPP_

#include <cassert>
#include <cstddef>    // for ptrdiff_t
#include <complex>
#include <algorithm>  // for fill
#include <memory>     // for unique_ptr
#include <numeric>    // for accumulate
#include <type_traits>
#include <vector>
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "fail.hpp"  // for fail

namespace gemmi {

/// Order of grid axis. Some Grid functionality works only with the XYZ order.
/// The values XYZ and ZYX are used only when the grid covers whole unit cell.
enum class AxisOrder : unsigned char {
  Unknown,
  XYZ,  // default, corresponds to CCP4 map with axis order XYZ,
        // i.e. index X (H in reciprocal space) is fast and Z (or L) is slow
  ZYX   // fast Z (or L), may not be fully supported everywhere
};

enum class GridSizeRounding {
  Nearest,
  Up,
  Down
};

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

inline int round_with_small_factorization(double exact, GridSizeRounding rounding) {
  int n;
  if (rounding == GridSizeRounding::Up) {
    n = int(std::ceil(exact));
    while (!has_small_factorization(n))
      ++n;
  } else if (rounding == GridSizeRounding::Down) {
    n = std::max(int(std::floor(exact)), 1);
    while (!has_small_factorization(n))
      --n;
  } else {  // GridSizeRounding::Nearest
    n = int(std::round(exact));
    int sign = n > exact ? -1 : 1;
    int k = 1;
    while (n <= 0 || !has_small_factorization(n)) {
      // for sign=1 we want sequence n+1, n-1, n+2, n-2, n+3, ...
      n += sign * k;
      sign = -sign;
      ++k;
    }
  }
  return n;
}

inline std::array<int, 3> good_grid_size(const std::array<double, 3>& limit,
                                         GridSizeRounding rounding,
                                         const SpaceGroup* sg) {
  std::array<int, 3> m = {{0, 0, 0}};
  GroupOps gops;
  if (sg)
    gops = sg->operations();
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
      int n = round_with_small_factorization(limit[i] / f, rounding);
      m[i] = n * f;
    }
  }
  for (int i = 1; i != 3; ++i)
    for (int j = 0; j != i; ++j)
      if (gops.are_directions_symmetry_related(i, j) && m[i] != m[j]) {
        bool denser = rounding != GridSizeRounding::Down;
        m[i] = m[j] = (denser ? std::max(m[i], m[j]) : std::min(m[i], m[j]));
      }
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

inline void check_grid_factors(const SpaceGroup* sg, std::array<int,3> size) {
  if (sg) {
    GroupOps gops = sg->operations();
    auto factors = gops.find_grid_factors();
    for (int i = 0; i != 3; ++i)
      if (size[i] % factors[i] != 0)
        fail("Grid not compatible with the space group " + sg->xhm());
    for (int i = 1; i != 3; ++i)
      for (int j = 0; j != i; ++j)
        if (gops.are_directions_symmetry_related(i, j) && size[i] != size[j])
          fail("Grid must have the same size in symmetry-related directions");
  }
}

inline double lerp_(double a, double b, double t) {
  return a + (b - a) * t;
}
template<typename T>
std::complex<T> lerp_(std::complex<T> a, std::complex<T> b, double t) {
  return a + (b - a) * (T) t;
}

/// Catmull–Rom spline interpolation. CINT_u from:
/// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
/// The same as (24) in https://journals.iucr.org/d/issues/2018/06/00/ic5103/
inline double cubic_interpolation(double u, double a, double b, double c, double d) {
  //return 0.5 * u * (u * (u * (3*b - 3*c + d - a) + (2*a - 5*b + 4*c - d)) + (c - a)) + b;
  // equivalent form that is faster on my computer:
  return -0.5 * (a * u * ((u-2)*u + 1) - b * ((3*u - 5) * u*u + 2) +
                 u * (c * ((3*u - 4) * u - 1) - d * (u-1) * u));
}

/// df/du (from Wolfram Alpha)
inline double cubic_interpolation_der(double u, double a, double b, double c, double d) {
  return a * (-1.5*u*u + 2*u - 0.5) + c * (-4.5*u*u + 4*u + 0.5)
         + u * (4.5*b*u - 5*b + 1.5*d*u - d);
}


/// The base of Grid classes that does not depend on stored data type.
struct GridMeta {
  UnitCell unit_cell;
  const SpaceGroup* spacegroup = nullptr;
  int nu = 0, nv = 0, nw = 0;
  AxisOrder axis_order = AxisOrder::Unknown;

  size_t point_count() const { return (size_t)nu * nv * nw; }
  /// u,v,w are not normalized here
  Fractional get_fractional(int u, int v, int w) const {
    return {u * (1.0 / nu), v * (1.0 / nv), w * (1.0 / nw)};
  }
  Position get_position(int u, int v, int w) const {
    return unit_cell.orthogonalize(get_fractional(u, v, w));
  }

  // operations re-scaled for faster later calculations; identity not included
  std::vector<GridOp> get_scaled_ops_except_id() const {
    std::vector<GridOp> grid_ops;
    if (!spacegroup || spacegroup->number == 1)
      return grid_ops;
    if (axis_order != AxisOrder::XYZ)
      fail("grid can use symmetries only if it is setup in the XYZ order");
    GroupOps gops = spacegroup->operations();
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

  /// Quick(est) index function, but works only if `0 <= u < nu`, etc.
  size_t index_q(int u, int v, int w) const {
    return size_t(w * nv + v) * nu + u;
  }

  /// Faster than index_s(), but works only if `-nu <= u < 2*nu`, etc.
  size_t index_n(int u, int v, int w) const { return index_n_ref(u, v, w); }

  /// The same as index_n(), but modifies arguments.
  size_t index_n_ref(int& u, int& v, int& w) const {
    if (u >= nu) u -= nu; else if (u < 0) u += nu;
    if (v >= nv) v -= nv; else if (v < 0) v += nv;
    if (w >= nw) w -= nw; else if (w < 0) w += nw;
    return this->index_q(u, v, w);
  }

  /// Faster than index_n(), but works only if -nu <= u < nu, etc.
  size_t index_near_zero(int u, int v, int w) const {
    return this->index_q(u >= 0 ? u : u + nu,
                         v >= 0 ? v : v + nv,
                         w >= 0 ? w : w + nw);
  }
};

/// A common subset of Grid and ReciprocalGrid.
template<typename T>
struct GridBase : GridMeta {
  /// grid coordinates (modulo size) and a pointer to value
  struct Point {
    int u, v, w;
    T* value;
  };

  std::vector<T> data;

  void check_not_empty() const {
    if (data.empty())
      fail("grid is empty");
  }

  void set_size_without_checking(int nu_, int nv_, int nw_) {
    nu = nu_, nv = nv_, nw = nw_;
    data.resize((size_t)nu_ * nv_ * nw_);
  }

  T get_value_q(int u, int v, int w) const { return data[index_q(u, v, w)]; }

  size_t point_to_index(const Point& p) const { return p.value - data.data(); }

  Point index_to_point(size_t idx) {
    auto d1 = std::div((ptrdiff_t)idx, (ptrdiff_t)nu);
    auto d2 = std::div(d1.quot, (ptrdiff_t)nv);
    int u = (int) d1.rem;
    int v = (int) d2.rem;
    int w = (int) d2.quot;
    assert(index_q(u, v, w) == idx);
    return {u, v, w, &data.at(idx)};
  }

  void fill(T value) {
    data.resize(point_count());
    std::fill(data.begin(), data.end(), value);
  }

  using Tsum = typename std::conditional<std::is_integral<T>::value,
                                         std::ptrdiff_t, T>::type;
  Tsum sum() const { return std::accumulate(data.begin(), data.end(), Tsum()); }


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

/// Real-space grid.
/// For simplicity, some operations work only if the grid covers whole unit cell
/// and axes u,v,w correspond to a,b,c in the unit cell.
template<typename T=float>
struct Grid : GridBase<T> {
  using Point = typename GridBase<T>::Point;
  using GridBase<T>::nu;
  using GridBase<T>::nv;
  using GridBase<T>::nw;
  using GridBase<T>::unit_cell;
  using GridBase<T>::spacegroup;
  using GridBase<T>::data;

  /// spacing between virtual planes, not between points
  double spacing[3] = {0., 0., 0.};

  /// copy unit_cell, spacegroup, nu, nv, nw, axis_order and set spacing
  void copy_metadata_from(const GridMeta& g) {
    unit_cell = g.unit_cell;
    spacegroup = g.spacegroup;
    nu = g.nu;
    nv = g.nv;
    nw = g.nw;
    this->axis_order = g.axis_order;
    calculate_spacing();
  }

  /// set #spacing
  void calculate_spacing() {
    spacing[0] = 1.0 / (nu * unit_cell.ar);
    spacing[1] = 1.0 / (nv * unit_cell.br);
    spacing[2] = 1.0 / (nw * unit_cell.cr);
  }

  double min_spacing() const {
    return std::min(std::min(spacing[0], spacing[1]), spacing[2]);
  }

  void set_size_without_checking(int nu_, int nv_, int nw_) {
    GridBase<T>::set_size_without_checking(nu_, nv_, nw_);
    calculate_spacing();
    this->axis_order = AxisOrder::XYZ;
  }

  void set_size(int nu_, int nv_, int nw_) {
    check_grid_factors(spacegroup, {{nu_, nv_, nw_}});
    set_size_without_checking(nu_, nv_, nw_);
  }

  // The resulting spacing can be smaller (if denser=true) or greater than arg.
  void set_size_from_spacing(double approx_spacing, GridSizeRounding rounding) {
    std::array<double, 3> limit = {{1. / (unit_cell.ar * approx_spacing),
                                    1. / (unit_cell.br * approx_spacing),
                                    1. / (unit_cell.cr * approx_spacing)}};
    auto m = good_grid_size(limit, rounding, spacegroup);
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

  template<typename S>
  void setup_from(const S& st, double approx_spacing) {
    spacegroup = st.find_spacegroup();
    set_unit_cell(st.cell);
    set_size_from_spacing(approx_spacing, GridSizeRounding::Up);
  }

  /// Returns index in data array for (u,v,w). Safe but slower than index_q().
  size_t index_s(int u, int v, int w) const {
    this->check_not_empty();
    return this->index_q(modulo(u, nu), modulo(v, nv), modulo(w, nw));
  }

  /// returns `data[index_s(u, v, w)]`
  T get_value(int u, int v, int w) const {
    return data[index_s(u, v, w)];
  }

  void set_value(int u, int v, int w, T x) {
    data[index_s(u, v, w)] = x;
  }

  /// Point stores normalizes indices (not the original u,v,w).
  Point get_point(int u, int v, int w) {
    u = modulo(u, nu);
    v = modulo(v, nv);
    w = modulo(w, nw);
    return {u, v, w, &data[this->index_q(u, v, w)]};
  }

  Point get_nearest_point(const Fractional& f) {
    if (this->axis_order != AxisOrder::XYZ)
      fail("grid is not fully setup");
    return get_point(iround(f.x * nu), iround(f.y * nv), iround(f.z * nw));
  }
  Point get_nearest_point(const Position& pos) {
    return get_nearest_point(unit_cell.fractionalize(pos));
  }

  size_t get_nearest_index(const Fractional& f) {
    return index_s(iround(f.x * nu), iround(f.y * nv), iround(f.z * nw));
  }

  /// Point stores normalized indices, so fractional coordinates are in [0,1).
  Fractional point_to_fractional(const Point& p) const {
    return this->get_fractional(p.u, p.v, p.w);
  }
  Position point_to_position(const Point& p) const {
    return this->get_position(p.u, p.v, p.w);
  }

  static double grid_modulo(double x, int n, int* iptr) {
    double f = std::floor(x);
    *iptr = modulo((int)f, n);
    return x - f;
  }

  /// https://en.wikipedia.org/wiki/Trilinear_interpolation
  /// x,y,z are grid coordinates (x=1.5 is between 2nd and 3rd grid point).
  T interpolate_value(double x, double y, double z) const {
    this->check_not_empty();
    int u, v, w;
    double xd = grid_modulo(x, nu, &u);
    double yd = grid_modulo(y, nv, &v);
    double zd = grid_modulo(z, nw, &w);
    assert(u >= 0 && v >= 0 && w >= 0);
    assert(u < nu && v < nv && w < nw);
    T avg[2];
    for (int i = 0; i < 2; ++i) {
      int wi = (i == 0 || w + 1 != nw ? w + i : 0);
      size_t idx1 = this->index_q(u, v, wi);
      int v2 = v + 1 != nv ? v + 1 : 0;
      size_t idx2 = this->index_q(u, v2, wi);
      int u_add = u + 1 != nu ? 1 : -u;
      avg[i] = (T) lerp_(lerp_(data[idx1], data[idx1 + u_add], xd),
                         lerp_(data[idx2], data[idx2 + u_add], xd),
                         yd);
    }
    return (T) lerp_(avg[0], avg[1], zd);
  }
  T interpolate_value(const Fractional& fctr) const {
    return interpolate_value(fctr.x * nu, fctr.y * nv, fctr.z * nw);
  }
  T interpolate_value(const Position& ctr) const {
    return interpolate_value(unit_cell.fractionalize(ctr));
  }

  /// https://en.wikipedia.org/wiki/Tricubic_interpolation
  /// x,y,z are grid coordinates (x=1.5 is between 2nd and 3rd grid point).
  double tricubic_interpolation(double x, double y, double z) const {
    std::array<std::array<std::array<T,4>,4>,4> copy;
    copy_4x4x4(x, y, z, copy);
    auto s = [&copy](int i, int j, int k) { return copy[i][j][k]; };
    double a[4], b[4];
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j)
        a[j] = cubic_interpolation(z, s(i,j,0), s(i,j,1), s(i,j,2), s(i,j,3));
      b[i] = cubic_interpolation(y, a[0], a[1], a[2], a[3]);
    }
    return cubic_interpolation(x, b[0], b[1], b[2], b[3]);
  }
  double tricubic_interpolation(const Fractional& fctr) const {
    return tricubic_interpolation(fctr.x * nu, fctr.y * nv, fctr.z * nw);
  }
  double tricubic_interpolation(const Position& ctr) const {
    return tricubic_interpolation(unit_cell.fractionalize(ctr));
  }
  /// returns the same as above + derivatives df/dx, df/dy, df/dz
  std::array<double,4> tricubic_interpolation_der(double x, double y, double z) const {
    std::array<std::array<std::array<T,4>,4>,4> copy;
    copy_4x4x4(x, y, z, copy);
    auto s = [&copy](int i, int j, int k) { return copy[i][j][k]; };
    double a[4][4];
    double b[4];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        a[i][j] = cubic_interpolation(z, s(i,j,0), s(i,j,1), s(i,j,2), s(i,j,3));
    for (int i = 0; i < 4; ++i)
      b[i] = cubic_interpolation(y, a[i][0], a[i][1], a[i][2], a[i][3]);
    std::array<double, 4> ret;
    ret[0] = cubic_interpolation(x, b[0], b[1], b[2], b[3]);
    ret[1] = cubic_interpolation_der(x, b[0], b[1], b[2], b[3]);
    for (int i = 0; i < 4; ++i)
      b[i] = cubic_interpolation(x, a[0][i], a[1][i], a[2][i], a[3][i]);
    ret[2] = cubic_interpolation_der(y, b[0], b[1], b[2], b[3]);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        a[i][j] = cubic_interpolation(y, s(i,0,j), s(i,1,j), s(i,2,j), s(i,3,j));
    for (int i = 0; i < 4; ++i)
      b[i] = cubic_interpolation(x, a[0][i], a[1][i], a[2][i], a[3][i]);
    ret[3] = cubic_interpolation_der(z, b[0], b[1], b[2], b[3]);
    return ret;
  }
  std::array<double,4> tricubic_interpolation_der(const Fractional& fctr) const {
    auto r = tricubic_interpolation_der(fctr.x * nu, fctr.y * nv, fctr.z * nw);
    return {r[0], r[1] * nu, r[2] * nv, r[3] * nw};
  }
  /// @private
  void copy_4x4x4(double& x, double& y, double& z,
                  std::array<std::array<std::array<T,4>,4>,4>& copy) const {
    this->check_not_empty();
    auto prepare_indices = [this](double& r, int nt, int (&indices)[4]) {
      int t;
      r = this->grid_modulo(r, nt, &t);
      indices[0] = (t != 0 ? t : nt) - 1;
      indices[1] = t;
      if (t + 2 < nt) {
        indices[2] = t + 1;
        indices[3] = t + 2;
      } else {
        indices[2] = t + 2 == nt ? t + 1 : 0;
        indices[3] = t + 2 == nt ? 0 : 1;
      }
    };
    int u_indices[4], v_indices[4], w_indices[4];
    prepare_indices(x, nu, u_indices);
    prepare_indices(y, nv, v_indices);
    prepare_indices(z, nw, w_indices);
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k)
          copy[i][j][k] = this->get_value_q(u_indices[i], v_indices[j], w_indices[k]);
  }

  /// @param order 1=nearest, 2=linear, 3=cubic interpolation
  T interpolate(const Fractional& f, int order) const {
    switch (order) {
      case 1: return *const_cast<Grid<T>*>(this)->get_nearest_point(f).value;
      case 2: return interpolate_value(f);
      case 3: return (T) tricubic_interpolation(f);
    }
    throw std::invalid_argument("interpolation \"order\" must 1, 2 or 3");
  }

  void get_subarray(T* dest, std::array<int,3> start, std::array<int,3> shape) const {
    this->check_not_empty();
    if (this->axis_order != AxisOrder::XYZ)
      fail("get_subarray() is for Grids in XYZ order");
    const int u_start0 = modulo(start[0], nu);
    for (int w = 0; w < shape[2]; w++) {
      const int w0 = modulo(start[2] + w, nw);
      for (int v = 0; v < shape[1]; v++) {
        const int v0 = modulo(start[1] + v, nv);
        int u_start = u_start0;
        const T* src0 = &data[this->index_q(u_start, v0, w0)];
        int len = shape[0];
        while (len > nu - u_start) {
          int elem = nu - u_start;
          std::copy(src0, src0 + elem, dest);
          src0 -= u_start;
          dest += elem;
          len -= elem;
          u_start = 0;
        }
        std::copy(src0, src0 + len, dest);
        dest += len;
      }
    }
  }

  void set_subarray(const T* src, std::array<int,3> start, std::array<int,3> shape) {
    this->check_not_empty();
    if (this->axis_order != AxisOrder::XYZ)
      fail("set_subarray() is for Grids in XYZ order");
    const int u_start0 = modulo(start[0], nu);
    for (int w = 0; w < shape[2]; w++) {
      const int w0 = modulo(start[2] + w, nw);
      for (int v = 0; v < shape[1]; v++) {
        const int v0 = modulo(start[1] + v, nv);
        int u_start = u_start0;
        T* dst0 = &data[this->index_q(u_start, v0, w0)];
        int len = shape[0];
        while (len > nu - u_start) {
          int elem = nu - u_start;
          std::copy(src, src + elem, dst0);
          dst0 -= u_start;
          src += elem;
          len -= elem;
          u_start = 0;
        }
        std::copy(src, src + len, dst0);
        src += len;
      }
    }
  }

  template <bool UsePbc>
  void check_size_for_points_in_box(int& du, int& dv, int& dw,
                                    bool fail_on_too_large_radius) const {
    if (fail_on_too_large_radius) {
      if (2 * du >= nu || 2 * dv >= nv || 2 * dw >= nw)
        fail("grid operation failed: radius bigger than half the unit cell?");
    }
    if (UsePbc && !fail_on_too_large_radius) {
      // If we'd use the minimum image convention the max would be (nu-1)/2.
      // The limits set here are necessary for index_n() that is used below.
      du = std::min(du, nu - 1);
      dv = std::min(dv, nv - 1);
      dw = std::min(dw, nw - 1);
    }
  }

  template <bool UsePbc, typename Func>
  void do_use_points_in_box(Fractional fctr, int du, int dv, int dw, Func&& func) {
    int u0 = iround(fctr.x * nu);
    int v0 = iround(fctr.y * nv);
    int w0 = iround(fctr.z * nw);
    int u_lo = u0 - du;
    int u_hi = u0 + du;
    int v_lo = v0 - dv;
    int v_hi = v0 + dv;
    int w_lo = w0 - dw;
    int w_hi = w0 + dw;
    if (!UsePbc) {
      u_lo = std::max(u_lo, 0);
      u_hi = std::min(u_hi, nu - 1);
      v_lo = std::max(v_lo, 0);
      v_hi = std::min(v_hi, nv - 1);
      w_lo = std::max(w_lo, 0);
      w_hi = std::min(w_hi, nw - 1);
    }
    const Position orth0(unit_cell.orth.mat.column_copy(0));
    for (int w = w_lo; w <= w_hi; ++w) {
      int w_ = UsePbc ? modulo(w, nw) : w;
      double fw = w * (1.0 / nw);
      for (int v = v_lo; v <= v_hi; ++v) {
        int v_ = UsePbc ? modulo(v, nv) : v;
        double fv = v * (1.0 / nv);
        size_t idx0 = this->index_q(0, v_, w_);
        Position delta0 = unit_cell.orthogonalize_difference(fctr - Fractional(0., fv, fw));
        for (int u = u_lo; u <= u_hi; ++u) {
          int u_ = UsePbc ? modulo(u, nu) : u;
          double fu = u * (1.0 / nu);
          Position delta = delta0 - orth0 * fu;
          func(data[idx0 + u_], delta, u, v, w);
        }
      }
    }
  }

  template <bool UsePbc, typename Func>
  void use_points_in_box(Fractional fctr, int du, int dv, int dw,
                         Func&& func, bool fail_on_too_large_radius=true) {
    check_size_for_points_in_box<UsePbc>(du, dv, dw, fail_on_too_large_radius);
    do_use_points_in_box<UsePbc>(fctr, du, dv, dw, func);
  }

  template <bool UsePbc, typename Func>
  void use_points_around(const Fractional& fctr_, double radius, Func&& func,
                         bool fail_on_too_large_radius=true) {
    int du = (int) std::ceil(radius / spacing[0]);
    int dv = (int) std::ceil(radius / spacing[1]);
    int dw = (int) std::ceil(radius / spacing[2]);
    use_points_in_box<UsePbc>(fctr_, du, dv, dw,
                      [&](T& ref, const Position& delta, int, int, int) {
                        double d2 = delta.length_sq();
                        if (d2 < radius * radius)
                          func(ref, d2);
                      },
                      fail_on_too_large_radius);
  }

  void set_points_around(const Position& ctr, double radius, T value, bool use_pbc=true) {
    Fractional fctr = unit_cell.fractionalize(ctr);
    if (use_pbc)
      use_points_around<true>(fctr, radius, [&](T& ref, double) { ref = value; });
    else
      use_points_around<false>(fctr, radius, [&](T& ref, double) { ref = value; });
  }

  void change_values(T old_value, T new_value) {
    for (auto& d : data)
      if (impl::is_same(d, old_value))
        d = new_value;
  }

  /// Use \par func to reduce values of all symmetry mates of each
  /// grid point, then assign the result to all the points.
  /// \par func takes two values and returns a value.
  template<typename Func>
  void symmetrize(Func func) {
    symmetrize_using_ops(this->get_scaled_ops_except_id(), func);
  }

  template<typename Func>
  void symmetrize_using_ops(const std::vector<GridOp>& ops, Func func) {
    if (ops.empty())
      return;
    std::vector<size_t> mates(ops.size(), 0);
    std::vector<bool> visited(data.size(), false);
    size_t idx = 0;
    for (int w = 0; w != nw; ++w)
      for (int v = 0; v != nv; ++v)
        for (int u = 0; u != nu; ++u, ++idx) {
          assert(idx == this->index_q(u, v, w));
          if (visited[idx])
            continue;
          for (size_t k = 0; k < ops.size(); ++k) {
            std::array<int,3> t = ops[k].apply(u, v, w);
            mates[k] = this->index_n(t[0], t[1], t[2]);
          }
          T value = data[idx];
          for (size_t k : mates) {
            if (visited[k])
              fail("grid size is not compatible with space group");
            value = func(value, data[k]);
          }
          data[idx] = value;
          visited[idx] = true;
          for (size_t k : mates) {
            data[k] = value;
            visited[k] = true;
          }
        }
    assert(idx == data.size());
  }

  // most common symmetrize functions
  void symmetrize_min() {
    symmetrize([](T a, T b) { return (a < b || !(b == b)) ? a : b; });
  }
  void symmetrize_max() {
    symmetrize([](T a, T b) { return (a > b || !(b == b)) ? a : b; });
  }
  void symmetrize_abs_max() {
    symmetrize([](T a, T b) { return (std::abs(a) > std::abs(b) || !(b == b)) ? a : b; });
  }
  /// multiplies grid points on special position
  void symmetrize_sum() {
    symmetrize([](T a, T b) { return a + b; });
  }
  void symmetrize_nondefault(T default_) {
    symmetrize([default_](T a, T b) { return impl::is_same(a, default_) ? b : a; });
  }

  /// scale the data to get mean == 0 and rmsd == 1 (doesn't work for T=complex)
  void normalize() {
    DataStats stats = calculate_data_statistics(data);
    for (T& x : data)
      x = static_cast<T>((x - stats.dmean) / stats.rms);
  }

  // TODO: can it be replaced with interpolate_grid(dest, src, Transform(), order)?
  void resample_to(Grid<T>& dest, int order) const {
    dest.check_not_empty();
    int idx = 0;
    for (int w = 0; w < dest.nw; ++w)
      for (int v = 0; v < dest.nv; ++v)
        for (int u = 0; u < dest.nu; ++u, ++idx) {
          const Fractional f = dest.get_fractional(u, v, w);
          dest.data[idx] = interpolate(f, order);
        }
  }
};


template<typename T>
Correlation calculate_correlation(const GridBase<T>& a, const GridBase<T>& b) {
  if (a.data.size() != b.data.size() || a.nu != b.nu || a.nv != b.nv || a.nw != b.nw)
    fail("calculate_correlation(): grids have different sizes");
  Correlation corr;
  for (size_t i = 0; i != a.data.size(); ++i)
    if (!std::isnan(a.data[i]) && !std::isnan(b.data[i]))
      corr.add_point(a.data[i], b.data[i]);
  return corr;
}

} // namespace gemmi
#endif
