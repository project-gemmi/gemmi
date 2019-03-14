// Copyright 2017 Global Phasing Ltd.
//
// 3d grid used by CCP4 maps, cell-method search and hkl data.

#ifndef GEMMI_GRID_HPP_
#define GEMMI_GRID_HPP_

#include <cassert>
#include <complex>  // for std::conj
#include <functional> // for function
#include <vector>
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "util.hpp"  // for fail

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

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
template<typename T=float>
struct Grid {
  int nu = 0, nv = 0, nw = 0;
  UnitCell unit_cell;
  bool full_canonical; // grid for the whole unit cell with X,Y,Z order
  const SpaceGroup* space_group = nullptr;
  double spacing[3];

  std::vector<T> data;

  void calculate_spacing() {
    spacing[0] = 1.0 / (nu * unit_cell.ar);
    spacing[1] = 1.0 / (nv * unit_cell.br);
    spacing[2] = 1.0 / (nw * unit_cell.cr);
  }

  void set_size_without_checking(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
    calculate_spacing();
    full_canonical = true;
  }

  void set_size(int u, int v, int w) {
    if (space_group) {
      auto factors = space_group->operations().find_grid_factors();
      if (u % factors[0] != 0 || v % factors[1] != 0 || w % factors[2] != 0)
        fail("Grid not compatible with the space group " + space_group->xhm());
    }
    set_size_without_checking(u, v, w);
  }

  std::array<int, 3> pick_good_size(const std::array<double, 3>& limit,
                                    bool denser) {
    std::array<int, 3> m = {{0, 0, 0}};
    const SpaceGroup& sg = space_group ? *space_group : get_spacegroup_p1();
    GroupOps gops = sg.operations();
    std::array<int, 3> sg_fac = gops.find_grid_factors();
    for (int i = 0; i != 3; ++i) {
      for (int j = 0; j < i; ++j)
        if (fabs(limit[i] - limit[j]) < 0.5 && sg_fac[i] == sg_fac[j]) {
          m[i] = m[j];
          break;
        }
      if (m[i] == 0) {
        int f = std::max(2, sg_fac[i]);
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

  void set_size_from(std::array<double, 3> limit, bool denser) {
    auto m = pick_good_size(limit, denser);
    set_size_without_checking(m[0], m[1], m[2]);
  }

  // The resulting spacing can be smaller (if denser=true) or greater than arg.
  void set_size_from_spacing(double approx_spacing, bool denser) {
    std::array<double, 3> limit = {{1. / (unit_cell.ar * approx_spacing),
                                    1. / (unit_cell.br * approx_spacing),
                                    1. / (unit_cell.cr * approx_spacing)}};
    set_size_from(limit, denser);
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

  double get_voxel_size() const { return unit_cell.volume / (nu * nv * nw); }

  // Quick but unsafe. assumes (for efficiency) that 0 <= u < nu, etc.
  int index_q(int u, int v, int w) const { return w * nu * nv + v * nu + u; }

  // Assumes (for efficiency) that -nu <= u < 2*nu, etc.
  int index_n(int u, int v, int w) const {
    if (u >= nu) u -= nu; else if (u < 0) u += nu;
    if (v >= nv) v -= nv; else if (v < 0) v += nv;
    if (w >= nw) w -= nw; else if (w < 0) w += nw;
    return index_q(u, v, w);
  }

  // Safe but slower.
  int index_s(int u, int v, int w) const {
    return index_q(modulo(u, nu), modulo(v, nv), modulo(w, nw));
  }

  T get_value_q(int u, int v, int w) const { return data[index_q(u, v, w)]; }

  T get_value(int u, int v, int w) const { return data[index_s(u, v, w)]; }

  void set_value(int u, int v, int w, T x) { data[index_s(u, v, w)] = x; }

  void fill(T value) { std::fill(data.begin(), data.end(), value); }

  void set_points_around(const Position& ctr, double radius, T value) {
    int du = (int) std::ceil(radius / spacing[0]);
    int dv = (int) std::ceil(radius / spacing[1]);
    int dw = (int) std::ceil(radius / spacing[2]);
    if (du > nu || dv > nv || dw > nw)
      fail("Masking radius bigger than the unit cell?");
    Fractional fctr = unit_cell.fractionalize(ctr).wrap_to_unit();
    int u0 = iround(fctr.x * nu);
    int v0 = iround(fctr.y * nv);
    int w0 = iround(fctr.z * nw);
    for (int w = w0-dw; w <= w0+dw; ++w)
      for (int v = v0-dv; v <= v0+dv; ++v)
        for (int u = u0-du; u <= u0+du; ++u) {
          Fractional fdelta{fctr.x - u * (1.0 / nu),
                            fctr.y - v * (1.0 / nv),
                            fctr.z - w * (1.0 / nw)};
          fdelta.move_toward_zero_by_one();
          Position d = unit_cell.orthogonalize(fdelta);
          if (d.x*d.x + d.y*d.y + d.z*d.z < radius*radius) {
            data[index_n(u, v, w)] = value;
          }
        }
  }

  void mask_atom(double x, double y, double z, double radius) {
    set_points_around(Position(x, y, z), radius, 1);
  }

  void make_zeros_and_ones(double threshold) {
    for (auto& d : data)
      d = d > threshold ? 1 : 0;
  }

  void symmetrize_using_ops(std::vector<Op> ops, std::function<T(T, T)> func) {
    auto id = std::find(ops.begin(), ops.end(), Op::identity());
    if (id != ops.end())
      ops.erase(id);
    for (Op& op : ops) {
      op.tran[0] = op.tran[0] * nu / Op::TDEN;
      op.tran[1] = op.tran[1] * nv / Op::TDEN;
      op.tran[2] = op.tran[2] * nw / Op::TDEN;
    }
    std::vector<int> mates(ops.size(), 0);
    std::vector<bool> visited(data.size(), false);
    int idx = 0;
    for (int w = 0; w != nw; ++w)
      for (int v = 0; v != nv; ++v)
        for (int u = 0; u != nu; ++u, ++idx) {
          assert(idx == index_q(u, v, w));
          if (visited[idx])
            continue;
          for (size_t k = 0; k < ops.size(); ++k) {
            int tu = u, tv = v, tw = w;
            ops[k].apply_in_place_mult(tu, tv, tw, 1);
            mates[k] = index_n(tu, tv, tw);
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
  void symmetrize(std::function<T(T, T)> func) {
    if (!space_group || space_group->number == 1 || !full_canonical)
      return;
    symmetrize_using_ops(space_group->operations().all_ops_sorted(), func);
  }

  // two most common symmetrize functions
  void symmetrize_min() {
    symmetrize([](T a, T b) { return (a < b || !(b == b)) ? a : b; });
  }
  void symmetrize_max() {
    symmetrize([](T a, T b) { return (a > b || !(b == b)) ? a : b; });
  }

  // makes sense only for hkl data
  void add_friedel_mates() {
    const T default_val{};
    for (int u = 0; u != nu; ++u) {
      int u_ = u == 0 ? 0 : nu - u;
      for (int v = 0; v != nv; ++v) {
        int v_ = v == 0 ? 0 : nv - v;
        for (int w = 0; w != nw; ++w) {
          int idx = index_q(u, v, w);
          if (data[idx] == default_val) {
            int w_ = w == 0 ? 0 : nw - w;
            int inv_idx = index_q(u_, v_, w_);
            data[idx] = std::conj(data[inv_idx]);
          }
        }
      }
    }
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
