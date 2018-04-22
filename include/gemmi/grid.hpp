// Copyright 2017 Global Phasing Ltd.
//
// 3d grid used by CCP4 maps and cell-method search.

#ifndef GEMMI_GRID_HPP_
#define GEMMI_GRID_HPP_

#include <cassert>
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
  for (int k : {2, 3, 5})
    while (n % k == 0)
      n /= k;
  return n == 1 || n == -1;
}

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
template<typename T=float>
struct Grid {
  int nu, nv, nw;
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

  // The resulting spacing can be smaller (if denser=true) or greater than arg.
  void set_size_from_spacing(double approx_spacing, bool denser) {
    const SpaceGroup& sg = space_group ? *space_group : get_spacegroup_p1();
    std::array<int, 3> sg_fac = sg.operations().find_grid_factors();
    int m[3];
    double abc[3] = {1./unit_cell.ar, 1./unit_cell.br, 1./unit_cell.cr};
    for (int i = 0; i != 3; ++i) {
      if (i > 0 && fabs(abc[i] - abc[0]) < 0.1 && sg_fac[i] == sg_fac[0]) {
        m[i] = m[0];
        continue;
      }
      if (i > 1 && fabs(abc[i] - abc[1]) < 0.1 && sg_fac[i] == sg_fac[1]) {
        m[i] = m[1];
        continue;
      }
      int f = std::max(2, sg_fac[i]);
      int n;
      if (denser) {
        n = int(std::ceil(abc[i] / (approx_spacing * f)));
        while (!has_small_factorization(n))
          ++n;
      } else {
        n = int(std::floor(abc[i] / (approx_spacing * f)));
        if (n > 1)
          while (!has_small_factorization(n))
            --n;
        else
          n = 1;
      }
      m[i] = n * f;
    }
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

  // _s stands for safe
  T get_value_s(int u, int v, int w) const { return data[index_s(u, v, w)]; }

  void set_value_s(int u, int v, int w, T x) { data[index_s(u, v, w)] = x; }

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

  // Use provided function to reduce values of all symmetry mates of each
  // grid points, then assign the result to all the points.
  void symmetrize(std::function<T(T, T)> func) {
    if (!space_group || space_group->number == 1 || !full_canonical)
      return;
    std::vector<Op> ops = space_group->operations().all_ops_sorted();
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

  // two most common symmetrize functions
  void symmetrize_min() { symmetrize([](T a, T b) { return a < b ? a : b; }); }
  void symmetrize_max() { symmetrize([](T a, T b) { return a > b ? a : b; }); }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
