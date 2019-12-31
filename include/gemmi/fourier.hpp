// Copyright 2019 Global Phasing Ltd.
//
// Fourier transform applied to map coefficients.

#ifndef GEMMI_FOURIER_HPP_
#define GEMMI_FOURIER_HPP_

#include <array>
#include <complex>       // for std::conj
#include "grid.hpp"      // for Grid
#include "math.hpp"      // for rad
#include "symmetry.hpp"  // for GroupOps, Op
#include "fail.hpp"      // for fail

#ifdef  __INTEL_COMPILER
// warning #2196: routine is both "inline" and "noinline"
# pragma warning disable 2196
#endif

#ifdef __MINGW32__  // MinGW may have problem with std::mutex etc
# define POCKETFFT_CACHE_SIZE 0
#endif
#define POCKETFFT_NO_MULTITHREADING
#include "third_party/pocketfft_hdronly.h"

namespace gemmi {

template<typename T>
double phase_in_angles(const std::complex<T>& v) {
  double angle = gemmi::deg(std::arg(v));
  return angle >= 0. ? angle : angle + 360.;
}


template<typename DataProxy>
std::array<int, 3> get_size_for_hkl(const DataProxy& data,
                                    std::array<int, 3> min_size,
                                    double sample_rate) {
  // adjust min_size by checking Miller indices in the data
  auto hkl_col = data.hkl_col();
  for (size_t i = 0; i < data.size(); i += data.stride())
    for (int j = 0; j != 3; ++j) {
      int v = 2 * std::abs(data.get_int(i + hkl_col[j])) + 1;
      if (v > min_size[j])
        min_size[j] = v;
    }
  std::array<double, 3> dsize{{(double)min_size[0],
                               (double)min_size[1],
                               (double)min_size[2]}};
  if (sample_rate > 0) {
    const UnitCell& cell = data.unit_cell();
    double max_1_d2 = 0;
    for (size_t i = 0; i < data.size(); i += data.stride()) {
      int h = data.get_int(i + hkl_col[0]);
      int k = data.get_int(i + hkl_col[1]);
      int l = data.get_int(i + hkl_col[2]);
      max_1_d2 = std::max(max_1_d2, cell.calculate_1_d2(h, k, l));
    }
    double inv_d_min = std::sqrt(max_1_d2);
    std::array<double, 3> cellr = {{cell.ar, cell.br, cell.cr}};
    for (int i = 0; i < 3; ++i)
      dsize[i] = std::max(dsize[i], sample_rate * inv_d_min / cellr[i]);
  }
  return good_grid_size(dsize, /*denser=*/true, data.spacegroup());
}

template<typename DataProxy>
void check_if_hkl_fits_in(const DataProxy& data, std::array<int, 3> size) {
  auto hkl_col = data.hkl_col();
  for (size_t i = 0; i < data.size(); i += data.stride())
    for (int j = 0; j != 3; ++j) {
      int index = data.get_int(i + hkl_col[j]);
      if (2 * std::abs(index) >= size[j])
        fail("grid size is too small for hkl data");
    }
}

inline float friedel_mate_value(float v) { return v; }
inline double friedel_mate_value(double v) { return v; }

template<typename T>
std::complex<T> friedel_mate_value(const std::complex<T>& v) {
  return std::conj(v);
}

template<typename T>
void add_friedel_mates(Grid<T>& grid) {
  const T default_val = T(); // initialized to 0 or 0+0i
  if (grid.hkl_orient == HklOrient::HKL) {
    for (int w = 0; w != (grid.half_l ? 1 : grid.nw); ++w) {
      int w_ = w == 0 ? 0 : grid.nw - w;
      for (int v = 0; v != grid.nv; ++v) {
        int v_ = v == 0 ? 0 : grid.nv - v;
        for (int u = 0; u != grid.nu; ++u) {
          int idx = grid.index_q(u, v, w);
          if (grid.data[idx] == default_val) {
            int u_ = u == 0 ? 0 : grid.nu - u;
            int inv_idx = grid.index_q(u_, v_, w_);
            grid.data[idx] = friedel_mate_value(grid.data[inv_idx]);
          }
        }
      }
    }
  } else { // grid.hkl_orient == HklOrient::LKH
    for (int w = 0; w != grid.nw; ++w) {
      int w_ = w == 0 ? 0 : grid.nw - w;
      for (int v = 0; v != grid.nv; ++v) {
        int v_ = v == 0 ? 0 : grid.nv - v;
        if (grid.half_l) {
          int idx = grid.index_q(0, v, w);
          if (grid.data[idx] == default_val) {
            int inv_idx = grid.index_q(0, v_, w_);
            grid.data[idx] = friedel_mate_value(grid.data[inv_idx]);
          }
        } else {
          for (int u = 0; u != grid.nu; ++u) {
            int idx = grid.index_q(u, v, w);
            if (grid.data[idx] == default_val) {
              int u_ = u == 0 ? 0 : grid.nu - u;
              int inv_idx = grid.index_q(u_, v_, w_);
              grid.data[idx] = friedel_mate_value(grid.data[inv_idx]);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DataProxy>
void initialize_hkl_grid(Grid<T>& grid, const DataProxy& data,
                         std::array<int, 3> size, bool half_l,
                         HklOrient hkl_orient) {
  if (!data.ok() || data.stride() < 5)
    fail("No data.");
  if (!data.spacegroup())
    fail("No spacegroup.");
  check_grid_factors(data.spacegroup(), size[0], size[1], size[2]);
  grid.unit_cell = data.unit_cell();
  grid.half_l = half_l;
  grid.hkl_orient = hkl_orient;
  grid.spacegroup = data.spacegroup();
  if (half_l)
    size[2] = size[2] / 2 + 1;
  if (hkl_orient == HklOrient::LKH)
    std::swap(size[0], size[2]);
  grid.set_size_without_checking(size[0], size[1], size[2]);
  grid.full_canonical = false; // disable some real-space functionality
}

// If half_l is true, grid has only data with l>=0.
// Parameter size can be obtained from get_size_for_hkl().
template<typename T, typename DataProxy>
Grid<std::complex<T>> get_f_phi_on_grid(const DataProxy& data,
                                        size_t f_col, size_t phi_col,
                                        std::array<int, 3> size, bool half_l,
                                        HklOrient hkl_orient=HklOrient::HKL) {
  Grid<std::complex<T>> grid;
  initialize_hkl_grid(grid, data, size, half_l, hkl_orient);

  if (f_col >= data.stride() || phi_col >= data.stride())
    fail("Map coefficients not found.");
  const std::complex<T> default_val; // initialized to 0+0i
  GroupOps ops = grid.spacegroup->operations();
  auto hkl_col = data.hkl_col();
  for (size_t i = 0; i < data.size(); i += data.stride()) {
    Miller hkl = data.get_hkl(i, hkl_col);
    T f = (T) data.get_num(i + f_col);
    if (f > 0.f) {
      double phi = rad(data.get_num(i + phi_col));
      for (const Op& op : ops.sym_ops) {
        auto hklp = op.apply_to_hkl(hkl);
        double shifted_phi = phi + op.phase_shift(hkl);
        int lp = hklp[2];
        if (hkl_orient == HklOrient::LKH)
          std::swap(hklp[0], hklp[2]);
        if (!half_l || lp >= 0) {
          int idx = grid.index_n(hklp[0], hklp[1], hklp[2]);
          if (grid.data[idx] == default_val)
            grid.data[idx] = std::polar(f, (T) shifted_phi);
        } else {
          int idx = grid.index_n(-hklp[0], -hklp[1], -hklp[2]);
          if (grid.data[idx] == default_val)
            grid.data[idx] = std::polar(f, (T) -shifted_phi);
        }
      }
    }
  }
  if (!ops.is_centric())
    add_friedel_mates(grid);
  return grid;
}

template<typename T, typename DataProxy>
Grid<T> get_value_on_grid(const DataProxy& data, size_t column,
                          std::array<int, 3> size, bool half_l,
                          HklOrient hkl_orient=HklOrient::HKL) {
  Grid<T> grid;
  initialize_hkl_grid(grid, data, size, half_l, hkl_orient);

  if (column >= data.stride())
    fail("Map coefficients not found.");
  GroupOps ops = grid.spacegroup->operations();
  auto hkl_col = data.hkl_col();
  for (size_t i = 0; i < data.size(); i += data.stride()) {
    Miller hkl = data.get_hkl(i, hkl_col);
    T val = (T) data.get_num(i + column);
    if (val != 0.) {
      for (const Op& op : ops.sym_ops) {
        auto hklp = op.apply_to_hkl(hkl);
        int lp = hklp[2];
        if (hkl_orient == HklOrient::LKH)
          std::swap(hklp[0], hklp[2]);
        if (!half_l || lp >= 0) {
          int idx = grid.index_n(hklp[0], hklp[1], hklp[2]);
          if (grid.data[idx] == 0.)  // 0 is the default value
            grid.data[idx] = val;
        } else {
          int idx = grid.index_n(-hklp[0], -hklp[1], -hklp[2]);
          if (grid.data[idx] == 0.)
            grid.data[idx] = val;
        }
      }
    }
  }
  if (!ops.is_centric())
    add_friedel_mates(grid);
  return grid;
}


template<typename T>
void transform_f_phi_grid_to_map_(Grid<std::complex<T>>&& hkl, Grid<T>& map) {
  // x -> conj(x) is equivalent to changing axis direction before FFT
  for (std::complex<T>& x : hkl.data)
    x.imag(-x.imag());
  map.spacegroup = hkl.spacegroup;
  map.unit_cell = hkl.unit_cell;
  map.hkl_orient = hkl.hkl_orient;
  if (hkl.hkl_orient == HklOrient::HKL) {
    int nw = hkl.half_l ? 2 * (hkl.nw - 1) : hkl.nw;
    map.set_size(hkl.nu, hkl.nv, nw);
  } else { // hkl.hkl_orient == HklOrient::LKH
    int nu = hkl.half_l ? 2 * (hkl.nu - 1) : hkl.nu;
    check_grid_factors(map.spacegroup, hkl.nw, hkl.nv, nu);
    map.set_size_without_checking(nu, hkl.nv, hkl.nw);
  }
  map.full_canonical = hkl.hkl_orient == HklOrient::HKL;
  pocketfft::shape_t shape{(size_t)hkl.nw, (size_t)hkl.nv, (size_t)hkl.nu};
  std::ptrdiff_t s = sizeof(T);
  pocketfft::stride_t stride{2*s * hkl.nv * hkl.nu, 2*s * hkl.nu, 2*s};
  pocketfft::shape_t axes{2, 1, 0};
  if (hkl.hkl_orient == HklOrient::LKH)
    std::swap(axes[0], axes[2]);
  T norm = T(1.0 / hkl.unit_cell.volume);
  if (hkl.half_l) {
    size_t last_axis = axes.back();
    axes.pop_back();
    pocketfft::c2c<T>(shape, stride, stride, axes, pocketfft::BACKWARD,
                      &hkl.data[0], &hkl.data[0], norm);
    pocketfft::stride_t stride_out{map.nv * map.nu * s, map.nu * s, s};
    shape[0] = (size_t) map.nw;
    shape[2] = (size_t) map.nu;
    pocketfft::c2r<T>(shape, stride, stride_out, last_axis, pocketfft::BACKWARD,
                      &hkl.data[0], &map.data[0], 1.0f);
  } else {
    pocketfft::c2c<T>(shape, stride, stride, axes, pocketfft::BACKWARD,
                      &hkl.data[0], &hkl.data[0], norm);
    assert(map.data.size() == hkl.data.size());
    for (size_t i = 0; i != map.data.size(); ++i)
      map.data[i] = hkl.data[i].real();
  }
}

template<typename T>
Grid<T> transform_f_phi_grid_to_map(Grid<std::complex<T>>&& hkl) {
  Grid<T> map;
  transform_f_phi_grid_to_map_(std::forward<Grid<std::complex<T>>>(hkl), map);
  return map;
}

template<typename T, typename DataProxy>
Grid<T> transform_f_phi_to_map(const DataProxy& data,
                               size_t f_col, size_t phi_col,
                               std::array<int, 3> size,
                               double sample_rate,
                               bool exact_size=false) {
  if (exact_size) {
    gemmi::check_if_hkl_fits_in(data, size);
    gemmi::check_grid_factors(data.spacegroup(), size[0], size[1], size[2]);
  } else {
    size = get_size_for_hkl(data, size, sample_rate);
  }
  return transform_f_phi_grid_to_map(get_f_phi_on_grid<T>(data, f_col, phi_col,
                                                          size, true));
}

template<typename T>
Grid<std::complex<T>> transform_map_to_f_phi(const Grid<T>& map, bool half_l) {
  Grid<std::complex<T>> hkl;
  hkl.unit_cell = map.unit_cell;
  hkl.half_l = half_l;
  hkl.spacegroup = map.spacegroup;
  int half_nw = map.nw / 2 + 1;
  hkl.set_size_without_checking(map.nu, map.nv, half_l ? half_nw : map.nw);
  hkl.full_canonical = false; // disable some real-space functionality
  T norm = T(map.unit_cell.volume / map.point_count());
  pocketfft::shape_t shape{(size_t)map.nw, (size_t)map.nv, (size_t)map.nu};
  std::ptrdiff_t s = sizeof(T);
  pocketfft::stride_t stride_in{s * hkl.nv * hkl.nu, s * hkl.nu, s};
  pocketfft::stride_t stride{2*s * hkl.nv * hkl.nu, 2*s * hkl.nu, 2*s};
  pocketfft::r2c<T>(shape, stride_in, stride, /*axis=*/0, pocketfft::FORWARD,
                    &map.data[0], &hkl.data[0], norm);
  shape[0] = half_nw;
  pocketfft::c2c<T>(shape, stride, stride, {1, 2}, pocketfft::FORWARD,
                    &hkl.data[0], &hkl.data[0], 1.0f);
  if (!half_l)  // add Friedel pairs
    for (int w = half_nw; w != hkl.nw; ++w) {
      int w_ = hkl.nw - w;
      for (int v = 0; v != hkl.nv; ++v) {
        int v_ = v == 0 ? 0 : hkl.nv - v;
        for (int u = 0; u != hkl.nu; ++u) {
          int u_ = u == 0 ? 0 : hkl.nu - u;
          int idx = hkl.index_q(u, v, w);
          int inv_idx = hkl.index_q(u_, v_, w_);
          hkl.data[idx] = hkl.data[inv_idx];  // conj() is called later
        }
      }
    }
  for (int i = 0; i != hkl.nu * hkl.nv * half_nw; ++i)
    hkl.data[i].imag(-hkl.data[i].imag());
  return hkl;
}

} // namespace gemmi
#endif
