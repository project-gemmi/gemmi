// Copyright 2019 Global Phasing Ltd.
//
// Fourier transform applied to map coefficients.

#ifndef GEMMI_FOURIER_HPP_
#define GEMMI_FOURIER_HPP_

#include <array>
#include "grid.hpp"      // for Grid
#include "math.hpp"      // for rad
#include "symmetry.hpp"  // for GroupOps, Op
#include "util.hpp"      // for fail

#if !(__GNUC__+0 >= 5 || __clang__+0 >= 5)
#define POCKETFFT_NO_VECTORS
#endif
#include <pocketfft_hdronly.h>

namespace gemmi {

// If half_l is true, grid has only data with l>=0.
template<typename T, typename DataProxy>
Grid<std::complex<T>> get_f_phi_on_grid(const DataProxy& data,
                                        size_t f_col, size_t phi_col,
                                        bool half_l,
                                        std::array<int, 3> min_size) {
  if (!data.ok() || data.stride() < 5)
    fail("No data.");
  if (!data.spacegroup())
    fail("No spacegroup.");
  Grid<std::complex<T>> grid;
  grid.unit_cell = data.unit_cell();
  grid.spacegroup = data.spacegroup();
  auto hkl_col = data.hkl_col();

  { // set grid size
    for (size_t i = 0; i < data.size(); i += data.stride())
      for (int j = 0; j != 3; ++j) {
        int v = 2 * std::abs(data.get_int(i + hkl_col[j])) + 1;
        if (v > min_size[j])
          min_size[j] = v;
      }
    std::array<int, 3> size = grid.pick_good_size(
            {{(double)min_size[0], (double)min_size[1], (double)min_size[2]}},
            /*denser=*/true);
    if (half_l)
      size[2] = size[2] / 2 + 1;
    // index H is fast and L is slow here -- not ideal
    grid.set_size_without_checking(size[0], size[1], size[2]);
    grid.full_canonical = false; // disable some real-space functionality
  }

  if (f_col >= data.stride() || phi_col >= data.stride())
    fail("Map coefficients not found.");
  const std::complex<T> default_val; // initialized to 0+0i
  GroupOps ops = grid.spacegroup->operations();
  for (size_t i = 0; i < data.size(); i += data.stride()) {
    int h = data.get_int(i + hkl_col[0]);
    int k = data.get_int(i + hkl_col[1]);
    int l = data.get_int(i + hkl_col[2]);
    T f = (T) data.get_num(i + f_col);
    if (f > 0.f) {
      double phi = rad(data.get_num(i + phi_col));
      for (const Op& op : ops.sym_ops) {
        auto hklp = op.apply_to_hkl({{h, k, l}});
        double shifted_phi = phi + op.phase_shift(h, k, l);
        if (!half_l || hklp[2] >= 0) {
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
  // add Friedel pairs
  if (!ops.is_centric()) {
    for (int w = 0; w != (half_l ? 1 : grid.nw); ++w) {
      int w_ = w == 0 ? 0 : grid.nw - w;
      for (int v = 0; v != grid.nv; ++v) {
        int v_ = v == 0 ? 0 : grid.nv - v;
        for (int u = 0; u != grid.nu; ++u) {
          int idx = grid.index_q(u, v, w);
          if (grid.data[idx] == default_val) {
            int u_ = u == 0 ? 0 : grid.nu - u;
            int inv_idx = grid.index_q(u_, v_, w_);
            grid.data[idx] = std::conj(grid.data[inv_idx]);
          }
        }
      }
    }
  }
  return grid;
}


template<typename T>
Grid<T> transform_f_phi_half_to_map(Grid<std::complex<T>>&& hkl) {
  Grid<T> map;
  // x -> conj(x) is equivalent to changing axis direction before FFT
  for (std::complex<T>& x : hkl.data)
    x.imag(-x.imag());
  map.spacegroup = hkl.spacegroup;
  map.unit_cell = hkl.unit_cell;
  int full_nw = 2 * (hkl.nw - 1);
  map.set_size(hkl.nu, hkl.nv, full_nw);
  map.full_canonical = true;
  pocketfft::shape_t shape{(size_t)hkl.nw, (size_t)hkl.nv, (size_t)hkl.nu};
  std::ptrdiff_t s = sizeof(T);
  pocketfft::stride_t stride{hkl.nv * hkl.nu * 2*s, hkl.nu * 2*s, 2*s};
  pocketfft::c2c<T>(shape, stride, stride, {1, 2}, pocketfft::BACKWARD,
                    &hkl.data[0], &hkl.data[0], 1.0f);
  pocketfft::stride_t stride_out{hkl.nv * hkl.nu * s, hkl.nu * s, s};
  shape[0] = full_nw;
  T norm = T(1.0 / hkl.unit_cell.volume);
  pocketfft::c2r<T>(shape, stride, stride_out, /*axis=*/0,
                    &hkl.data[0], &map.data[0], norm);
  return map;
}


template<typename T, typename DataProxy>
Grid<T> transform_f_phi_to_map(const DataProxy& data,
                               size_t f_col, size_t phi_col,
                               std::array<int, 3> min_size) {
  return transform_f_phi_half_to_map(
                  get_f_phi_on_grid<T>(data, f_col, phi_col, true, min_size));
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
