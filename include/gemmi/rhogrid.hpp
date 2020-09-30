// Copyright 2019 Global Phasing Ltd.
//
// Tools to make a grid with:
// - values of electron density of a model,
// - bulk solvent mask.

#ifndef GEMMI_RHOGRID_HPP_
#define GEMMI_RHOGRID_HPP_

#include <cassert>
#include "grid.hpp"    // for Grid
#include "model.hpp"   // for Structure, ...

namespace gemmi {

// mask utilities
template<typename Real>
void mask_points_in_constant_radius(Grid<Real>& mask, const Model& model,
                                    double radius, Real value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, radius, value);
}

template<typename Real>
void mask_points_in_vdw_radius(Grid<Real>& mask, const Model& model,
                               double r_probe, Real value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, atom.element.vdw_r() + r_probe, value);
}

template<typename Real>
void set_margin_around_mask(Grid<Real>& mask, double r, Real old_value, Real new_value) {
  int du = (int) std::floor(r / mask.spacing[0]);
  int dv = (int) std::floor(r / mask.spacing[1]);
  int dw = (int) std::floor(r / mask.spacing[2]);
  if (2 * du >= mask.nu || 2 * dv >= mask.nv || 2 * dw >= mask.nw)
    fail("grid operation failed: radius bigger than half the unit cell?");
  std::vector<std::array<int,3>> neighbour_indices;
  for (int w = -dw; w <= dw; ++w)
    for (int v = -dv; v <= dv; ++v)
      for (int u = -du; u <= du; ++u) {
        Fractional fdelta{u * (1.0 / mask.nu), v * (1.0 / mask.nv), w * (1.0 / mask.nw)};
        Position delta = mask.unit_cell.orthogonalize_difference(fdelta);
        if (delta.length_sq() <= r * r)
          neighbour_indices.push_back({{u, v, w}});
      }
  for (int w = 0; w < mask.nw; ++w)
    for (int v = 0; v < mask.nv; ++v)
      for (int u = 0; u < mask.nu; ++u)
        if (mask.data[mask.index_q(u, v, w)] == old_value) {
          for (const std::array<int,3>& d : neighbour_indices) {
            Real& point = mask.data[mask.index_n(u+d[0], v+d[1], w+d[2])];
            if (point != old_value)
              point = new_value;
          }
        }
}


template <typename F>
double determine_cutoff_radius(const F& func, float cutoff_level) {
  float x1 = 3.5f;
  float y1 = func(x1);
  float x2 = x1;
  float y2 = y1;
  if (y1 < cutoff_level)
    while (y1 < cutoff_level) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5f;
      y1 = func(x1);
    }
  else
    while (y2 > cutoff_level) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5f;
      y2 = func(x2);
    }
  while (x2 - x1 > 0.02f) {
    float new_x = 0.5f * (x2 + x1);
    float new_y = func(new_x);
    if (new_y < cutoff_level) {
      x2 = new_x;
      y2 = new_y;
    } else {
      x1 = new_x;
      y1 = new_y;
    }
  }
  return x2;
}

// Usual usage:
// - set d_min and optionally also other parameters,
// - set fprimes to f' values for your wavelength (see fprime.hpp)
// - use set_grid_cell_and_spacegroup() to set grid's unit cell and space group
// - check that Table has SF coefficients for all elements that are to be used
// - call put_model_density_on_grid()
// - do FFT using transform_map_to_f_phi()
// - if blur is used, multiply the SF by reciprocal_space_multiplier()
template <typename Table, typename Real>
struct DensityCalculator {
  Grid<Real> grid;
  double d_min = 0.;
  double rate = 1.5;
  double blur = 0.;
  float r_cut = 5e-5f;
  std::vector<float> fprimes = std::vector<float>((int)El::END, 0.f);
  // parameters for used only in put_solvent_mask_on_grid()
  double rprobe = 1.0;
  double rshrink = 1.1;

  // pre: check if Table::has(atom.element)
  void add_atom_density_to_grid(const Atom& atom) {
    constexpr double UtoB = 8 * sq(pi());
    auto& scat = Table::get(atom.element);
    float fprime = fprimes[atom.element.ordinal()];
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    if (!atom.aniso.nonzero()) {
      // isotropic
      double b = atom.b_iso + blur;
      auto precal = scat.precalculate_density_iso(b, fprime);
      double radius = determine_cutoff_radius(
          [&](float r) { return (float)precal.calculate(r*r); },
          r_cut);
      grid.use_points_around(fpos, radius, [&](Real& point, double r2) {
          point += Real(atom.occ * precal.calculate((Real)r2));
      }, /*fail_on_too_large_radius=*/false);
    } else {
      // anisotropic
      SMat33<double> aniso_b = atom.aniso.scaled(UtoB).added_kI(blur);
      // rough estimate, so we don't calculate eigenvalues
      double b_max = std::max(std::max(aniso_b.u11, aniso_b.u22), aniso_b.u33);
      auto precal_iso = scat.precalculate_density_iso(b_max, fprime);
      double radius = determine_cutoff_radius(
          [&](float r) { return (float)precal_iso.calculate(r*r); },
          r_cut);
      auto precal = scat.precalculate_density_aniso_b(aniso_b, fprime);
      int du = (int) std::ceil(radius / grid.spacing[0]);
      int dv = (int) std::ceil(radius / grid.spacing[1]);
      int dw = (int) std::ceil(radius / grid.spacing[2]);
      grid.use_points_in_box(fpos, du, dv, dw,
                             [&](Real& point, const Position& delta) {
        if (delta.length_sq() < radius * radius)
          point += Real(atom.occ * precal.calculate(delta));
      }, false);
    }
  }

  void add_model_density_to_grid(const Model& model) {
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          add_atom_density_to_grid(atom);
  }

  void put_model_density_on_grid(const Model& model) {
    grid.data.clear();
    grid.set_size_from_spacing(d_min / (2 * rate), true);
    add_model_density_to_grid(model);
    grid.symmetrize([](Real a, Real b) { return a + b; });
  }

  void set_grid_cell_and_spacegroup(const Structure& st) {
    grid.unit_cell = st.cell;
    grid.spacegroup = st.find_spacegroup();
  }

  // The argument is 1/d^2 - as outputted by unit_cell.calculate_1_d2(hkl).
  double reciprocal_space_multiplier(double inv_d2) {
    return std::exp(blur * 0.25 * inv_d2);
  }

  void put_solvent_mask_on_grid(const Model& model) {
    assert(!grid.data.empty());
    std::fill(grid.data.begin(), grid.data.end(), 1);
    mask_points_in_vdw_radius<Real>(grid, model, rprobe, 0);
    set_margin_around_mask<Real>(grid, rshrink, 1, -1);
    grid.change_values(-1, 1);
  }
};

} // namespace gemmi
#endif
