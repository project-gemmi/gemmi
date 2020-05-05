// Copyright 2019 Global Phasing Ltd.
//
// Prepares a grid with values of electron density of a model.

#ifndef GEMMI_RHOGRID_HPP_
#define GEMMI_RHOGRID_HPP_

#include <complex>
#include "grid.hpp"    // for Grid
#include "model.hpp"   // for Structure, ...

namespace gemmi {

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
};

} // namespace gemmi
#endif
