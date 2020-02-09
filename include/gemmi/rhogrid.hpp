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
// - set d_min and optionally also other parameters
// - set grid's unit cell and space group using set_grid_cell_and_spacegroup()
// - check that Table has SF coefficients for all elements that are to be used
// - call put_model_density_on_grid()
// - if blur is used, the SF must be then multiplied by the value returned by
//   reciprocal_space_multiplier()
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
    auto& scat = Table::get(atom.element);
    float fprime = fprimes[(int)atom.element.elem];
    double b = atom.b_iso + blur;
    auto precal = scat.precalculate_density(b, fprime);
    double radius = determine_cutoff_radius(
                              [&](float r) { return precal.calculate(r*r); },
                              r_cut);
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    grid.use_points_around(fpos, radius, [&](Real& point, double r2) {
        point += Real(atom.occ * precal.calculate((Real)r2));
    }, /*fail_on_too_large_radius=*/false);
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
    grid.symmetrize([](double a, double b) { return a + b; });
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
