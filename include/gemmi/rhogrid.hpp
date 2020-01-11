// Copyright 2019 Global Phasing Ltd.
//
// Prepares a grid with values of electron density of a model.

#ifndef GEMMI_RHOGRID_HPP_
#define GEMMI_RHOGRID_HPP_

#include <complex>
#include "grid.hpp"    // for Grid
#include "model.hpp"   // for Structure, ...

namespace gemmi {

template <typename Coef>
double determine_effective_radius(const Coef& coef, double b, double cutoff) {
  double x1 = 3.5;
  double y1 = coef.calculate_density(x1*x1, b);
  double x2 = x1;
  double y2 = y1;
  if (y1 < cutoff)
    while (y1 < cutoff) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5;
      y1 = coef.calculate_density(x1*x1, b);
    }
  else
    while (y2 > cutoff) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5;
      y2 = coef.calculate_density(x2*x2, b);
    }
  while (x2 - x1 > 0.02) {
    double new_x = 0.5 * (x2 + x1);
    double new_y = coef.calculate_density(new_x*new_x, b);
    if (new_y < cutoff) {
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

  // pre: check if Table::has(atom.element)
  void add_atom_density_to_grid(const Atom& atom) {
    auto& scat = Table::get(atom.element);
    double b = atom.b_iso + blur;
    double radius = determine_effective_radius(scat, (float) b, r_cut);
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    grid.use_points_around(fpos, radius, [&](Real& point, double r2) {
        point += Real(atom.occ * scat.calculate_density((Real)r2, (Real)b));
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
    grid.spacegroup = find_spacegroup_by_name(st.spacegroup_hm);
  }

  // The argument is 1/d^2 - as outputted by unit_cell.calculate_1_d2(hkl).
  double reciprocal_space_multiplier(double inv_d2) {
    return std::exp(blur * 0.25 * inv_d2);
  }
};

} // namespace gemmi
#endif
