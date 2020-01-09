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

template <typename T>
void set_grid_cell_and_spacegroup(Grid<T>& grid, const Structure& st) {
  grid.unit_cell = st.cell;
  grid.spacegroup = find_spacegroup_by_name(st.spacegroup_hm);
}

struct RhoGridOptions {
  double d_min;
  double rate = 1.5;
  double blur = 0.;
  float r_cut = 5e-5f;
};

// pre: check if Table::has(atom.element)
template <typename Table, typename T>
void add_atom_density_to_grid(const Atom& atom, Grid<T>& grid,
                              const RhoGridOptions& opt) {
  auto& scat = Table::get(atom.element);
  double b = atom.b_iso + opt.blur;
  double radius = determine_effective_radius(scat, (float) b, opt.r_cut);
  Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
  grid.use_points_around(fpos, radius, [&](T& point, double r2) {
      point += T(atom.occ * scat.calculate_density((T)r2, (T)b));
  }, /*fail_on_too_large_radius=*/false);
}

template <typename Table, typename T>
void add_model_density_to_grid(const Model& model, Grid<T>& grid,
                               const RhoGridOptions& opt) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        add_atom_density_to_grid<Table>(atom, grid, opt);
}

template <typename Table, typename T>
void put_first_model_density_on_grid(const Structure& st, Grid<T>& grid,
                                     const RhoGridOptions& opt) {
  grid.data.clear();
  set_grid_cell_and_spacegroup(grid, st);
  grid.set_size_from_spacing(opt.d_min / (2 * opt.rate), true);
  add_model_density_to_grid<Table>(st.models.at(0), grid, opt);
  grid.symmetrize([](double a, double b) { return a + b; });
}

} // namespace gemmi
#endif
