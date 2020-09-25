// Copyright 2020 Global Phasing Ltd.

#ifndef GEMMI_MASK_HPP_
#define GEMMI_MASK_HPP_

#include <gemmi/grid.hpp>
#include <gemmi/model.hpp>

namespace gemmi {

template<typename T>
void mask_points_in_constant_radius(Grid<T>& mask, const Model& model,
                                    double radius, T value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, radius, value);
}

template<typename T>
void mask_points_in_vdw_radius(Grid<T>& mask, const Model& model, double r_probe, T value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, atom.element.vdw_r() + r_probe, value);
}

template<typename T>
void set_margin_around_mask(Grid<T>& mask, double r, T old_value, T new_value) {
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
            T& point = mask.data[mask.index_n(u+d[0], v+d[1], w+d[2])];
            if (point != old_value)
              point = new_value;
          }
        }
}

} // namespace gemmi
#endif
