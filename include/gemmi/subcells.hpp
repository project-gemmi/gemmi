// Copyright 2018 Global Phasing Ltd.
//
// Cell-linked lists method for atom searching (a.k.a. grid search, binning,
// bucketing, cell techinque for neighbours search, etc).

#ifndef GEMMI_SUBCELLS_HPP_
#define GEMMI_SUBCELLS_HPP_

#include <vector>
#include "grid.hpp"
#include "model.hpp"

namespace gemmi {

struct SubCells {
  struct AtomImage {
    float pos[3];
    char altloc;
    El element;
    int image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;

    AtomImage(const Position& p, char alt, El el, int im,
              int chain, int res, int atom)
    : pos{(float)p.x, (float)p.y, (float)p.z}, altloc(alt), element(el),
      image_idx(im), chain_idx(chain), residue_idx(res), atom_idx(atom) {}
  };

  using item_type = std::vector<AtomImage>;
  Grid<item_type> grid;

  SubCells(const Model& model, const UnitCell& cell, double max_radius);

  // assumes data in [0, 1), but uses index_n to handle numeric deviations
  item_type& get_subcell(const Fractional& fr) {
    return grid.data[grid.index_n(int(fr.x * grid.nu),
                                  int(fr.y * grid.nv),
                                  int(fr.z * grid.nw))];
  }

  template<typename T>
  void for_each(const Position& pos, char alt, float radius, const T& func);

  std::vector<AtomImage*> find(const Position& pos, char alt, float radius) {
    std::vector<AtomImage*> out;
    for_each(pos, alt, radius,
             [&out](AtomImage& a, float) { out.push_back(&a); });
    return out;
  }
};


inline SubCells::SubCells(const Model& model, const UnitCell& cell,
                          double max_radius) {
  if (cell.is_crystal()) {
    grid.set_unit_cell(cell);
  } else {
    // TODO: determine boundaries and add 2 empty cells as a margin
    fail("not a crystal");
  }
  grid.set_size_from_spacing(max_radius, false);
  if (grid.nu < 3 || grid.nv < 3 || grid.nw < 3)
    grid.set_size_without_checking(std::max(grid.nu, 3), std::max(grid.nv, 3),
                                   std::max(grid.nw, 3));
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    const Chain& chain = model.chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Atom& atom = res.atoms[n_atom];
        Fractional frac0 = cell.fractionalize(atom.pos);
        {
          Fractional frac = frac0.wrap_to_unit();
          Position pos = cell.orthogonalize(frac);
          get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                         0, n_ch, n_res, n_atom);
        }
        for (int n_im = 0; n_im != (int) cell.images.size(); ++n_im) {
          Fractional frac = cell.images[n_im].apply(frac0).wrap_to_unit();
          Position pos = cell.orthogonalize(frac);
          get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                         n_im + 1, n_ch, n_res, n_atom);
        }
      }
    }
  }
}

template<typename T>
void SubCells::for_each(const Position& pos, char alt, float radius,
                        const T& func) {
  Fractional fr = grid.unit_cell.fractionalize(pos).wrap_to_unit();
  int u0 = int(fr.x * grid.nu);
  int v0 = int(fr.y * grid.nv);
  int w0 = int(fr.z * grid.nw);
  for (int w = w0 - 1; w < w0 + 2; ++w) {
    int dw = w >= grid.nw ? -1 : w < 0 ? 1 : 0;
    for (int v = v0 - 1; v < v0 + 2; ++v) {
      int dv = v >= grid.nv ? -1 : v < 0 ? 1 : 0;
      for (int u = u0 - 1; u < u0 + 2; ++u) {
        int du = u >= grid.nu ? -1 : u < 0 ? 1 : 0;
        int idx = grid.index_q(u + du * grid.nu,
                               v + dv * grid.nv,
                               w + dw * grid.nw);
        Position p = grid.unit_cell.orthogonalize(Fractional(fr.x + du,
                                                             fr.y + dv,
                                                             fr.z + dw));
        for (AtomImage& a : grid.data[idx]) {
          float dist_sq = sq((float) p.x - a.pos[0]) +
                          sq((float) p.y - a.pos[1]) +
                          sq((float) p.z - a.pos[2]);
          if (dist_sq < sq(radius) && is_same_conformer(alt, a.altloc))
            func(a, dist_sq);
        }
      }
    }
  }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
