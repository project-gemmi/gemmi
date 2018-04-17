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
    int image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;

    AtomImage(const Position& p, int im, int chain, int residue, int atom)
    : pos{(float)p.x, (float)p.y, (float)p.z},
      image_idx(im), chain_idx(chain), residue_idx(residue), atom_idx(atom) {}
  };

  using item_type = std::vector<AtomImage>;
  Grid<item_type> grid;

  SubCells(const Structure& st, double max_radius);

  // assumes data in [0, 1), but uses index_n to handle numeric deviations
  item_type& get_subcell(const Fractional& fr) {
    return grid.data[grid.index_n(int(fr.x * grid.nu),
                                  int(fr.y * grid.nv),
                                  int(fr.z * grid.nw))];
  }

  template<typename T>
  void for_each(const Position& pos, float radius, const T& func);

  std::vector<AtomImage*> find(const Position& pos, float radius) {
    std::vector<AtomImage*> out;
    for_each(pos, radius, [&out](AtomImage& a) { out.push_back(&a); });
    return out;
  }
};


inline SubCells::SubCells(const Structure& st, double max_radius) {
  grid.set_unit_cell(st.cell);
  grid.set_size_from_spacing(max_radius, false);
  const Model& model = st.models.at(0);
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    const Chain& chain = model.chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Position& orig_pos = res.atoms[n_atom].pos;
        Fractional frac = st.cell.fractionalize(orig_pos);
        get_subcell(frac).emplace_back(orig_pos, 0, n_ch, n_res, n_atom);
        for (int n_im = 0; n_im != (int) st.cell.images.size(); ++n_im) {
          Fractional ifrac = st.cell.images[n_im].apply(frac).wrap_to_unit();
          Position pos = st.cell.orthogonalize(ifrac);
          get_subcell(ifrac).emplace_back(pos, n_im + 1, n_ch, n_res, n_atom);
        }
      }
    }
  }
}

template<typename T>
void SubCells::for_each(const Position& pos, float radius, const T& func) {
  Fractional fr = grid.unit_cell.fractionalize(pos).wrap_to_unit();
  int u0 = int(fr.x * grid.nu);
  int v0 = int(fr.y * grid.nv);
  int w0 = int(fr.z * grid.nw);
  for (int w = w0 - 1; w <= w0 + 1; ++w) {
    int dw = w >= grid.nw ? -1 : w < 0 ? 1 : 0;
    for (int v = v0 - 1; v <= v0 + 1; ++v) {
      int dv = v >= grid.nv ? -1 : v < 0 ? 1 : 0;
      for (int u = u0 - 1; u <= u0 + 1; ++u) {
        int du = u >= grid.nu ? -1 : u < 0 ? 1 : 0;
        int idx = grid.index_q(u + du * grid.nu,
                               v + dv * grid.nv,
                               w + dw * grid.nw);
        Position p = grid.unit_cell.orthogonalize(Fractional(fr.x + du,
                                                             fr.y + dv,
                                                             fr.z + dw));
        for (AtomImage& a : grid.data[idx]) {
          float dx = (float) p.x - a.pos[0];
          float dy = (float) p.y - a.pos[1];
          float dz = (float) p.z - a.pos[2];
          if (sq(dx) + sq(dy) + sq(dz) < sq(radius))
            func(a);
        }
      }
    }
  }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
