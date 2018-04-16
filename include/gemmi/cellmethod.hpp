// Copyright 2018 Global Phasing Ltd.
//
// Cell-method search.

#ifndef GEMMI_CELLMETHOD_HPP_
#define GEMMI_CELLMETHOD_HPP_

#include <vector>
#include "grid.hpp"
#include "model.hpp"

namespace gemmi {

struct CellMethod {
  struct AtomImage {
    float pos[3];
    int image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;
  };
  struct ChainImage {
    int image_idx;
    int chain_idx;
  };
  Grid<std::vector<AtomImage>> grid;

  CellMethod(const Structure& st, double max_radius);
  std::vector<AtomImage*> find(const Position& pos, double radius);
  std::vector<ChainImage> find_chains(const Position& pos, double radius);
};


inline CellMethod::CellMethod(const Structure& st, double max_radius) {
  grid.set_unit_cell(st.cell);
  grid.set_size_from_spacing(max_radius, false);
  const Model& model = st.models.at(0);
  for (int ch_idx = 0; ch_idx != (int) model.chains.size(); ++ch_idx) {
    const Chain& chain = model.chains[ch_idx];
    for (int res_idx = 0; res_idx != (int) chain.residues.size(); ++res_idx) {
      const Residue& res = chain.residues[res_idx];
      for (int a_idx = 0; a_idx != (int) res.atoms.size(); ++a_idx) {
        Fractional frac = st.cell.fractionalize(res.atoms[a_idx].pos);
        for (int im_idx = 0; im_idx != (int) st.cell.images.size(); ++im_idx) {
          Fractional ifrac = st.cell.images[im_idx].apply(frac).wrap_to_unit();
          int idx = grid.index_n((int) (ifrac.x * grid.nu),
                                 (int) (ifrac.y * grid.nv),
                                 (int) (ifrac.z * grid.nw));
          Position pos = st.cell.orthogonalize(ifrac);
          grid.data[idx].push_back({{(float)pos.x, (float)pos.y, (float)pos.z},
                                     im_idx, ch_idx, res_idx, a_idx});
        }
      }
    }
  }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
