// Copyright 2018 Global Phasing Ltd.
//
// Cell-linked lists method for atom searching (a.k.a. grid search, binning,
// bucketing, cell technique for neighbours search, etc).

#ifndef GEMMI_SUBCELLS_HPP_
#define GEMMI_SUBCELLS_HPP_

#include <vector>
#include <cmath>  // for INFINITY, sqrt

#include "fail.hpp"       // for fail
#include "grid.hpp"
#include "model.hpp"

namespace gemmi {

struct SubCells {

  struct Mark {
    float x, y, z;
    char altloc;
    El element;
    int image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;

    Mark(const Position& p, char alt, El el, int im, int ch, int res, int atom)
    : x(float(p.x)), y(float(p.y)), z(float(p.z)), altloc(alt), element(el),
      image_idx(im), chain_idx(ch), residue_idx(res), atom_idx(atom) {}

    Position pos() const { return {x, y, z}; }

    CRA to_cra(Model& mdl) const {
      Chain& c = mdl.chains.at(chain_idx);
      Residue& r = c.residues.at(residue_idx);
      Atom& a = r.atoms.at(atom_idx);
      return {&c, &r, &a};
    }
    const_CRA to_cra(const Model& mdl) const {
      const Chain& c = mdl.chains.at(chain_idx);
      const Residue& r = c.residues.at(residue_idx);
      const Atom& a = r.atoms.at(atom_idx);
      return {&c, &r, &a};
    }

    float dist_sq(const Position& p) const {
      return sq((float)p.x - x) + sq((float)p.y - y) + sq((float)p.z - z);
    }
  };

  using item_type = std::vector<Mark>;
  Grid<item_type> grid;
  double radius_specified = 0.;
  const Model* model = nullptr;

  SubCells() = default;
  SubCells(const Model& model_, const UnitCell& cell, double max_radius) {
    initialize(model_, cell, max_radius);
  }
  void initialize(const Model& model, const UnitCell& cell, double max_radius);
  void populate(bool include_h=true);
  void add_atom(const Atom& atom, int n_ch, int n_res, int n_atom);

  // assumes data in [0, 1), but uses index_n to handle numeric deviations
  item_type& get_subcell(const Fractional& fr) {
    return grid.data[grid.index_n(int(fr.x * grid.nu),
                                  int(fr.y * grid.nv),
                                  int(fr.z * grid.nw))];
  }

  template<typename T>
  void for_each(const Position& pos, char alt, float radius, const T& func);

  std::vector<Mark*> find_atoms(const Position& pos, char alt, float radius) {
    std::vector<Mark*> out;
    for_each(pos, alt, radius, [&out](Mark& a, float) { out.push_back(&a); });
    return out;
  }

  std::vector<Mark*> find_neighbors(const Atom& atom,
                                    float min_dist, float max_dist) {
    std::vector<Mark*> out;
    for_each(atom.pos, atom.altloc, max_dist, [&](Mark& a, float dist_sq) {
        if (dist_sq > sq(min_dist))
          out.push_back(&a);
    });
    return out;
  }

  Mark* find_nearest_atom(const Position& pos) {
    Mark* mark = nullptr;
    float nearest_dist_sq = float(radius_specified * radius_specified);
    for_each(pos, '\0', nearest_dist_sq, [&](Mark& a, float dist_sq) {
        if (dist_sq < nearest_dist_sq) {
          mark = &a;
          nearest_dist_sq = dist_sq;
        }
    });
    return mark;
  }

  float dist_sq(const Position& pos1, const Position& pos2) const {
    return (float) grid.unit_cell.distance_sq(pos1, pos2);
  }
  float dist(const Position& pos1, const Position& pos2) const {
    return std::sqrt(dist_sq(pos1, pos2));
  }
};


inline void SubCells::initialize(const Model& model_, const UnitCell& cell,
                                 double max_radius) {
  model = &model_;
  radius_specified = max_radius;
  if (cell.is_crystal()) {
    grid.unit_cell = cell;
  } else {
    BoundingBox box;
    for (const Chain& chain : model->chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          box.add(atom.pos);
    double margin = 4 * max_radius;
    Vec3 size = box.get_size() + Vec3(margin, margin, margin);
    grid.unit_cell.set(size.x, size.y, size.z, 90, 90, 90);
  }
  grid.set_size_from_spacing(max_radius, false);
  if (grid.nu < 3 || grid.nv < 3 || grid.nw < 3)
    grid.set_size_without_checking(std::max(grid.nu, 3), std::max(grid.nv, 3),
                                   std::max(grid.nw, 3));
}

inline void SubCells::populate(bool include_h) {
  if (!model)
    fail("SubCells not initialized");
  for (int n_ch = 0; n_ch != (int) model->chains.size(); ++n_ch) {
    const Chain& chain = model->chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Atom& atom = res.atoms[n_atom];
        if (include_h || !atom.is_hydrogen())
          add_atom(atom, n_ch, n_res, n_atom);
      }
    }
  }
}

inline void SubCells::add_atom(const Atom& atom,
                               int n_ch, int n_res, int n_atom) {
  const UnitCell& gcell = grid.unit_cell;
  Fractional frac0 = gcell.fractionalize(atom.pos);
  {
    Fractional frac = frac0.wrap_to_unit();
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   0, n_ch, n_res, n_atom);
  }
  for (int n_im = 0; n_im != (int) gcell.images.size(); ++n_im) {
    Fractional frac = gcell.images[n_im].apply(frac0).wrap_to_unit();
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   n_im + 1, n_ch, n_res, n_atom);
  }
}

template<typename T>
void SubCells::for_each(const Position& pos, char alt, float radius,
                        const T& func) {
  if (radius <= 0.f)
    return;
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
        for (Mark& a : grid.data[idx]) {
          float dist_sq = a.dist_sq(p);
          if (a.dist_sq(p) < sq(radius) && is_same_conformer(alt, a.altloc))
            func(a, dist_sq);
        }
      }
    }
  }
}

} // namespace gemmi
#endif
