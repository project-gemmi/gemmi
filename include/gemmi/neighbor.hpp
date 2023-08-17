// Copyright 2018 Global Phasing Ltd.
//
// Cell-linked lists method for atom searching (a.k.a. grid search, binning,
// bucketing, cell technique for neighbor search, etc).

#ifndef GEMMI_NEIGHBOR_HPP_
#define GEMMI_NEIGHBOR_HPP_

#include <vector>
#include <cmath>  // for INFINITY, sqrt

#include "fail.hpp"      // for fail
#include "grid.hpp"
#include "model.hpp"
#include "small.hpp"

namespace gemmi {

struct NeighborSearch {

  struct Mark {
    Position pos;
    char altloc;
    Element element;
    short image_idx;
    int chain_idx;
    int residue_idx;
    int atom_idx;

    Mark(const Position& p, char alt, El el, short im, int ch, int res, int atom)
    : pos(p), altloc(alt), element(el),
      image_idx(im), chain_idx(ch), residue_idx(res), atom_idx(atom) {}

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
    SmallStructure::Site& to_site(SmallStructure& small_st) const {
      return small_st.sites.at(atom_idx);
    }
    const SmallStructure::Site& to_site(const SmallStructure& small_st) const {
      return small_st.sites.at(atom_idx);
    }
  };

  Grid<std::vector<Mark>> grid;
  double radius_specified = 0.;
  Model* model = nullptr;
  SmallStructure* small_structure = nullptr;
  bool use_pbc = true;
  bool include_h = true;

  NeighborSearch() = default;
  // Model is not const so it can be modified in for_each_contact()
  NeighborSearch(Model& model_, const UnitCell& cell, double radius) {
    model = &model_;
    radius_specified = radius;
    set_bounding_cell(cell);
    set_grid_size();
  }
  NeighborSearch(SmallStructure& small_st, double radius) {
    small_structure = &small_st;
    radius_specified = radius;
    grid.unit_cell = small_st.cell;
    set_grid_size();
  }

  NeighborSearch& populate(bool include_h_=true);
  void add_chain(const Chain& chain, bool include_h_=true);
  void add_chain_n(const Chain& chain, int n_ch);
  void add_atom(const Atom& atom, int n_ch, int n_res, int n_atom);
  void add_site(const SmallStructure::Site& site, int n);

  // assumes data in [0, 1), but uses index_n to account for numerical errors
  std::vector<Mark>& get_subcell(const Fractional& fr) {
    return grid.data[grid.index_n(int(fr.x * grid.nu),
                                  int(fr.y * grid.nv),
                                  int(fr.z * grid.nw))];
  }

  template<typename Func>
  void for_each_cell(const Position& pos, const Func& func, int k=1);

  template<typename Func>
  void for_each(const Position& pos, char alt, double radius, const Func& func, int k=1) {
    if (radius <= 0)
      return;
    for_each_cell(pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
        Position p = use_pbc ? grid.unit_cell.orthogonalize(fr) : pos;
        for (Mark& m : marks) {
          double dist_sq = m.pos.dist_sq(p);
          if (dist_sq < sq(radius) && is_same_conformer(alt, m.altloc))
            func(m, dist_sq);
        }
    }, k);
  }

  int sufficient_k(double r) const {
    // .00001 is added to account for possible numeric error in r
    return r <= radius_specified ? 1 : int(r / radius_specified + 1.00001);
  }

  // with radius==0 it uses radius_specified
  std::vector<Mark*> find_atoms(const Position& pos, char alt,
                                double min_dist, double radius) {
    int k = sufficient_k(radius);
    if (radius == 0)
      radius = radius_specified;
    std::vector<Mark*> out;
    for_each(pos, alt, radius, [&](Mark& a, double dist_sq) {
        if (dist_sq >= sq(min_dist))
          out.push_back(&a);
    }, k);
    return out;
  }

  std::vector<Mark*> find_neighbors(const Atom& atom, double min_dist, double max_dist) {
    return find_atoms(atom.pos, atom.altloc, min_dist, max_dist);
  }
  std::vector<Mark*> find_site_neighbors(const SmallStructure::Site& site,
                                         double min_dist, double max_dist) {
    Position pos = grid.unit_cell.orthogonalize(site.fract);
    return find_atoms(pos, '\0', min_dist, max_dist);
  }

  std::pair<Mark*, double>
  find_nearest_atom_within_k(const Position& pos, int k, double radius) {
    Mark* mark = nullptr;
    double nearest_dist_sq = radius * radius;
    for_each_cell(pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
        Position p = use_pbc ? grid.unit_cell.orthogonalize(fr) : pos;
        for (Mark& m : marks) {
          double dist_sq = m.pos.dist_sq(p);
          if (dist_sq < nearest_dist_sq) {
            mark = &m;
            nearest_dist_sq = dist_sq;
          }
        }
    }, k);
    return {mark, nearest_dist_sq};
  }

  // it would be good to return also NearestImage
  Mark* find_nearest_atom(const Position& pos, double radius=INFINITY) {
    double r_spec = radius_specified;
    if (radius == 0.f)
      radius = r_spec;
    int max_k = std::max(std::max(std::max(grid.nu, grid.nv), grid.nw), 2);
    for (int k = 1; k < max_k; k *= 2) {
      auto result = find_nearest_atom_within_k(pos, k, radius);
      // if Mark was not found, result.second is set to radius^2.
      if (result.second < sq(k * r_spec))
        return result.first;
      if (result.first != nullptr) {
        // We found an atom, but because it was further away than k*r_spec,
        // so now it's sufficient to find the nearest atom in dist:
        double dist = std::sqrt(result.second);
        return find_nearest_atom_within_k(pos, sufficient_k(dist), radius).first;
      }
    }
    if (!use_pbc)
      // pos can be outside of bounding box. In such case, although it's slow,
      // search in all cells. Using large number that will be clipped.
      return find_nearest_atom_within_k(pos, INT_MAX/4, radius).first;
    return nullptr;
  }

  double dist_sq(const Position& pos1, const Position& pos2) const {
    return grid.unit_cell.distance_sq(pos1, pos2);
  }
  double dist(const Position& pos1, const Position& pos2) const {
    return std::sqrt(dist_sq(pos1, pos2));
  }

  FTransform get_image_transformation(int image_idx) const {
    // 0 is for identity, other indices are shifted by one.
    if (image_idx == 0)
      return Transform{};
    if ((size_t)image_idx <= grid.unit_cell.images.size())
      return grid.unit_cell.images[image_idx-1];
    fail("No such image index: " + std::to_string(image_idx));
  }

private:
  void set_grid_size() {
    // We don't use set_size_from_spacing() etc because we don't need
    // FFT-friendly size nor symmetry.
    double inv_radius = 1 / radius_specified;
    const UnitCell& uc = grid.unit_cell;
    grid.set_size_without_checking(std::max(int(inv_radius / uc.ar), 1),
                                   std::max(int(inv_radius / uc.br), 1),
                                   std::max(int(inv_radius / uc.cr), 1));
  }

  void set_bounding_cell(const UnitCell& cell) {
    use_pbc = cell.is_crystal();
    if (use_pbc) {
      grid.unit_cell = cell;
    } else {
      // cf. calculate_box()
      Box<Position> box;
      for (CRA cra : model->all())
        box.extend(cra.atom->pos);
      // The box needs to include all NCS images (strict NCS from MTRIXn).
      // To avoid additional function parameter that would pass Structure::ncs,
      // here we obtain NCS transformations from UnitCell::images.
      std::vector<FTransform> ncs = cell.get_ncs_transforms();
      if (!ncs.empty()) {
        for (CRA cra : model->all())
          // images store fractional transforms, but for non-crystal
          // it should be the same as Cartesian transform.
          for (const Transform& tr : ncs)
            box.extend(Position(tr.apply(cra.atom->pos)));
      }
      box.add_margin(0.01);
      Position size = box.get_size();
      grid.unit_cell.set(size.x, size.y, size.z, 90, 90, 90);
      grid.unit_cell.frac.vec -= grid.unit_cell.fractionalize(box.minimum);
      grid.unit_cell.orth.vec += box.minimum;
      for (const Transform& tr : ncs) {
        UnitCell& c = grid.unit_cell;
        // cf. add_ncs_images_to_cs_images()
        c.images.push_back(c.frac.combine(tr.combine(c.orth)));
      }
    }
  }
};

inline NeighborSearch& NeighborSearch::populate(bool include_h_) {
  include_h = include_h_;
  if (model) {
    for (int n_ch = 0; n_ch != (int) model->chains.size(); ++n_ch)
      add_chain_n(model->chains[n_ch], n_ch);
  } else if (small_structure) {
    for (int n = 0; n != (int) small_structure->sites.size(); ++n) {
      SmallStructure::Site& site = small_structure->sites[n];
      if (include_h || !site.element.is_hydrogen())
        add_site(site, n);
    }
  } else {
    fail("NeighborSearch not initialized");
  }
  return *this;
}

inline void NeighborSearch::add_chain(const Chain& chain, bool include_h_) {
  if (!model)
    fail("NeighborSearch.add_chain(): model not initialized yet");
  // to be safe avoid (&chain - model.chains[0]) which could be UB
  for (int n_ch = 0; n_ch != (int) model->chains.size(); ++n_ch)
    if (&model->chains[n_ch] == &chain) {
      include_h = include_h_;
      add_chain_n(chain, n_ch);
      return;
    }
  fail("NeighborSearch.add_chain(): chain not in this model");
}

inline void NeighborSearch::add_chain_n(const Chain& chain, int n_ch) {
  for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
    const Residue& res = chain.residues[n_res];
    for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
      const Atom& atom = res.atoms[n_atom];
      if (include_h || !atom.is_hydrogen())
        add_atom(atom, n_ch, n_res, n_atom);
    }
  }
}

inline void NeighborSearch::add_atom(const Atom& atom,
                                     int n_ch, int n_res, int n_atom) {
  const UnitCell& gcell = grid.unit_cell;
  Fractional frac0 = gcell.fractionalize(atom.pos);
  {
    Fractional frac = frac0.wrap_to_unit();
    // for non-crystals, frac==frac0 => pos = atom.pos
    Position pos = use_pbc ? gcell.orthogonalize(frac) : atom.pos;
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   0, n_ch, n_res, n_atom);
  }
  for (int n_im = 0; n_im != (int) gcell.images.size(); ++n_im) {
    Fractional frac = gcell.images[n_im].apply(frac0).wrap_to_unit();
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, atom.altloc, atom.element.elem,
                                   short(n_im + 1), n_ch, n_res, n_atom);
  }
}

// We exclude special position images of atoms here, but not in add_atom.
// This choice is somewhat arbitrary, but it also reflects the fact that
// in MX files occupances of atoms on special positions are (almost always)
// fractional and all images are to be taken into account.
inline void NeighborSearch::add_site(const SmallStructure::Site& site, int n) {
  const double SPECIAL_POS_TOL = 0.4;
  const UnitCell& gcell = grid.unit_cell;
  std::vector<Fractional> others;
  others.reserve(gcell.images.size());
  Fractional frac0 = site.fract.wrap_to_unit();
  {
    Position pos = gcell.orthogonalize(frac0);
    get_subcell(frac0).emplace_back(pos, '\0', site.element.elem, 0, -1, -1, n);
  }
  for (int n_im = 0; n_im != (int) gcell.images.size(); ++n_im) {
    Fractional frac = gcell.images[n_im].apply(site.fract).wrap_to_unit();
    if (gcell.distance_sq(frac, frac0) < sq(SPECIAL_POS_TOL) ||
        std::any_of(others.begin(), others.end(), [&](const Fractional& f) {
          return gcell.distance_sq(frac, f) < sq(SPECIAL_POS_TOL);
        }))
      continue;
    Position pos = gcell.orthogonalize(frac);
    get_subcell(frac).emplace_back(pos, '\0', site.element.elem,
                                   short(n_im + 1), -1, -1, n);
    others.push_back(frac);
  }
}

template<typename Func>
void NeighborSearch::for_each_cell(const Position& pos, const Func& func, int k) {
  Fractional fr = grid.unit_cell.fractionalize(pos);
  if (use_pbc)
    fr = fr.wrap_to_unit();
  int u0 = int(fr.x * grid.nu) - k;
  int v0 = int(fr.y * grid.nv) - k;
  int w0 = int(fr.z * grid.nw) - k;
  int uend = u0 + 2 * k + 1;
  int vend = v0 + 2 * k + 1;
  int wend = w0 + 2 * k + 1;
  if (use_pbc) {
    auto shift = [](int j, int n) {
      if (j < 0)
        return (j + 1) / n - 1;
      if (j >= n)
        return j / n;
      return 0;
    };
    for (int w = w0; w < wend; ++w) {
      int dw = shift(w, grid.nw);
      for (int v = v0; v < vend; ++v) {
        int dv = shift(v, grid.nv);
        size_t idx0 = grid.index_q(0, v - dv * grid.nv, w - dw * grid.nw);
        for (int u = u0; u < uend; ++u) {
          int du = shift(u, grid.nu);
          size_t idx = idx0 + (u - du * grid.nu);
          func(grid.data[idx], Fractional(fr.x - du, fr.y - dv, fr.z - dw));
        }
      }
    }
  } else {
    u0 = std::max(0, u0);
    v0 = std::max(0, v0);
    w0 = std::max(0, w0);
    uend = std::min(uend, grid.nu);
    vend = std::min(vend, grid.nv);
    wend = std::min(wend, grid.nw);
    for (int w = w0; w < wend; ++w)
      for (int v = v0; v < vend; ++v)
        for (int u = u0; u < uend; ++u) {
          size_t idx = grid.index_q(u, v, w);
          func(grid.data[idx], fr);
        }
  }
}

} // namespace gemmi
#endif
