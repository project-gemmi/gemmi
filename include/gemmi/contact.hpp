// Copyright 2020 Global Phasing Ltd.
//
// Contact search, based on SubCells from subcells.hpp.

#ifndef GEMMI_CONTACT_HPP_
#define GEMMI_CONTACT_HPP_

#include "model.hpp"
#include "subcells.hpp"
#include "polyheur.hpp"  // for check_polymer_type, are_connected

namespace gemmi {

struct ContactSearch {
  // parameters used to configure the search
  float search_radius;
  bool skip_intra_residue = true;  // ignore "contacts" in the same residue
  bool skip_adjacent_residue = false;  // ignore contacts with prev/next res.
  bool twice = false;  // report both A-B and B-A
  float special_pos_cutoff_sq = 0.8f * 0.8f;
  float min_occupancy = 0.f;
  std::vector<float> radii;

  ContactSearch(float radius) noexcept : search_radius(radius) {}

  // a helper function that sets per-atom radii basing on covalent_radius()
  void setup_atomic_radii(double multiplier, double tolerance) {
    radii.resize((size_t)El::END);
    for (int i = 0; i != (int) El::END; ++i)
      radii[i] = float(multiplier * Element(i).covalent_r() + tolerance / 2);
  }
  float get_radius(El el) const { return radii.empty() ? 0.f : radii[(int)el]; }
  void set_radius(El el, float r) {
    if (!radii.empty())
      radii[(int)el] = r;
  }

  template<typename Func>
  void for_each_contact(SubCells& sc, const Func& func);
};

template<typename Func>
void ContactSearch::for_each_contact(SubCells& sc, const Func& func) {
  if (!sc.model)
    fail("SubCells not initialized");
  for (int n_ch = 0; n_ch != (int) sc.model->chains.size(); ++n_ch) {
    Chain& chain = sc.model->chains[n_ch];
    PolymerType pt = PolymerType::Unknown;
    if (skip_adjacent_residue)
      pt = check_polymer_type(chain.get_polymer());
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        Atom& atom = res.atoms[n_atom];
        if (!sc.include_h && is_hydrogen(atom.element))
          continue;
        if (atom.occ < min_occupancy)
          continue;
        sc.for_each(atom.pos, atom.altloc, search_radius,
                    [&](SubCells::Mark& m, float dist_sq) {
            // do not consider connections inside a residue
            if (skip_intra_residue && m.image_idx == 0 &&
                m.chain_idx == n_ch && m.residue_idx == n_res)
              return;
            // additionally, we may have per-element distances
            if (!radii.empty()) {
              float d = radii[atom.element.ordinal()] + radii[(int)m.element];
              if (d < 0 || dist_sq > d * d)
                return;
            }
            // do not consider connections between adjacent residues
            if (skip_adjacent_residue && m.image_idx == 0 &&
                m.chain_idx == n_ch)
              if (are_connected(res, chain.residues[m.residue_idx], pt) ||
                  are_connected(chain.residues[m.residue_idx], res, pt))
                return;
            // avoid reporting connections twice (A-B and B-A)
            if (!twice)
              if (m.chain_idx < n_ch || (m.chain_idx == n_ch &&
                    (m.residue_idx < n_res || (m.residue_idx == n_res &&
                                               m.atom_idx < n_atom))))
                return;
            // atom can be linked with its image, but if the image
            // is too close the atom is likely on special position.
            if (m.chain_idx == n_ch && m.residue_idx == n_res &&
                m.atom_idx == n_atom && dist_sq < special_pos_cutoff_sq)
              return;
            CRA cra2 = m.to_cra(*sc.model);
            // ignore atoms with occupancy below the specified value
            if (cra2.atom->occ < min_occupancy)
              return;
            func(CRA{&chain, &res, &atom}, cra2, m.image_idx, dist_sq);
        });
      }
    }
  }
}

} // namespace gemmi
#endif
