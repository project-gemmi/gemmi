// Copyright 2020 Global Phasing Ltd.
//
// Contact search, based on NeighborSearch from neighbor.hpp.

#ifndef GEMMI_CONTACT_HPP_
#define GEMMI_CONTACT_HPP_

#include "model.hpp"
#include "neighbor.hpp"
#include "polyheur.hpp"  // for check_polymer_type, are_connected

namespace gemmi {

//! Find contacts (close approaches) between atoms.
//!
//! Uses NeighborSearch for efficient spatial queries with configurable
//! filtering to ignore intra-residue, adjacent, or same-chain contacts.
struct ContactSearch {
  //! Filtering options for contact search.
  enum class Ignore {
    Nothing=0,           //!< Report all contacts
    SameResidue,         //!< Ignore contacts within same residue
    AdjacentResidues,    //!< Ignore contacts in adjacent residues
    SameChain,           //!< Ignore contacts within same chain
    SameAsu              //!< Ignore contacts within same asymmetric unit
  };

  // Configuration parameters
  double search_radius;                    //!< Maximum search distance
  Ignore ignore = Ignore::SameResidue;     //!< Contact filtering mode
  bool twice = false;                      //!< Report both A-B and B-A contacts
  float min_occupancy = 0.f;               //!< Minimum occupancy threshold
  double special_pos_cutoff_sq = 0.8 * 0.8; //!< Cutoff for special positions
  std::vector<float> radii;                //!< Per-element distance cutoffs

  //! @brief Construct contact search with specified radius.
  //! @param radius Maximum search distance
  ContactSearch(double radius) noexcept : search_radius(radius) {}

  //! @brief Set per-element radii based on covalent radii.
  //! @param multiplier Multiplier for covalent radius
  //! @param tolerance Additional tolerance to add
  //!
  //! Initializes element-specific distance cutoffs using covalent radii.
  void setup_atomic_radii(double multiplier, double tolerance) {
    radii.resize((size_t)El::END);
    for (int i = 0; i != (int) El::END; ++i)
      radii[i] = float(multiplier * Element(i).covalent_r() + tolerance / 2);
  }

  //! @brief Get radius for element.
  //! @param el Element
  //! @return Radius or 0 if not set
  float get_radius(El el) const { return radii.empty() ? 0.f : radii[(int)el]; }

  //! @brief Set radius for element.
  //! @param el Element
  //! @param r Radius value
  void set_radius(El el, float r) {
    if (!radii.empty())
      radii[(int)el] = r;
  }

  //! @brief Apply function to each contact found.
  //! @tparam Func Callable type
  //! @param ns Populated NeighborSearch object
  //! @param func Function called for each contact
  template<typename Func>
  void for_each_contact(NeighborSearch& ns, const Func& func);

  //! Contact result containing both partners and distance.
  struct Result {
    CRA partner1;      //!< First contact partner
    CRA partner2;      //!< Second contact partner
    int image_idx;     //!< Symmetry image index
    double dist_sq;    //!< Squared distance
  };

  //! @brief Find all contacts and return as vector.
  //! @param ns Populated NeighborSearch object
  //! @return Vector of contact results
  std::vector<Result> find_contacts(NeighborSearch& ns) {
    std::vector<Result> out;
    for_each_contact(ns, [&out](const CRA& cra1, const CRA& cra2,
                                int image_idx, double dist_sq) {
        out.push_back({cra1, cra2, image_idx, dist_sq});
    });
    return out;
  }
};

template<typename Func>
void ContactSearch::for_each_contact(NeighborSearch& ns, const Func& func) {
  if (!ns.model)
    fail(ns.small_structure ? "ContactSearch does not work with SmallStructure"
                            : "NeighborSearch not initialized");
  for (int n_ch = 0; n_ch != (int) ns.model->chains.size(); ++n_ch) {
    Chain& chain = ns.model->chains[n_ch];
    PolymerType pt = PolymerType::Unknown;
    if (ignore == Ignore::AdjacentResidues)
      pt = check_polymer_type(chain.get_polymer());
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        Atom& atom = res.atoms[n_atom];
        if (!ns.include_h && is_hydrogen(atom.element))
          continue;
        if (atom.occ < min_occupancy)
          continue;
        ns.for_each(atom.pos, atom.altloc, search_radius,
                    [&](NeighborSearch::Mark& m, double dist_sq) {
            // do not consider connections inside a residue
            if (ignore != Ignore::Nothing && m.image_idx == 0 &&
                m.chain_idx == n_ch && m.residue_idx == n_res)
              return;
            switch (ignore) {
              case Ignore::Nothing:
                break;
              case Ignore::SameResidue:
                if (m.image_idx == 0 && m.chain_idx == n_ch)
                  if (m.residue_idx == n_res)
                    return;
                break;
              case Ignore::AdjacentResidues:
                if (m.image_idx == 0 && m.chain_idx == n_ch)
                  if (m.residue_idx == n_res ||
                      are_connected(res, chain.residues[m.residue_idx], pt) ||
                      are_connected(chain.residues[m.residue_idx], res, pt))
                    return;
                break;
              case Ignore::SameChain:
                if (m.image_idx == 0 && m.chain_idx == n_ch)
                  return;
                break;
              case Ignore::SameAsu:
                if (m.image_idx == 0)
                  return;
                break;
            }
            // additionally, we may have per-element distances
            if (!radii.empty()) {
              double d = radii[atom.element.ordinal()] + radii[m.element.ordinal()];
              if (d < 0 || dist_sq > d * d)
                return;
            }
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
            CRA cra2 = m.to_cra(*ns.model);
            // ignore atoms with occupancy below the specified value
            if (cra2.atom->occ < min_occupancy)
              return;
            func(CRA{&chain, &res, &atom}, cra2, m.image_idx, dist_sq);
        }, ns.sufficient_k(search_radius));
      }
    }
  }
}

} // namespace gemmi
#endif
