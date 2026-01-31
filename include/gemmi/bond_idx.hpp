//! @file
//! @brief BondIndex for atom connectivity and graph distance.
//!
//! BondIndex: for checking which atoms are bonded, calculating graph distance.

// Copyright 2018 Global Phasing Ltd.
//
// BondIndex: for checking which atoms are bonded, calculating graph distance.

#ifndef GEMMI_BOND_IDX_HPP_
#define GEMMI_BOND_IDX_HPP_

#include <map>
#include "model.hpp"    // for Residue, Atom
#include "monlib.hpp"   // for MonLib

namespace gemmi {

//! @brief Index of bonds for checking connectivity and graph distance.
struct BondIndex {
  const Model& model;  //!< Reference to model

  //! @brief Atom reference (serial number and symmetry image flag).
  struct AtomImage {
    int atom_serial;  //!< Atom serial number
    bool same_image;  //!< True if same symmetry image

    //! @brief Equality comparison.
    //! @param o Other AtomImage
    //! @return True if both serial and same_image match
    bool operator==(const AtomImage& o) const {
      return atom_serial == o.atom_serial && same_image == o.same_image;
    }
  };

  std::map<int, std::vector<AtomImage>> index;  //!< Bond index map

  //! @brief Constructor.
  //! @param model_ Model reference
  //! @throws std::runtime_error if duplicate serial numbers found
  BondIndex(const Model& model_) : model(model_) {
    for (const_CRA cra : model.all())
      if (!index.emplace(cra.atom->serial, std::vector<AtomImage>()).second)
        fail("duplicated serial numbers");
  }

  //! @brief Add one-way bond link.
  //! @param a First atom
  //! @param b Second atom
  //! @param same_image True if same symmetry image
  void add_oneway_link(const Atom& a, const Atom& b, bool same_image) {
    std::vector<AtomImage>& list_a = index.at(a.serial);
    AtomImage ai{b.serial, same_image};
    if (!in_vector(ai, list_a))
      list_a.push_back(ai);
  }

  //! @brief Add bidirectional bond link.
  //! @param a First atom
  //! @param b Second atom
  //! @param same_image True if same symmetry image
  void add_link(const Atom& a, const Atom& b, bool same_image) {
    add_oneway_link(a, b, same_image);
    add_oneway_link(b, a, same_image);
  }

  //! @brief Add monomer bonds from monomer library.
  //! @param monlib Monomer library
  //!
  //! add_monomer_bonds() is not aware of modifications associated with links.
  //! Modifications that add bonds are rare, but to be more correct, use bonds
  //! from topology (Topo::bonds).
  void add_monomer_bonds(MonLib& monlib) {
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues) {
        std::string altlocs;
        add_distinct_altlocs(res, altlocs);
        if (altlocs.empty())
          altlocs += '*';
        auto monomer = monlib.monomers.find(res.name);
        if (monomer == monlib.monomers.end())
          fail("Monomer description not found: " + res.name);
        for (const Restraints::Bond& bond : monomer->second.rt.bonds)
          for (char alt : altlocs)
            if (const Atom* at1 = res.find_atom(bond.id1.atom, alt))
              if (const Atom* at2 = res.find_atom(bond.id2.atom, alt)) {
                add_link(*at1, *at2, true);
                if (!at1->altloc && !at2->altloc)
                  break;
              }
      }
  }

  //! @brief Check if two atoms are directly bonded.
  //! @param a First atom
  //! @param b Second atom
  //! @param same_image True if same symmetry image
  //! @return True if atoms are linked
  bool are_linked(const Atom& a, const Atom& b, bool same_image) const {
    return in_vector({b.serial, same_image}, index.at(a.serial));
  }

  //! @brief Calculate graph distance between atoms.
  //! @param a First atom
  //! @param b Second atom
  //! @param same_image True if same symmetry image
  //! @param max_distance Maximum distance to search
  //! @return Graph distance (number of bonds), or max_distance+1 if not connected
  int graph_distance(const Atom& a, const Atom& b, bool same_image,
                     int max_distance=4) const {
    std::vector<AtomImage> neighbors(1, {a.serial, true});
    for (int distance = 1; distance <= max_distance; ++distance) {
      for (size_t n = neighbors.size(); n--; ) {
        for (AtomImage ai : index.at(neighbors[n].atom_serial)) {
          if (!neighbors[n].same_image)
            ai.same_image = !ai.same_image;
          if (ai.atom_serial == b.serial && ai.same_image == same_image)
            return distance;
          if (!in_vector(ai, neighbors))
            neighbors.push_back(ai);
        }
      }
    }
    return max_distance + 1;
  }
};

} // namespace gemmi
#endif
