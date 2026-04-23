// Copyright 2018 Global Phasing Ltd.
//
// BondIndex: for checking which atoms are bonded, calculating graph distance.

#ifndef GEMMI_BOND_IDX_HPP_
#define GEMMI_BOND_IDX_HPP_

#include <map>
#include "model.hpp"    // for Residue, Atom
#include "monlib.hpp"   // for MonLib

namespace gemmi {

/// @brief Index for efficient bond topology queries in a crystal structure
/// Enables checking atom connectivity and calculating graph distances, including
/// handling of atoms in different unit cell images.
struct BondIndex {
  const Model& model;

  /// @brief Represents an atom and whether it's in the same unit cell image
  struct AtomImage {
    int atom_serial;     ///< Serial number of the atom
    bool same_image;     ///< True if atom is in the same unit cell image as reference
    /// @brief Equality comparison
    /// @param o the other AtomImage to compare
    /// @return true if both serial and image flag match
    bool operator==(const AtomImage& o) const {
      return atom_serial == o.atom_serial && same_image == o.same_image;
    }
  };
  std::map<int, std::vector<AtomImage>> index;

  /// @brief Construct a BondIndex for the given model
  /// @details Initializes the index with all atoms from the model.
  ///          Fails if duplicate atom serial numbers are found.
  /// @param model_ the crystallographic model to index
  BondIndex(const Model& model_) : model(model_) {
    for (const_CRA cra : model.all())
      if (!index.emplace(cra.atom->serial, std::vector<AtomImage>()).second)
        fail("duplicated serial numbers");
  }

  /// @brief Add a unidirectional bond link between two atoms
  /// @details Does not add the reverse link (a->b without b->a).
  ///          Does not add duplicate links.
  /// @param a the first atom
  /// @param b the second atom
  /// @param same_image whether both atoms are in the same unit cell image
  void add_oneway_link(const Atom& a, const Atom& b, bool same_image) {
    std::vector<AtomImage>& list_a = index.at(a.serial);
    AtomImage ai{b.serial, same_image};
    if (!in_vector(ai, list_a))
      list_a.push_back(ai);
  }

  /// @brief Add a bidirectional bond link between two atoms
  /// @param a the first atom
  /// @param b the second atom
  /// @param same_image whether both atoms are in the same unit cell image
  void add_link(const Atom& a, const Atom& b, bool same_image) {
    add_oneway_link(a, b, same_image);
    add_oneway_link(b, a, same_image);
  }

  /// @brief Add bonds from monomer library restraints to the index
  /// @details Populates the bond index with standard bonds defined for each
  ///          residue type in the monomer library. Does not handle custom
  ///          bond modifications; for more accurate results, use bonds from
  ///          topology (Topo::bonds) which accounts for modifications.
  /// @param monlib the monomer library containing bond definitions
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

  /// @brief Check if two atoms are directly bonded
  /// @param a the first atom
  /// @param b the second atom
  /// @param same_image whether both atoms should be in the same unit cell image
  /// @return true if a direct bond exists between the atoms
  bool are_linked(const Atom& a, const Atom& b, bool same_image) const {
    return in_vector({b.serial, same_image}, index.at(a.serial));
  }

  /// @brief Calculate the minimum graph distance between two atoms
  /// @details Uses breadth-first search to find the shortest path through bonds.
  ///          Automatically handles transitions between unit cell images.
  /// @param a the starting atom
  /// @param b the target atom
  /// @param same_image whether both atoms should be in the same unit cell image
  /// @param max_distance maximum distance to search (default 4)
  /// @return the graph distance in bonds, or (max_distance + 1) if no path exists
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
