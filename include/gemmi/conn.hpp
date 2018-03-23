// Copyright 2018 Global Phasing Ltd.
//
// Connection search.

#ifndef GEMMI_CONN_HPP_
#define GEMMI_CONN_HPP_

#include "unitcell.hpp"
#include "model.hpp"

namespace gemmi {

// Number of SG atoms is relatively small; checking all pairs should be fast.
std::vector<Connection> find_disulfide_bonds(const Model& model,
                                             const UnitCell& cell) {
  const double max_dist = 3.0;
  // Find all SG sulfur atoms.
  std::vector<const_CRA> atoms;
  const std::string sg = "SG";
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        if (atom.element == El::S && atom.name == sg)
          atoms.push_back(const_CRA{&chain, &res, &atom});
  // Check distances.
  std::vector<Connection> ret;
  for (size_t i = 0; i < atoms.size(); ++i) {
    const Atom* a1 = atoms[i].atom;
    for (size_t j = 0; j < i; ++j) {
      const Atom* a2 = atoms[j].atom;
      if (a1->same_conformer(*a2)) {
        NearbyImage im = cell.find_nearest_image(a1->pos, a2->pos,
                                                 SymmetryImage::Unspecified);
        if (im.dist_sq < max_dist * max_dist) {
          Connection c;
          c.name = "disulf" + std::to_string(ret.size() + 1);
          c.type = Connection::Disulf;
          c.image = im.same_image() ? SymmetryImage::Same
                                    : SymmetryImage::Different;
          c.atom[0] = AtomAddress(*atoms[i].chain, *atoms[i].residue, *a1);
          c.atom[1] = AtomAddress(*atoms[j].chain, *atoms[j].residue, *a2);
          ret.push_back(c);
        }
      }
    }
  }
  return ret;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
