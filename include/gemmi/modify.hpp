// Copyright 2018 Global Phasing Ltd.
//
// Modify the model.

#ifndef GEMMI_MODIFY_HPP_
#define GEMMI_MODIFY_HPP_

#include <algorithm>  // for std::remove_if
#include "model.hpp"
#include "resinfo.hpp"

namespace gemmi {

namespace impl {
  template<class Item, class UnaryPredicate>
  void remove(std::vector<Item>& v, UnaryPredicate p) {
    v.erase(std::remove_if(v.begin(), v.end(), p), v.end());
  }
} // namespace impl

// Remove hydrogens.
template<class T> void remove_hydrogens(T& obj) {
  for (auto& child : obj.children())
    remove_hydrogens(child);
}
template<> inline void remove_hydrogens(Residue& res) {
  impl::remove(res.atoms, [](const Atom& a) {
    return a.element == El::H || a.element == El::D;
  });
}

// Remove waters. It may leave empty chains.
template<class T> void remove_waters(T& obj) {
  for (auto& child : obj.children())
    remove_waters(child);
}
template<> inline void remove_waters(Chain& ch) {
  impl::remove(ch.residues, [](const Residue& res) {
      return find_tabulated_residue(res.name).is_water();
  });
}

// Remove ligands and waters. It may leave empty chains.
inline void remove_ligands_and_waters(Chain& ch,
                                      EntityType etype=EntityType::Unknown) {
  if (etype == EntityType::NonPolymer || etype == EntityType::Water) {
    ch.residues.clear();
  } else if (etype == EntityType::Unknown) {
    impl::remove(ch.residues, [](const Residue& res) {
        ResidueInfo info = find_tabulated_residue(res.name);
        return !info.is_nucleic() && !info.is_amino();
    });
  }
}
inline void remove_ligands_and_waters(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains) {
      EntityType etype = EntityType::Unknown;
      if (const Entity* ent = st.get_entity_of(chain))
        etype = ent->entity_type;
      remove_ligands_and_waters(chain, etype);
    }
}

// Remove empty chains.
inline void remove_empty_chains(Model& m) {
  m.chains.erase(std::remove_if(m.chains.begin(), m.chains.end(),
                   [](const Chain& chain) { return chain.residues.empty(); }),
                 m.chains.end());
}
inline void remove_empty_chains(Structure& st) {
  for (Model& model : st.models)
    remove_empty_chains(model);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
