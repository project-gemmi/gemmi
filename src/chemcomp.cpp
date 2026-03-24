// Copyright 2026 Global Phasing Ltd.

#include <gemmi/chemcomp.hpp>
#include <gemmi/model.hpp>

namespace gemmi {

Atom* Restraints::AtomId::get_from(Residue& res1, Residue* res2,
                                   char alt, char altloc2) const {
  Residue* residue = &res1;
  if (comp == 2 && res2 != nullptr) {
    residue = res2;
    if (altloc2 != '\0')
      alt = altloc2;
  }
  Atom* a = residue->find_atom(atom, alt, El::X, false);
  // Special case: microheterogeneity may have shared atoms only in
  // the first residue. Example: in 1ejg N is shared between PRO and SER.
  if (a == nullptr && alt != '\0' && residue->group_idx > 0)
    a = (residue - residue->group_idx)->find_atom(atom, alt, El::X, false);
  return a;
}

const Atom* Restraints::AtomId::get_from(const Residue& res1,
                                         const Residue* res2,
                                         char alt, char alt2) const {
  return get_from(const_cast<Residue&>(res1), const_cast<Residue*>(res2),
                  alt, alt2);
}

} // namespace gemmi
