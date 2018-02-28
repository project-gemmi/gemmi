// Copyright 2018 Global Phasing Ltd.
//
// Modify the model.

#ifndef GEMMI_MODIFY_HPP_
#define GEMMI_MODIFY_HPP_

#include <algorithm>  // for std::remove_if
#include "model.hpp"

namespace gemmi {

template<class T> void remove_hydrogens(T& obj) {
  for (auto& child : obj.children())
    remove_hydrogens(child);
}
template<> inline void remove_hydrogens(Residue& res) {
  res.atoms.erase(std::remove_if(res.atoms.begin(), res.atoms.end(),
                  [](const Atom& a) { return a.element == El::H ||
                                             a.element == El::D; }),
                  res.atoms.end());
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
