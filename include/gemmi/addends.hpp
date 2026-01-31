//! @file
//! @brief Addends to scattering form factors.
//!
//! Addends to scattering form factors used in DensityCalculator
//! and StructureFactorCalculator.

// Copyright 2020 Global Phasing Ltd.
//
// Addends to scattering form factors used in DensityCalculator
// and StructureFactorCalculator.

#ifndef GEMMI_ADDENDS_HPP_
#define GEMMI_ADDENDS_HPP_

#include <array>
#include "elem.hpp"  // for El, Element

namespace gemmi {

//! @brief Addends to scattering form factors per element.
struct Addends {
  std::array<float, (int)El::END> values = {};  //!< Addend values per element

  //! @brief Set addend value for element.
  //! @param el Element
  //! @param val Addend value
  void set(Element el, float val) { values[el.ordinal()] = val; }

  //! @brief Get addend value for element.
  //! @param el Element
  //! @return Addend value
  float get(Element el) const { return values[el.ordinal()]; }

  //! @brief Get number of elements.
  //! @return Size of values array
  size_t size() const { return values.size(); }

  //! @brief Clear all addend values to zero.
  void clear() {
    for (size_t i = 0; i != size(); ++i)
      values[i] = 0.;
  }

  //! @brief Subtract atomic number from each element's addend.
  //! @param except_hydrogen If true, don't subtract for hydrogen/deuterium
  void subtract_z(bool except_hydrogen=false) {
    for (int z = 2; z < (int)El::D; ++z)
      values[z] -= z;
    if (!except_hydrogen) {
      values[(int)El::H] -= 1;
      values[(int)El::D] -= 1;
    }
  }
};

} // namespace gemmi
#endif
