// Copyright 2020 Global Phasing Ltd.
//
// Addends to scattering form factors used in DensityCalculator
// and StructureFactorCalculator.

#ifndef GEMMI_ADDENDS_HPP_
#define GEMMI_ADDENDS_HPP_

#include <array>
#include "elem.hpp"  // for El, Element

namespace gemmi {

/// @brief Container for anomalous scattering correction addends
/// Stores addend values for each element used in density and structure factor calculations.
struct Addends {
  std::array<float, (int)El::END> values = {};

  /// @brief Set the addend value for a given element
  /// @param el the chemical element
  /// @param val the addend value to set
  void set(Element el, float val) { values[el.ordinal()] = val; }

  /// @brief Get the addend value for a given element
  /// @param el the chemical element
  /// @return the addend value for the element
  float get(Element el) const { return values[el.ordinal()]; }

  /// @brief Get the total number of elements in the array
  /// @return the size of the addends array
  size_t size() const { return values.size(); }

  /// @brief Clear all addend values to zero
  void clear() {
    for (size_t i = 0; i != size(); ++i)
      values[i] = 0.;
  }

  /// @brief Subtract atomic number Z from each element's addend value
  /// Optionally preserves hydrogen and deuterium values.
  /// @param except_hydrogen if true, skip subtracting from hydrogen and deuterium
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
