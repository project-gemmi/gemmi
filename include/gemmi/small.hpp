// Copyright 2018 Global Phasing Ltd.
//
// Representation of small molecule or inorganic crystal.
// Flat list of atom sites. Minimal functionality.

#ifndef GEMMI_SMALL_HPP_
#define GEMMI_SMALL_HPP_

#include <cctype>        // for isalpha
#include <algorithm>     // for any_of
#include <bitset>
#include <string>
#include <vector>
#include "elem.hpp"      // Element
#include "unitcell.hpp"  // UnitCell, Fractional

namespace gemmi {

struct SmallStructure {
  struct Site {
    std::string label;
    std::string type_symbol;
    Fractional fract;
    double occ = 1.0;
    double u_iso = 0.;
    Element element = El::X;
    signed char charge = 0;  // [-8, +8]
  };

  struct AtomType {
    std::string symbol;
    Element element = El::X;
    signed char charge = 0;  // [-8, +8]
    double dispersion_real;
    double dispersion_imag;
  };

  std::string name;
  UnitCell cell;
  std::string spacegroup_hm;
  std::vector<Site> sites;
  std::vector<AtomType> atom_types;
  double wavelength = 0.; // the first wavelength if multiple

  std::vector<Site> get_all_unit_cell_sites() const;

  const AtomType* get_atom_type(const std::string& symbol) const {
    for (const AtomType& at : atom_types)
      if (at.symbol == symbol)
        return &at;
    return nullptr;
  }

  // similar to Model::present_elements() from model.hpp
  std::bitset<(size_t)El::END> present_elements() const {
    std::bitset<(size_t)El::END> table;
    for (const Site& atom : sites)
      table.set((size_t)atom.element.elem);
    return table;
  }
};

template<typename T>
inline void split_element_and_charge(const std::string& label, T* dest) {
  int len = label.size() > 1 && std::isalpha(label[1]) ? 2 : 1;
  dest->element = len == 1 ? impl::find_single_letter_element(label[0] & ~0x20)
                           : find_element(label.c_str());
  if (dest->element != El::X && label[len] >= '0' && label[len] <= '9')
    dest->charge = label[len] - '0';
  if (label[len+1] == '-')
    dest->charge = -dest->charge;
}

inline std::vector<SmallStructure::Site>
SmallStructure::get_all_unit_cell_sites() const {
  std::vector<Site> all;
  for (const Site& site : sites) {
    size_t start = all.size();
    all.push_back(site);
    for (const FTransform& image : cell.images) {
      Fractional fpos = image.apply(site.fract);
      if (std::any_of(all.begin() + start, all.end(), [&](const Site& other) {
            return cell.distance_sq(fpos, other.fract) < 0.5 * 0.5;
          }))
        continue;
      all.push_back(site);
      all.back().fract = fpos;
    }
  }
  return all;
}

} // namespace gemmi
#endif
