// Copyright 2018 Global Phasing Ltd.
//
// Representation of small molecule or inorganic crystal.
// Flat list of atom sites. Minimal functionality.

#ifndef GEMMI_SMODEL_HPP_
#define GEMMI_SMODEL_HPP_

#include <cctype>        // for isalpha
#include <algorithm>     // for any_of
#include <string>
#include <vector>
#include "elem.hpp"      // Element
#include "unitcell.hpp"  // UnitCell, Fractional

namespace gemmi {

struct AtomicStructure {
  struct Site {
    std::string label;
    std::string type_symbol;
    Fractional fract;
    double occ = 1.0;
    double u_iso = 0.;
    Element element = El::X;
    signed char charge = 0;  // [-8, +8]

    void fill_in_element_and_charge();
  };

  std::string name;
  UnitCell cell;
  std::string spacegroup_hm;
  std::vector<Site> sites;

  std::vector<Site> get_all_unit_cell_sites() const;
};

inline void AtomicStructure::Site::fill_in_element_and_charge() {
  const std::string& s = type_symbol.empty() ? label : type_symbol;
  int len = s.size() > 1 && std::isalpha(s[1]) ? 2 : 1;
  element = len == 1 ? impl::find_single_letter_element(s[0])
                     : find_element(s.c_str());
  if (element != El::X && s[len] >= '0' && s[len] <= '9')
    charge = (s[len] - '0') * (s[len+1] == '-' ? -1 : 1);
}

inline std::vector<AtomicStructure::Site>
AtomicStructure::get_all_unit_cell_sites() const {
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
