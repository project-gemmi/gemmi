// Copyright 2018 Global Phasing Ltd.
//
// Read small molecule CIF file into AtomicStructure.
// Minimal functionality.

#ifndef GEMMI_SMCIF_HPP_
#define GEMMI_SMCIF_HPP_

#include <algorithm>     // for any_of
#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "numb.hpp"      // for as_number
#include "symmetry.hpp"  // SpaceGroup
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

inline
AtomicStructure make_atomic_structure_from_block(const cif::Block& block_) {
  using cif::as_number;
  using cif::as_string;
  cif::Block& block = const_cast<cif::Block&>(block_);
  AtomicStructure st;
  st.name = block.name;

  // unit cell and symmetry
  cif::Table cell = block.find("_cell_",
                               {"length_a", "length_b", "length_c",
                                "angle_alpha", "angle_beta", "angle_gamma"});
  if (cell.ok()) {
    auto c = cell.one();
    if (!cif::is_null(c[0]) && !cif::is_null(c[1]) && !cif::is_null(c[2]))
      st.cell.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
                  as_number(c[3]), as_number(c[4]), as_number(c[5]));
  }
  st.spacegroup_hm =
                as_string(block.find_value("_symmetry_space_group_name_H-M"));

  enum { kLabel, kSymbol, kX, kY, kZ, kUiso, kOcc };
  cif::Table atom_table = block.find("_atom_site_",
                                     {"label",
                                      "?type_symbol",
                                      "?fract_x",
                                      "?fract_y",
                                      "?fract_z",
                                      "?U_iso_or_equiv",
                                      "?occupancy"});
  for (auto row : atom_table) {
    AtomicStructure::Site site;
    site.label = as_string(row[kLabel]);
    if (row.has(kSymbol))
      site.type_symbol = as_string(row[kSymbol]);
    if (row.has(kX))
      site.fract.x = as_number(row[kX]);
    if (row.has(kY))
      site.fract.y = as_number(row[kY]);
    if (row.has(kZ))
      site.fract.z = as_number(row[kZ]);
    if (row.has(kUiso))
      site.occ = as_number(row[kUiso], 0.0);
    if (row.has(kOcc))
      site.occ = as_number(row[kOcc], 1.0);
    site.fill_in_element_and_charge();
    st.sites.push_back(site);
  }
  const SpaceGroup* sg = find_spacegroup_by_name(st.spacegroup_hm);
  st.cell.set_cell_images_from_spacegroup(sg);
  return st;
}

inline void AtomicStructure::Site::fill_in_element_and_charge() {
  const std::string& s = type_symbol.empty() ? label : type_symbol;
  int len = s.size() > 1 && std::isalpha(s[1]) ? 2 : 1;
  element = len == 1 ? impl::find_single_letter_element(s[0])
                     : find_element(s.c_str());
  if (element != El::X && std::isdigit(s[len]))
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
