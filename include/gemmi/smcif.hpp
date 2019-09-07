// Copyright 2018 Global Phasing Ltd.
//
// Read small molecule CIF file into AtomicStructure.
// Minimal functionality.

#ifndef GEMMI_SMCIF_HPP_
#define GEMMI_SMCIF_HPP_

#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "math.hpp" // Fractional
#include "numb.hpp"  // for as_number
#include "symmetry.hpp" // SpaceGroup
#include "unitcell.hpp" // UnitCell

namespace gemmi {

struct AtomicStructure {
  struct Site {
    std::string label;
    std::string type_symbol;
    Fractional fract;
    float occ = 1.0f;
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

  enum { kLabel, kSymbol, kX, kY, kZ, kOcc };
  cif::Table atom_table = block.find("_atom_site_",
                                     {"label",
                                      "?type_symbol",
                                      "?fract_x",
                                      "?fract_y",
                                      "?fract_z",
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
    if (row.has(kOcc))
      site.occ = (float) as_number(row[kOcc], 1.0);
    st.sites.push_back(site);
  }
  const SpaceGroup* sg = find_spacegroup_by_name(st.spacegroup_hm);
  st.cell.set_cell_images_from_spacegroup(sg);
  return st;
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
