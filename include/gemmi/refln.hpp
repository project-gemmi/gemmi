// Copyright 2019 Global Phasing Ltd.
//
// Reads reflection data from the mmCIF format.

#ifndef GEMMI_REFLN_HPP_
#define GEMMI_REFLN_HPP_

#include <array>
#include "cifdoc.hpp"
#include "numb.hpp"      // for as_number
#include "unitcell.hpp"  // for UnitCell
#include "symmetry.hpp"  // for find_spacegroup_by_name

namespace gemmi {

namespace impl {
inline void set_cell_from_mmcif(cif::Block& block, UnitCell& cell) {
  cif::Table tab = block.find("_cell.",
                              {"length_a", "length_b", "length_c",
                               "angle_alpha", "angle_beta", "angle_gamma"});
  if (tab.ok()) {
    auto c = tab.one();
    using namespace gemmi::cif;
    if (!is_null(c[0]) && !is_null(c[1]) && !is_null(c[2]))
      cell.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
               as_number(c[3]), as_number(c[4]), as_number(c[5]));
  }
}
} // namespace impl


struct ReflnBlock {
  cif::Block block;
  std::string entry_id;
  UnitCell cell;
  const SpaceGroup* spacegroup = nullptr;
  double wavelength;
  cif::Loop* refln_loop = nullptr;

  ReflnBlock(cif::Block&& block_) : block(block_) {
    entry_id = cif::as_string(block.find_value("_entry.id"));
    impl::set_cell_from_mmcif(block, cell);
    const char* hm_tag = "_symmetry.space_group_name_H-M";
    if (const std::string* hm = block.find_value(hm_tag))
      spacegroup = find_spacegroup_by_name(cif::as_string(*hm));
    const char* wave_tag = "_diffrn_radiation_wavelength.wavelength";
    cif::Column wave_col = block.find_values(wave_tag);
    wavelength = wave_col.length() == 1 ? cif::as_number(wave_col[0]) : 0.;
    refln_loop = block.find_loop("_refln.index_h").get_loop();
  }

  size_t get_column_index(const std::string& tag) const {
    for (size_t i = 0; i != refln_loop->tags.size(); ++i)
      if (refln_loop->tags[i].compare(7, std::string::npos, tag))
        return i;
    throw std::runtime_error("Column not found: " + tag);
  }

  template <typename T>
  std::vector<T> refln_column(const std::string& tag, T null) {
    size_t i = get_column_index(tag);
    std::vector<T> v(refln_loop->length());
    for (size_t j = 0; j != v.size(); i += refln_loop->width(), ++j)
      v[j] = cif::as_any(refln_loop->values[i], null);
    return v;
  }

  std::vector<std::array<int,3>> refln_indices() {
    size_t h = get_column_index("index_h");
    size_t k = get_column_index("index_k");
    size_t l = get_column_index("index_l");
    std::vector<std::array<int,3>> v(refln_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += refln_loop->width()) {
      v[j][0] = cif::as_int(refln_loop->values[n + h]);
      v[j][1] = cif::as_int(refln_loop->values[n + k]);
      v[j][2] = cif::as_int(refln_loop->values[n + l]);
    }
    return v;
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
