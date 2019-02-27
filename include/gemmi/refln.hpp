// Copyright 2019 Global Phasing Ltd.
//
// Reads reflection data from the mmCIF format.

#ifndef GEMMI_REFLN_HPP_
#define GEMMI_REFLN_HPP_

#include <array>
#include "cifdoc.hpp"
#include "numb.hpp"       // for as_number
#include "unitcell.hpp"   // for UnitCell
#include "symmetry.hpp"   // for find_spacegroup_by_name
#include "mmcif_impl.hpp" // for set_cell_from_mmcif

namespace gemmi {

struct ReflnBlock {
  cif::Block block;
  std::string entry_id;
  UnitCell cell;
  const SpaceGroup* spacegroup = nullptr;
  double wavelength;
  cif::Loop* refln_loop = nullptr;
  cif::Loop* diffrn_refln_loop = nullptr;
  cif::Loop* default_loop = nullptr;

  ReflnBlock(cif::Block&& block_) : block(block_) {
    entry_id = cif::as_string(block.find_value("_entry.id"));
    impl::set_cell_from_mmcif(block, cell);
    if (const std::string* hm = impl::find_spacegroup_hm_value(block))
      spacegroup = find_spacegroup_by_name(cif::as_string(*hm));
    const char* wave_tag = "_diffrn_radiation_wavelength.wavelength";
    cif::Column wave_col = block.find_values(wave_tag);
    wavelength = wave_col.length() == 1 ? cif::as_number(wave_col[0]) : 0.;
    refln_loop = block.find_loop("_refln.index_h").get_loop();
    diffrn_refln_loop = block.find_loop("_diffrn_refln.index_h").get_loop();
    default_loop = refln_loop ? refln_loop : diffrn_refln_loop;
  }

  size_t get_column_index(const std::string& tag) const {
    int name_pos = refln_loop ? 7 : 13;
    for (size_t i = 0; i != default_loop->tags.size(); ++i)
      if (default_loop->tags[i].compare(name_pos, std::string::npos, tag) == 0)
        return i;
    throw std::runtime_error("Column not found: " + tag);
  }

  template <typename T>
  std::vector<T> make_vector(const std::string& tag, T null) {
    size_t n = get_column_index(tag);
    std::vector<T> v(default_loop->length());
    for (size_t j = 0; j != v.size(); n += default_loop->width(), ++j)
      v[j] = cif::as_any(default_loop->values[n], null);
    return v;
  }

  std::vector<std::array<int,3>> make_index_vector() const {
    size_t h_idx = get_column_index("index_h");
    size_t k_idx = get_column_index("index_k");
    size_t l_idx = get_column_index("index_l");
    std::vector<std::array<int,3>> v(default_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += default_loop->width()) {
      v[j][0] = cif::as_int(default_loop->values[n + h_idx]);
      v[j][1] = cif::as_int(default_loop->values[n + k_idx]);
      v[j][2] = cif::as_int(default_loop->values[n + l_idx]);
    }
    return v;
  }

  std::vector<double> make_1_d2_vector() const {
    if (!cell.is_crystal() || cell.a <= 0)
      fail("Unit cell is not known");
    size_t h_idx = get_column_index("index_h");
    size_t k_idx = get_column_index("index_k");
    size_t l_idx = get_column_index("index_l");
    std::vector<double> r(default_loop->length());
    for (size_t j = 0, n = 0; j != r.size(); j++, n += default_loop->width()) {
      int h = cif::as_int(default_loop->values[n + h_idx]);
      int k = cif::as_int(default_loop->values[n + k_idx]);
      int l = cif::as_int(default_loop->values[n + l_idx]);
      r[j] = cell.calculate_1_d2(h, k, l);
    }
    return r;
  }
};

// moves blocks from the argument to the return value
inline
std::vector<ReflnBlock> as_refln_blocks(std::vector<cif::Block>&& blocks) {
  std::vector<ReflnBlock> r;
  r.reserve(blocks.size());
  for (cif::Block& block : blocks)
    r.emplace_back(std::move(block));
  blocks.clear();
  return r;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
