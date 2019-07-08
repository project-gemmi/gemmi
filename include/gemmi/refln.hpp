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

  bool ok() const { return default_loop != nullptr; }
  void check_ok() const { if (!ok()) fail("Invalid ReflnBlock"); }

  // position after "_refln." or "_diffrn_refln".
  int tag_offset() const { return refln_loop ? 7 : 13; }

  void use_unmerged(bool unmerged) {
    default_loop = unmerged ? diffrn_refln_loop : refln_loop;
  }
  bool is_unmerged() const { return ok() && default_loop == diffrn_refln_loop; }

  std::vector<std::string> column_labels() const {
    check_ok();
    std::vector<std::string> labels(default_loop->tags.size());
    for (size_t i = 0; i != labels.size(); ++i)
      labels[i].assign(default_loop->tags[i], tag_offset(), std::string::npos);
    return labels;
  }

  int find_column_index(const std::string& tag) const {
    if (!ok())
      return -1;
    int name_pos = tag_offset();
    for (int i = 0; i != (int) default_loop->tags.size(); ++i)
      if (default_loop->tags[i].compare(name_pos, std::string::npos, tag) == 0)
        return i;
    return -1;
  }

  size_t get_column_index(const std::string& tag) const {
    int idx = find_column_index(tag);
    if (idx == -1)
      throw std::runtime_error("Column not found: " + tag);
    return idx;
  }

  template<typename T>
  std::vector<T> make_vector(const std::string& tag, T null) const {
    size_t n = get_column_index(tag);
    std::vector<T> v(default_loop->length());
    for (size_t j = 0; j != v.size(); n += default_loop->width(), ++j)
      v[j] = cif::as_any(default_loop->values[n], null);
    return v;
  }

  std::array<size_t,3> get_hkl_column_indices() const {
    return {{get_column_index("index_h"),
             get_column_index("index_k"),
             get_column_index("index_l")}};
  }

  std::vector<std::array<int,3>> make_index_vector() const {
    auto hkl_idx = get_hkl_column_indices();
    std::vector<std::array<int,3>> v(default_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += default_loop->width())
      for (int i = 0; i != 3; ++i)
        v[j][i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
    return v;
  }

  std::vector<double> make_1_d2_vector() const {
    if (!cell.is_crystal() || cell.a <= 0)
      fail("Unit cell is not known");
    auto hkl_idx = get_hkl_column_indices();
    std::vector<double> r(default_loop->length());
    for (size_t j = 0, n = 0; j != r.size(); j++, n += default_loop->width()) {
      int h = cif::as_int(default_loop->values[n + hkl_idx[0]]);
      int k = cif::as_int(default_loop->values[n + hkl_idx[1]]);
      int l = cif::as_int(default_loop->values[n + hkl_idx[2]]);
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
  // Some blocks miss space group tag, try to fill it in.
  const gemmi::SpaceGroup* first_sg = nullptr;
  for (gemmi::ReflnBlock& rblock : r)
    if (!first_sg)
      first_sg = rblock.spacegroup;
    else if (!rblock.spacegroup)
      rblock.spacegroup = first_sg;
  return r;
}

// Abstraction of data source, cf. MtzDataProxy.
struct ReflnDataProxy {
  const ReflnBlock& rb_;
  const cif::Loop& loop() const { rb_.check_ok(); return *rb_.default_loop; }
  bool ok() const { return rb_.ok(); }
  std::array<size_t,3> hkl_col() const { return rb_.get_hkl_column_indices(); }
  size_t stride() const { return loop().tags.size(); }
  size_t size() const { return loop().values.size(); }
  int get_int(size_t n) const { return cif::as_int(loop().values[n]); }
  double get_num(size_t n) const { return cif::as_number(loop().values[n]); }
  const UnitCell& unit_cell() const { return rb_.cell; }
  const SpaceGroup* spacegroup() const { return rb_.spacegroup; }
};

} // namespace gemmi
#endif
