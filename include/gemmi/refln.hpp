//! @file
//! @brief Reads reflection data from the mmCIF format.
//!
//! Provides classes for parsing and accessing reflection data from SF-mmCIF files,
//! including both merged (_refln) and unmerged (_diffrn_refln) data.

// Copyright 2019 Global Phasing Ltd.

#ifndef GEMMI_REFLN_HPP_
#define GEMMI_REFLN_HPP_

#include <array>
#include "cifdoc.hpp"
#include "fail.hpp"       // for fail
#include "mmcif_impl.hpp" // for set_cell_from_mmcif, read_spacegroup_from_block
#include "numb.hpp"       // for as_number
#include "symmetry.hpp"   // for SpaceGroup
#include "unitcell.hpp"   // for UnitCell

namespace gemmi {

//! @brief Wrapper for mmCIF block containing reflection data.
//!
//! Provides access to both merged (_refln) and unmerged (_diffrn_refln)
//! reflection data from SF-mmCIF files. Stores crystallographic metadata
//! (unit cell, space group, wavelength) and provides methods to extract
//! reflection indices and data columns.
struct ReflnBlock {
  cif::Block block;  //!< Underlying CIF data block
  std::string entry_id;  //!< Entry identifier from _entry.id
  UnitCell cell;  //!< Unit cell parameters
  const SpaceGroup* spacegroup = nullptr;  //!< Space group symmetry
  double wavelength;  //!< X-ray wavelength (single wavelength case)
  int wavelength_count;  //!< Number of wavelengths in the dataset
  cif::Loop* refln_loop = nullptr;  //!< Pointer to merged reflection loop (_refln)
  cif::Loop* diffrn_refln_loop = nullptr;  //!< Pointer to unmerged reflection loop (_diffrn_refln)
  cif::Loop* default_loop = nullptr;  //!< Active loop (merged or unmerged)

  ReflnBlock() = default;
  ReflnBlock(ReflnBlock&& rblock_) = default;

  //! @brief Construct from CIF block, extracting reflection metadata.
  //! @param block_ CIF block (moved, not copied)
  //!
  //! Extracts entry ID, unit cell, space group, wavelength, and locates
  //! reflection data loops (_refln and/or _diffrn_refln).
  ReflnBlock(cif::Block&& block_) : block(std::move(block_)) {
    entry_id = cif::as_string(block.find_value("_entry.id"));
    impl::set_cell_from_mmcif(block, cell);
    if (const std::string* hm = impl::find_spacegroup_hm_value(block))
      spacegroup = find_spacegroup_by_name(cif::as_string(*hm),
                                           cell.alpha, cell.gamma);
    cell.set_cell_images_from_spacegroup(spacegroup);
    const char* wave_tag = "_diffrn_radiation_wavelength.wavelength";
    cif::Column wave_col = block.find_values(wave_tag);
    wavelength_count = wave_col.length();
    wavelength = wavelength_count == 1 ? cif::as_number(wave_col[0]) : 0.;
    refln_loop = block.find_loop("_refln.index_h").get_loop();
    diffrn_refln_loop = block.find_loop("_diffrn_refln.index_h").get_loop();
    default_loop = refln_loop ? refln_loop : diffrn_refln_loop;
  }
  ReflnBlock& operator=(ReflnBlock&&) = default;
  ReflnBlock& operator=(const ReflnBlock& o) {
    if (this == &o)
      return *this;
    block = o.block;
    entry_id = o.entry_id;
    cell = o.cell;
    spacegroup = o.spacegroup;
    wavelength = o.wavelength;
    wavelength_count = o.wavelength_count;
    if (o.refln_loop)
      refln_loop = block.find_loop("_refln.index_h").get_loop();
    if (o.diffrn_refln_loop)
      diffrn_refln_loop = block.find_loop("_diffrn_refln.index_h").get_loop();
    default_loop = refln_loop ? refln_loop : diffrn_refln_loop;
    return *this;
  }

  //! @brief Check if block contains valid reflection data.
  //! @return True if a reflection loop is available
  bool ok() const { return default_loop != nullptr; }

  //! @brief Throw error if block is invalid.
  //! @throws std::runtime_error if no reflection loop is available
  void check_ok() const { if (!ok()) fail("Invalid ReflnBlock"); }

  //! @brief Get position after "_refln." or "_diffrn_refln." in tag names.
  //! @return Offset (7 for merged, 14 for unmerged)
  size_t tag_offset() const { return default_loop == refln_loop ? 7 : 14; }

  //! @brief Switch between merged and unmerged reflection data.
  //! @param unmerged If true, use _diffrn_refln; if false, use _refln
  void use_unmerged(bool unmerged) {
    default_loop = unmerged ? diffrn_refln_loop : refln_loop;
  }

  //! @brief Check if currently using merged reflection data.
  //! @return True if using _refln loop
  bool is_merged() const { return ok() && default_loop == refln_loop; }

  //! @deprecated Use !is_merged() instead
  //! @brief Check if currently using unmerged reflection data.
  //! @return True if using _diffrn_refln loop
  bool is_unmerged() const { return ok() && default_loop == diffrn_refln_loop; }

  //! @brief Get list of column labels (without category prefix).
  //! @return Vector of column names (e.g., "index_h", "F_meas")
  //!
  //! Strips the "_refln." or "_diffrn_refln." prefix from tag names.
  std::vector<std::string> column_labels() const {
    check_ok();
    std::vector<std::string> labels(default_loop->tags.size());
    for (size_t i = 0; i != labels.size(); ++i)
      labels[i].assign(default_loop->tags[i], tag_offset(), std::string::npos);
    return labels;
  }

  //! @brief Find column index by tag name.
  //! @param tag Column name without prefix (e.g., "index_h", "F_meas")
  //! @return Column index, or -1 if not found
  int find_column_index(const std::string& tag) const {
    if (!ok())
      return -1;
    size_t name_pos = tag_offset();
    for (int i = 0; i != (int) default_loop->tags.size(); ++i)
      if (default_loop->tags[i].compare(name_pos, std::string::npos, tag) == 0)
        return i;
    return -1;
  }

  //! @brief Get column index by tag name, throwing if not found.
  //! @param tag Column name without prefix (e.g., "index_h", "F_meas")
  //! @return Column index
  //! @throws std::runtime_error if column not found or block invalid
  size_t get_column_index(const std::string& tag) const {
    int idx = find_column_index(tag);
    if (idx == -1) {
      if (!ok())
        fail("No reflection data in block ", block.name);
      const char* category = default_loop == refln_loop ? "_refln." : "_diffrn_refln.";
      fail("Column not found in block ", block.name, ": ", category, tag);
    }
    return idx;
  }

  //! @brief Extract a column as a vector of values.
  //! @tparam T Value type (int, double, std::string, etc.)
  //! @param tag Column name without prefix
  //! @param null Value to use for missing/null CIF values
  //! @return Vector of column values
  template<typename T>
  std::vector<T> make_vector(const std::string& tag, T null) const {
    size_t n = get_column_index(tag);
    std::vector<T> v(default_loop->length());
    for (size_t j = 0; j != v.size(); n += default_loop->width(), ++j)
      v[j] = cif::as_any(default_loop->values[n], null);
    return v;
  }

  //! @brief Get column indices for Miller indices (h, k, l).
  //! @return Array of three column indices [h_idx, k_idx, l_idx]
  //! @throws std::runtime_error if any index column not found
  std::array<size_t,3> get_hkl_column_indices() const {
    return {{get_column_index("index_h"),
             get_column_index("index_k"),
             get_column_index("index_l")}};
  }

  //! @brief Extract all Miller indices as a vector.
  //! @return Vector of Miller indices [h,k,l] for all reflections
  std::vector<Miller> make_miller_vector() const {
    auto hkl_idx = get_hkl_column_indices();
    std::vector<Miller> v(default_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += default_loop->width())
      for (int i = 0; i != 3; ++i)
        v[j][i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
    return v;
  }

  //! @brief Calculate 1/d² for all reflections.
  //! @return Vector of 1/d² values in Ų
  //! @throws std::runtime_error if unit cell is not set
  std::vector<double> make_1_d2_vector() const {
    if (!cell.is_crystal() || cell.a <= 0)
      fail("Unit cell is not known");
    auto hkl_idx = get_hkl_column_indices();
    std::vector<double> r(default_loop->length());
    for (size_t j = 0, n = 0; j != r.size(); j++, n += default_loop->width()) {
      Miller hkl;
      for (int i = 0; i != 3; ++i)
        hkl[i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
      r[j] = cell.calculate_1_d2(hkl);
    }
    return r;
  }

  //! @brief Calculate resolution (d-spacing) for all reflections.
  //! @return Vector of d values in Angstroms
  //! @throws std::runtime_error if unit cell is not set
  std::vector<double> make_d_vector() const {
    std::vector<double> vec = make_1_d2_vector();
    for (double& d : vec)
      d = 1.0 / std::sqrt(d);
    return vec;
  }
};

//! @brief Convert vector of CIF blocks to ReflnBlocks.
//! @param blocks Vector of CIF blocks (moved, not copied)
//! @return Vector of ReflnBlock wrappers
//!
//! Moves blocks from the argument to the return value. Some blocks may miss
//! space group or unit cell information; this function attempts to fill in
//! missing metadata by propagating from the first valid block.
inline
std::vector<ReflnBlock> as_refln_blocks(std::vector<cif::Block>&& blocks) {
  std::vector<ReflnBlock> rvec;
  rvec.reserve(blocks.size());
  for (cif::Block& block : blocks)
    rvec.emplace_back(std::move(block));
  blocks.clear();
  // Some blocks miss space group or unit cell, try to fill it in.
  const SpaceGroup* first_sg = nullptr;
  const UnitCell* first_cell = nullptr;
  for (ReflnBlock& rblock : rvec) {
    if (!first_sg)
      first_sg = rblock.spacegroup;
    else if (!rblock.spacegroup)
      rblock.spacegroup = first_sg;
    if (rblock.cell.is_crystal()) {
      if (!first_cell)
        first_cell = &rblock.cell;
    } else if (first_cell) {
      rblock.cell = *first_cell;
    }
  }
  return rvec;
}

//! @brief Get the first merged block with required column labels.
//! @param blocks Vector of CIF blocks to search (moved, not copied)
//! @param labels Required column names (without _refln. prefix)
//! @param block_name Optional block name filter (nullptr for any block)
//! @return ReflnBlock containing all required columns
//! @throws std::runtime_error if required columns or block not found
//!
//! Searches for the first (merged) block with all required labels.
//! Optionally, block name can be specified to filter by block.
inline ReflnBlock get_refln_block(std::vector<cif::Block>&& blocks,
                                  const std::vector<std::string>& labels,
                                  const char* block_name=nullptr) {
  const SpaceGroup* first_sg = nullptr;
  bool has_block = false;
  for (cif::Block& block : blocks) {
    if (!first_sg)
      if (const std::string* hm = impl::find_spacegroup_hm_value(block)) {
        first_sg = find_spacegroup_by_name(cif::as_string(*hm));
        if (first_sg && first_sg->ext == 'H') {
          UnitCell cell;
          impl::set_cell_from_mmcif(block, cell);
          first_sg = find_spacegroup_by_name(cif::as_string(*hm), cell.alpha, cell.gamma);
        }
      }
    if (!block_name || block.name == block_name) {
      has_block = true;
      if (cif::Loop* loop = block.find_loop("_refln.index_h").get_loop())
        if (std::all_of(labels.begin(), labels.end(),
              [&](const std::string& s) { return s.empty() || loop->has_tag("_refln."+s); })) {
          ReflnBlock rblock(std::move(block));
          if (!rblock.spacegroup && first_sg)
            rblock.spacegroup = first_sg;
          return rblock;
        }
    }
  }
  if (block_name && !has_block)
    fail("Required block not found in SF-mmCIF file.");
  fail("Tags not found in SF-mmCIF file: _refln.", join_str(labels, ", _refln."));
}

//! @brief Convert HKL-format CIF block to ReflnBlock.
//! @param block CIF block in HKL format (swapped, not copied)
//! @return ReflnBlock wrapper
//!
//! Handles the older HKL-style CIF format (non-mmCIF) with tags like
//! _refln_index_h instead of _refln.index_h. Extracts unit cell and
//! space group from legacy tags.
inline ReflnBlock hkl_cif_as_refln_block(cif::Block& block) {
  ReflnBlock rb;
  rb.block.swap(block);
  rb.entry_id = rb.block.name;
  impl::set_cell_from_mmcif(rb.block, rb.cell, /*mmcif=*/false);
  const char* hm_tag = "_symmetry_space_group_name_H-M";
  if (const std::string* hm = rb.block.find_value(hm_tag))
    rb.spacegroup = find_spacegroup_by_name(cif::as_string(*hm),
                                            rb.cell.alpha, rb.cell.gamma);
  rb.cell.set_cell_images_from_spacegroup(rb.spacegroup);
  rb.refln_loop = rb.block.find_loop("_refln_index_h").get_loop();
  rb.default_loop = rb.refln_loop;
  return rb;
}

//! @brief Abstraction of reflection data source (analogous to MtzDataProxy).
//!
//! Provides uniform access to reflection data for algorithms that work with
//! multiple data formats. Allows generic code to iterate over reflections
//! and access Miller indices, intensities, and other columns.
struct ReflnDataProxy {
  const ReflnBlock& rb_;  //!< Reference to underlying ReflnBlock
  std::array<size_t,3> hkl_cols_;  //!< Cached column indices for h, k, l

  //! @brief Construct proxy for a ReflnBlock.
  //! @param rb ReflnBlock to wrap
  explicit ReflnDataProxy(const ReflnBlock& rb)
    : rb_(rb), hkl_cols_(rb_.get_hkl_column_indices()) {}

  //! @brief Get row stride (number of columns).
  //! @return Number of tags in the reflection loop
  size_t stride() const { return loop().tags.size(); }

  //! @brief Get total number of values (rows × columns).
  //! @return Total size of the value array
  size_t size() const { return loop().values.size(); }

  using num_type = double;  //!< Numeric type for reflection data

  //! @brief Get numeric value at given position.
  //! @param n Linear offset into values array
  //! @return Parsed double value
  double get_num(size_t n) const { return cif::as_number(loop().values[n]); }

  //! @brief Get unit cell.
  //! @return Reference to unit cell parameters
  const UnitCell& unit_cell() const { return rb_.cell; }

  //! @brief Get space group.
  //! @return Pointer to space group (may be nullptr)
  const SpaceGroup* spacegroup() const { return rb_.spacegroup; }

  //! @brief Get Miller indices at given row offset.
  //! @param offset Row offset (stride × row_number)
  //! @return Miller indices [h, k, l]
  Miller get_hkl(size_t offset) const {
    return {{get_int(offset + hkl_cols_[0]),
             get_int(offset + hkl_cols_[1]),
             get_int(offset + hkl_cols_[2])}};
  }

  //! @brief Get column index by label.
  //! @param label Column name without prefix
  //! @return Column index
  size_t column_index(const std::string& label) const { return rb_.get_column_index(label); }

private:
  const cif::Loop& loop() const { rb_.check_ok(); return *rb_.default_loop; }
  int get_int(size_t n) const { return cif::as_int(loop().values[n]); }
};

//! @brief Create data proxy for a ReflnBlock.
//! @param rb ReflnBlock to wrap
//! @return ReflnDataProxy for the block
inline ReflnDataProxy data_proxy(const ReflnBlock& rb) { return ReflnDataProxy(rb); }

} // namespace gemmi
#endif
