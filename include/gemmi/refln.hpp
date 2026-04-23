// Copyright 2019 Global Phasing Ltd.
//
// Reads reflection data from the mmCIF format.

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

/// @file
/// @brief Structure and accessors for reflection data from SF-mmCIF files.

/// @brief Wrapper for reflection data block from mmCIF file.
///
/// Provides column access and HKL index management for reflection data
/// stored in CIF loops (either merged _refln or unmerged _diffrn_refln categories).
struct ReflnBlock {
  /// CIF data block containing reflection data.
  cif::Block block;
  /// Entry identifier from _entry.id tag.
  std::string entry_id;
  /// Unit cell parameters.
  UnitCell cell;
  /// Pointer to space group; nullptr if unknown.
  const SpaceGroup* spacegroup = nullptr;
  /// X-ray wavelength in Angstroms (0 if multiple wavelengths).
  double wavelength;
  /// Number of wavelengths in the data block.
  int wavelength_count;
  /// Pointer to _refln loop (merged data) or nullptr if absent.
  cif::Loop* refln_loop = nullptr;
  /// Pointer to _diffrn_refln loop (unmerged data) or nullptr if absent.
  cif::Loop* diffrn_refln_loop = nullptr;
  /// Points to active loop (refln_loop or diffrn_refln_loop).
  cif::Loop* default_loop = nullptr;

  /// Default constructor (empty block).
  ReflnBlock() = default;
  /// Move constructor.
  ReflnBlock(ReflnBlock&& rblock_) = default;
  /// Construct from CIF block; extracts cell, spacegroup, wavelength, and reflection loops.
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
  /// Move assignment.
  ReflnBlock& operator=(ReflnBlock&&) = default;
  /// Copy assignment (deep copy of block and pointers).
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

  /// @brief Check if block contains valid reflection data.
  /// @return True if default_loop is set and not null.
  bool ok() const { return default_loop != nullptr; }
  /// @brief Throw exception if block is not valid.
  void check_ok() const { if (!ok()) fail("Invalid ReflnBlock"); }

  /// @brief Get offset to tag name (after "_refln." or "_diffrn_refln.").
  /// @return 7 for merged, 14 for unmerged.
  size_t tag_offset() const { return default_loop == refln_loop ? 7 : 14; }

  /// @brief Switch between merged and unmerged reflection loops.
  /// @param unmerged If true, use _diffrn_refln loop; otherwise use _refln.
  void use_unmerged(bool unmerged) {
    default_loop = unmerged ? diffrn_refln_loop : refln_loop;
  }
  /// @brief Check if active loop is merged reflection data.
  bool is_merged() const { return ok() && default_loop == refln_loop; }
  /// @brief Check if active loop is unmerged reflection data (deprecated).
  bool is_unmerged() const { return ok() && default_loop == diffrn_refln_loop; }

  /// @brief Get list of column labels (without category prefix).
  /// @return Vector of tag names from the active loop.
  std::vector<std::string> column_labels() const {
    check_ok();
    std::vector<std::string> labels(default_loop->tags.size());
    for (size_t i = 0; i != labels.size(); ++i)
      labels[i].assign(default_loop->tags[i], tag_offset(), std::string::npos);
    return labels;
  }

  /// @brief Find column index by tag name (without category prefix).
  /// @param tag Column name (e.g., "index_h", "F_meas_sigma_au").
  /// @return Column index, or -1 if not found.
  int find_column_index(const std::string& tag) const {
    if (!ok())
      return -1;
    size_t name_pos = tag_offset();
    for (int i = 0; i != (int) default_loop->tags.size(); ++i)
      if (default_loop->tags[i].compare(name_pos, std::string::npos, tag) == 0)
        return i;
    return -1;
  }

  /// @brief Get column index by tag name, throwing exception if not found.
  /// @param tag Column name.
  /// @return Column index.
  /// @throws gemmi::fail if column does not exist.
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

  /// @brief Extract column values as typed vector.
  /// @tparam T Value type (e.g., int, double).
  /// @param tag Column name.
  /// @param null Default value for missing data.
  /// @return Vector of converted values.
  template<typename T>
  std::vector<T> make_vector(const std::string& tag, T null) const {
    size_t n = get_column_index(tag);
    std::vector<T> v(default_loop->length());
    for (size_t j = 0; j != v.size(); n += default_loop->width(), ++j)
      v[j] = cif::as_any(default_loop->values[n], null);
    return v;
  }

  /// @brief Get column indices for h, k, l indices.
  /// @return Array of 3 column indices for index_h, index_k, index_l.
  std::array<size_t,3> get_hkl_column_indices() const {
    return {{get_column_index("index_h"),
             get_column_index("index_k"),
             get_column_index("index_l")}};
  }

  /// @brief Extract Miller indices from all reflections.
  /// @return Vector of Miller indices.
  std::vector<Miller> make_miller_vector() const {
    auto hkl_idx = get_hkl_column_indices();
    std::vector<Miller> v(default_loop->length());
    for (size_t j = 0, n = 0; j != v.size(); j++, n += default_loop->width())
      for (int i = 0; i != 3; ++i)
        v[j][i] = cif::as_int(default_loop->values[n + hkl_idx[i]]);
    return v;
  }

  /// @brief Calculate 1/d^2 for all reflections using unit cell parameters.
  /// @return Vector of 1/d^2 values.
  /// @throws gemmi::fail if unit cell is not set.
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

  /// @brief Calculate d-spacing for all reflections using unit cell parameters.
  /// @return Vector of d-spacing values in Angstroms.
  std::vector<double> make_d_vector() const {
    std::vector<double> vec = make_1_d2_vector();
    for (double& d : vec)
      d = 1.0 / std::sqrt(d);
    return vec;
  }
};

/// @brief Convert CIF blocks to ReflnBlocks, propagating cell and spacegroup info.
///
/// Fills in missing cell or spacegroup data by copying from the first block
/// that contains it. Moves blocks from input to output.
/// @param blocks Input CIF blocks (consumed).
/// @return Vector of ReflnBlocks.
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

/// @brief Find the first merged reflection block containing specified columns.
///
/// Searches blocks for one with _refln loop and required column labels.
/// Propagates spacegroup from first block if needed.
/// @param blocks Input CIF blocks (consumed).
/// @param labels Required column names (empty string means optional).
/// @param block_name Optional: if provided, only this named block is considered.
/// @return The matching ReflnBlock.
/// @throws gemmi::fail if no matching block or required columns not found.
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

/// @brief Convert CIF block from non-mmCIF format to ReflnBlock.
/// @param block Input CIF block (e.g., from hkl file).
/// @return ReflnBlock with data from _refln_index_h loop.
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

/// @brief Generic data source abstraction over a ReflnBlock for row iteration.
///
/// Provides uniform interface for accessing reflection columns (similar to MtzDataProxy).
struct ReflnDataProxy {
  /// Reference to underlying ReflnBlock.
  const ReflnBlock& rb_;
  /// Cached indices for h, k, l columns.
  std::array<size_t,3> hkl_cols_;
  /// Initialize proxy from a ReflnBlock.
  explicit ReflnDataProxy(const ReflnBlock& rb)
    : rb_(rb), hkl_cols_(rb_.get_hkl_column_indices()) {}
  /// Number of columns (values per reflection).
  size_t stride() const { return loop().tags.size(); }
  /// Total number of values (stride * reflection count).
  size_t size() const { return loop().values.size(); }
  /// Numeric value type.
  using num_type = double;
  /// Get numeric value at flattened loop index.
  double get_num(size_t n) const { return cif::as_number(loop().values[n]); }
  /// Get unit cell.
  const UnitCell& unit_cell() const { return rb_.cell; }
  /// Get spacegroup.
  const SpaceGroup* spacegroup() const { return rb_.spacegroup; }
  /// Get Miller indices at given offset in loop.
  Miller get_hkl(size_t offset) const {
    return {{get_int(offset + hkl_cols_[0]),
             get_int(offset + hkl_cols_[1]),
             get_int(offset + hkl_cols_[2])}};
  }
  /// Get column index by label.
  size_t column_index(const std::string& label) const { return rb_.get_column_index(label); }
private:
  /// Get active loop (with bounds check).
  const cif::Loop& loop() const { rb_.check_ok(); return *rb_.default_loop; }
  /// Get integer value at flattened loop index.
  int get_int(size_t n) const { return cif::as_int(loop().values[n]); }
};

/// @brief Create a data proxy over a ReflnBlock.
/// @param rb ReflnBlock to wrap.
/// @return ReflnDataProxy for generic access.
inline ReflnDataProxy data_proxy(const ReflnBlock& rb) { return ReflnDataProxy(rb); }

} // namespace gemmi
#endif
