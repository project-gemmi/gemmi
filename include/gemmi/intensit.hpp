// Copyright 2020 Global Phasing Ltd.
//
/// @file
/// @brief Intensities class for reading and merging intensity data from various formats.
//
// Class Intensities that reads multi-record data from MTZ, mmCIF or XDS_ASCII
// and merges it into mean or anomalous intensities.
// It can also read merged data.

#ifndef GEMMI_INTENSIT_HPP_
#define GEMMI_INTENSIT_HPP_

#include <cstdint>      // for int8_t
#include <unordered_map>
#include "symmetry.hpp"
#include "unitcell.hpp"
#include "util.hpp"     // for vector_remove_if
#include "stats.hpp"    // for Correlation

namespace gemmi {

struct Binner;
struct Mtz;
struct XdsAscii;
struct ReflnBlock;
namespace cif { struct Block; }
using std::int8_t;

/// Data type of intensity data: unmerged, mean intensity, or anomalous intensities.
///
/// When requesting a particular data type, the MergedMA/MergedAM/UAM variants
/// allow fallback to secondary options (e.g., MergedMA = Mean if available, else Anomalous).
enum class DataType {
  Unknown,     ///< Unknown or unspecified intensity type
  Unmerged,    ///< Unmerged (multi-record) intensity data
  Mean,        ///< Mean intensity <I>
  Anomalous,   ///< Anomalous intensities (I+/I-)
  MergedMA,    ///< Mean if available, otherwise Anomalous (fallback type)
  MergedAM,    ///< Anomalous if available, otherwise Mean (fallback type)
  UAM          ///< Unmerged if available, otherwise MergedAM (fallback type)
};

/// @brief Statistics calculated for a resolution shell (bin) of merged intensities.
///
/// Accumulates numerators and denominators for R-merge, R-meas, and R-pim.
/// Can be summed across shells to compute overall statistics.
struct GEMMI_DLL MergingStats {
  int all_refl = 0;        ///< Total number of observations (all reflections)
  int unique_refl = 0;     ///< Number of unique reflections
  int stats_refl = 0;      ///< Unique reflections with 2+ observations (for statistics)
  double r_merge_num = 0;  ///< Numerator for R-merge calculation
  double r_meas_num = 0;   ///< Numerator for R-meas (redundancy-weighted) calculation
  double r_pim_num = 0;    ///< Numerator for R-pim (precision-indicating merge) calculation
  double r_denom = 0;      ///< Denominator for R-merge/meas/pim calculations
  double sum_ibar = 0;     ///< Sum of mean intensities for CC1/2
  double sum_ibar2 = 0;    ///< Sum of squared mean intensities for CC1/2
  double sum_sig2_eps = 0; ///< Sum of variance terms for CC1/2

  /// Accumulate statistics from another shell.
  /// @param o Statistics from another resolution shell.
  /// Adding two MergingStats gives the same result as calculating statistics
  /// for the combined shells from the start.
  void add_other(const MergingStats& o) {
    all_refl += o.all_refl;
    unique_refl += o.unique_refl;
    stats_refl += o.stats_refl;
    r_merge_num += o.r_merge_num;
    r_meas_num += o.r_meas_num;
    r_pim_num += o.r_pim_num;
    r_denom += o.r_denom;
    sum_ibar += o.sum_ibar;
    sum_ibar2 += o.sum_ibar2;
    sum_sig2_eps += o.sum_sig2_eps;
  }

  /// Compute R-merge for this shell.
  double r_merge() const { return r_merge_num / r_denom; }
  /// Compute redundancy-weighted R-meas for this shell.
  double r_meas() const { return r_meas_num / r_denom; }
  /// Compute precision-indicating R-pim for this shell.
  double r_pim() const { return r_pim_num / r_denom; }
  /// Compute CC1/2 (correlation coefficient of half-datasets).
  double cc_half() const; // calculated using sigma-tau method
  /// Compute CC* using Spearman-Brown prediction formula from CC1/2.
  double cc_full() const {
    double cc = cc_half();
    return 2 * cc / (1 + cc);
  }
  /// Estimate of overall correlation coefficient from CC1/2.
  double cc_star() const { return std::sqrt(cc_full()); }
};

/// @brief Extract STARANISO anisotropy B-tensor from MTZ file.
/// @param mtz MTZ file object to read from.
/// @param output Anisotropic B-tensor (3x3 symmetric matrix in Voigt notation).
/// @return STARANISO version string if found, empty string otherwise.
GEMMI_DLL std::string read_staraniso_b_from_mtz(const Mtz& mtz, SMat33<double>& output);

/// @brief Container for intensity data from reflection measurements.
///
/// Stores multi-record (unmerged) or merged intensities with metadata such as
/// unit cell, space group, and wavelength. Supports merging operations and
/// import from MTZ, mmCIF, and XDS_ASCII formats.
struct GEMMI_DLL Intensities {
  /// @brief A single reflection record with intensity, sigma, and metadata.
  struct Refl {
    Miller hkl;       ///< Miller indices (h, k, l)
    int8_t isign;     ///< Intensity component: 1=I(+), -1=I(-), 0=mean/unmerged
    int8_t isym;      ///< Symmetry operator encoding (ISYM in MTZ for unmerged data)
    short nobs;       ///< Number of observations (used during merging)
    double value;     ///< Intensity value
    double sigma;     ///< Standard deviation of intensity

    /// Compare reflections by (h,k,l, isign) for sorting.
    bool operator<(const Refl& o) const {
      return std::tie(hkl[0], hkl[1], hkl[2], isign) <
             std::tie(o.hkl[0], o.hkl[1], o.hkl[2], o.isign);
    }
    /// @brief Get intensity label ("&lt;I&gt;", "I(+)", or "I(-)").
    const char* intensity_label() const {
      if (isign == 0)
        return "<I>";
      return isign > 0 ? "I(+)" : "I(-)";
    }
    /// @brief Format HKL and intensity label as a string.
    std::string hkl_label() const {
      return cat(intensity_label(), " (", hkl[0], ' ', hkl[1], ' ', hkl[2], ')');
    }
  };

  /// @brief Anisotropic scaling tensor for STARANISO B-factor correction.
  struct AnisoScaling {
    SMat33<double> b = {0., 0., 0., 0., 0., 0.}; ///< Symmetric B-tensor in Voigt notation

    /// Check if anisotropic tensor is set (non-zero).
    bool ok() const { return !b.all_zero(); }
    /// @brief Compute scaling factor for a reflection at given HKL.
    /// @param hkl Miller indices
    /// @param cell Unit cell to convert HKL to reciprocal-space vector
    /// @return Exponential scaling factor exp(0.5 * B * s * s)
    double scale(const Miller& hkl, const UnitCell& cell) const {
      Vec3 s = cell.frac.mat.left_multiply(Vec3(hkl[0], hkl[1], hkl[2]));
      return std::exp(0.5 * b.r_u_r(s));
    }
  };

  std::vector<Refl> data;                  ///< Reflection records
  const SpaceGroup* spacegroup = nullptr;  ///< Space group (not owned by this object)
  UnitCell unit_cell;                      ///< Crystal unit cell parameters
  double unit_cell_rmsd[6] = {0., 0., 0., 0., 0., 0.}; ///< RMSDs of unit cell parameters
  double wavelength;                       ///< Diffraction wavelength in Angstroms
  DataType type = DataType::Unknown;       ///< Type of intensity data stored
  std::vector<Op> isym_ops;                ///< Symmetry operators (for unmerged data)
  AnisoScaling staraniso_b;                ///< STARANISO anisotropy correction tensor

  /// @brief Get string representation of intensity data type.
  /// @param data_type Type to convert.
  /// @return Human-readable label ("I", "&lt;I&gt;", "I+/I-", "n/a", etc.).
  static const char* type_str(DataType data_type) {
    switch (data_type) {
      case DataType::Unmerged: return "I";
      case DataType::Mean: return "<I>";
      case DataType::Anomalous: return "I+/I-";
      case DataType::MergedAM:
      case DataType::MergedMA:
      case DataType::UAM:
      case DataType::Unknown: return "n/a";
    }
    unreachable();
  }

  /// @brief Get string representation of this object's intensity data type.
  const char* type_str() const { return Intensities::type_str(type); }

  /// @brief Get space group name (Hermann-Mauguin symbol) or "none".
  std::string spacegroup_str() const { return spacegroup ? spacegroup->xhm() : "none"; }

  /// @brief Get minimum and maximum resolution of reflections.
  /// @return {d_max, d_min} as array of two doubles.
  std::array<double,2> resolution_range() const;

  /// @brief Calculate correlation of intensity values between two sorted lists.
  /// @param other Another Intensities object to correlate with.
  /// @return Correlation object with matching reflections analyzed.
  Correlation calculate_correlation(const Intensities& other) const;

  /// @brief Add a single reflection if its data is valid.
  /// Skips reflections with NaN or non-positive sigma (rejected by XDS, etc.).
  /// @param hkl Miller indices
  /// @param isign Intensity sign (1 for I+, -1 for I-, 0 for mean)
  /// @param isym Symmetry operator encoding
  /// @param value Intensity value
  /// @param sigma Standard deviation of intensity
  void add_if_valid(const Miller& hkl, int8_t isign, int8_t isym, double value, double sigma) {
    // XDS marks rejected reflections with negative sigma.
    // Sigma 0.0 rarely happens (e.g. 5tkn), but is also problematic.
    if (!std::isnan(value) && sigma > 0)
      data.push_back({hkl, isign, isym, /*nobs=*/0, value, sigma});
  }

  /// @brief Remove reflections that are forbidden by space group symmetry.
  void remove_systematic_absences() {
    if (!spacegroup)
      return;
    GroupOps gops = spacegroup->operations();
    vector_remove_if(data, [&](Refl& x) { return gops.is_systematically_absent(x.hkl); });
  }

  /// @brief Sort reflections by (h, k, l, isign) in ascending order.
  void sort() { std::sort(data.begin(), data.end()); }

  /// @brief Merge reflections in-place to a specified data type (mean or anomalous).
  /// @param new_type Target data type for merged intensities.
  void merge_in_place(DataType new_type);

  /// @brief Create a merged copy without modifying this object.
  /// @param new_type Target data type for merged intensities.
  /// @return New Intensities object with merged data.
  Intensities merged(DataType new_type) {
    Intensities m(*this);
    m.merge_in_place(new_type);
    return m;
  }

  /// @brief Calculate R-merge and related statistics for each resolution shell.
  /// @param binner Resolution shell binning (or nullptr for all data in one shell).
  /// @param use_weights Weighting scheme: 'Y'=Aimless style, 'U'=unweighted, 'X'=XDS style.
  /// @return Vector of MergingStats, one per shell.
  std::vector<MergingStats> calculate_merging_stats(const Binner* binner,
                                                    char use_weights='Y') const;

  /// @brief Prepare data for merging and classify anomalous/mean component.
  /// Call with DataType::Anomalous before calculate_merging_stats() to get I+/I- stats.
  /// @param new_type Target data type.
  /// @return Classified data type after preparation.
  DataType prepare_for_merging(DataType new_type);

  /// @brief Convert unmerged ISYM indices to ASU indices using stored isym_ops.
  void switch_to_asu_indices();

  /// @brief Load unmerged intensities from MTZ file.
  /// @param mtz MTZ object to read from.
  void import_unmerged_intensities_from_mtz(const Mtz& mtz);
  /// @brief Load mean intensities from MTZ file.
  /// @param mtz MTZ object to read from.
  void import_mean_intensities_from_mtz(const Mtz& mtz);
  /// @brief Load anomalous intensities (I+/I-) from MTZ file.
  /// @param mtz MTZ object to read from.
  /// @param check_complete If true, throw if anomalous data is null where expected.
  void import_anomalous_intensities_from_mtz(const Mtz& mtz, bool check_complete=false);

  /// @brief Load intensities from MTZ file, auto-detecting data type.
  /// @param mtz MTZ object to read from.
  /// @param data_type Requested data type; DataType::Unknown auto-detects.
  void import_mtz(const Mtz& mtz, DataType data_type=DataType::Unknown);

  /// @brief Load unmerged intensities from mmCIF reflection block.
  /// @param rb Reflection block to read from.
  void import_unmerged_intensities_from_mmcif(const ReflnBlock& rb);
  /// @brief Load mean intensities from mmCIF reflection block.
  /// @param rb Reflection block to read from.
  void import_mean_intensities_from_mmcif(const ReflnBlock& rb);
  /// @brief Load anomalous intensities (I+/I-) from mmCIF reflection block.
  /// @param rb Reflection block to read from.
  /// @param check_complete If true, throw if anomalous data is null where expected.
  void import_anomalous_intensities_from_mmcif(const ReflnBlock& rb, bool check_complete=false);

  /// @brief Load structure factor squared (F^2) from mmCIF reflection block.
  /// @param rb Reflection block to read from.
  void import_f_squared_from_mmcif(const ReflnBlock& rb);

  /// @brief Load intensities from mmCIF reflection block, auto-detecting data type.
  /// @param rb Reflection block to read from.
  /// @param data_type Requested data type; DataType::Unknown auto-detects.
  void import_refln_block(const ReflnBlock& rb, DataType data_type=DataType::Unknown);

  /// @brief Load intensities from XDS_ASCII file.
  /// @param xds XDS_ASCII object to read from.
  void import_xds(const XdsAscii& xds);

  /// @brief Extract and store STARANISO B-tensor from MTZ file.
  /// @param mtz MTZ object to read from.
  /// @return STARANISO version string if found, empty string otherwise.
  std::string take_staraniso_b_from_mtz(const Mtz& mtz);

  /// @brief Extract and store STARANISO B-tensor from mmCIF block.
  /// @param block CIF block to read from.
  /// @return True if tensor was found and loaded, false otherwise.
  bool take_staraniso_b_from_mmcif(const cif::Block& block);

  /// @brief Create a merged MTZ file from these intensities.
  /// @param with_nobs If true, include NOBS (observation count) column.
  /// @return Mtz object with merged data.
  Mtz prepare_merged_mtz(bool with_nobs);
};

/// @brief Adapter providing DataProxy interface to Intensities data.
/// Enables use of Intensities with generic algorithms expecting standard data proxy interface.
struct IntensitiesDataProxy {
  const Intensities& intensities_; ///< Reference to underlying Intensities object

  size_t stride() const { return 1; }
  size_t size() const { return intensities_.data.size(); }
  const SpaceGroup* spacegroup() const { return intensities_.spacegroup; }
  const UnitCell& unit_cell() const { return intensities_.unit_cell; }
  Miller get_hkl(size_t offset) const { return intensities_.data[offset].hkl; }
  double get_num(size_t n) const { return intensities_.data[n].value; }
};

/// @brief Infer intensity data type from reflection data under symmetry.
///
/// Examines unique reflections and detects whether data is unmerged (multiple
/// copies of same HKL), mean intensity (single copy), or anomalous (both I+/I-).
///
/// @tparam DataProxy Type with spacegroup(), unit_cell(), size(), stride(),
///                    get_hkl() interface.
/// @param proxy Data proxy object to analyze.
/// @return Pair of (inferred DataType, number of unique HKLs in ASU).
template<typename DataProxy>
std::pair<DataType, size_t> check_data_type_under_symmetry(const DataProxy& proxy) {
  const SpaceGroup* sg = proxy.spacegroup();
  if (!sg)
    return {DataType::Unknown, 0};
  std::unordered_map<Op::Miller, int, MillerHash> seen;
  ReciprocalAsu asu(sg);
  GroupOps gops = sg->operations();
  bool centric = gops.is_centrosymmetric();
  DataType data_type = DataType::Mean;
  for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
    auto hkl_sign = asu.to_asu_sign(proxy.get_hkl(i), gops);
    int sign = hkl_sign.second ? 2 : 1;  // 2=positive, 1=negative
    auto r = seen.emplace(hkl_sign.first, sign);
    if (data_type != DataType::Unmerged && !r.second) {
      if ((r.first->second & sign) != 0 || centric) {
        data_type = DataType::Unmerged;
      } else {
        r.first->second |= sign;
        data_type = DataType::Anomalous;
      }
    }
  }
  return {data_type, seen.size()};
}

} // namespace gemmi
#endif
