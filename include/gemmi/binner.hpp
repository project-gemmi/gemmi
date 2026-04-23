// Copyright 2022 Global Phasing Ltd.
//
/// @file
/// @brief Binner class for organizing reflections into resolution shells.

#ifndef GEMMI_BINNER_HPP_
#define GEMMI_BINNER_HPP_

#include <vector>
#include <limits>        // for numeric_limits
#include "unitcell.hpp"  // for UnitCell

namespace gemmi {

/// @brief Divide reflections into resolution shells (bins) for statistics.
///
/// Organizes reflection data into bins by resolution. Supports multiple
/// binning schemes (equal count, linear/quadratic/cubic in d or 1/d^2).
/// Provides fast bin lookup for both sorted and unsorted data.
struct Binner {
  /// @brief Strategy for dividing reflections into resolution shells.
  enum class Method {
    EqualCount, ///< Bins with approximately equal number of reflections
    Dstar,      ///< Linear spacing in 1/d (resolution)
    Dstar2,     ///< Linear spacing in 1/d^2 (squared reciprocal spacing)
    Dstar3,     ///< Cubic spacing (linear in (1/d^2)^(1/3))
  };

  /// @brief Set up bins using pre-calculated 1/d^2 values.
  /// @param nbins Number of resolution bins to create.
  /// @param method Binning scheme (EqualCount, Dstar, Dstar2, Dstar3).
  /// @param inv_d2 Vector of 1/d^2 values for all reflections (moved from caller).
  /// @param cell_ Unit cell for resolution calculations (nullptr uses already-set cell).
  void setup_from_1_d2(int nbins, Method method, std::vector<double>&& inv_d2,
                       const UnitCell* cell_) {
    if (nbins < 1)
      fail("Binner: nbins argument must be positive");
    if (inv_d2.empty())
      fail("Binner: no data");
    if (cell_)
      cell = *cell_;
    if (!cell.is_crystal())
      fail("Binner: unknown unit cell");
    // first setup 2N bins to get both bin limits and middle points
    limits.resize(2 * nbins);
    if (method == Method::EqualCount) {
      std::sort(inv_d2.begin(), inv_d2.end());
      min_1_d2 = inv_d2.front();
      max_1_d2 = inv_d2.back();
    } else {
      min_1_d2 = max_1_d2 = inv_d2.front();
      for (double x : inv_d2) {
        if (x < min_1_d2)
          min_1_d2 = x;
        if (x > max_1_d2)
          max_1_d2 = x;
      }
    }
    switch (method) {
      case Method::EqualCount: {
        double avg_count = double(inv_d2.size()) / limits.size();
        for (size_t i = 1; i < limits.size(); ++i)
          limits[i-1] = inv_d2[int(avg_count * i)];
        break;
      }
      case Method::Dstar2: {
        double step = (max_1_d2 - min_1_d2) / limits.size();
        for (size_t i = 1; i < limits.size(); ++i)
          limits[i-1] = min_1_d2 + i * step;
        break;
      }
      case Method::Dstar: {
        double min_1_d = std::sqrt(min_1_d2);
        double max_1_d = std::sqrt(max_1_d2);
        double step = (max_1_d - min_1_d) / limits.size();
        for (size_t i = 1; i < limits.size(); ++i)
          limits[i-1] = sq(min_1_d + i * step);
        break;
      }
      case Method::Dstar3: {
        double min_1_d3 = min_1_d2 * std::sqrt(min_1_d2);
        double max_1_d3 = max_1_d2 * std::sqrt(max_1_d2);
        double step = (max_1_d3 - min_1_d3) / limits.size();
        for (size_t i = 1; i < limits.size(); ++i)
          limits[i-1] = sq(std::cbrt(min_1_d3 + i * step));
        break;
      }
    }
    limits.back() = std::numeric_limits<double>::infinity();

    mids.resize(nbins);
    for (int i = 0; i < nbins; ++i) {
      mids[i] = limits[2*i];
      limits[i] = limits[2*i+1];
    }
    limits.resize(nbins);
  }

  /// @brief Set up bins from reflection data using DataProxy interface.
  /// Automatically calculates 1/d^2 values from HKLs and unit cell.
  /// @tparam DataProxy Type with spacegroup(), unit_cell(), size(), stride(),
  ///                    get_hkl(), get_num() interface.
  /// @param nbins Number of resolution bins to create.
  /// @param method Binning scheme (EqualCount, Dstar, Dstar2, Dstar3).
  /// @param proxy Data proxy object (typically Intensities or MTZ).
  /// @param cell_ Unit cell for resolution calculations (nullptr uses proxy.unit_cell()).
  /// @param col_idx Column index for filtering NaN values (0 = skip filtering).
  template<typename DataProxy>
  void setup(int nbins, Method method, const DataProxy& proxy,
             const UnitCell* cell_=nullptr, size_t col_idx=0) {
    if (col_idx >= proxy.stride())
      fail("wrong col_idx in Binner::setup()");
    cell = cell_ ? *cell_ : proxy.unit_cell();
    std::vector<double> inv_d2;
    inv_d2.reserve(proxy.size() / proxy.stride());
    for (size_t offset = 0; offset < proxy.size(); offset += proxy.stride())
      if (col_idx == 0 || !std::isnan(proxy.get_num(offset + col_idx)))
        inv_d2.push_back(cell.calculate_1_d2(proxy.get_hkl(offset)));
    setup_from_1_d2(nbins, method, std::move(inv_d2), nullptr);
  }

  /// @brief Check that bins have been set up; throw if not.
  void ensure_limits_are_set() const {
    if (limits.empty())
      fail("Binner not set up");
  }

  /// @brief Find bin index for a given 1/d^2 value (binary search, generic).
  /// @param inv_d2 Squared reciprocal resolution (1/d^2).
  /// @return Bin index (0 to size()-1).
  int get_bin_from_1_d2(double inv_d2) {
    ensure_limits_are_set();
    auto it = std::lower_bound(limits.begin(), limits.end(), inv_d2);
    // it can't be limits.end() b/c limits.back() is +inf
    return int(it - limits.begin());
  }

  /// @brief Find bin index for a reflection by Miller indices.
  /// @param hkl Miller indices (h,k,l).
  /// @return Bin index (0 to size()-1).
  int get_bin(const Miller& hkl) {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_from_1_d2(inv_d2);
  }

  /// @brief Find bin index using a hint (fast path for sorted reflections).
  /// Updates hint to the found bin index for the next call.
  /// Assumes sorted reflections: each bin is same or next to previous.
  /// @param inv_d2 Squared reciprocal resolution (1/d^2).
  /// @param hint In/out: previous bin index (input), current bin index (output).
  /// @return Bin index (0 to size()-1).
  int get_bin_from_1_d2_hinted(double inv_d2, int& hint) const {
    if (inv_d2 <= limits[hint]) {
      while (hint != 0 && limits[hint-1] > inv_d2)
        --hint;
    } else {
      // limits.back() is +inf, so we won't overrun
      while (limits[hint] < inv_d2)
        ++hint;
    }
    return hint;
  }

  /// @brief Find bin index for a reflection using a hint (fast path).
  /// @param hkl Miller indices (h,k,l).
  /// @param hint In/out: previous bin index (input), current bin index (output).
  /// @return Bin index (0 to size()-1).
  int get_bin_hinted(const Miller& hkl, int& hint) const {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_from_1_d2_hinted(inv_d2, hint);
  }

  /// @brief Get bin indices for all reflections in a DataProxy.
  /// Uses hinting for fast processing of sorted data.
  /// @tparam DataProxy Type with size(), stride(), get_hkl() interface.
  /// @param proxy Data proxy object.
  /// @return Vector of bin indices (length = proxy.size() / proxy.stride()).
  template<typename DataProxy>
  std::vector<int> get_bins(const DataProxy& proxy) const {
    ensure_limits_are_set();
    int hint = 0;
    std::vector<int> nums(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < nums.size(); ++i, offset += proxy.stride())
      nums[i] = get_bin_hinted(proxy.get_hkl(offset), hint);
    return nums;
  }

  /// @brief Get bin indices for an array of 1/d^2 values.
  /// Uses hinting for fast processing of sorted data.
  /// @param inv_d2 Pointer to array of squared reciprocal resolutions.
  /// @param size Number of values in array.
  /// @return Vector of bin indices (length = size).
  std::vector<int> get_bins_from_1_d2(const double* inv_d2, size_t size) const {
    ensure_limits_are_set();
    int hint = 0;
    std::vector<int> nums(size);
    for (size_t i = 0; i < size; ++i)
      nums[i] = get_bin_from_1_d2_hinted(inv_d2[i], hint);
    return nums;
  }

  /// @brief Get bin indices for a vector of 1/d^2 values.
  /// @param inv_d2 Vector of squared reciprocal resolutions.
  /// @return Vector of bin indices (length = inv_d2.size()).
  std::vector<int> get_bins_from_1_d2(const std::vector<double>& inv_d2) const {
    return get_bins_from_1_d2(inv_d2.data(), inv_d2.size());
  }

  /// @brief Get minimum resolution (highest 1/d^2) of a bin.
  /// @param n Bin index (0 to size()-1).
  /// @return Minimum resolution d in Angstroms (highest angle, tightest spacing).
  double dmin_of_bin(int n) const {
    return 1. / std::sqrt(n == (int) size() - 1 ? max_1_d2 : limits.at(n));
  }
  /// @brief Get maximum resolution (lowest 1/d^2) of a bin.
  /// @param n Bin index (0 to size()-1).
  /// @return Maximum resolution d in Angstroms (lowest angle, loosest spacing).
  double dmax_of_bin(int n) const {
    return 1. / std::sqrt(n == 0 ? min_1_d2 : limits.at(n-1));
  }

  /// @brief Number of bins (resolution shells).
  size_t size() const { return limits.size(); }

  UnitCell cell;                  ///< Unit cell for calculating 1/d^2 from HKL
  double min_1_d2;               ///< Minimum 1/d^2 in data (lowest resolution)
  double max_1_d2;               ///< Maximum 1/d^2 in data (highest resolution)
  std::vector<double> limits;    ///< Upper limit (1/d^2) of each bin
  std::vector<double> mids;      ///< Midpoint (1/d^2) of each bin
};

} // namespace gemmi
#endif
