// Copyright 2022 Global Phasing Ltd.
//
// Binning - resolution shells for reflections.

#ifndef GEMMI_BINNER_HPP_
#define GEMMI_BINNER_HPP_

#include <vector>
#include <limits>        // for numeric_limits
#include "unitcell.hpp"  // for UnitCell

namespace gemmi {

struct Binner {
  enum class Method {
    EqualCount,
    Dstar,
    Dstar2,
    Dstar3,
  };

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

  void ensure_limits_are_set() const {
    if (limits.empty())
      fail("Binner not set up");
  }

  // Generic. Method-specific versions could be faster.
  int get_bin_from_1_d2(double inv_d2) {
    ensure_limits_are_set();
    auto it = std::lower_bound(limits.begin(), limits.end(), inv_d2);
    // it can't be limits.end() b/c limits.back() is +inf
    return int(it - limits.begin());
  }

  int get_bin(const Miller& hkl) {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_from_1_d2(inv_d2);
  }

  // We assume that bins are seeked mostly for sorted reflections,
  // so it's usually either the same bin as previously, or the next one.
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

  int get_bin_hinted(const Miller& hkl, int& hint) const {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_from_1_d2_hinted(inv_d2, hint);
  }

  template<typename DataProxy>
  std::vector<int> get_bins(const DataProxy& proxy) const {
    ensure_limits_are_set();
    int hint = 0;
    std::vector<int> nums(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < nums.size(); ++i, offset += proxy.stride())
      nums[i] = get_bin_hinted(proxy.get_hkl(offset), hint);
    return nums;
  }

  std::vector<int> get_bins_from_1_d2(const double* inv_d2, size_t size) const {
    ensure_limits_are_set();
    int hint = 0;
    std::vector<int> nums(size);
    for (size_t i = 0; i < size; ++i)
      nums[i] = get_bin_from_1_d2_hinted(inv_d2[i], hint);
    return nums;
  }

  std::vector<int> get_bins_from_1_d2(const std::vector<double>& inv_d2) const {
    return get_bins_from_1_d2(inv_d2.data(), inv_d2.size());
  }

  double dmin_of_bin(int n) const {
    return 1. / std::sqrt(n == (int) size() - 1 ? max_1_d2 : limits.at(n));
  }
  double dmax_of_bin(int n) const {
    return 1. / std::sqrt(n == 0 ? min_1_d2 : limits.at(n-1));
  }

  size_t size() const { return limits.size(); }

  UnitCell cell;
  double min_1_d2;
  double max_1_d2;
  std::vector<double> limits;  // upper limit of each bin
  std::vector<double> mids;    // the middle of each bin
};

} // namespace gemmi
#endif
