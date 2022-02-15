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

  int setup_from_1_d2(int nbins, Method method, std::vector<double>&& inv_d2,
                      const UnitCell* cell_) {
    if (nbins < 1)
      fail("Binner: nbins argument must be positive");
    if (inv_d2.empty())
      fail("Binner: no data");
    if (cell_)
      cell = *cell_;
    if (!cell.is_crystal())
      fail("Binner: unknown unit cell");
    bin_limits.resize(nbins);
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
        double avg_count = double(inv_d2.size()) / nbins;
        for (int i = 1; i < nbins; ++i)
          bin_limits[i-1] = inv_d2[int(avg_count * i)];
        break;
      }
      case Method::Dstar2: {
        double step = (max_1_d2 - min_1_d2) / nbins;
        for (int i = 1; i < nbins; ++i)
          bin_limits[i-1] = min_1_d2 + i * step;
        break;
      }
      case Method::Dstar: {
        double min_1_d = std::sqrt(min_1_d2);
        double max_1_d = std::sqrt(max_1_d2);
        double step = (max_1_d - min_1_d) / nbins;
        for (int i = 1; i < nbins; ++i)
          bin_limits[i-1] = sq(min_1_d + i * step);
        break;
      }
      case Method::Dstar3: {
        double min_1_d3 = min_1_d2 * std::sqrt(min_1_d2);
        double max_1_d3 = max_1_d2 * std::sqrt(max_1_d2);
        double step = (max_1_d3 - min_1_d3) / nbins;
        for (int i = 1; i < nbins; ++i)
          bin_limits[i-1] = sq(std::cbrt(min_1_d3 + i * step));
        break;
      }
    }
    bin_limits.back() = std::numeric_limits<double>::infinity();
    return (int) bin_limits.size();
  }

  template<typename DataProxy>
  int setup(int nbins, Method method, const DataProxy& proxy, const UnitCell* cell_) {
    cell = cell_ ? *cell_ : proxy.unit_cell();
    std::vector<double> inv_d2(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < inv_d2.size(); ++i, offset += proxy.stride())
      inv_d2[i] = cell.calculate_1_d2(proxy.get_hkl(offset));
    return setup_from_1_d2(nbins, method, std::move(inv_d2), nullptr);
  }

  void ensure_limits_are_set() const {
    if (bin_limits.empty())
      fail("Binner not set up");
  }

  // Generic. Method-specific versions could be faster.
  int get_bin_number_from_1_d2(double inv_d2) {
    ensure_limits_are_set();
    auto it = std::lower_bound(bin_limits.begin(), bin_limits.end(), inv_d2);
    // it can't be bin_limits.end() b/c bin_limits.back() is +inf
    return int(it - bin_limits.begin());
  }

  int get_bin_number(const Miller& hkl) {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_number_from_1_d2(inv_d2);
  }

  // We assume that the bin number is seeked mostly for sorted reflections,
  // so it's usually either the same bin as previously, or the next one.
  int get_bin_number_from_1_d2_hinted(double inv_d2, int& hint) const {
    if (inv_d2 <= bin_limits[hint]) {
      while (hint != 0 && bin_limits[hint-1] > inv_d2)
        --hint;
    } else {
      // bin_limits.back() is +inf, so we won't overrun
      while (bin_limits[hint] < inv_d2)
        ++hint;
    }
    return hint;
  }

  int get_bin_number_hinted(const Miller& hkl, int& hint) const {
    double inv_d2 = cell.calculate_1_d2(hkl);
    return get_bin_number_from_1_d2_hinted(inv_d2, hint);
  }

  template<typename DataProxy>
  std::vector<int> get_bin_numbers(const DataProxy& proxy) const {
    ensure_limits_are_set();
    int hint = 0;
    std::vector<int> nums(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < nums.size(); ++i, offset += proxy.stride())
      nums[i] = get_bin_number_hinted(proxy.get_hkl(offset), hint);
    return nums;
  }

  double get_bin_dmin(int n) const {
    return 1. / std::sqrt(n == 0 ? min_1_d2 : bin_limits.at(n-1));
  }
  double get_bin_dmax(int n) const {
    return 1. / std::sqrt(bin_limits.at(n));
  }

  size_t bin_count() const { return bin_limits.size(); }

  UnitCell cell;
  double min_1_d2;
  double max_1_d2;
  std::vector<double> bin_limits;  // upper limit of each bin
};

} // namespace gemmi
#endif
