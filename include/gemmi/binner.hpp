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
    Dstar3,
    Dstar2,
    Dstar,
    LogDstar,
    EqualCount,
    Refmac,
  };

  int setup_from_1_d2(int nbins, Method method, std::vector<double>&& dm2) {
    if (nbins < 2)
      fail("Binner: at least two bins are needed");
    if (dm2.empty())
      fail("Binner: no data");
    bin_limits.resize(nbins);
    if (method == Method::EqualCount) {
      std::sort(dm2.begin(), dm2.end());
      min_1_d2 = dm2.front();
      max_1_d2 = dm2.back();
    } else {
      min_1_d2 = max_1_d2 = dm2.front();
      for (double x : dm2) {
        if (x < min_1_d2)
          min_1_d2 = x;
        if (x > max_1_d2)
          max_1_d2 = x;
      }
    }
    if (method == Method::EqualCount) {
      double avg_count = double(dm2.size()) / nbins;
      for (int i = 1; i < nbins; ++i)
        bin_limits[i-1] = dm2[int(avg_count * i)];
    } else if (method == Method::Dstar2) {
      double step = (max_1_d2 - min_1_d2) / nbins;
      for (int i = 1; i < nbins; ++i)
        bin_limits[i-1] = min_1_d2 + i * step;
    } else if (method == Method::Dstar) {
      double min_1_d = std::sqrt(min_1_d2);
      double max_1_d = std::sqrt(max_1_d2);
      double step = (max_1_d - min_1_d) / nbins;
      for (int i = 1; i < nbins; ++i)
        bin_limits[i-1] = sq(min_1_d + i * step);
    } else if (method == Method::Dstar3) {
      double min_1_d3 = min_1_d2 * std::sqrt(min_1_d2);
      double max_1_d3 = max_1_d2 * std::sqrt(max_1_d2);
      double step = (max_1_d3 - min_1_d3) / nbins;
      for (int i = 1; i < nbins; ++i)
        bin_limits[i-1] = sq(std::cbrt(min_1_d3 + i * step));
    } else if (method == Method::LogDstar) {
      // TODO
    } else if (method == Method::Refmac) {
      // TODO
    }
    bin_limits.back() = std::numeric_limits<double>::infinity();
    return (int) bin_limits.size();
  }

  template<typename DataProxy>
  int setup(int nbins, Method method, const DataProxy& proxy) {
    cell = proxy.unit_cell();
    std::vector<double> dm2(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < dm2.size(); ++i, offset += proxy.stride())
      dm2[i] = cell.calculate_1_d2(proxy.get_hkl(offset));
    return setup_from_1_d2(nbins, method, std::move(dm2));
  }

  // Generic. Method-specific versions could be faster.
  int get_bin_number_from_1_d2(double inv_d2) {
    if (bin_limits.empty())
      fail("Binner not set up");
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
    if (bin_limits.empty())
      fail("Binner not set up");
    int hint = 0;
    std::vector<int> nums(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < nums.size(); ++i, offset += proxy.stride())
      nums[i] = get_bin_number_hinted(proxy.get_hkl(offset), hint);
    return nums;
  }

  size_t bin_count() const { return bin_limits.size(); }

  UnitCell cell;
  double min_1_d2;
  double max_1_d2;
  std::vector<double> bin_limits;  // upper limit of each bin
};

} // namespace gemmi
#endif
