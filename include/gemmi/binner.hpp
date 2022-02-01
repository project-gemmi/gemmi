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

  int setup_from_1_d2(int nbins, Method method, std::vector<double> dm2) {
    if (nbins < 2)
      fail("Binner: at least two bins are needed");
    if (dm2.empty())
      fail("Binner: no data");
    bins.resize(nbins);
    std::sort(dm2.begin(), dm2.end());
    min_1_d2 = dm2.front();
    max_1_d2 = dm2.back();
    if (method == Method::EqualCount) {
      double avg_count = double(dm2.size()) / nbins;
      for (int i = 1; i < nbins; ++i)
        bins[i-1] = dm2[int(avg_count * i)];
    } else if (method == Method::Dstar2) {
      double step = (dm2.back() - dm2.front()) / nbins;
      for (int i = 1; i < nbins; ++i)
        bins[i-1] = i * step;
    } else if (method == Method::Dstar3) {
    } else if (method == Method::Dstar) {
    } else if (method == Method::LogDstar) {
    } else if (method == Method::Refmac) {
    }
    bins.back() = std::numeric_limits<double>::infinity();
    return (int) bins.size();
  }

  template<typename DataProxy>
  int setup(int nbins, Method method, const DataProxy& proxy) {
    cell = proxy.unit_cell();
    std::vector<double> dm2(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < dm2.size(); ++i, offset += proxy.stride())
      dm2[i] = cell.calculate_1_d2(proxy.get_hkl(offset));
    return setup_from_1_d2(nbins, method, dm2);
  }

  // Generic. Method-specific versions could be faster.
  int get_bin_number(const Miller& hkl) {
    if (bins.empty())
      fail("Binner not set up");
    double inv_d2 = cell.calculate_1_d2(hkl);
    auto it = std::lower_bound(bins.begin(), bins.end(), inv_d2);
    // it can't be bins.end() b/c bins.back() is +inf
    return int(it - bins.begin());
  }

  // We assume that the bin number is seeked mostly for sorted reflections,
  // so it's usually either the same bin as previously, or the next one.
  int get_bin_number_hinted(const Miller& hkl, int& hint) const {
    double inv_d2 = cell.calculate_1_d2(hkl);
    if (inv_d2 <= bins[hint]) {
      while (hint != 0 && bins[hint-1] > inv_d2)
        --hint;
    } else {
      while (hint + 1 < (int)bins.size() && bins[++hint] < inv_d2) {}
    }
    return hint;
  }

  template<typename DataProxy>
  std::vector<int> get_bin_numbers(const DataProxy& proxy) const {
    if (bins.empty())
      fail("Binner not set up");
    int hint = (int)bins.size() - 1;
    std::vector<int> nums(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < nums.size(); ++i, offset += proxy.stride())
      nums[i] = get_bin_number_hinted(proxy.get_hkl(offset), hint);
    return nums;
  }

  size_t bin_count() const { return bins.size(); }

  UnitCell cell;
  double min_1_d2;
  double max_1_d2;
  std::vector<double> bins;  // upper limit of each bin
};

} // namespace gemmi
#endif
