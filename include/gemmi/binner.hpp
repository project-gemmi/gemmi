// Copyright 2022 Global Phasing Ltd.
//
// Binning - resolution shells for reflections.

#ifndef GEMMI_BINNER_HPP_
#define GEMMI_BINNER_HPP_

#include <vector>
#include <limits>        // for numeric_limits
#include <unordered_map> // for unordered_map
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
    limits.resize(nbins);
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
          limits[i-1] = inv_d2[int(avg_count * i)];
        break;
      }
      case Method::Dstar2: {
        double step = (max_1_d2 - min_1_d2) / nbins;
        for (int i = 1; i < nbins; ++i)
          limits[i-1] = min_1_d2 + i * step;
        break;
      }
      case Method::Dstar: {
        double min_1_d = std::sqrt(min_1_d2);
        double max_1_d = std::sqrt(max_1_d2);
        double step = (max_1_d - min_1_d) / nbins;
        for (int i = 1; i < nbins; ++i)
          limits[i-1] = sq(min_1_d + i * step);
        break;
      }
      case Method::Dstar3: {
        double min_1_d3 = min_1_d2 * std::sqrt(min_1_d2);
        double max_1_d3 = max_1_d2 * std::sqrt(max_1_d2);
        double step = (max_1_d3 - min_1_d3) / nbins;
        for (int i = 1; i < nbins; ++i)
          limits[i-1] = sq(std::cbrt(min_1_d3 + i * step));
        break;
      }
    }
    limits.back() = std::numeric_limits<double>::infinity();
    return (int) limits.size();
  }

  template<typename DataProxy>
  int setup(int nbins, Method method, const DataProxy& proxy,
            const UnitCell* cell_=nullptr) {
    cell = cell_ ? *cell_ : proxy.unit_cell();
    std::vector<double> inv_d2(proxy.size() / proxy.stride());
    for (size_t i = 0, offset = 0; i < inv_d2.size(); ++i, offset += proxy.stride())
      inv_d2[i] = cell.calculate_1_d2(proxy.get_hkl(offset));
    return setup_from_1_d2(nbins, method, std::move(inv_d2), nullptr);
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

  double dmin_of_bin(int n) const {
    return 1. / std::sqrt(limits.at(n));
  }
  double dmax_of_bin(int n) const {
    return 1. / std::sqrt(n == 0 ? min_1_d2 : limits.at(n-1));
  }

  size_t size() const { return limits.size(); }

  UnitCell cell;
  double min_1_d2;
  double max_1_d2;
  std::vector<double> limits;  // upper limit of each bin
};

// the result is in a
inline Correlation combine_two_correlations(const Correlation& a, const Correlation& b) {
  Correlation r;
  r.n = a.n + b.n;
  r.mean_x = (a.n * a.mean_x + b.n * b.mean_x) / r.n;
  r.mean_y = (a.n * a.mean_y + b.n * b.mean_y) / r.n;
  r.sum_xx = a.sum_xx + a.n * sq(a.mean_x - r.mean_x)
           + b.sum_xx + b.n * sq(b.mean_x - r.mean_x);
  r.sum_yy = a.sum_yy + a.n * sq(a.mean_y - r.mean_y)
           + b.sum_yy + b.n * sq(b.mean_y - r.mean_y);
  r.sum_xy = a.sum_xy + a.n * (a.mean_x - r.mean_x) * (a.mean_y - r.mean_y)
           + b.sum_xy + b.n * (b.mean_x - r.mean_x) * (b.mean_y - r.mean_y);
  return r;
}

inline Correlation combine_correlations(const std::vector<Correlation>& cors) {
  Correlation result;
  for (const Correlation& cor : cors)
    result = combine_two_correlations(result, cor);
  return result;
}

struct HklMatch {
  std::vector<int> pos;
  size_t hkl_size;

  HklMatch(const Miller* hkl, size_t hkl_size_, const Miller* ref, size_t ref_size)
      : pos(ref_size, -1), hkl_size(hkl_size_) {
    // Usually, both datasets are sorted. This make things faster.
    if (std::is_sorted(hkl, hkl + hkl_size) &&
        std::is_sorted(ref, ref + ref_size)) {
      // cf. for_matching_reflections()
      auto a = hkl;
      auto b = ref;
      while (a != hkl + hkl_size && b != ref + ref_size) {
        if (*a == *b)
          pos[b++ - ref] = static_cast<int>(a++ - hkl);
        else if (*a < *b)
          ++a;
        else
          ++b;
      }
    } else {
      std::unordered_map<Miller, int, MillerHash> hkl_index;
      for (int i = 0; i != (int)hkl_size; ++i)
        hkl_index.emplace(hkl[i], i);
      for (size_t i = 0; i != ref_size; ++i) {
        auto it = hkl_index.find(ref[i]);
        if (it != hkl_index.end())
          pos[i] = it->second;
      }
    }
  }

  HklMatch(const std::vector<Miller>& hkl, const std::vector<Miller>& ref)
    : HklMatch(hkl.data(), hkl.size(), ref.data(), ref.size()) {}

  template <typename T> std::vector<T> aligned(const std::vector<T>& v, T nan) {
    if (v.size() != hkl_size)
      fail("HklMatch.aligned(): wrong data, size differs");
    std::vector<T> result(pos.size());
    for (size_t i = 0; i != pos.size(); ++i)
      result[i] = pos[i] >= 0 ? v[pos[i]] : nan;
    return result;
  }
};

} // namespace gemmi
#endif
