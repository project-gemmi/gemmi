//! @file
//! @brief AsuData for storing reflection data.

// Copyright 2020 Global Phasing Ltd.
//
// AsuData for storing reflection data.

#ifndef GEMMI_ASUDATA_HPP_
#define GEMMI_ASUDATA_HPP_

#include <complex>       // for arg, abs
#include <algorithm>     // for sort, is_sorted
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "stats.hpp"     // for Correlation
#include "util.hpp"      // for vector_remove_if

namespace gemmi {

//! @brief Complex correlation calculator.
struct ComplexCorrelation {
  int n = 0;  //!< Number of points
  double sum_xx = 0.;  //!< Sum of |x-mean_x|^2
  double sum_yy = 0.;  //!< Sum of |y-mean_y|^2
  std::complex<double> sum_xy = 0.;  //!< Sum of (x-mean_x)*conj(y-mean_y)
  std::complex<double> mean_x = 0.;  //!< Mean of x
  std::complex<double> mean_y = 0.;  //!< Mean of y

  //! @brief Add data point.
  //! @param x X value (complex)
  //! @param y Y value (complex)
  void add_point(std::complex<double> x, std::complex<double> y) {
    ++n;
    double inv_n = 1.0 / n;
    double weight = (n - 1.0) * inv_n;
    std::complex<double> dx = x - mean_x;
    std::complex<double> dy = y - mean_y;
    sum_xx += weight * std::norm(dx);
    sum_yy += weight * std::norm(dy);
    sum_xy += weight * (dx * std::conj(dy));
    mean_x += dx * inv_n;
    mean_y += dy * inv_n;
  }
  //! @brief Add data point (float version).
  //! @param x X value (complex float)
  //! @param y Y value (complex float)
  void add_point(std::complex<float> x, std::complex<float> y) {
    add_point(std::complex<double>(x), std::complex<double>(y));
  }

  //! @brief Get correlation coefficient.
  //! @return Complex correlation coefficient
  std::complex<double> coefficient() const { return sum_xy / std::sqrt(sum_xx * sum_yy); }

  //! @brief Get ratio of means.
  //! @return |mean_y| / |mean_x|
  double mean_ratio() const { return std::abs(mean_y) / std::abs(mean_x); }
};


//! @brief Apply function to matching reflections in two sorted lists.
//! @tparam Func Function type
//! @tparam T HklValue type
//! @param a First sorted vector
//! @param b Second sorted vector
//! @param func Function to apply to matching pairs
//!
//! Precondition: a and b are sorted.
template<typename Func, typename T>
void for_matching_reflections(const std::vector<T>& a,
                              const std::vector<T>& b,
                              const Func& func) {
  auto r1 = a.begin();
  auto r2 = b.begin();
  while (r1 != a.end() && r2 != b.end()) {
    if (r1->hkl == r2->hkl) {
      func(*r1, *r2);
      ++r1;
      ++r2;
    } else if (r1->hkl < r2->hkl) {
      ++r1;
    } else {
      ++r2;
    }
  }
}

//! @brief Calculate correlation between reflection values.
//! @tparam T HklValue type
//! @param a First sorted vector
//! @param b Second sorted vector
//! @return Correlation object
//!
//! Precondition: a and b are sorted.
template<typename T>
Correlation calculate_hkl_value_correlation(const std::vector<T>& a,
                                            const std::vector<T>& b) {
  Correlation cor;
  for_matching_reflections(a, b, [&cor](const T& x, const T& y) {
      cor.add_point(x.value, y.value);
  });
  return cor;
}

//! @brief Calculate complex correlation between reflections.
//! @tparam T HklValue type with complex values
//! @param a First sorted vector
//! @param b Second sorted vector
//! @return ComplexCorrelation object
//!
//! Precondition: a and b are sorted.
template<typename T>
ComplexCorrelation calculate_hkl_complex_correlation(const std::vector<T>& a,
                                                     const std::vector<T>& b) {
  ComplexCorrelation cor;
  for_matching_reflections(a, b, [&cor](const T& x, const T& y) {
      cor.add_point(x.value, y.value);
  });
  return cor;
}

//! @brief Count matching reflections with equal values.
//! @tparam T HklValue type
//! @param a First sorted vector
//! @param b Second sorted vector
//! @return Number of matching reflections with equal values
//!
//! Precondition: a and b are sorted.
template<typename T>
int count_equal_values(const std::vector<T>& a, const std::vector<T>& b) {
  int count = 0;
  for_matching_reflections(a, b, [&count](const T& x, const T& y) {
      if (x.value == y.value)
        ++count;
  });
  return count;
}

//! @brief HKL indices with associated value.
//! @tparam T Value type
template<typename T>
struct HklValue {
  alignas(8) Miller hkl;  //!< Miller indices
  T value;  //!< Associated value

  bool operator<(const Miller& m) const { return hkl < m; }  //!< Compare with Miller indices
  bool operator<(const HklValue& o) const { return operator<(o.hkl); }  //!< Compare with another HklValue
};

//! @brief Value with standard uncertainty.
//! @tparam T Value type
template<typename T>
struct ValueSigma {
  using value_type = T;  //!< Value type
  T value, sigma;  //!< Value and standard uncertainty

  //! @brief Equality comparison.
  //! @param o Other ValueSigma
  //! @return True if both value and sigma are equal
  bool operator==(const ValueSigma& o) const {
    return value == o.value && sigma == o.sigma;
  }
};

namespace impl {
template<typename T>
void move_to_asu(const GroupOps&, const Miller& hkl, int, HklValue<T>& hkl_value) {
  hkl_value.hkl = hkl;
}

template<typename R>
void move_to_asu(const GroupOps& gops, const Miller& hkl, int isym,
                 HklValue<std::complex<R>>& v) {
  v.hkl = hkl;
  // cf. Mtz::ensure_asu()
  const Op& op = gops.sym_ops[(isym - 1) / 2];
  double shift = op.phase_shift(hkl);
  if (shift != 0) {
    double phase = std::arg(v.value) + shift;
    v.value = std::polar(std::abs(v.value), (R)phase);
  }
  if (isym % 2 == 0)
    v.value.imag(-v.value.imag());
}
} // namespace impl

//! @brief Reflection data in asymmetric unit.
//! @tparam T Value type
template<typename T>
struct AsuData {
  std::vector<HklValue<T>> v;  //!< HKL values
  UnitCell unit_cell_;  //!< Unit cell
  const SpaceGroup* spacegroup_ = nullptr;  //!< Space group pointer

  // Functions defining FPhiProxy interface
  size_t stride() const { return 1; }  //!< Get stride (always 1)
  size_t size() const { return v.size(); }  //!< Get number of reflections
  Miller get_hkl(size_t n) const { return v[n].hkl; }  //!< Get Miller indices
  double get_f(size_t n) const { return std::abs(v[n].value); }  //!< Get amplitude
  double get_phi(size_t n) const { return std::arg(v[n].value); }  //!< Get phase
  const UnitCell& unit_cell() const { return unit_cell_; }  //!< Get unit cell
  const SpaceGroup* spacegroup() const { return spacegroup_; }  //!< Get space group

  //! @brief Ensure data is sorted by HKL.
  void ensure_sorted() {
    if (!std::is_sorted(v.begin(), v.end()))
      std::sort(v.begin(), v.end());
  }

  //! @brief Move all reflections to asymmetric unit.
  //! @param tnt_asu If true, use TNT ASU convention
  void ensure_asu(bool tnt_asu=false) {
    if (!spacegroup_)
      fail("AsuData::ensure_asu(): space group not set");
    GroupOps gops = spacegroup_->operations();
    ReciprocalAsu asu(spacegroup_, tnt_asu);
    for (HklValue<T>& hkl_value : v) {
      const Miller& hkl = hkl_value.hkl;
      if (asu.is_in(hkl))
        continue;
      auto result = asu.to_asu(hkl, gops);
      impl::move_to_asu(gops, result.first, result.second, hkl_value);
    }
  }

  //! @brief Load values from one column.
  //! @tparam DataProxy Data proxy type
  //! @param proxy Data proxy
  //! @param label Column label
  //! @param as_is If false, ensure ASU and sort
  template<typename DataProxy>
  void load_values(const DataProxy& proxy, const std::string& label,
                   bool as_is=false) {
    std::size_t col = proxy.column_index(label);
    unit_cell_ = proxy.unit_cell();
    spacegroup_ = proxy.spacegroup();
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      auto num = proxy.get_num(i + col);
      if (!std::isnan(num))
        v.push_back({proxy.get_hkl(i), (T)num});
    }
    if (!as_is) {
      ensure_asu();
      ensure_sorted();
    }
  }

  //! @brief Load values from multiple columns.
  //! @tparam N Number of columns
  //! @tparam DataProxy Data proxy type
  //! @param proxy Data proxy
  //! @param labels Column labels
  //! @param as_is If false, ensure ASU and sort
  template<int N, typename DataProxy>
  void load_values(const DataProxy& proxy, const std::array<std::string,N>& labels,
                   bool as_is=false) {
    std::array<std::size_t, N> cols;
    for (int i = 0; i < N; ++i)
      cols[i] = proxy.column_index(labels[i]);
    unit_cell_ = proxy.unit_cell();
    spacegroup_ = proxy.spacegroup();
    using Val = typename T::value_type;
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      std::array<Val, N> nums;
      for (int j = 0; j < N; ++j)
        nums[j] = (Val) proxy.get_num(i + cols[j]);
      if (!std::any_of(nums.begin(), nums.end(), [](Val f) { return std::isnan(f); })) {
        v.emplace_back();
        v.back().hkl = proxy.get_hkl(i);
        set_value_from_array(v.back().value, nums);
      }
    }
    if (!as_is) {
      ensure_asu();
      ensure_sorted();
    }
  }

private:
  // for T being a number, std::array and std::complex, respectively:
  static void set_value_from_array(T& val, const std::array<T,1>& nums) { val = nums[0]; }
  static void set_value_from_array(T& val, const T& nums) { val = nums; }
  template<typename R>
  static void set_value_from_array(std::complex<R>& val, const std::array<R,2>& nums) {
    R theta = (R)rad(nums[1]);
    val = {nums[0] * std::cos(theta), nums[0] * std::sin(theta)};
  }
  template<typename R>
  static void set_value_from_array(ValueSigma<R>& val, const std::array<R,2>& nums) {
    val.value = nums[0];
    val.sigma = nums[1];
  }
};

//! @brief Create AsuData from multiple columns.
//! @tparam T Value type
//! @tparam N Number of columns
//! @tparam Data Data source type
//! @param data Data source
//! @param labels Column labels
//! @param as_is If false, ensure ASU and sort
//! @return AsuData object
template<typename T, int N, typename Data>
AsuData<T> make_asu_data(const Data& data, const std::array<std::string,N>& labels,
                         bool as_is=false) {
  AsuData<T> asu_data;
  asu_data.template load_values<N>(data_proxy(data), labels, as_is);
  return asu_data;
}

//! @brief Create AsuData from one column.
//! @tparam T Value type
//! @tparam Data Data source type
//! @param data Data source
//! @param label Column label
//! @param as_is If false, ensure ASU and sort
//! @return AsuData object
template<typename T, typename Data>
AsuData<T> make_asu_data(const Data& data, const std::string& label, bool as_is) {
  AsuData<T> asu_data;
  asu_data.load_values(data_proxy(data), label, as_is);
  return asu_data;
}

//! @brief Discard weak reflections by sigma ratio.
//! @tparam T Value type
//! @param asu_data ASU data (modified in place)
//! @param cutoff F/sigma cutoff
//!
//! Retains only points with positive SIGF and F/SIGF > cutoff.
template<typename T>
void discard_by_sigma_ratio(AsuData<ValueSigma<T>>& asu_data, double cutoff) {
  vector_remove_if(asu_data.v, [cutoff](const HklValue<ValueSigma<T>>& p) {
      auto& v = p.value;
      return v.sigma <= 0 || v.value <= cutoff * v.sigma;
  });
}

} // namespace gemmi
#endif
