// Copyright 2020 Global Phasing Ltd.
//
// AsuData for storing reflection data

#ifndef GEMMI_ASUDATA_HPP_
#define GEMMI_ASUDATA_HPP_

#include <complex>       // for arg, abs
#include <tuple>         // for tie
#include <algorithm>     // for sort, is_sorted
#include "unitcell.hpp"
#include "symmetry.hpp"

namespace gemmi {

template<typename T>
struct HklValue {
  Miller hkl;
  T value;

  bool operator<(const Miller& m) const {
    return std::tie(hkl[0], hkl[1], hkl[2]) < std::tie(m[0], m[1], m[2]);
  }
  bool operator<(const HklValue& o) const { return operator<(o.hkl); }
};


template<typename T>
struct AsuData {
  std::vector<HklValue<T>> v;
  UnitCell unit_cell_;
  const SpaceGroup* spacegroup_ = nullptr;
  // function defining FPhiProxy interface
  size_t stride() const { return 1; }
  size_t size() const { return v.size(); }
  Miller get_hkl(size_t n) const { return v[n].hkl; }
  double get_f(size_t n) const { return std::abs(v[n].value); }
  double get_phi(size_t n) const { return std::arg(v[n].value); }
  const UnitCell& unit_cell() const { return unit_cell_; }
  const SpaceGroup* spacegroup() const { return spacegroup_; }
  void ensure_sorted() {
    if (!std::is_sorted(v.begin(), v.end()))
      std::sort(v.begin(), v.end());
  }

  template<int N, typename DataProxy>
  void load_values(const DataProxy& proxy, const std::array<std::string,N>& labels) {
    std::array<std::size_t, N> cols;
    for (int i = 0; i < N; ++i)
      cols[i] = proxy.column_index(labels[i]);
    unit_cell_ = proxy.unit_cell();
    spacegroup_ = proxy.spacegroup();
    using Val = typename T::value_type;
    for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
      std::array<Val, N> nums;
      for (int j = 0; j < N; ++j)
        nums[j] = proxy.get_num(i + cols[j]);
      if (!std::any_of(nums.begin(), nums.end(), [](Val f) { return std::isnan(f); })) {
        v.emplace_back();
        v.back().hkl = proxy.get_hkl(i);
        set_value_from_array(v.back().value, nums);
      }
    }
  }

private:
  // for T being a number, std::array and std::complex, respectively:
  static void set_value_from_array(T& val, const std::array<T,1>& nums) { val = nums[0]; }
  static void set_value_from_array(T& val, const T& nums) { val = nums; }
  template<typename R>
  static void set_value_from_array(std::complex<R>& val, const std::array<R,2>& nums) {
    val = std::polar(nums[0], gemmi::rad(nums[1]));
  }
};

} // namespace gemmi
#endif
