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
};

} // namespace gemmi
#endif
