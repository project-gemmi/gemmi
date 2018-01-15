// Copyright 2017 Global Phasing Ltd.
//
// Macromolecular structure + symmetry operations

#ifndef GEMMI_CRYSTAL_HPP_
#define GEMMI_CRYSTAL_HPP_

#include <cmath>  // for INFINITY
#include "model.hpp"
#include "symmetry.hpp"

namespace gemmi {

struct Bounds {
  Position low = {INFINITY, INFINITY, INFINITY};
  Position high = {-INFINITY, -INFINITY, -INFINITY};

  void include(const Position& p) {
    if (p.x < low.x) low.x = p.x;
    if (p.y < low.y) low.y = p.y;
    if (p.z < low.z) low.z = p.z;
    if (p.x > high.x) high.x = p.x;
    if (p.y > high.y) high.y = p.y;
    if (p.z > high.z) high.z = p.z;
  }
  void include(const Atom& atom) { include(atom.pos); }
  template<typename T> void include(const T& obj) {
    for (const auto& child : obj.children())
      include(child);
  }

  void add(const Bounds& o) {
    include(o.low);
    include(o.high);
  }

  bool has_overlap(const Bounds& o) const {
    return (low.x < o.low.x) != (high.x < o.high.x) &&
           (low.y < o.low.y) != (high.y < o.high.y) &&
           (low.z < o.low.z) != (high.z < o.high.z);
  }

  void apply(const Op& op) {
    // TODO
    // apply to low and high, swap if low._ > high._
  }
};

struct Crystal : Structure {
  SpaceGroup space_group;

  Bounds calculate_xyz_bounds() const {
    Bounds b;
    b.include(*this);
    //TODO: NCS
    return b;
  }
  std::vector<Op> find_images_in_box(const Bounds& b) const {
    std::vector<Op> ops;
    //TODO
    return ops;
  }
  std::vector<Op> find_images_in_sphere(const Position& ctr, double r) const {
    std::vector<Op> ops;
    //TODO
    return ops;
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
