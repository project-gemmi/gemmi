// Copyright 2017 Global Phasing Ltd.
//
// Unit cell.

#ifndef GEMMI_UNITCELL_HPP_
#define GEMMI_UNITCELL_HPP_

#include <cmath>      // for cos, sin, sqrt, floor
#include <linalg.h>
#include "util.hpp"

namespace gemmi {

struct Position {
  double x, y, z;
  constexpr double operator[](int i) const { return (&x)[i]; }
  double& operator[](int i) { return (&x)[i]; }
  Position& wrap_to_unit() {
    x -= std::floor(x);
    y -= std::floor(y);
    z -= std::floor(z);
    return *this;
  }
};

struct Matrix33 {
  double a11, a12, a13,
         a21, a22, a23,
         a31, a32, a33;

  Position multiply(const Position& p) const {
    return {a11 * p.x  + a12 * p.y  + a13 * p.z,
          /*a21 * p.x*/+ a22 * p.y  + a23 * p.z,
          /*a31 * p.x  + a32 * p.y*/+ a33 * p.z};
  }
};

// discussion: https://stackoverflow.com/questions/20305272/
inline double calculate_dihedral(const Position& p0, const Position& p1,
                                 const Position& p2, const Position& p3) {
  typedef linalg::vec<double,3> double3;
  double3 b0(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
  double3 b1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
  double3 b2(p3.x - p2.x, p3.y - p2.y, p3.z - p2.z);
  double3 u = linalg::cross(b1, b0);
  double3 w = linalg::cross(b2, b1);
  double y = linalg::dot(linalg::cross(u, w), b1);
  double x = linalg::dot(u, w) * linalg::length(b1);
  return std::atan2(y, x);
}


struct UnitCell {
  double a = 1.0, b = 1.0, c = 1.0;
  double alpha = 90.0, beta = 90.0, gamma = 90.0;
  constexpr double operator[](int i) const { return (&a)[i]; }
  Matrix33 orth = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  Matrix33 frac = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  Position shift = {0., 0., 0.};
  // volume and reciprocal parameters a*, b*, c*, alpha*, beta*, gamma*
  double volume = 1.0;
  double ar = 1.0, br = 1.0, cr = 1.0;
  double cos_alphar = 0.0, cos_betar = 0.0, cos_gammar = 0.0;
  bool explicit_matrices = false;

  void calculate_properties() {
    constexpr double deg2rad = 3.1415926535897932384626433832795029 / 180.0;
    // ensure exact values for right angles
    double cos_alpha = alpha == 90. ? 0. : std::cos(deg2rad * alpha);
    double cos_beta  = beta  == 90. ? 0. : std::cos(deg2rad * beta);
    double cos_gamma = gamma == 90. ? 0. : std::cos(deg2rad * gamma);
    double sin_alpha = alpha == 90. ? 1. : std::sin(deg2rad * alpha);
    double sin_beta  = beta  == 90. ? 1. : std::sin(deg2rad * beta);
    double sin_gamma = gamma == 90. ? 1. : std::sin(deg2rad * gamma);
    if (sin_alpha == 0 || sin_beta == 0 || sin_gamma == 0)
      gemmi::fail("Impossible angle - N*180deg.");

		// volume - formula from Giacovazzo p.62
    volume = a * b * c * sqrt(1 - cos_alpha * cos_alpha - cos_beta * cos_beta
                              - cos_gamma * cos_gamma
                              + 2 * cos_alpha * cos_beta * cos_gamma);

    // reciprocal parameters a*, b*, ... (Giacovazzo, p. 64)
    ar = b * c * sin_alpha / volume;
    br = a * c * sin_beta / volume;
    cr = a * b * sin_gamma / volume;
    double cos_alphar_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
    cos_alphar = cos_alphar_sin_beta / sin_beta;
    //cos_alphar = (cos_beta * cos_gamma - cos_alpha) / (sin_beta * sin_gamma);
    cos_betar = (cos_alpha * cos_gamma - cos_beta) / (sin_alpha * sin_gamma);
    cos_gammar = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta);

    if (explicit_matrices)
      return;

    // The orthogonalization matrix we use is described in ITfC B p.262:
    // "An alternative mode of orthogonalization, used by the Protein
    // Data Bank and most programs, is to align the a1 axis of the unit
    // cell with the Cartesian X_1 axis, and to align the a*_3 axis with the
    // Cartesian X_3 axis."
    double sin_alphar = std::sqrt(1.0 - cos_alphar * cos_alphar);
    orth = {a,  b * cos_gamma,  c * cos_beta,
            0., b * sin_gamma, -c * cos_alphar_sin_beta,
            0., 0.           ,  c * sin_beta * sin_alphar};

    double o12 = -cos_gamma / (sin_gamma * a);
    double o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
                  / (sin_alphar * sin_beta * sin_gamma * a);
    double o23 = cos_alphar / (sin_alphar * sin_gamma * b);
    frac = {1 / a,  o12,           o13,
            0.,     1 / orth.a22,  o23,
            0.,     0.,            1 / orth.a33};
  }

  void set_matrices_from_fract(const linalg::mat<double,4,4>& fract) {
    frac = {fract.x.x, fract.y.x, fract.z.x,
            fract.x.y, fract.y.y, fract.z.y,
            fract.x.z, fract.y.z, fract.z.z};
    shift = {fract.w.x, fract.w.y, fract.w.z};
    auto ortho = linalg::inverse(fract);
    orth = {ortho.x.x, ortho.y.x, ortho.z.x,
            ortho.x.y, ortho.y.y, ortho.z.y,
            ortho.x.z, ortho.y.z, ortho.z.z};
    explicit_matrices = true;
  }

  void set(double a_, double b_, double c_,
           double alpha_, double beta_, double gamma_) {
    if (gamma_ == 0.0)  // ignore empty/partial CRYST1 (example: 3iyp)
      return;
    a = a_;
    b = b_;
    c = c_;
    alpha = alpha_;
    beta = beta_;
    gamma = gamma_;
    calculate_properties();
  }

  Position orthogonalize(const Position& f) const { return orth.multiply(f); }
  Position fractionalize(const Position& o) const { return frac.multiply(o); }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
