// Copyright 2017 Global Phasing Ltd.
//
// Unit cell.

#ifndef GEMMI_UNITCELL_HH_
#define GEMMI_UNITCELL_HH_

#include <cmath>

namespace gemmi {
namespace mol {

struct Position {
  double x, y, z;
};

// if we need more involved math we may switch to using a third-party lib
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

struct UnitCell {
  double a = 1.0, b = 1.0, c = 1.0;
  double alpha = 90.0, beta = 90.0, gamma = 90.0;
  Matrix33 orth = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  Matrix33 frac = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

  void calculate_matrices() {
    double deg2rad = 3.1415926535897932384626433832795029 / 180.0;
    double cos_alpha = std::cos(deg2rad * alpha);
    double cos_beta = std::cos(deg2rad * beta);
    double cos_gamma = std::cos(deg2rad * gamma);
    double sin_alpha = std::sin(deg2rad * alpha);
    double sin_beta = std::sin(deg2rad * beta);
    double sin_gamma = std::sin(deg2rad * gamma);
    if (sin_alpha == 0 || sin_beta == 0 || sin_gamma == 0)
      throw std::runtime_error("Impossible angle - N*180deg.");
    double cos_alpha_star_sin_beta = (cos_beta * cos_gamma - cos_alpha) /
                                      sin_gamma;
    double cos_alpha_star = cos_alpha_star_sin_beta / sin_beta;
    double s1rca2 = std::sqrt(1.0 - cos_alpha_star * cos_alpha_star);
    double sb_s1rca2 = sin_beta * s1rca2;

    // The orthogonalization matrix we use is described in ITfC B p.262:
    // "An alternative mode of orthogonalization, used by the Protein
    // Data Bank and most programs, is to align the a1 axis of the unit
    // cell with the Cartesian X_1 axis, and to align the a*_3 axis with the
    // Cartesian X_3 axis."
    orth = {a,  b * cos_gamma,  c * cos_beta,
            0., b * sin_gamma, -c * cos_alpha_star_sin_beta,
            0., 0.           ,  c * sb_s1rca2};

    double o13 = -(cos_gamma * cos_alpha_star_sin_beta + cos_beta * sin_gamma)
                  / (sb_s1rca2 * sin_gamma * a);
    double o23 = cos_alpha_star / (s1rca2 * sin_gamma * b);
    frac = {1 / a,  -cos_gamma / (sin_gamma * a), o13,
            0.,     1 / (sin_gamma * b),          o23,
            0.,     0.,                           1 / (sb_s1rca2 * c)};
  }

  void set(double a_, double b_, double c_,
           double alpha_, double beta_, double gamma_) {
    a = a_;
    b = b_;
    c = c_;
    alpha = alpha_;
    beta = beta_;
    gamma = gamma_;
    calculate_matrices();
  }

  Position orthogonalize(const Position& f) const { return orth.multiply(f); }
  Position fractionalize(const Position& o) const { return frac.multiply(o); }
};

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
