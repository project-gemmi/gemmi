// Copyright 2021 Global Phasing Ltd.
//
// Unit cell reduction to Buerger or Niggli cell.

#ifndef GEMMI_NIGGLI_HPP_
#define GEMMI_NIGGLI_HPP_

#include <cmath>
#include <array>

namespace gemmi {

// G6 vector ("G" for Gruber). Used in cell reduction algorithms.
// Originally, in B. Gruber, Acta Cryst. A29, 433 (1973), it was called
// "characteristic" of a lattice/cell.
struct GruberVector {
  //    a.a  b.b c.c 2b.c 2a.c 2a.b
  double A, B, C, xi, eta, zeta;  // the 1973 paper uses names A B C ξ η ζ

  GruberVector(const std::array<double,6>& g6)
    : A(g6[0]), B(g6[1]), C(g6[2]), xi(g6[3]), eta(g6[4]), zeta(g6[5]) {}

  bool is_normalized() const {
    // eq(3) from Gruber 1973
    return A <= B && B && C &&
           (A != B || std::abs(xi) <= std::abs(eta)) &&
           (B != C || std::abs(eta) <= std::abs(zeta)) &&
           (xi > 0) == (eta > 0) && (xi > 0) == (zeta > 0);
  }

  bool is_buerger() const {
    return is_normalized() &&
      // eq (4) from Gruber 1973
      std::abs(xi) <= B && std::abs(eta) <= A && std::abs(zeta) <= A;
  }

  // Algorithm N from Gruber (1973).
  // Returns branch taken in N3.
  bool normalize() {
    if (A > B || (A == B && std::abs(xi) > std::abs(eta))) { // N1
      std::swap(A, B);
      std::swap(xi, eta);
    }
    if (B > C || (B == C && std::abs(eta) > std::abs(zeta))) { // N2
      std::swap(B, C);
      std::swap(eta, zeta);
      if (A > B || (A == B && std::abs(xi) > std::abs(eta))) { // N1 again
        std::swap(A, B);
        std::swap(xi, eta);
      }
    }
    // N3
    bool cond = (xi * eta * zeta > 0);
    double sgn =  cond ? 1 : -1;
    xi = std::copysign(xi, sgn);
    eta = std::copysign(eta, sgn);
    zeta = std::copysign(zeta, sgn);
    return cond;
  }

  // Algorithm B from Gruber (1973).
  // Returns number of iterations.
  int buerger_reduce() {
    int n = 0;
    while (++n < 100) {
      normalize();
      // B2
      if (xi > B) {
        double j = std::floor(0.5*xi/B + 0.5);
        C += j * (j*B - xi);
        xi -= 2 * j * B;
        eta -= j * zeta;
        continue;
      }
      // B3
      if (eta > A) {
        double j = std::floor(0.5*eta/A + 0.5);
        C += j * (j*A - eta);
        xi -= j * zeta;
        eta -= 2 * j * A;
        continue;
      }
      // B4
      if (zeta > A) {
        double j = std::floor(0.5*zeta/A + 0.5);
        B += j * (j*A - zeta);
        xi -= j * eta;
        zeta -= 2 * j * A;
        continue;
      }
      // B5
      if (xi + eta + zeta + A + B < 0) {
        double j = std::floor(0.5 * (xi + eta) / (A + B + zeta) + 0.5);
        C += j * (j * (A + B + zeta) - (xi + eta));
        xi -= j * (2*B + zeta);
        eta -= j * (2*A + zeta);
        continue;
      }
      break;
    }
    return n;
  }
};

} // namespace gemmi
#endif
