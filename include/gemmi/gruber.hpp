// Copyright 2021 Global Phasing Ltd.
//
// Unit cell reduction to Buerger or Niggli cell. Uses G6 (Gruber) vector.

#ifndef GEMMI_GRUBER_HPP_
#define GEMMI_GRUBER_HPP_

#include <cmath>
#include <array>
#include "math.hpp"  // for deg

namespace gemmi {

// G6 vector ("G" for Gruber). Used in cell reduction algorithms.
// Originally, in B. Gruber, Acta Cryst. A29, 433 (1973), it was called
// "characteristic" of a lattice/cell.
struct GruberVector {
  //    a.a  b.b c.c 2b.c 2a.c 2a.b
  double A, B, C, xi, eta, zeta;  // the 1973 paper uses names A B C ξ η ζ

  GruberVector(const std::array<double,6>& g6)
    : A(g6[0]), B(g6[1]), C(g6[2]), xi(g6[3]), eta(g6[4]), zeta(g6[5]) {}

  std::array<double,6> parameters() { return {A, B, C, xi, eta, zeta}; }
  std::array<double,6> cell_parameters() const {
    // inverse of UnitCell::g6()
    double a = std::sqrt(A);
    double b = std::sqrt(B);
    double c = std::sqrt(C);
    return {a, b, c,
            deg(std::acos(xi/(2*b*c))),
            deg(std::acos(eta/(2*a*c))),
            deg(std::acos(zeta/(2*a*b)))};
  }
  std::array<double,6> niggli_parameters() { return {A, B, C, xi/2, eta/2, zeta/2}; }

  bool is_normalized() const {
    // eq(3) from Gruber 1973
    return A <= B && B && C &&
           (A != B || std::abs(xi) <= std::abs(eta)) &&
           (B != C || std::abs(eta) <= std::abs(zeta)) &&
           (xi > 0) == (eta > 0) && (xi > 0) == (zeta > 0);
  }

  bool is_buerger(double epsilon=0) const {
    return is_normalized() &&
           // eq (4) from Gruber 1973
           std::abs(xi) <= B + epsilon &&
           std::abs(eta) <= A + epsilon &&
           std::abs(zeta) <= A + epsilon;
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
    double prev_sum = -1;
    int stall_count = 0;
    for (;;) {
      ++n;
      normalize();
      if (std::abs(xi) > B) { // B2
        double j = std::floor(0.5*xi/B + 0.5);
        C += j * (j*B - xi);
        xi -= 2 * j * B;
        eta -= j * zeta;
      } else if (std::abs(eta) > A) { // B3
        double j = std::floor(0.5*eta/A + 0.5);
        C += j * (j*A - eta);
        xi -= j * zeta;
        eta -= 2 * j * A;
      } else if (std::abs(zeta) > A) { // B4
        double j = std::floor(0.5*zeta/A + 0.5);
        B += j * (j*A - zeta);
        xi -= j * eta;
        zeta -= 2 * j * A;
      } else if (xi + eta + zeta + A + B < 0) { // B5
        double j = std::floor(0.5 * (xi + eta) / (A + B + zeta) + 0.5);
        C += j * (j * (A + B + zeta) - (xi + eta));
        xi -= j * (2*B + zeta);
        eta -= j * (2*A + zeta);
      } else {
        break;
      }
      // In rare cases numerical errors push the algorithm into infinite loop,
      // as described in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
      // Ad-hoc solution: stop if a+b+c is stalled for 5 iterations.
      if (n > 8) {  // don't waste time during the first few iterations
        double sum = std::sqrt(A) + std::sqrt(B) + std::sqrt(C);
        if (std::abs(sum - prev_sum) < sum * 1e-6) {
          if (++stall_count == 5) {
            normalize();
            break;
          }
        } else {
          stall_count = 0;
        }
        prev_sum = sum;
      }
    }
    return n;
  }

  // Algorithm from Krivy & Gruber, Acta Cryst. (1976) A32, 297.
  int niggli_reduce() {
    int n = 0;
    for (;;) {
      ++n;
      normalize();
      if (n == 100)
        break;
      if (std::abs(xi) > B || (xi == B && 2 * eta < zeta) ||
                              (xi == -B && zeta < 0)) { // 5.
        double sign_xi = xi >= 0 ? 1 : -1;
        C += B - xi * sign_xi;
        eta -= zeta * sign_xi;
        xi -= 2 * B * sign_xi;
      } else if (std::abs(eta) > A || (eta == A && 2 * xi < zeta) ||
                                      (eta == -A && zeta < 0)) { // 6.
        double sign_eta = eta >= 0 ? 1 : -1;
        C += A - eta * sign_eta;
        xi -= zeta * sign_eta;
        eta -= 2 * A * sign_eta;
      } else if (std::abs(zeta) > A || (zeta == A && 2 * xi < eta) ||
                                       (zeta == -A && eta < 0)) { // 7.
        double sign_zeta = zeta >= 0 ? 1 : -1;
        B += A - zeta * sign_zeta;
        xi -= eta * sign_zeta;
        zeta -= 2 * A * sign_zeta;
      } else if (xi + eta + zeta + A + B < 0 ||
                 (xi + eta + zeta + A + B == 0 && 2 * (A + eta) + zeta > 0)) { // 8.
        C += A + B + xi + eta + zeta;
        xi += 2 * B + zeta;
        eta += 2 * A + zeta;
      } else {
        break;
      }
    }
    return n;
  }
};

} // namespace gemmi
#endif
