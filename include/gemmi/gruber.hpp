// Copyright 2021 Global Phasing Ltd.
//
// Unit cell reduction to Buerger or Niggli cell. Uses G6 (Gruber) vector.

#ifndef GEMMI_GRUBER_HPP_
#define GEMMI_GRUBER_HPP_

#include <cmath>
#include <array>
#include "math.hpp"  // for deg

namespace gemmi {

// GruberVector contains G6 vector (G for Gruber) and cell reduction algorithms.
// Originally, in B. Gruber, Acta Cryst. A29, 433 (1973), the vector was called
// "characteristic" of a lattice/cell.
// Functions that take epsilon as a parameter use it for comparisons,
// as proposed in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
struct GruberVector {
  //    a.a  b.b c.c 2b.c 2a.c 2a.b
  double A, B, C, xi, eta, zeta;  // the 1973 paper uses names A B C ξ η ζ

  GruberVector(const std::array<double,6>& g6)
    : A(g6[0]), B(g6[1]), C(g6[2]), xi(g6[3]), eta(g6[4]), zeta(g6[5]) {}

  std::array<double,6> parameters() const { return {A, B, C, xi, eta, zeta}; }
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

  bool is_normalized() const {
    // eq(3) from Gruber 1973
    return A <= B && B <= C &&
           (A != B || std::abs(xi) <= std::abs(eta)) &&
           (B != C || std::abs(eta) <= std::abs(zeta)) &&
           (xi > 0) == (eta > 0) && (xi > 0) == (zeta > 0);
  }

  bool is_buerger(double epsilon=1e-9) const {
    return is_normalized() &&
           // eq (4) from Gruber 1973
           std::abs(xi) <= B + epsilon &&
           std::abs(eta) <= A + epsilon &&
           std::abs(zeta) <= A + epsilon;
  }

  // Algorithm N from Gruber (1973).
  // Returns branch taken in N3.
  void normalize(double eps=1e-9) {
    if (A - B > eps || (A - B >= -eps && std::abs(xi) > std::abs(eta) + eps)) { // N1
      std::swap(A, B);
      std::swap(xi, eta);
    }
    if (B - C > eps || (B - C >= -eps && std::abs(eta) > std::abs(zeta) + eps)) { // N2
      std::swap(B, C);
      std::swap(eta, zeta);
      // To make it faster, instead of "go to the point N1" we repeat N1 once
      // (which is equivalent - three swaps are sufficient to reorder ABC).
      if (A - B > eps || (A - B >= -eps && std::abs(xi) > std::abs(eta) + eps)) { // N1
        std::swap(A, B);
        std::swap(xi, eta);
      }
    }
    // N3
    // xi * eta * zeta > 0 <=> positive count is 1 or 3 and no zeros
    int pos_count = (xi > eps) + (eta > eps) + (zeta > eps);
    int nonneg_count = (xi >= -eps) + (eta >= -eps) + (zeta >= -eps);
    double sgn = (pos_count == nonneg_count && pos_count % 2 == 1) ? 1 : -1;
    xi = std::copysign(xi, sgn);
    eta = std::copysign(eta, sgn);
    zeta = std::copysign(zeta, sgn);
  }

  // Algorithm B from Gruber (1973).
  // Returns true if no change was needed.
  bool buerger_step() {
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
      return true;
    }
    return false;
  }

  // Returns number of iterations.
  int buerger_reduce() {
    int n = 0;
    double prev_sum = -1;
    int stall_count = 0;
    for (;;) {
      normalize();
      // In rare cases numerical errors push the algorithm into infinite loop,
      // as described in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
      // Ad-hoc solution: stop if a+b+c is stalled for 5 iterations.
      if (++n > 8) {  // don't waste time during the first few iterations
        double sum = std::sqrt(A) + std::sqrt(B) + std::sqrt(C);
        if (std::abs(sum - prev_sum) < sum * 1e-6) {
          if (++stall_count == 5)
            break;
        } else {
          stall_count = 0;
        }
        prev_sum = sum;
      }
      if (buerger_step())
        break;
    }
    return n;
  }

  // To be called after normalize() or is_normalized().
  // Returns true if it already was Niggli cell.
  // Algorithm from Krivy & Gruber, Acta Cryst. (1976) A32, 297.
  bool niggli_step(double epsilon=1e-9) {
    if (std::abs(xi) > B + epsilon ||  // step 5. from Krivy-Gruber (1976)
        (xi >= B - epsilon && 2 * eta < zeta - epsilon) ||
        (xi <= -(B - epsilon) && zeta < -epsilon)) {
      double sign_xi = xi >= 0 ? 1 : -1;
      C += B - xi * sign_xi;
      eta -= zeta * sign_xi;
      xi -= 2 * B * sign_xi;
    } else if (std::abs(eta) > A + epsilon ||  // step 6.
               (eta >= A - epsilon && 2 * xi < zeta - epsilon) ||
               (eta <= -(A - epsilon) && zeta < -epsilon)) {
      double sign_eta = eta >= 0 ? 1 : -1;
      C += A - eta * sign_eta;
      xi -= zeta * sign_eta;
      eta -= 2 * A * sign_eta;
    } else if (std::abs(zeta) > A + epsilon ||  // step 7.
               (zeta >= A - epsilon && 2 * xi < eta - epsilon) ||
               (zeta <= -(A - epsilon) && eta < -epsilon)) {
      double sign_zeta = zeta >= 0 ? 1 : -1;
      B += A - zeta * sign_zeta;
      xi -= eta * sign_zeta;
      zeta -= 2 * A * sign_zeta;
    } else if (xi + eta + zeta + A + B < -epsilon || // step 8.
               (xi + eta + zeta + A + B <= epsilon && 2 * (A + eta) + zeta > epsilon)) {
      C += A + B + xi + eta + zeta;
      xi += 2 * B + zeta;
      eta += 2 * A + zeta;
    } else {
      return true;
    }
    return false;
  }

  // Returns number of iterations.
  int niggli_reduce(double epsilon=1e-9, int iteration_limit=100) {
    int n = 0;
    for (;;) {
      normalize(epsilon);
      if (++n == iteration_limit || niggli_step(epsilon))
        break;
    }
    return n;
  }

  bool is_niggli(double epsilon=1e-9) const {
    return is_normalized() && GruberVector(*this).niggli_step(epsilon);
  }
};

} // namespace gemmi
#endif
