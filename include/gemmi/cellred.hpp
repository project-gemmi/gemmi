// Copyright 2021 Global Phasing Ltd.
//
// Unit cell reductions: Buerger, Niggli, Selling-Delaunay.

#ifndef GEMMI_CELLRED_HPP_
#define GEMMI_CELLRED_HPP_

#include <cmath>
#include <array>
#include <memory>    // for unique_ptr
#include "math.hpp"  // for deg
#include "symmetry.hpp"  // for Op
#include "unitcell.hpp"  // for UnitCell

namespace gemmi {

struct SellingVector;

/// @brief G6 Gruber vector representing a lattice cell
/// @details Contains Buerger/Niggli reduction algorithms
/// Originally, in B. Gruber, Acta Cryst. A29, 433 (1973), the vector was called
/// "characteristic" of a lattice/cell.
/// Functions that take epsilon as a parameter use it for comparisons,
/// as proposed in Grosse-Kunstleve et al, Acta Cryst. (2004) A60, 1.
struct GruberVector {
  //    a.a  b.b c.c 2b.c 2a.c 2a.b
  /// @brief G6 vector elements (A, B, C, ξ, η, ζ) from Gruber 1973
  double A, B, C, xi, eta, zeta;
  /// @brief Change of basis transformation (Op::Rot only)
  std::unique_ptr<Op> change_of_basis;

  /// @brief Construct from orthogonalization matrix of primitive cell
  /// @param m orthogonalization matrix
  explicit GruberVector(const Mat33& m)
    : A(m.column_dot(0,0)),
      B(m.column_dot(1,1)),
      C(m.column_dot(2,2)),
      xi(2 * m.column_dot(1,2)),
      eta(2 * m.column_dot(0,2)),
      zeta(2 * m.column_dot(0,1)) {}

  /// @brief Construct from G6 vector array
  /// @param g6 array {A, B, C, ξ, η, ζ}
  explicit GruberVector(const std::array<double,6>& g6)
    : A(g6[0]), B(g6[1]), C(g6[2]), xi(g6[3]), eta(g6[4]), zeta(g6[5]) {}

  /// @brief Construct from UnitCell with centring
  /// @param u unit cell
  /// @param centring centring type character
  /// @param track_change_of_basis if true, track change of basis
  GruberVector(const UnitCell& u, char centring, bool track_change_of_basis=false)
    : GruberVector(u.primitive_orth_matrix(centring)) {
    if (track_change_of_basis)
      set_change_of_basis(Op{centred_to_primitive(centring), {0,0,0}, 'x'});
  }

  /// @brief Construct from UnitCell with SpaceGroup
  /// @param u unit cell
  /// @param sg space group (may be null for P)
  /// @param track_change_of_basis if true, track change of basis
  GruberVector(const UnitCell& u, const SpaceGroup* sg, bool track_change_of_basis=false)
  : GruberVector(u, sg ? sg->centring_type() : 'P', track_change_of_basis) {}

  /// @brief Set change of basis transformation
  /// @param op transformation operation
  void set_change_of_basis(const Op& op) { change_of_basis.reset(new Op(op)); }

  /// @brief Get G6 vector as array
  /// @return array {A, B, C, ξ, η, ζ}
  std::array<double,6> parameters() const { return {A, B, C, xi, eta, zeta}; }
  /// @brief Get unit cell parameters
  /// @details Convert from G6 vector to cell parameters
  /// @return array {a, b, c, α, β, γ} with lengths in Å and angles in degrees
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
  /// @brief Construct UnitCell from G6 vector
  /// @return unit cell
  UnitCell get_cell() const { return UnitCell(cell_parameters()); }

  /// @brief Convert to Selling-Delaunay vector
  /// @return Selling vector
  SellingVector selling() const;

  /// @brief Check if G6 satisfies Gruber 1973 normalization conditions
  /// @return true if normalized (eq(3) from Gruber 1973)
  bool is_normalized() const {
    // eq(3) from Gruber 1973
    return A <= B && B <= C &&
           (A != B || std::abs(xi) <= std::abs(eta)) &&
           (B != C || std::abs(eta) <= std::abs(zeta)) &&
           (xi > 0) == (eta > 0) && (xi > 0) == (zeta > 0);
  }

  /// @brief Check if this is a Buerger-reduced cell
  /// @param epsilon tolerance for comparisons
  /// @return true if Buerger-reduced
  bool is_buerger(double epsilon=1e-9) const {
    return is_normalized() &&
           // eq (4) from Gruber 1973
           std::abs(xi) <= B + epsilon &&
           std::abs(eta) <= A + epsilon &&
           std::abs(zeta) <= A + epsilon;
  }

  /// @brief Normalize using Gruber Algorithm N
  /// @details Apply Gruber normalization to the G6 vector (Algorithm N from Gruber 1973)
  /// @param eps tolerance for comparisons
  void normalize(double eps=1e-9) {
    auto step_N1 = [&]() {
      if (A - B > eps || (A - B >= -eps && std::abs(xi) > std::abs(eta) + eps)) { // N1
        std::swap(A, B);
        std::swap(xi, eta);
        if (change_of_basis)
          swap_columns_and_negate(0, 1);
      }
    };
    step_N1();
    if (B - C > eps || (B - C >= -eps && std::abs(eta) > std::abs(zeta) + eps)) { // N2
      std::swap(B, C);
      std::swap(eta, zeta);
      if (change_of_basis)
        swap_columns_and_negate(1, 2);
      // To make it faster, instead of "go to the point N1" we repeat N1 once
      // (which is equivalent - three swaps are sufficient to reorder ABC).
      step_N1();
    }
    // N3
    // xi * eta * zeta > 0 <=> positive count is 1 or 3 and no zeros
    int pos_count = (xi > eps) + (eta > eps) + (zeta > eps);
    int nonneg_count = (xi >= -eps) + (eta >= -eps) + (zeta >= -eps);
    double sgn = (pos_count == nonneg_count && pos_count % 2 == 1) ? 1 : -1;
    if (change_of_basis) {
      if (sgn * xi < -eps)   negate_column(0);
      if (sgn * eta < -eps)  negate_column(1);
      if (sgn * zeta < -eps) negate_column(2);
      if (pos_count != nonneg_count && pos_count % 2 == 1)
        negate_column(std::fabs(zeta) <= eps ? 2 :
                      std::fabs(eta) <= eps ? 1 : 0);
    }
    xi = std::copysign(xi, sgn);
    eta = std::copysign(eta, sgn);
    zeta = std::copysign(zeta, sgn);
  }

  /// @brief Perform one step of Gruber Algorithm B
  /// @details Execute one iteration of the Buerger reduction algorithm
  /// @return true if already Buerger-reduced (no change made)
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

  /// @brief Reduce to Buerger cell
  /// @details Apply normalize() and buerger_step() repeatedly until convergence
  /// @return iteration count
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

  /// @brief Perform one step of Krivy-Gruber Niggli reduction
  /// @details To be called after normalize() or is_normalized()
  /// Algorithm from Krivy & Gruber, Acta Cryst. (1976) A32, 297
  /// @param epsilon tolerance for comparisons
  /// @return true if already Niggli cell (no change made)
  bool niggli_step(double epsilon=1e-9) {
    if (std::abs(xi) > B + epsilon ||  // step 5. from Krivy-Gruber (1976)
        (xi >= B - epsilon && 2 * eta < zeta - epsilon) ||
        (xi <= -(B - epsilon) && zeta < -epsilon)) {
      double sign_xi = xi >= 0 ? 1 : -1;
      C += B - xi * sign_xi;
      eta -= zeta * sign_xi;
      xi -= 2 * B * sign_xi;
      if (change_of_basis)
        add_column(1, 2, -int(sign_xi));
    } else if (std::abs(eta) > A + epsilon ||  // step 6.
               (eta >= A - epsilon && 2 * xi < zeta - epsilon) ||
               (eta <= -(A - epsilon) && zeta < -epsilon)) {
      double sign_eta = eta >= 0 ? 1 : -1;
      C += A - eta * sign_eta;
      xi -= zeta * sign_eta;
      eta -= 2 * A * sign_eta;
      if (change_of_basis)
        add_column(0, 2, -int(sign_eta));
    } else if (std::abs(zeta) > A + epsilon ||  // step 7.
               (zeta >= A - epsilon && 2 * xi < eta - epsilon) ||
               (zeta <= -(A - epsilon) && eta < -epsilon)) {
      double sign_zeta = zeta >= 0 ? 1 : -1;
      B += A - zeta * sign_zeta;
      xi -= eta * sign_zeta;
      zeta -= 2 * A * sign_zeta;
      if (change_of_basis)
        add_column(0, 1, -int(sign_zeta));
    } else if (xi + eta + zeta + A + B < -epsilon || // step 8.
               (xi + eta + zeta + A + B <= epsilon && 2 * (A + eta) + zeta > epsilon)) {
      C += A + B + xi + eta + zeta;
      xi += 2 * B + zeta;
      eta += 2 * A + zeta;
      if (change_of_basis) {
        add_column(0, 2, 1);
        add_column(1, 2, 1);
      }
    } else {
      return true;
    }
    return false;
  }

  /// @brief Reduce to Niggli cell
  /// @details Apply normalize() and niggli_step() repeatedly until convergence
  /// @param epsilon tolerance for comparisons
  /// @param iteration_limit maximum iterations
  /// @return iteration count
  int niggli_reduce(double epsilon=1e-9, int iteration_limit=100) {
    int n = 0;
    for (;;) {
      normalize(epsilon);
      if (++n == iteration_limit || niggli_step(epsilon))
        break;
    }
    return n;
  }

  /// @brief Check if this is a Niggli-reduced cell
  /// @param epsilon tolerance for comparisons
  /// @return true if Niggli-reduced
  bool is_niggli(double epsilon=1e-9) const {
    return is_normalized() && GruberVector(parameters()).niggli_step(epsilon);
  }

private:
  void swap_columns_and_negate(int i, int j) {
    for (auto& r : change_of_basis->rot)
      std::swap(r[i], r[j]);
    for (auto& r : change_of_basis->rot)
      for (auto& v : r)
        v = -v;
  }
  void negate_column(int i) {
    for (auto& r : change_of_basis->rot)
      r[i] = -r[i];
  }
  void add_column(int pos, int dest, int sign) {
    for (auto& r : change_of_basis->rot)
      r[dest] += sign * r[pos];
  }
};


/// @brief Selling-Delaunay vector for lattice reduction
/// @details Represents a lattice in terms of 6 dot products among four basis vectors
/// Used for Selling reduction. Based on:
/// - Chapter "Delaunay reduction and standardization" in
///   International Tables for Crystallography vol. A (2016), sec. 3.1.2.3
/// - Patterson & Love (1957), Acta Cryst. 10, 111
/// - Andrews et al (2019), Acta Cryst. A75, 115
struct SellingVector {
  /// @brief Selling vector elements s (6 scalar products)
  /// Order: b.c, a.c, a.b, a.d, b.d, c.d
  std::array<double,6> s;

  /// @brief Construct from Selling vector array
  /// @param s_ array of 6 dot products
  explicit SellingVector(const std::array<double,6>& s_) : s(s_) {}

  /// @brief Construct from orthogonalization matrix
  /// @param orth orthogonalization matrix
  explicit SellingVector(const Mat33& orth) {
    Vec3 b[4];
    for (int i = 0; i < 3; ++i)
      b[i] = orth.column_copy(i);
    b[3]= -b[0] - b[1] - b[2];
    s[0] = b[1].dot(b[2]);
    s[1] = b[0].dot(b[2]);
    s[2] = b[0].dot(b[1]);
    s[3] = b[0].dot(b[3]);
    s[4] = b[1].dot(b[3]);
    s[5] = b[2].dot(b[3]);
  }

  /// @brief Construct from UnitCell with centring
  /// @param u unit cell
  /// @param centring centring type
  SellingVector(const UnitCell& u, char centring)
    : SellingVector(u.primitive_orth_matrix(centring)) {}
  /// @brief Construct from UnitCell with SpaceGroup
  /// @param u unit cell
  /// @param sg space group (may be null for P)
  SellingVector(const UnitCell& u, const SpaceGroup* sg)
    : SellingVector(u, sg ? sg->centring_type() : 'P') {}

  /// @brief Sum of squared basis vector lengths
  /// @details The reduction minimizes sum(b_i²) which equals -2·sum(s_i)
  /// @return sum of squared lengths
  double sum_b_squared() const {
    return -2 * (s[0] + s[1] + s[2] + s[3] + s[4] + s[5]);
  }

  /// @brief Check if reduced
  /// @details A Selling vector is reduced if all s[i] ≤ 0 within tolerance
  /// @param eps tolerance
  /// @return true if reduced
  bool is_reduced(double eps=1e-9) const {
    return std::all_of(s.begin(), s.end(), [eps](double x) { return x <= eps; });
  }

  /// @brief Perform one reduction step
  /// @details Apply one iteration of Selling reduction
  /// @param eps tolerance
  /// @return true if a step was applied
  bool reduce_step(double eps=1e-9) {
    //printf(" s = %g %g %g %g %g %g  sum=%g\n",
    //       s[0], s[1], s[2], s[3], s[4], s[5], sum_b_squared());
    const int table[6][5] = {
      // When negating s[n] we need to apply operations from table[n]:
      // 2 x add, subtract, 2 x swap&add
      {2, 4, 3, 1, 5},  // 0
      {2, 3, 4, 0, 5},  // 1
      {1, 3, 5, 0, 4},  // 2
      {1, 2, 0, 4, 5},  // 3
      {0, 2, 1, 3, 5},  // 4
      {0, 1, 2, 3, 4},  // 5
    };

    double max_s = eps;
    int max_s_pos = -1;
    for (int i = 0; i < 6; ++i)
      if (s[i] > max_s) {
        max_s = s[i];
        max_s_pos = i;
      }
    if (max_s_pos < 0)
      return false;
    const int (&indices)[5] = table[max_s_pos];
    s[max_s_pos] = -max_s;
    s[indices[0]] += max_s;
    s[indices[1]] += max_s;
    s[indices[2]] -= max_s;
    std::swap(s[indices[3]], s[indices[4]]);
    s[indices[3]] += max_s;
    s[indices[4]] += max_s;
    //printf("  s[%d]=%g  sum: %g\n", max_s_pos, max_s, sum_b_squared());
    return true;
  }

  /// @brief Reduce to Selling form
  /// @details Apply reduce_step() repeatedly until convergence
  /// @param eps tolerance
  /// @param iteration_limit maximum iterations
  /// @return iteration count
  int reduce(double eps=1e-9, int iteration_limit=100) {
    int n = 0;
    while (++n != iteration_limit)
      if (!reduce_step(eps))
        break;
    return n;
  }

  /// @brief Convert to G6 parameters
  /// @return GruberVector parameters
  std::array<double,6> g6_parameters() const {
    return {-s[1]-s[2]-s[3], -s[0]-s[2]-s[4], -s[0]-s[1]-s[5], 2*s[0], 2*s[1], 2*s[2]};
  }

  /// @brief Convert to GruberVector
  /// @return Gruber vector equivalent
  GruberVector gruber() const { return GruberVector(g6_parameters()); }

  /// @brief Sort basis vectors by squared length
  /// @details Swap values to make a² ≤ b² ≤ c² ≤ d²
  /// @param eps tolerance for comparisons
  void sort(double eps=1e-9) {
    double abcd_sq_neg[4] = {
      // -a^2, -b^2, -c^2, -d^2 (negated - to be sorted in descending order)
      s[1]+s[2]+s[3], s[0]+s[2]+s[4], s[0]+s[1]+s[5], s[3]+s[4]+s[5]
    };
    // First, make sure that d >= a,b,c (therefore -d^2 <= -a^2,...).
    int min_idx = 3;
    for (int i = 0; i < 3; ++i)
      if (abcd_sq_neg[i] < abcd_sq_neg[min_idx] - eps)
        min_idx = i;
    switch (min_idx) {
      case 0:  // a <-> d
        std::swap(s[1], s[5]);
        std::swap(s[2], s[4]);
        break;
      case 1:  // b <-> d
        std::swap(s[0], s[5]);
        std::swap(s[2], s[3]);
        break;
      case 2:  // c <-> d
        std::swap(s[0], s[4]);
        std::swap(s[1], s[3]);
        break;
    }
    // we could stop here and not care about the order of a,b,c.
    std::swap(abcd_sq_neg[min_idx], abcd_sq_neg[3]);
    if (abcd_sq_neg[0] < abcd_sq_neg[1] - eps) {  // a <-> b
      std::swap(s[0], s[1]);
      std::swap(s[3], s[4]);
      std::swap(abcd_sq_neg[0], abcd_sq_neg[1]);
    }
    if (abcd_sq_neg[1] < abcd_sq_neg[2] - eps) {  // b <-> c
      std::swap(s[1], s[2]);
      std::swap(s[4], s[5]);
      std::swap(abcd_sq_neg[1], abcd_sq_neg[2]);
    }
    if (abcd_sq_neg[0] < abcd_sq_neg[1] - eps) {  // a <-> b
      std::swap(s[0], s[1]);
      std::swap(s[3], s[4]);
      //std::swap(abcd_sq_neg[0], abcd_sq_neg[1]);
    }
  }

  /// @brief Get unit cell parameters
  /// @return array {a, b, c, α, β, γ} via Gruber conversion
  std::array<double,6> cell_parameters() const {
    return gruber().cell_parameters();
  }
  /// @brief Construct UnitCell from Selling vector
  /// @return unit cell
  UnitCell get_cell() const { return UnitCell(cell_parameters()); }
};

/// @brief Convert GruberVector to Selling-Delaunay vector
/// @details Compute Selling vector from G6 Gruber parameters
/// @return Selling vector
inline SellingVector GruberVector::selling() const {
  double s0 = 0.5 * xi;
  double s1 = 0.5 * eta;
  double s2 = 0.5 * zeta;
  return SellingVector({s0, s1, s2, -A - s1 - s2, -B - s0 - s2, -C - s0 - s1});
}

} // namespace gemmi
#endif
