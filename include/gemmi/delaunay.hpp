// Copyright 2021 Global Phasing Ltd.
//
// Selling-Delaunay reduction of unit cell / lattice basis.

#ifndef GEMMI_DELAUNAY_HPP_
#define GEMMI_DELAUNAY_HPP_

#include <cmath>
#include <array>
#include <algorithm>  // for all_of
#include "math.hpp"  // for deg

namespace gemmi {

void assert_almost_eq(const char* s, double a, double b) {
  if (std::fabs(a-b) > 1e-9) {
    printf("FAIL %s: %g %g  diff %g\n", s, a, b, std::fabs(a-b));
    exit(1);
  }
}

// See:
// - chapter "Delaunay reduction and standardization" in
//   International Tables for Crystallography vol. A (2016), sec. 3.1.2.3.
//   https://onlinelibrary.wiley.com/iucr/itc/Ac/ch3o1v0001/
// - Patterson & Love (1957), Acta Cryst. 10, 111,
//   "Remarks on the Delaunay reduction", doi:10.1107/s0365110x57000328
// - Andrews et al (2019), Acta Cryst. A75, 115,
//   "Selling reduction versus Niggli reduction for crystallographic lattices".
struct Delaunay {
  // b.c a.c a.b a.d b.d c.d
  double s[6];
  // for now we also store:
  Vec3 b[4];

  Delaunay() = default;

  Delaunay(const Mat33& orth) {
    for (int i = 0; i < 3; ++i)
      b[i] = orth.column_copy(i);
    b[3]= -b[0] - b[1] - b[2];
    s[0] = dot(1,2);
    s[1] = dot(0,2);
    s[2] = dot(0,1);
    s[3] = dot(0,3);
    s[4] = dot(1,3);
    s[5] = dot(2,3);
  }

  double dot(int i, int j) const { return b[i].dot(b[j]); }

  // The reduction minimizes the sum b_i^2.
  double sum_b_squared() const {
    double sum = -(s[0] + s[1] + s[2] + s[3] + s[4] + s[5]);
    assert_almost_eq("b_squared", 2 * sum,
        b[0].length_sq() + b[1].length_sq() + b[2].length_sq() + b[3].length_sq());
    return sum;
  }

  bool is_reduced(double eps=1e-9) const {
    return std::all_of(s, s+6, [eps](double x) { return x <= eps; });
  }

  bool reduce_step(double eps=1e-9) {
    //printf("  a, b, c, d = %g %g %g %g   sum: %g\n",
    //       b[0].length(), b[1].length(), b[2].length(), b[3].length(),
    //       sum_b_squared());
    //printf(" s = %g %g %g %g %g %g\n", s[0], s[1], s[2], s[3], s[4], s[5]);
    double max_prod = eps;
    int sel_i = 0, sel_j = 0;
    for (int i = 0; i < 3; ++i)
      for (int j = i+1; j < 4; ++j) {
        double prod = dot(i, j);
        if (prod > max_prod) {
          sel_i = i;
          sel_j = j;
          max_prod = prod;
        }
      }
    if (sel_j == 0)
      return false;
    for (int i = 0; i < 4; ++i)
      if (i != sel_i && i != sel_j)
        b[i] += b[sel_j];
    b[sel_j] = -b[sel_j];

    const int table[6][5] = {
      // 2 x add, subtract, 2 x swap&add
      {2, 4, 3, 1, 5},  // negate s[0]
      {2, 3, 4, 0, 5},  // negate s[1]
      {1, 3, 5, 0, 4},  // negate s[2]
      {1, 2, 0, 4, 5},  // negate s[3]
      {0, 2, 1, 3, 5},  // negate s[4]
      {0, 1, 2, 3, 4},  // negate s[5]
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
    double v = s[max_s_pos];
    const int (&indices)[5] = table[max_s_pos];
    s[max_s_pos] = -v;
    s[indices[0]] += v;
    s[indices[1]] += v;
    s[indices[2]] -= v;
    std::swap(s[indices[3]], s[indices[4]]);
    s[indices[3]] += v;
    s[indices[4]] += v;

    //printf("[debug]  (%d %d) product:%g\n", sel_i, sel_j, max_prod);
    assert_almost_eq("s[0]", s[0], dot(1,2));
    assert_almost_eq("s[1]", s[1], dot(0,2));
    assert_almost_eq("s[2]", s[2], dot(0,1));
    assert_almost_eq("s[3]", s[3], dot(0,3));
    assert_almost_eq("s[4]", s[4], dot(1,3));
    assert_almost_eq("s[5]", s[5], dot(2,3));
    //printf("[debug]  sum: %g\n", sum_b_squared());
    return true;
  }

  // Returns number of iterations.
  int reduce(double eps=1e-9, int iteration_limit=100) {
    int n = 0;
    while (++n != iteration_limit)
      if (!reduce_step(eps))
        break;
    return n;
  }

  std::array<double,6> g6_parameters() const {
    return {-s[1]-s[2]-s[3], -s[0]-s[2]-s[4], -s[0]-s[1]-s[5], 2*s[0], 2*s[1], 2*s[2]};
  }

  void sort_last() {
    double abcd_sq_neg[4] = {
      s[1]+s[2]+s[3], s[0]+s[2]+s[4], s[0]+s[1]+s[5], s[3]+s[4]+s[5]
    };
    int min_idx = 0;
    for (int i = 1; i < 4; ++i)
      if (abcd_sq_neg[i] < abcd_sq_neg[min_idx])
        min_idx = i;
    if (min_idx == 0) {         // a <-> d
      //printf("a <-> d\n");
      std::swap(s[1], s[5]);
      std::swap(s[2], s[4]);
    } else if (min_idx == 1) {  // b <-> d
      //printf("b <-> d\n");
      std::swap(s[0], s[5]);
      std::swap(s[2], s[3]);
    } else if (min_idx == 2) {  // c <-> d
      //printf("c <-> d\n");
      std::swap(s[0], s[4]);
      std::swap(s[1], s[3]);
    }
  }


  std::array<double,6> cell_parameters() const {
    double A = std::sqrt(-s[1] - s[2] - s[3]);
    double B = std::sqrt(-s[0] - s[2] - s[4]);
    double C = std::sqrt(-s[0] - s[1] - s[5]);
    //double D = std::sqrt(-s[3] - s[4] - s[5]);
    return {A, B, C,
            deg(std::acos(s[0]/(B*C))),
            deg(std::acos(s[1]/(A*C))),
            deg(std::acos(s[2]/(A*B)))};
  }
};

} // namespace gemmi
#endif
