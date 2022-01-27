// Functions derived from modified Bessel functions I1(x) and I0(x).
//
// Crystallographic codes (including Refmac and cctbx) often use polynomial
// approximation of I0 and I1 from p. 378 of Abramowitz and Stegun.
// Gemmi uses approximation based on polynomial coefficients from bessel_i0
// and bessel_i1(float) from Boost.Math:
// https://www.boost.org/doc/libs/1_76_0/libs/math/doc/html/math_toolkit/bessel/mbessel.html
// This approximation was derived in 2017 by John Maddock,
// building on the work of Pavel Holoborodko:
// https://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision/
// The efficiency is similar to that of scitbx.math.bessel_i1_over_i0.
// Using std::cyl_bessel_if was not considered, because it requires C++17,

#ifndef GEMMI_BESSEL_HPP_
#define GEMMI_BESSEL_HPP_

#include <cmath>

namespace gemmi {

template<int N>
inline double evaluate_polynomial(const double(&poly)[N], double x) {
  static_assert(N > 1, "");
  double result = poly[N-1];
  for (int i = N-2; i >= 0; --i)
    result = result * x + poly[i];
  return result;
}

template<class Dummy>
struct BesselTables_
{
  static const double P1[8];
  static const double Q1[9];
  static const double P2[5];
  static const double Q2[5];
  static const double Q3[3];
};
template<class Dummy> const double BesselTables_<Dummy>::P1[8] = {
  8.333333221e-02,
  6.944453712e-03,
  3.472097211e-04,
  1.158047174e-05,
  2.739745142e-07,
  5.135884609e-09,
  5.262251502e-11,
  1.331933703e-12
};
template<class Dummy> const double BesselTables_<Dummy>::Q1[9] = {
  1.00000003928615375e+00,
  2.49999576572179639e-01,
  2.77785268558399407e-02,
  1.73560257755821695e-03,
  6.96166518788906424e-05,
  1.89645733877137904e-06,
  4.29455004657565361e-08,
  3.90565476357034480e-10,
  1.48095934745267240e-11
};
template<class Dummy> const double BesselTables_<Dummy>::P2[5] = {
  3.98942115977513013e-01,
  -1.49581264836620262e-01,
  -4.76475741878486795e-02,
  -2.65157315524784407e-02,
  -1.47148600683672014e-01
};
template<class Dummy> const double BesselTables_<Dummy>::Q2[5] = {
  3.98942651588301770e-01,
  4.98327234176892844e-02,
  2.91866904423115499e-02,
  1.35614940793742178e-02,
  1.31409251787866793e-01
};
template<class Dummy> const double BesselTables_<Dummy>::Q3[3] = {
  3.98942391532752700e-01,
  4.98455950638200020e-02,
  2.94835666900682535e-02
};


inline double bessel_i1_over_i0(double x) {
  using B = BesselTables_<void>;
  if (x < 0)
    return -bessel_i1_over_i0(-x);
  if (x < 7.75) {
     double a = x * x / 4;
     double bessel0 = a * evaluate_polynomial(B::Q1, a) + 1;
     double R[3] = { 1, 0.5f, evaluate_polynomial(B::P1, a) };
     double bessel1 = x * evaluate_polynomial(R, a) / 2;
     return bessel1 / bessel0;
  }
  double p = evaluate_polynomial(B::P2, 1 / x);
  double q = x < 50 ? evaluate_polynomial(B::Q2, 1 / x)
                    : evaluate_polynomial(B::Q3, 1 / x);
  return p / q;
}

// Simplified function from Boost.Math.
// Similar to std::cyl_bessel_i(0, x), but much faster, less exact and doesn't
// throw out_of_range on negative argument. Relative error < 5.02e-08.
inline double bessel_i0(double x) {
  using B = BesselTables_<void>;
  x = std::fabs(x);
  if (x < 7.75) {
    double a = x * x / 4;
    return a * evaluate_polynomial(B::Q1, a) + 1;
  }
  if (x < 50)
    return std::exp(x) * evaluate_polynomial(B::Q2, 1 / x) / std::sqrt(x);
  double ex = std::exp(x / 2);
  return ex * evaluate_polynomial(B::Q3, 1 / x) / std::sqrt(x) * ex;
}

// Relative error < 4e-08.
inline double log_bessel_i0(double x) {
  using B = BesselTables_<void>;
  x = std::fabs(x);
  if (x < 7.75) {
    double a = x * x / 4;
    return std::log1p(a * evaluate_polynomial(B::Q1, a));
  }
  double q = x < 50 ? evaluate_polynomial(B::Q2, 1 / x)
                    : evaluate_polynomial(B::Q3, 1 / x);
  return x + std::log(q / std::sqrt(x));
}

} // namespace gemmi
#endif
