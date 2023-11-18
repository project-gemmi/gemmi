// Copyright 2019 Global Phasing Ltd.

// Calculation of atomic form factors approximated by a sum of Gaussians.
// Tables with numeric coefficient are in it92.hpp and c4322.hpp.

#ifndef GEMMI_FORMFACT_HPP_
#define GEMMI_FORMFACT_HPP_

#include <cmath>     // for exp, sqrt
#include <cstring>   // for memcpy
#include <limits>    // for numeric_limits
#include <utility>   // for pair
#include "math.hpp"  // for pi()

namespace gemmi {

// NOTE: the argument x must be between -88 and 88.
// It is based on expapprox() from
// https://github.com/jhjourdan/SIMD-math-prims/blob/master/simd_math_prims.h
// Relative error is below 1e-5.
inline float unsafe_expapprox(float x) {
  static_assert(std::numeric_limits<float>::is_iec559, "float is not IEEE 754?");
  //static float zero = 0.f;  // non-const to disable optimization
  float val = 12102203.1615614f * x + 1065353216.f;
  //val = std::max(val, zero);  // check if x < -88.02969
  std::int32_t vali = static_cast<std::int32_t>(val);
  std::int32_t xu1 = vali & 0x7F800000;
  std::int32_t xu2 = (vali & 0x7FFFFF) | 0x3F800000;
  float a, b;
  std::memcpy(&a, &xu1, 4);
  std::memcpy(&b, &xu2, 4);
  return a * (0.509871020f + b * (0.312146713f + b * (0.166617139f + b *
          (-2.190619930e-3f + b * 1.3555747234e-2f))));
}

// precalculated density of an isotropic atom
template<int N, typename Real>
struct ExpSum {
  Real a[N], b[N];

  Real calculate(Real r2) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i] * r2);
    return density;
  }

  std::pair<Real,Real> calculate_with_derivative(Real r) const {
    Real density = 0;
    Real derivative = 0;
    for (int i = 0; i < N; ++i) {
      Real y = a[i] * std::exp(b[i] * (r * r));
      density += y;
      derivative += 2 * b[i] * r * y;
    }
    return std::make_pair(density, derivative);
  }
};

template<int N>
struct ExpSum<N, float> {
  float a[N], b[N];
  float calculate(float r2) const {
    float density = 0;
    float tmp[N];
    for (int i = 0; i < N; ++i)
      tmp[i] = std::max(b[i] * r2, -88.f);
    for (int i = 0; i < N; ++i)
      density += a[i] * unsafe_expapprox(tmp[i]);
    return density;
  }

  std::pair<float,float> calculate_with_derivative(float r) const {
    float density = 0;
    float derivative = 0;
    float tmp[N];
    for (int i = 0; i < N; ++i)
      tmp[i] = std::max(b[i] * (r * r), -88.f);
    for (int i = 0; i < N; ++i) {
      float y = a[i] * unsafe_expapprox(tmp[i]);
      density += y;
      derivative += y * b[i] * (2 * r);
    }
    return std::make_pair(density, derivative);
  }
};

// precalculated density of an anisotropic atom
template<int N, typename Real>
struct ExpAnisoSum {
  Real a[N];
  SMat33<Real> b[N];

  Real calculate(const Vec3& r) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i].r_u_r(r));
    return density;
  }
};

template<int N>
struct ExpAnisoSum<N, float> {
  float a[N];
  SMat33<float> b[N];

  float calculate(const Vec3& r_) const {
    Vec3f r((float)r_.x, (float)r_.y, (float)r_.z);
    float density = 0;
    float tmp[N];
    for (int i = 0; i < N; ++i)
      tmp[i] = std::max(b[i].r_u_r(r), -88.f);
    for (int i = 0; i < N; ++i)
      density += a[i] * unsafe_expapprox(tmp[i]);
    return density;
  }
};

template<typename Real>
Real pow15(Real x) { return x * std::sqrt(x); }

// Gaussian coefficients with functions to calculate sf and density.
template<int N, int WithC, typename Real>
struct GaussianCoef {
  using coef_type = Real;
  static const int ncoeffs = N;
  std::array<Real, 2*N+WithC> coefs;
  Real a(int n) const { return coefs[n]; }
  Real b(int n) const { return coefs[N+n]; }
  Real c() const { return WithC ? coefs[2*N] : 0; }

  void set_coefs(const std::array<Real, 2*N+WithC>& c) { coefs = c; }

  // argument: (sin(theta)/lambda)^2
  Real calculate_sf(Real stol2) const {
    Real sf = c();
    for (int i = 0; i < N; ++i)
      sf += a(i) * std::exp(-b(i)*stol2);
    return sf;
  }

  Real calculate_density_iso(Real r2, Real B) const {
    constexpr Real _4pi = Real(4 * pi());
    Real r2pi = Real(r2 * pi());
    Real density = c() * pow15(_4pi / B) * std::exp(-(_4pi / B) * r2pi);
    for (int i = 0; i < N; ++i) {
      Real t = _4pi / (b(i)+B);
      density += a(i) * pow15(t) * std::exp(-t*r2pi);
    }
    return density;
  }

  // note: addend is considered only if WithC (addend is usually dispersion f')
  ExpSum<N+WithC,Real> precalculate_density_iso(Real B, Real addend=0) const {
    ExpSum<N+WithC,Real> prec;
    constexpr Real _4pi = Real(4 * pi());
    for (int i = 0; i < N; ++i) {
      Real t = _4pi / (b(i)+B);
      prec.a[i] = a(i) * pow15(t);
      prec.b[i] = -t * Real(pi());
    }
    if (WithC) {
      Real t = _4pi / B;
      prec.a[N] = (c() + addend) * pow15(t);
      prec.b[N] = -t * Real(pi());
    }
    return prec;
  }

  Real calculate_density_aniso(const Vec3& r, const SMat33<float>& U) const {
    constexpr Real pi2 = sq(pi());
    const SMat33<Real> B = U.scaled(8 * pi2);
    Real density = c() * pow15(4 * pi()) / std::sqrt(B.determinant()) *
                   std::exp(-4 * pi2 * B.inverse().r_u_r(r));
    for (int i = 0; i < N; ++i) {
      SMat33<Real> Bb = B.added_kI(b(i));
      density += a(i) * pow15(4 * pi()) / std::sqrt(Bb.determinant()) *
                 std::exp(-4 * pi2 * Bb.inverse().r_u_r(r));
    }
    return density;
  }

  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_b(const SMat33<Real>& B,
                                                         Real addend=0) const {
    constexpr Real m4pi2 = Real(-4 * sq(pi()));
    constexpr Real pow_4pi_15 = (Real) 44.546623974653663; // pow15(4 * pi())
    ExpAnisoSum<N+WithC,Real> prec;
    for (int i = 0; i < N; ++i) {
      SMat33<Real> Bb = B.added_kI(b(i));
      prec.a[i] = a(i) * pow_4pi_15 / std::sqrt(Bb.determinant());
      prec.b[i] = Bb.inverse().scaled(m4pi2);
    }
    if (WithC) {
      prec.a[N] = (c() + addend) * pow_4pi_15 / std::sqrt(B.determinant());
      prec.b[N] = B.inverse().scaled(m4pi2);
    }
    return prec;
  }

  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_u(const SMat33<float>& U,
                                                         Real addend=0) const {
    constexpr Real UtoB = 8 * sq(pi());
    return precalculate_density_aniso_b(U.scaled(UtoB), addend);
  }
};

} // namespace gemmi
#endif
