//! @file
//! @brief Atomic form factors approximated by Gaussian sums.
//!
//! Calculation of atomic form factors approximated by a sum of Gaussians.
//! Tables with numerical coefficients are in it92.hpp and c4322.hpp.

// Copyright 2019 Global Phasing Ltd.

// Calculation of atomic form factors approximated by a sum of Gaussians.
// Tables with numerical coefficients are in it92.hpp and c4322.hpp.

#ifndef GEMMI_FORMFACT_HPP_
#define GEMMI_FORMFACT_HPP_

#include <cmath>     // for exp, sqrt
#include <cstdint>   // for int32_t
#include <cstring>   // for memcpy
#include <limits>    // for numeric_limits
#include <utility>   // for pair
#include "math.hpp"  // for pi()

namespace gemmi {

//! @brief Fast exponential approximation for float.
//! @param x Argument (must be in [-88, 88])
//! @return Approximation of exp(x) with relative error < 1e-5
//!
//! NOTE: the argument x must be between -88 and 88.
//! It is based on expapprox() from
//! https://github.com/jhjourdan/SIMD-math-prims/blob/master/simd_math_prims.h
//! Relative error is below 1e-5.
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

//! @brief Precalculated exponential sum for isotropic atomic density.
//! @tparam N Number of Gaussian terms
//! @tparam Real Floating-point type
//!
//! precalculated density of an isotropic atom
template<int N, typename Real>
struct ExpSum {
  Real a[N], b[N];  //!< Coefficients: density = sum(a[i] * exp(b[i] * r^2))

  //! @brief Calculate density at r^2.
  //! @param r2 Squared distance
  //! @return Electron density
  Real calculate(Real r2) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i] * r2);
    return density;
  }

  //! @brief Calculate density and derivative at distance r.
  //! @param r Distance (not squared)
  //! @return Pair of (density, derivative)
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

//! @brief Float specialization using fast exp approximation.
//! @tparam N Number of Gaussian terms
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

//! @brief Precalculated exponential sum for anisotropic atomic density.
//! @tparam N Number of Gaussian terms
//! @tparam Real Floating-point type
//!
//! precalculated density of an anisotropic atom
template<int N, typename Real>
struct ExpAnisoSum {
  Real a[N];  //!< Amplitude coefficients
  SMat33<Real> b[N];  //!< Symmetric 3x3 exponent matrices

  //! @brief Calculate density at position r.
  //! @param r Position vector
  //! @return Electron density
  Real calculate(const Vec3& r) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i].r_u_r(r));
    return density;
  }
};

//! @brief Float specialization using fast exp approximation.
//! @tparam N Number of Gaussian terms
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

//! @brief Calculate x^1.5 efficiently.
//! @tparam Real Floating-point type
//! @param x Value
//! @return x^1.5 = x * sqrt(x)
template<typename Real>
Real pow15(Real x) { return x * std::sqrt(x); }

//! @brief Gaussian coefficients for scattering factor and density calculations.
//! @tparam N Number of Gaussian terms
//! @tparam WithC Whether constant term c is included (0 or 1)
//! @tparam Real Floating-point type
//!
//! Gaussian coefficients with functions to calculate sf and density.
template<int N, int WithC, typename Real>
struct GaussianCoef {
  using coef_type = Real;  //!< Coefficient type
  static const int ncoeffs = N;  //!< Number of Gaussian terms
  std::array<Real, 2*N+WithC> coefs;  //!< Coefficient array [a1..aN, b1..bN, c]
  Real a(int n) const { return coefs[n]; }  //!< Get amplitude coefficient a_n
  Real b(int n) const { return coefs[N+n]; }  //!< Get exponent coefficient b_n
  Real c() const { return WithC ? coefs[2*N] : 0; }  //!< Get constant c (if WithC=1)

  //! @brief Set all coefficients.
  //! @param c Coefficient array
  void set_coefs(const std::array<Real, 2*N+WithC>& c) { coefs = c; }

  //! @brief Calculate structure factor.
  //! @param stol2 (sin(theta)/lambda)^2
  //! @return Structure factor f(s)
  //!
  //! argument: (sin(theta)/lambda)^2
  Real calculate_sf(Real stol2) const {
    Real sf = c();
    for (int i = 0; i < N; ++i)
      sf += a(i) * std::exp(-b(i)*stol2);
    return sf;
  }

  //! @brief Calculate electron density at distance r with isotropic B-factor.
  //! @param r2 Squared distance from atom center
  //! @param B Isotropic temperature factor
  //! @return Electron density at r
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

  //! @brief Precalculate exponential sum for fast isotropic density evaluation.
  //! @param B Isotropic temperature factor
  //! @param addend Additional constant to add to c (usually dispersion f')
  //! @return Precalculated ExpSum for efficient repeated evaluation
  //!
  //! note: addend is considered only if WithC (addend is usually dispersion f')
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

  //! @brief Calculate electron density at position r with anisotropic U tensor.
  //! @param r Position vector from atom center
  //! @param U Anisotropic displacement tensor
  //! @return Electron density at r
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

  //! @brief Precalculate exponential sum for fast anisotropic density (B tensor).
  //! @param B Anisotropic B-factor tensor (8π²U)
  //! @param addend Additional constant to add to c (usually dispersion f')
  //! @return Precalculated ExpAnisoSum for efficient repeated evaluation
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

  //! @brief Precalculate exponential sum for fast anisotropic density (U tensor).
  //! @param U Anisotropic displacement tensor
  //! @param addend Additional constant to add to c (usually dispersion f')
  //! @return Precalculated ExpAnisoSum for efficient repeated evaluation
  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_u(const SMat33<float>& U,
                                                         Real addend=0) const {
    constexpr Real UtoB = 8 * sq(pi());
    return precalculate_density_aniso_b(U.scaled(UtoB), addend);
  }
};

} // namespace gemmi
#endif
