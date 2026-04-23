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

/// @brief Fast approximate exp for float; relative error < 1e-5.
/// @note Input must be in [-88, 88].
/// @param x exponent value
/// @return approximate exp(x)
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

/// @brief Precalculated isotropic density as sum of N Gaussians in r².
/// Amplitude and exponent coefficients are stored for fast evaluation.
template<int N, typename Real>
struct ExpSum {
  Real a[N];  ///< Amplitude coefficients
  Real b[N];  ///< Exponent coefficients

  /// @brief Calculate density at squared distance.
  /// @param r2 squared distance
  /// @return density at r²
  Real calculate(Real r2) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i] * r2);
    return density;
  }

  /// @brief Calculate density and its derivative.
  /// @param r distance
  /// @return pair of (density, derivative)
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

/// @brief Float specialisation of ExpSum using unsafe_expapprox for speed.
template<int N>
struct ExpSum<N, float> {
  float a[N];  ///< Amplitude coefficients
  float b[N];  ///< Exponent coefficients
  /// @brief Calculate density at squared distance.
  /// @param r2 squared distance
  /// @return density at r²
  float calculate(float r2) const {
    float density = 0;
    float tmp[N];
    for (int i = 0; i < N; ++i)
      tmp[i] = std::max(b[i] * r2, -88.f);
    for (int i = 0; i < N; ++i)
      density += a[i] * unsafe_expapprox(tmp[i]);
    return density;
  }

  /// @brief Calculate density and its derivative.
  /// @param r distance
  /// @return pair of (density, derivative)
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

/// @brief Precalculated anisotropic density as sum of N Gaussians.
/// Amplitude and tensor coefficients are stored for fast evaluation.
template<int N, typename Real>
struct ExpAnisoSum {
  Real a[N];           ///< Amplitude coefficients
  SMat33<Real> b[N];   ///< Anisotropic exponent tensor coefficients

  /// @brief Calculate density at vector position.
  /// @param r position vector
  /// @return density at r
  Real calculate(const Vec3& r) const {
    Real density = 0;
    for (int i = 0; i < N; ++i)
      density += a[i] * std::exp(b[i].r_u_r(r));
    return density;
  }
};

/// @brief Float specialisation of ExpAnisoSum using unsafe_expapprox for speed.
template<int N>
struct ExpAnisoSum<N, float> {
  float a[N];           ///< Amplitude coefficients
  SMat33<float> b[N];   ///< Anisotropic exponent tensor coefficients

  /// @brief Calculate density at vector position.
  /// @param r_ position vector
  /// @return density at r
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

/// @brief Compute x^1.5.
/// @tparam Real floating-point type
/// @param x input value
/// @return x^1.5
template<typename Real>
Real pow15(Real x) { return x * std::sqrt(x); }

/// @brief Gaussian approximation coefficients for atomic scattering.
/// N Gaussians plus optional constant term for efficient computation.
template<int N, int WithC, typename Real>
struct GaussianCoef {
  using coef_type = Real;          ///< Type alias for coefficient values
  static const int ncoeffs = N;    ///< Number of Gaussian terms (constant term separate)
  std::array<Real, 2*N+WithC> coefs;  ///< Coefficients: a[0..N-1], b[0..N-1], c (if WithC)
  /// @brief Get amplitude coefficient for Gaussian i.
  Real a(int n) const { return coefs[n]; }
  /// @brief Get exponent coefficient for Gaussian i.
  Real b(int n) const { return coefs[N+n]; }
  /// @brief Get constant term (0 if WithC=0).
  Real c() const { return WithC ? coefs[2*N] : 0; }

  /// @brief Set all coefficients.
  void set_coefs(const std::array<Real, 2*N+WithC>& c) { coefs = c; }

  /// @brief Calculate structure factor at specified reciprocal resolution.
  /// @param stol2 (sinθ/λ)²
  /// @return structure factor
  Real calculate_sf(Real stol2) const {
    Real sf = c();
    for (int i = 0; i < N; ++i)
      sf += a(i) * std::exp(-b(i)*stol2);
    return sf;
  }

  /// @brief Calculate isotropic density at given distance and B-factor.
  /// @param r2 squared distance
  /// @param B isotropic B-factor (Å²)
  /// @return density
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

  /// @brief Precalculate coefficients for fast isotropic density evaluation.
  /// @param B isotropic B-factor (Å²)
  /// @param addend added to constant term if WithC=1 (e.g., dispersion f')
  /// @return ExpSum ready for fast density calculation
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

  /// @brief Calculate anisotropic density at given position and U-tensor.
  /// @param r position vector
  /// @param U anisotropic displacement tensor
  /// @return density
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

  /// @brief Precalculate coefficients for fast anisotropic density evaluation with B-tensor.
  /// @param B anisotropic B-factor matrix
  /// @param addend added to constant term if WithC=1
  /// @return ExpAnisoSum ready for fast density calculation
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

  /// @brief Precalculate coefficients for fast anisotropic density evaluation with U-tensor.
  /// @param U anisotropic displacement tensor
  /// @param addend added to constant term if WithC=1
  /// @return ExpAnisoSum ready for fast density calculation
  ExpAnisoSum<N+WithC,Real> precalculate_density_aniso_u(const SMat33<float>& U,
                                                         Real addend=0) const {
    constexpr Real UtoB = 8 * sq(pi());
    return precalculate_density_aniso_b(U.scaled(UtoB), addend);
  }
};

} // namespace gemmi
#endif
