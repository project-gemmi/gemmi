// Copyright 2023 Global Phasing Ltd.

#include <gemmi/xds_ascii.hpp>
#include <gemmi/math.hpp>

namespace gemmi {

/// Based on Phil Evans' expertise and the literature, see:
/// https://github.com/project-gemmi/gemmi/discussions/248
/// \par p is defined as in XDS (p=0.5 for unpolarized beam).
void XdsAscii::apply_polarization_correction(double p, Vec3 normal) {
  if (!has_cell_axes())
    fail("unknown unit cell axes");
  Mat33 UB = cell_axes.inverse();
  Vec3 rot_axis = get_rotation_axis();
  Vec3 s0_dir = get_s0_direction();
  normal = normal.normalized();
  // The polarization normal is expected to be approx. orthogonal to the beam.
  // dot() is the same as cos_angle() for normalized vectors.
  if (normal.dot(s0_dir) > std::cos(rad(5.0)))
    fail("polarization normal is far from orthogonal to the incident beam");
  // make normal exactly orthogonal to the beam
  normal = s0_dir.cross(normal).cross(s0_dir).normalized();
  // wavevector
  Vec3 s0 = s0_dir / wavelength;
  double s0_m2 = 1. / s0.length_sq();  // s0^-2

  for (Refl& refl : data) {
    double phi = rad(rot_angle(refl));
    Vec3 h(refl.hkl[0], refl.hkl[1], refl.hkl[2]);
    Vec3 r0 = UB.multiply(h);
    Vec3 r = rotate_about_axis(r0, rot_axis, phi);
    Vec3 s = s0 + r;
#if 0
    double two_theta = s0.angle(s);
    // 2d sin(theta) = lambda
    double bragg_angle = std::asin(wavelength / (2 * unit_cell.calculate_d(refl.hkl)));
    printf("(%d %d %d) two-theta %g %g\n",
           refl.hkl[0], refl.hkl[1], refl.hkl[2], deg(two_theta), deg(2 * bragg_angle));
#endif
    // we should have |s| == |s0|, but just in case calculate it separately
    double s_m2 = 1. / s.length_sq();
    // 1 + cos^2(2theta) = 2 * correction for unpolarized beam
    double t = 1 + sq(s.dot(s0)) * s_m2 * s0_m2;
    double polariz_factor = (1 - 2*p) * (1 - sq(normal.dot(s)) * s_m2) + p * t;
    // We assume that the XDS files has polarization correction applied,
    // but for non-polarized beam. So we multiply intensities by P0=t/2
    // and divide by a hopefully more accurate polarization factor.
    double mult = 0.5 * t / polariz_factor;
    refl.iobs *= mult;
    refl.sigma *= mult;
    refl.rlp *= mult;
  }
}

}  // namespace gemmi
