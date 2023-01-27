// Copyright 2023 Global Phasing Ltd.

#include <gemmi/xds_ascii.hpp>
#include <gemmi/math.hpp>

namespace gemmi {

Mat33 XdsAscii::get_orientation() const {
  Vec3 a = get_unit_cell_axis(0);
  Vec3 b = get_unit_cell_axis(1);
  Vec3 c = get_unit_cell_axis(2);
  if (a.length_sq() == 0 || b.length_sq() == 0 || c.length_sq() == 0)
    fail("unknown unit cell axes");
  Vec3 ar = b.cross(c).normalized();
  Vec3 br = c.cross(a);
  Vec3 cr = ar.cross(br).normalized();
  br = cr.cross(ar);
  return Mat33::from_columns(ar, br, cr);
}

// Based on Phil Evans' expertise and also on
// R. Kahn et al, J. Appl. Cryst. (1982) 15, 330  doi:10.1107/S0021889882012060
void XdsAscii::apply_polarization_correction(double fraction, Vec3 normal) {
  normal = normal.normalized();
  Vec3 s0 = get_s0();
  // The polarization normal is expected to be approx. orthogonal to the beam.
  // dot() is the same as cos_angle() for normalized vectors.
  if (normal.dot(s0) > std::cos(rad(5.0)))
    fail("polarization normal is far from orthogonal to the incident beam");
  // make normal exactly orthogonal to the beam
  normal = s0.cross(normal).cross(s0).normalized();

  Vec3 rot_axis = get_rotation_axis();

  (void) fraction;
  // TODO
  /*
  for (Refl& refl : data) {
    refl.hkl;
    Vec3 s1 = ?;
    refl.iobs;
  }
  */
}

}  // namespace gemmi
