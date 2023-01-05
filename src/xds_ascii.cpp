// Copyright 2023 Global Phasing Ltd.

#include <gemmi/xds_ascii.hpp>
#include <gemmi/math.hpp>

namespace gemmi {

// Based on Phil Evans' expertise and also on
// R. Kahn et al, J. Appl. Cryst. (1982) 15, 330  doi:10.1107/S0021889882012060
void XdsAscii::apply_polarization_correction(double fraction, Vec3 normal) {
  normal = normal.normalized();
  Vec3 beam_dir(incident_beam_dir[0], incident_beam_dir[1], incident_beam_dir[2]);
  double beam_dir_length = beam_dir.length();
  if (beam_dir_length == 0)
    fail("unknown incident beam direction");
  beam_dir /= beam_dir_length;
  // The polarization normal is expected to be orthogonal to the beam.
  // dot() is the same as cos_angle() for normalized vectors.
  if (normal.cos_angle(beam_dir) > std::cos(rad(1.0)))
    fail("polarization normal is not orthogonal to the incident beam");
  /*
  Vec3 rot_axis(rotation_axis[0], rotation_axis[1], rotation_axis[2]);
  double rot_axis_length = rot_axis.length();
  if (rot_axis_length == 0)
    fail("unknown rotation axis");
  rot_axis /= rot_axis_length;
  */

  (void) fraction;
  // TODO
}

}  // namespace gemmi
