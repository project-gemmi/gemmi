// Copyright 2018 Global Phasing Ltd.

#include "gemmi/calculate.hpp"
#include "gemmi/eig3.hpp"

namespace gemmi {

std::array<double, 4> find_best_plane(const std::vector<Atom*>& atoms) {
  Vec3 mean;
  for (const Atom* atom : atoms)
    mean += atom->pos;
  mean /= (double) atoms.size();
  SMat33<double> m{0, 0, 0, 0, 0, 0};
  for (const Atom* atom : atoms) {
    Vec3 p = Vec3(atom->pos) - mean;
    m.u11 += p.x * p.x;
    m.u22 += p.y * p.y;
    m.u33 += p.z * p.z;
    m.u12 += p.x * p.y;
    m.u13 += p.x * p.z;
    m.u23 += p.y * p.z;
  }
  double eig[3] = {};
  Mat33 V = eigen_decomposition(m, eig);
  int smallest_idx = std::fabs(eig[0]) < std::fabs(eig[1]) ? 0 : 1;
  if (std::fabs(eig[2]) < std::fabs(eig[smallest_idx]))
    smallest_idx = 2;
  Vec3 eigvec = V.column_copy(smallest_idx);
  if (eigvec.x < 0)
    eigvec *= -1;
  return {{eigvec.x, eigvec.y, eigvec.z, -eigvec.dot(mean)}};
}

FTransform parse_triplet_as_ftransform(const std::string& s) {
  // c.f. parse_triplet()
  if (std::count(s.begin(), s.end(), ',') != 2)
    fail("expected exactly two commas in triplet");
  size_t comma1 = s.find(',');
  size_t comma2 = s.find(',', comma1 + 1);
  FTransform frac_tr;
  auto set_ftransform_row = [](FTransform& tr, int i, const std::string& part) {
    const double mult = 1. / Op::DEN;
    double decfr[4] = {};
    const auto op_row = parse_triplet_part(part, decfr);
    for (int j = 0; j < 3; ++j)
      tr.mat[i][j] = decfr[j] == 0. ? mult * op_row[j] : decfr[j];
    tr.vec.at(i) = decfr[3] == 0. ? mult * op_row[3] : decfr[3];
  };
  set_ftransform_row(frac_tr, 0, s.substr(0, comma1));
  set_ftransform_row(frac_tr, 1, s.substr(comma1 + 1, comma2 - (comma1 + 1)));
  set_ftransform_row(frac_tr, 2, s.substr(comma2 + 1));
  return frac_tr;
}

SMat33<double> calculate_u_from_tls(const TlsGroup& tls, const Position& pos) {
  Position r = (pos - tls.origin) * rad(1);
  SMat33<double> l_contrib = {
    r.y * r.y * tls.L.u33 + r.z * r.z * tls.L.u22 - 2 * r.z * r.y * tls.L.u23,
    r.x * r.x * tls.L.u33 + r.z * r.z * tls.L.u11 - 2 * r.z * r.x * tls.L.u13,
    r.x * r.x * tls.L.u22 + r.y * r.y * tls.L.u11 - 2 * r.y * r.x * tls.L.u12,
   -r.x * r.y * tls.L.u33 + r.z * ( r.x * tls.L.u23 + r.y * tls.L.u13 - r.z * tls.L.u12),
   -r.x * r.z * tls.L.u22 + r.y * ( r.x * tls.L.u23 - r.y * tls.L.u13 + r.z * tls.L.u12),
   -r.y * r.z * tls.L.u11 + r.x * (-r.x * tls.L.u23 + r.y * tls.L.u13 + r.z * tls.L.u12)
  };
  SMat33<double> s_contrib = {
    2 * (tls.S[1][0] * r.z - tls.S[2][0] * r.y),
    2 * (tls.S[2][1] * r.x - tls.S[0][1] * r.z),
    2 * (tls.S[0][2] * r.y - tls.S[1][2] * r.x),
    tls.S[2][0] * r.x - tls.S[2][1] * r.y + (tls.S[1][1] - tls.S[0][0]) * r.z,
    tls.S[1][2] * r.z - tls.S[1][0] * r.x + (tls.S[0][0] - tls.S[2][2]) * r.y,
    tls.S[0][1] * r.y - tls.S[0][2] * r.z + (tls.S[2][2] - tls.S[1][1]) * r.x
  };
  return tls.T + l_contrib + s_contrib;
}

} // namespace gemmi
