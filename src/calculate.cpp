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

} // namespace gemmi
