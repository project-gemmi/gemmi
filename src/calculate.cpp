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

} // namespace gemmi
