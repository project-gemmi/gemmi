// Copyright 2018 Global Phasing Ltd.

#include "gemmi/calculate.hpp"

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
  auto eig = m.calculate_eigenvalues();
  double smallest = eig[0];
  for (double d : {eig[1], eig[2]})
    if (std::fabs(d) < std::fabs(smallest))
      smallest = d;
  Vec3 eigvec = m.calculate_eigenvector(smallest);
  if (eigvec.x < 0)
    eigvec *= -1;
  return {{eigvec.x, eigvec.y, eigvec.z, -eigvec.dot(mean)}};
}

} // namespace gemmi
