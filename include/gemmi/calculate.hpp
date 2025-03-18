// Copyright 2017-2018 Global Phasing Ltd.
//
// Calculate various properties of the model.

#ifndef GEMMI_CALCULATE_HPP_
#define GEMMI_CALCULATE_HPP_

#include <array>
#include <algorithm>  // for std::min, std::minmax
#include "model.hpp"
#include "select.hpp"

namespace gemmi {

template<class T> bool has_hydrogen(const T& obj) {
  for (const auto& child : obj.children())
    if (has_hydrogen(child))
      return true;
  return false;
}
template<> inline bool has_hydrogen(const Atom& atom) {
  return atom.is_hydrogen();
}

template<class T> size_t count_atom_sites(const T& obj, const Selection* sel=nullptr) {
  size_t sum = 0;
  if (!sel || sel->matches(obj))
    for (const auto& child : obj.children())
      sum += count_atom_sites(child, sel);
  return sum;
}
template<> inline size_t count_atom_sites(const Atom& atom, const Selection* sel) {
  return (!sel || sel->matches(atom)) ? 1 : 0;
}

template<class T> double count_occupancies(const T& obj, const Selection* sel=nullptr) {
  double sum = 0;
  if (!sel || sel->matches(obj))
    for (const auto& child : obj.children())
        sum += count_occupancies(child, sel);
  return sum;
}
template<> inline double count_occupancies(const Atom& atom, const Selection* sel) {
  return (!sel || sel->matches(atom)) ? atom.occ : 0;
}


template<class T> double calculate_mass(const T& obj) {
  double sum = 0;
  for (const auto& child : obj.children())
    sum += calculate_mass(child);
  return sum;
}
template<> inline double calculate_mass(const Atom& atom) {
  return atom.occ * atom.element.weight();
}

struct CenterOfMass {
  Position weighted_sum;
  double mass;
  Position get() const { return Position(weighted_sum / mass); }
};

template<class T> CenterOfMass calculate_center_of_mass(const T& obj) {
  CenterOfMass total{{}, 0.};
  for (const auto& child : obj.children()) {
    CenterOfMass part = calculate_center_of_mass(child);
    total = {total.weighted_sum + part.weighted_sum, total.mass + part.mass};
  }
  return total;
}
template<> inline CenterOfMass calculate_center_of_mass(const Atom& atom) {
  double w_mass = atom.element.weight() * atom.occ;
  return CenterOfMass{Position(atom.pos * w_mass), w_mass};
}

template<class T> std::pair<float,float> calculate_b_iso_range(const T& obj) {
  std::pair<float, float> range{INFINITY, -INFINITY};
  for (const auto& child : obj.children()) {
    auto r = calculate_b_iso_range(child);
    range.first = std::min(range.first, r.first);
    range.second = std::max(range.second, r.second);
  }
  return range;
}
template<> inline std::pair<float,float> calculate_b_iso_range(const Atom& atom) {
  return {atom.b_iso, atom.b_iso};
}

/// uses min/max eigenvalues of Baniso, or Biso if B-factor is isotropic
inline std::pair<double, double> calculate_b_aniso_range(const Model& model) {
  std::pair<double, double> range{INFINITY, -INFINITY};
  for (const Chain& chain : model.chains)
    for (const Residue& residue : chain.residues)
      for (const Atom& atom : residue.atoms) {
        if (atom.occ == 0)
          continue;
        if (atom.aniso.nonzero()) {
          std::array<double,3> eig = atom.aniso.calculate_eigenvalues();
          auto u = std::minmax({eig[0], eig[1], eig[2]});
          range.first = std::min(range.first, u.first * u_to_b());
          range.second = std::max(range.second, u.second * u_to_b());
        } else {
          range.first = std::min(range.first, (double) atom.b_iso);
          range.second = std::max(range.second, (double) atom.b_iso);
        }
      }
  return range;
}


template<class T> void expand_box(const T& obj, Box<Position>& box) {
  for (const auto& child : obj.children())
    expand_box(child, box);
}
template<> inline void expand_box(const Atom& atom, Box<Position>& box) {
  box.extend(atom.pos);
}

// we don't take NCS into account here (cf. NeighborSearch::set_bounding_cell())
inline Box<Position> calculate_box(const Structure& st, double margin=0.) {
  Box<Position> box;
  expand_box(st, box);
  if (margin != 0.)
    box.add_margin(margin);
  return box;
}

inline Box<Fractional> calculate_fractional_box(const Structure& st, double margin=0.) {
  if (!st.cell.is_crystal())
    fail("calculate_fractional_box(): Structure has no unit cell for fractionalization");
  Box<Fractional> box;
  for (const Model& model : st.models)
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          box.extend(st.cell.fractionalize(atom.pos));
  if (margin != 0.)
    box.add_margins({margin * st.cell.ar, margin * st.cell.br, margin * st.cell.cr});
  return box;
}


// Calculate B_est from E. Merritt, Some B_eq are more equivalent than others,
// Acta Cryst. A67, 512 (2011)
// http://skuld.bmsc.washington.edu/parvati/ActaA_67_512.pdf
inline double calculate_b_est(const Atom& atom) {
  auto eig = atom.aniso.calculate_eigenvalues();
  return u_to_b() * std::sqrt((eig[0] + eig[1] + eig[2]) /
                              (1/eig[0] + 1/eig[1] + 1/eig[2]));
}

inline double calculate_angle(const Position& p0, const Position& p1,
                              const Position& p2) {
  return (p0 - p1).angle(p2 - p1);
}

// discussion: https://stackoverflow.com/questions/20305272/
inline double calculate_dihedral(const Position& p0, const Position& p1,
                                 const Position& p2, const Position& p3) {
  Vec3 b0 = p1 - p0;
  Vec3 b1 = p2 - p1;
  Vec3 b2 = p3 - p2;
  Vec3 u = b1.cross(b0);
  Vec3 w = b2.cross(b1);
  double y = u.cross(w).dot(b1);
  double x = u.dot(w) * b1.length();
  return std::atan2(y, x);
}

/// the return value is in the same range as that of atan2(), i.e. [-pi, pi]
inline double calculate_dihedral_from_atoms(const Atom* a, const Atom* b,
                                            const Atom* c, const Atom* d) {
  if (a && b && c && d)
    return calculate_dihedral(a->pos, b->pos, c->pos, d->pos);
  return NAN;
}

inline double calculate_omega(const Residue& res, const Residue& next) {
  return calculate_dihedral_from_atoms(res.get_ca(), res.get_c(),
                                       next.get_n(), next.get_ca());
}

inline bool is_peptide_bond_cis(const Atom* ca1, const Atom* c,
                                const Atom* n, const Atom* ca2,
                                double tolerance_deg) {
  double omega = calculate_dihedral_from_atoms(ca1, c, n, ca2);
  return std::fabs(omega) < rad(tolerance_deg);
}

inline double calculate_chiral_volume(const Position& actr, const Position& a1,
                                      const Position& a2, const Position& a3) {
  return (a1 - actr).dot((a2 - actr).cross(a3 - actr));
}

inline std::array<double, 2> calculate_phi_psi(const Residue* prev,
                                               const Residue& res,
                                               const Residue* next) {
  std::array<double, 2> phi_psi{{NAN, NAN}};
  if (prev || next) {
    const Atom* CA = res.get_ca();
    const Atom* C = res.get_c();
    const Atom* N = res.get_n();
    if (prev)
      phi_psi[0] = calculate_dihedral_from_atoms(prev->get_c(), N, CA, C);
    if (next)
      phi_psi[1] = calculate_dihedral_from_atoms(N, CA, C, next->get_n());
  }
  return phi_psi;
}

GEMMI_DLL std::array<double, 4> find_best_plane(const std::vector<Atom*>& atoms);

inline double get_distance_from_plane(const Position& pos,
                                      const std::array<double, 4>& coeff) {
  return coeff[0] * pos.x + coeff[1] * pos.y + coeff[2] * pos.z + coeff[3];
}

GEMMI_DLL FTransform parse_triplet_as_ftransform(const std::string& s);

GEMMI_DLL SMat33<double> calculate_u_from_tls(const TlsGroup& tls, const Position& pos);

} // namespace gemmi
#endif
