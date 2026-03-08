// Copyright 2026 Global Phasing Ltd.

#include "gemmi/chemcomp_xyz.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/levmar.hpp"
#include "gemmi/math.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>

namespace gemmi {

namespace {

Position position_from_angle_and_torsion_local(const Position& x1,
                                               const Position& x2,
                                               const Position& x3,
                                               double dist,
                                               double theta,
                                               double tau) {
  Vec3 u = x2 - x1;
  Vec3 v = x3 - x2;
  Vec3 e1 = v.normalized();
  Vec3 e2 = -(u - u.dot(e1) * e1).normalized();
  Vec3 e3 = e1.cross(e2);
  Vec3 e23 = std::cos(tau) * e2 + std::sin(tau) * e3;
  return x3 + Position(dist * (-std::cos(theta) * e1 + std::sin(theta) * e23));
}

Position arbitrary_position_from_angle_local(const Position& x2,
                                             const Position& x3,
                                             double dist,
                                             double theta) {
  Vec3 u(1, 0, 0);
  Vec3 v = x3 - x2;
  Vec3 e1 = v.normalized();
  Vec3 e2_ = -(u - u.dot(e1) * e1);
  if (e2_.length_sq() < 1e-6) {
    u = Vec3(0, 1, 0);
    e2_ = -(u - u.dot(e1) * e1);
  }
  Vec3 e2 = e2_.normalized();
  return x3 + Position(dist * (-std::cos(theta) * e1 + std::sin(theta) * e2));
}

Vec3 get_vector_to_line(const Position& point,
                        const Position& point_on_the_line,
                        const Vec3& unit_vector) {
  Vec3 ap = point_on_the_line - point;
  return ap - ap.dot(unit_vector) * unit_vector;
}

std::pair<Position, Position> trilaterate_local(const Position& p1, double r1sq,
                                                const Position& p2, double r2sq,
                                                const Position& p3, double r3sq) {
  Vec3 ex = (p2 - p1).normalized();
  double i = ex.dot(p3 - p1);
  Vec3 ey = (Vec3(p3) - p1 - i * ex).normalized();
  Vec3 ez = ex.cross(ey);
  double d = (p2 - p1).length();
  double j = ey.dot(p3 - p1);
  double x = (r1sq - r2sq + d * d) / (2 * d);
  double y = (r1sq - r3sq + i * i + j * j) / (2 * j) - x * i / j;
  double z2 = r1sq - x * x - y * y;
  double z = std::sqrt(z2);
  return std::make_pair(p1 + Position(x * ex + y * ey + z * ez),
                        p1 + Position(x * ex + y * ey - z * ez));
}

Position trilaterate_in_plane(const Position& p1, double r1sq,
                              const Position& p2, double r2sq,
                              const Position& p3, double r3sq) {
  Vec3 ex = (p2 - p1).normalized();
  double i = ex.dot(p3 - p1);
  Vec3 eyv = Vec3(p3) - p1 - i * ex;
  if (eyv.length_sq() < 1e-8)
    return Position(NAN, NAN, NAN);
  Vec3 ey = eyv.normalized();
  double d = (p2 - p1).length();
  double j = ey.dot(p3 - p1);
  if (std::fabs(d) < 1e-8 || std::fabs(j) < 1e-8)
    return Position(NAN, NAN, NAN);
  double x = (r1sq - r2sq + d * d) / (2 * d);
  double y = (r1sq - r3sq + i * i + j * j) / (2 * j) - x * i / j;
  return p1 + Position(x * ex + y * ey);
}

std::pair<Position, Position> position_from_two_angles_local(const Position& p1,
                                                             const Position& p2,
                                                             const Position& p3,
                                                             double dist14,
                                                             double theta214,
                                                             double theta314) {
  double d12sq = p1.dist_sq(p2);
  double d13sq = p1.dist_sq(p3);
  double d14sq = dist14 * dist14;
  double d24sq = d14sq + d12sq - 2 * std::sqrt(d14sq * d12sq) * std::cos(theta214);
  double d34sq = d14sq + d13sq - 2 * std::sqrt(d14sq * d13sq) * std::cos(theta314);
  return trilaterate_local(p1, d14sq, p2, d24sq, p3, d34sq);
}

double calculate_tetrahedral_delta(double theta0, double theta1, double theta2) {
  double x = std::cos(theta1);
  double y = (std::cos(theta2) - x * std::cos(theta0)) / std::sin(theta0);
  double z2 = 1 - x * x - y * y;
  double z = std::sqrt(std::max(0.0, z2));
  return std::asin(z);
}

double angle_in_triangle(double a, double b, double c) {
  return std::acos((a * a + c * c - b * b) / (2 * a * c));
}

double cone_theta_from_hh_angle(double hh_angle) {
  double c = (std::cos(hh_angle) + 0.5) / 1.5;
  c = std::max(0.0, std::min(1.0, c));
  return std::acos(-std::sqrt(c));
}

bool atom_has_finite_xyz(const ChemComp::Atom& atom) {
  return std::isfinite(atom.xyz.x) && std::isfinite(atom.xyz.y) && std::isfinite(atom.xyz.z);
}

double get_bond_dist_or_default(const ChemComp& cc,
                                const std::vector<std::vector<size_t>>& adjacency,
                                const std::vector<std::vector<double>>& bond_dist,
                                size_t a, size_t b);

double angle_rad_or_default(const ChemComp& cc,
                            const std::map<std::string, size_t>& atom_index,
                            const std::vector<std::vector<size_t>>& adjacency,
                            size_t a, size_t center, size_t b);

bool is_heavy_atom(const ChemComp::Atom& atom) {
  return !atom.is_hydrogen();
}

double default_angle_for_center(const std::vector<std::vector<size_t>>& adjacency,
                                size_t center) {
  size_t degree = adjacency[center].size();
  if (degree == 2)
    return 120.0;
  return 109.5;
}

const Restraints::Angle* find_angle_restraint(
    const ChemComp& cc, const std::map<std::string, size_t>& atom_index,
    size_t a, size_t center, size_t b) {
  auto match = [&](const Restraints::Angle& angle, size_t ia, size_t ic, size_t ib) {
    auto ita = atom_index.find(angle.id1.atom);
    auto itc = atom_index.find(angle.id2.atom);
    auto itb = atom_index.find(angle.id3.atom);
    return ita != atom_index.end() && itc != atom_index.end() && itb != atom_index.end() &&
           ita->second == ia && itc->second == ic && itb->second == ib;
  };
  for (const auto& angle : cc.rt.angles)
    if (match(angle, a, center, b) || match(angle, b, center, a))
      return &angle;
  return nullptr;
}

const Restraints::Torsion* find_torsion_restraint(
    const ChemComp& cc, const std::map<std::string, size_t>& atom_index,
    size_t a, size_t b, size_t c, size_t d) {
  auto match = [&](const Restraints::Torsion& tor, size_t i1, size_t i2, size_t i3, size_t i4) {
    auto it1 = atom_index.find(tor.id1.atom);
    auto it2 = atom_index.find(tor.id2.atom);
    auto it3 = atom_index.find(tor.id3.atom);
    auto it4 = atom_index.find(tor.id4.atom);
    return it1 != atom_index.end() && it2 != atom_index.end() &&
           it3 != atom_index.end() && it4 != atom_index.end() &&
           it1->second == i1 && it2->second == i2 &&
           it3->second == i3 && it4->second == i4;
  };
  for (const auto& tor : cc.rt.torsions)
    if (match(tor, a, b, c, d) || match(tor, d, c, b, a))
      return &tor;
  return nullptr;
}

double placement_contact_score(const ChemComp& cc, const Position& pos, size_t skip1, size_t skip2) {
  double score = 0.0;
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    if (i == skip1 || i == skip2 || !atom_has_finite_xyz(cc.atoms[i]))
      continue;
    score += 1.0 / std::max(1e-6, pos.dist_sq(cc.atoms[i].xyz));
  }
  return score;
}

double distance_sq_from_angle(double d12, double d13, double angle123) {
  return d12 * d12 + d13 * d13 - 2.0 * d12 * d13 * std::cos(angle123);
}

bool atom_in_plane_restraint(const ChemComp& cc,
                             const std::map<std::string, size_t>& atom_index,
                             size_t idx) {
  for (const Restraints::Plane& plane : cc.rt.planes) {
    for (const Restraints::AtomId& atom_id : plane.ids) {
      std::map<std::string, size_t>::const_iterator it = atom_index.find(atom_id.atom);
      if (it != atom_index.end() && it->second == idx)
        return true;
    }
  }
  return false;
}

bool plane_from_positions(const std::vector<Position>& positions,
                          Position& centroid, Vec3& normal) {
  if (positions.size() < 3)
    return false;
  centroid = Position();
  for (const Position& pos : positions)
    centroid += pos;
  centroid /= static_cast<double>(positions.size());

  Vec3 accum(0, 0, 0);
  for (size_t i = 0; i != positions.size(); ++i) {
    const Position& a = positions[i];
    const Position& b = positions[(i + 1) % positions.size()];
    accum.x += (a.y - b.y) * (a.z + b.z);
    accum.y += (a.z - b.z) * (a.x + b.x);
    accum.z += (a.x - b.x) * (a.y + b.y);
  }
  if (accum.length_sq() < 1e-8)
    return false;
  normal = accum.normalized();
  return true;
}

void enforce_plane_restraints(ChemComp& cc,
                              const std::map<std::string, size_t>& atom_index) {
  for (const Restraints::Plane& plane : cc.rt.planes) {
    std::vector<size_t> indices;
    std::vector<Position> positions;
    indices.reserve(plane.ids.size());
    positions.reserve(plane.ids.size());
    for (const Restraints::AtomId& atom_id : plane.ids) {
      auto it = atom_index.find(atom_id.atom);
      if (it == atom_index.end())
        continue;
      size_t idx = it->second;
      if (!atom_has_finite_xyz(cc.atoms[idx]))
        continue;
      indices.push_back(idx);
      positions.push_back(cc.atoms[idx].xyz);
    }
    Position centroid;
    Vec3 normal;
    if (!plane_from_positions(positions, centroid, normal))
      continue;
    for (size_t idx : indices) {
      Vec3 delta = cc.atoms[idx].xyz - centroid;
      cc.atoms[idx].xyz -= Position(normal * delta.dot(normal));
    }
  }
}

bool reflect_across_plane(Position& pos, const Position& p1, const Position& p2,
                          const Position& p3) {
  Vec3 n = (p2 - p1).cross(p3 - p1);
  if (n.length_sq() < 1e-8)
    return false;
  n = n.normalized();
  double signed_dist = (pos - p1).dot(n);
  pos -= Position(2.0 * signed_dist * n);
  return true;
}

void enforce_chirality_restraints(ChemComp& cc,
                                  const std::map<std::string, size_t>& atom_index,
                                  const std::vector<std::vector<size_t>>& adjacency,
                                  const std::vector<bool>& ring_atom) {
  auto candidate_score = [&](size_t idx) {
    if (cc.atoms[idx].is_hydrogen())
      return 4;
    if (!ring_atom[idx] && adjacency[idx].size() <= 2)
      return 3;
    if (!ring_atom[idx])
      return 2;
    return 0;
  };
  for (const Restraints::Chirality& chir : cc.rt.chirs) {
    if (chir.sign == ChiralityType::Both)
      continue;
    auto ctr_it = atom_index.find(chir.id_ctr.atom);
    auto a1_it = atom_index.find(chir.id1.atom);
    auto a2_it = atom_index.find(chir.id2.atom);
    auto a3_it = atom_index.find(chir.id3.atom);
    if (ctr_it == atom_index.end() || a1_it == atom_index.end() ||
        a2_it == atom_index.end() || a3_it == atom_index.end())
      continue;
    size_t ctr = ctr_it->second;
    size_t a1 = a1_it->second;
    size_t a2 = a2_it->second;
    size_t a3 = a3_it->second;
    if (!atom_has_finite_xyz(cc.atoms[ctr]) || !atom_has_finite_xyz(cc.atoms[a1]) ||
        !atom_has_finite_xyz(cc.atoms[a2]) || !atom_has_finite_xyz(cc.atoms[a3]))
      continue;
    double vol = calculate_chiral_volume(cc.atoms[ctr].xyz, cc.atoms[a1].xyz,
                                         cc.atoms[a2].xyz, cc.atoms[a3].xyz);
    if (!std::isfinite(vol) || !chir.is_wrong(vol))
      continue;
    struct Option { size_t move, p1, p2; int score; };
    Option options[3] = {
      {a1, a2, a3, candidate_score(a1)},
      {a2, a1, a3, candidate_score(a2)},
      {a3, a1, a2, candidate_score(a3)}
    };
    const Option* best = nullptr;
    for (const Option& opt : options)
      if (opt.score > 0 && (!best || opt.score > best->score))
        best = &opt;
    if (!best)
      continue;
    reflect_across_plane(cc.atoms[best->move].xyz, cc.atoms[ctr].xyz,
                         cc.atoms[best->p1].xyz, cc.atoms[best->p2].xyz);
  }
}

std::vector<size_t> shortest_path_excluding_edge(const std::vector<std::vector<size_t>>& adjacency,
                                                 size_t start, size_t goal,
                                                 size_t ban_a, size_t ban_b) {
  std::vector<int> parent((int) adjacency.size(), -1);
  std::queue<size_t> q;
  q.push(start);
  parent[start] = (int) start;
  while (!q.empty()) {
    size_t cur = q.front();
    q.pop();
    if (cur == goal)
      break;
    for (size_t nb : adjacency[cur]) {
      if ((cur == ban_a && nb == ban_b) || (cur == ban_b && nb == ban_a))
        continue;
      if (parent[nb] != -1)
        continue;
      parent[nb] = (int) cur;
      q.push(nb);
    }
  }
  if (parent[goal] == -1)
    return {};
  std::vector<size_t> path;
  for (size_t cur = goal; cur != start; cur = (size_t) parent[cur])
    path.push_back(cur);
  path.push_back(start);
  std::reverse(path.begin(), path.end());
  return path;
}

std::vector<std::vector<size_t>> detect_small_cycles(const std::vector<std::vector<size_t>>& adjacency) {
  std::set<std::vector<size_t>> unique;
  for (size_t a = 0; a != adjacency.size(); ++a) {
    for (size_t b : adjacency[a]) {
      if (a >= b)
        continue;
      std::vector<size_t> path = shortest_path_excluding_edge(adjacency, a, b, a, b);
      if (path.size() < 3 || path.size() > 6)
        continue;
      std::vector<size_t> cycle = path;
      std::sort(cycle.begin(), cycle.end());
      unique.insert(cycle);
    }
  }
  std::vector<std::vector<size_t>> cycles(unique.begin(), unique.end());
  std::sort(cycles.begin(), cycles.end(),
            [](const std::vector<size_t>& x, const std::vector<size_t>& y) {
              if (x.size() != y.size())
                return x.size() < y.size();
              return x < y;
            });
  return cycles;
}

bool order_cycle_nodes(const std::vector<std::vector<size_t>>& adjacency,
                       const std::vector<size_t>& cycle,
                       std::vector<size_t>& ordered) {
  std::unordered_set<size_t> ring_set(cycle.begin(), cycle.end());
  ordered.clear();
  ordered.push_back(*std::min_element(cycle.begin(), cycle.end()));
  while (ordered.size() < cycle.size()) {
    size_t cur = ordered.back();
    size_t prev = ordered.size() > 1 ? ordered[ordered.size() - 2] : SIZE_MAX;
    size_t next = SIZE_MAX;
    for (size_t nb : adjacency[cur])
      if (ring_set.count(nb) && nb != prev &&
          std::find(ordered.begin(), ordered.end(), nb) == ordered.end()) {
        next = nb;
        break;
      }
    if (next == SIZE_MAX)
      break;
    ordered.push_back(next);
  }
  if (ordered.size() != cycle.size())
    return false;
  return true;
}

bool seed_cycle_coordinates(ChemComp& cc, const std::vector<std::vector<size_t>>& adjacency,
                            const std::vector<size_t>& cycle, double x_shift) {
  if (cycle.size() < 5 || cycle.size() > 6)
    return false;
  for (size_t idx : cycle)
    if (atom_has_finite_xyz(cc.atoms[idx]))
      return false;

  std::vector<size_t> ordered;
  if (!order_cycle_nodes(adjacency, cycle, ordered))
    return false;

  double mean_bond = 0.0;
  int n_bonds = 0;
  for (size_t i = 0; i != ordered.size(); ++i) {
    auto bond = cc.rt.find_bond(cc.atoms[ordered[i]].id,
                                cc.atoms[ordered[(i + 1) % ordered.size()]].id);
    if (bond != cc.rt.bonds.end() && std::isfinite(bond->value)) {
      mean_bond += bond->value;
      ++n_bonds;
    }
  }
  if (n_bonds == 0)
    mean_bond = 1.4;
  else
    mean_bond /= n_bonds;
  double radius = mean_bond / (2.0 * std::sin(pi() / ordered.size()));
  for (size_t i = 0; i != ordered.size(); ++i) {
    double angle = 2.0 * pi() * i / ordered.size();
    cc.atoms[ordered[i]].xyz = Position(x_shift + radius * std::cos(angle),
                                        radius * std::sin(angle), 0.0);
  }
  return true;
}

std::vector<size_t> connected_fragment_from_seed(const std::vector<std::vector<size_t>>& adjacency,
                                                 const std::vector<bool>& in_fragment,
                                                 size_t seed) {
  std::vector<size_t> out;
  if (seed >= adjacency.size() || !in_fragment[seed])
    return out;
  std::vector<bool> seen(adjacency.size(), false);
  std::queue<size_t> q;
  q.push(seed);
  seen[seed] = true;
  while (!q.empty()) {
    size_t cur = q.front();
    q.pop();
    out.push_back(cur);
    for (size_t nb : adjacency[cur])
      if (in_fragment[nb] && !seen[nb]) {
        seen[nb] = true;
        q.push(nb);
      }
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<std::vector<size_t>> detect_plane_fragments(
    const ChemComp& cc,
    const std::map<std::string, size_t>& atom_index,
    const std::vector<std::vector<size_t>>& adjacency) {
  std::vector<std::vector<size_t> > raw;
  for (const Restraints::Plane& plane : cc.rt.planes) {
    std::vector<bool> in_fragment(cc.atoms.size(), false);
    std::vector<size_t> atoms;
    for (const Restraints::AtomId& atom_id : plane.ids) {
      std::map<std::string, size_t>::const_iterator it = atom_index.find(atom_id.atom);
      if (it == atom_index.end())
        continue;
      size_t idx = it->second;
      if (!is_heavy_atom(cc.atoms[idx]))
        continue;
      in_fragment[idx] = true;
      atoms.push_back(idx);
    }
    if (atoms.size() < 4 || atoms.size() > 12)
      continue;
    std::vector<size_t> frag = connected_fragment_from_seed(adjacency, in_fragment, atoms[0]);
    if (frag.size() != atoms.size())
      continue;
    raw.push_back(frag);
  }

  bool changed = true;
  while (changed) {
    changed = false;
    for (size_t i = 0; i != raw.size() && !changed; ++i) {
      std::set<size_t> aset(raw[i].begin(), raw[i].end());
      for (size_t j = i + 1; j != raw.size(); ++j) {
        int overlap = 0;
        for (size_t idx : raw[j])
          if (aset.count(idx))
            ++overlap;
        if (overlap < 2)
          continue;
        aset.insert(raw[j].begin(), raw[j].end());
        raw[i].assign(aset.begin(), aset.end());
        raw.erase(raw.begin() + j);
        changed = true;
        break;
      }
    }
  }

  std::set<std::vector<size_t> > unique;
  for (size_t i = 0; i != raw.size(); ++i) {
    const std::vector<size_t>& frag = raw[i];
    if (frag.size() < 4 || frag.size() > 12)
      continue;
    std::vector<bool> in_fragment(cc.atoms.size(), false);
    for (size_t idx : frag)
      in_fragment[idx] = true;
    int edges = 0;
    for (size_t idx : frag) {
      int local_degree = 0;
      for (size_t nb : adjacency[idx])
        if (in_fragment[nb])
          ++local_degree;
      if (local_degree == 0) {
        edges = -1;
        break;
      }
      edges += local_degree;
    }
    if (edges < 2 * ((int) frag.size() - 1))
      continue;
    unique.insert(frag);
  }
  return std::vector<std::vector<size_t> >(unique.begin(), unique.end());
}

bool seed_plane_fragment_coordinates(ChemComp& cc,
                                     const std::map<std::string, size_t>& atom_index,
                                     const std::vector<std::vector<size_t>>& adjacency,
                                     const std::vector<std::vector<double>>& bond_dist,
                                     const std::vector<size_t>& fragment,
                                     double x_shift) {
  if (fragment.size() < 4 || fragment.size() > 12)
    return false;
  for (size_t idx : fragment)
    if (atom_has_finite_xyz(cc.atoms[idx]))
      return false;

  std::set<size_t> frag_set(fragment.begin(), fragment.end());
  size_t seed = fragment[0];
  size_t second = SIZE_MAX;
  for (size_t nb : adjacency[seed])
    if (frag_set.count(nb)) {
      second = nb;
      break;
    }
  if (second == SIZE_MAX)
    return false;

  cc.atoms[seed].xyz = Position(x_shift, 0, 0);
  cc.atoms[second].xyz = Position(x_shift + get_bond_dist_or_default(cc, adjacency, bond_dist, seed, second),
                                  0, 0);

  size_t third = SIZE_MAX;
  for (size_t nb : adjacency[second])
    if (nb != seed && frag_set.count(nb)) {
      third = nb;
      break;
    }
  if (third != SIZE_MAX) {
    Position pos = arbitrary_position_from_angle_local(
        cc.atoms[seed].xyz, cc.atoms[second].xyz,
        get_bond_dist_or_default(cc, adjacency, bond_dist, second, third),
        angle_rad_or_default(cc, atom_index, adjacency, seed, second, third));
    cc.atoms[third].xyz = Position(pos.x, pos.y, 0.0);
  }

  bool progress = true;
  while (progress) {
    progress = false;
    for (size_t target : fragment) {
      if (atom_has_finite_xyz(cc.atoms[target]))
        continue;
      std::vector<size_t> placed_neighbors;
      for (size_t nb : adjacency[target])
        if (frag_set.count(nb) && atom_has_finite_xyz(cc.atoms[nb]))
          placed_neighbors.push_back(nb);
      if (placed_neighbors.empty())
        continue;

      bool placed = false;
      for (size_t i = 0; i != placed_neighbors.size() && !placed; ++i) {
        size_t center = placed_neighbors[i];
        std::vector<size_t> anchors;
        for (size_t nb : adjacency[center])
          if (nb != target && frag_set.count(nb) && atom_has_finite_xyz(cc.atoms[nb]))
            anchors.push_back(nb);
        if (anchors.size() >= 2) {
          size_t a = anchors[0];
          size_t b = anchors[1];
          std::pair<Position, Position> pair = position_from_two_angles_local(
              cc.atoms[center].xyz, cc.atoms[a].xyz, cc.atoms[b].xyz,
              get_bond_dist_or_default(cc, adjacency, bond_dist, center, target),
              angle_rad_or_default(cc, atom_index, adjacency, a, center, target),
              angle_rad_or_default(cc, atom_index, adjacency, b, center, target));
          Position options[3] = {pair.first, pair.second,
                                 trilaterate_in_plane(cc.atoms[center].xyz,
                                                      std::pow(get_bond_dist_or_default(
                                                                   cc, adjacency, bond_dist, center, target), 2),
                                                      cc.atoms[a].xyz,
                                                      distance_sq_from_angle(
                                                          get_bond_dist_or_default(cc, adjacency, bond_dist, center, target),
                                                          cc.atoms[center].xyz.dist(cc.atoms[a].xyz),
                                                          angle_rad_or_default(cc, atom_index, adjacency, a, center, target)),
                                                      cc.atoms[b].xyz,
                                                      distance_sq_from_angle(
                                                          get_bond_dist_or_default(cc, adjacency, bond_dist, center, target),
                                                          cc.atoms[center].xyz.dist(cc.atoms[b].xyz),
                                                          angle_rad_or_default(cc, atom_index, adjacency, b, center, target)))};
          double best_score = INFINITY;
          Position best(NAN, NAN, NAN);
          for (int k = 0; k != 3; ++k) {
            Position pos = Position(options[k].x, options[k].y, 0.0);
            if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z))
              continue;
            double score = placement_contact_score(cc, pos, target, center);
            if (score < best_score) {
              best_score = score;
              best = pos;
            }
          }
          if (std::isfinite(best.x)) {
            cc.atoms[target].xyz = best;
            placed = true;
          }
        } else if (!anchors.empty()) {
          size_t anchor = anchors[0];
          Position pos = arbitrary_position_from_angle_local(
              cc.atoms[anchor].xyz, cc.atoms[center].xyz,
              get_bond_dist_or_default(cc, adjacency, bond_dist, center, target),
              angle_rad_or_default(cc, atom_index, adjacency, anchor, center, target));
          cc.atoms[target].xyz = Position(pos.x, pos.y, 0.0);
          placed = atom_has_finite_xyz(cc.atoms[target]);
        }
      }
      if (!placed) {
        size_t center = placed_neighbors[0];
        cc.atoms[target].xyz = cc.atoms[center].xyz +
                               Position(get_bond_dist_or_default(cc, adjacency, bond_dist, center, target), 0, 0);
        placed = atom_has_finite_xyz(cc.atoms[target]);
      }
      progress = placed || progress;
    }
  }

  for (size_t idx : fragment)
    if (!atom_has_finite_xyz(cc.atoms[idx]))
      return false;
  return true;
}

void regularize_small_cycles(ChemComp& cc,
                             const std::vector<std::vector<size_t>>& adjacency,
                             const std::vector<std::vector<size_t>>& small_cycles) {
  for (const std::vector<size_t>& cycle : small_cycles) {
    if (cycle.size() < 5 || cycle.size() > 6)
      continue;
    std::vector<size_t> ordered;
    if (!order_cycle_nodes(adjacency, cycle, ordered))
      continue;
    std::vector<Position> positions;
    positions.reserve(ordered.size());
    for (size_t idx : ordered) {
      if (!atom_has_finite_xyz(cc.atoms[idx])) {
        positions.clear();
        break;
      }
      positions.push_back(cc.atoms[idx].xyz);
    }
    if (positions.size() != ordered.size())
      continue;

    Position centroid;
    Vec3 normal;
    if (!plane_from_positions(positions, centroid, normal))
      continue;

    double mean_bond = 0.0;
    int n_bonds = 0;
    for (size_t i = 0; i != ordered.size(); ++i) {
      auto bond = cc.rt.find_bond(cc.atoms[ordered[i]].id,
                                  cc.atoms[ordered[(i + 1) % ordered.size()]].id);
      if (bond != cc.rt.bonds.end()) {
        double value = std::isfinite(bond->value) ? bond->value : bond->value_nucleus;
        if (std::isfinite(value)) {
          mean_bond += value;
          ++n_bonds;
        }
      }
    }
    if (n_bonds == 0)
      continue;
    mean_bond /= n_bonds;
    double radius = mean_bond / (2.0 * std::sin(pi() / ordered.size()));

    Vec3 u = positions[0] - centroid;
    u -= u.dot(normal) * normal;
    if (u.length_sq() < 1e-8)
      continue;
    u = u.normalized();
    Vec3 v = normal.cross(u).normalized();

    double phase = 0.0;
    for (size_t i = 0; i != ordered.size(); ++i) {
      double angle = phase + 2.0 * pi() * i / ordered.size();
      cc.atoms[ordered[i]].xyz = centroid + Position(radius * (std::cos(angle) * u +
                                                                std::sin(angle) * v));
    }
  }
}

Position position_from_axis_angle(const Position& anchor, const Position& center,
                                  double dist, double theta, double phi) {
  Vec3 e1 = (center - anchor).normalized();
  Vec3 trial(1, 0, 0);
  if ((trial - trial.dot(e1) * e1).length_sq() < 1e-6)
    trial = Vec3(0, 1, 0);
  Vec3 e2 = (trial - trial.dot(e1) * e1).normalized();
  Vec3 e3 = e1.cross(e2);
  Vec3 radial = std::cos(phi) * e2 + std::sin(phi) * e3;
  return center + Position(dist * (-std::cos(theta) * e1 + std::sin(theta) * radial));
}

void spread_terminal_children(ChemComp& cc,
                              const std::vector<std::vector<size_t>>& adjacency,
                              const std::map<std::string, size_t>& atom_index) {
  for (size_t center = 0; center != cc.atoms.size(); ++center) {
    if (!atom_has_finite_xyz(cc.atoms[center]))
      continue;
    std::vector<size_t> anchors;
    std::vector<size_t> terminals;
    for (size_t nb : adjacency[center]) {
      if (!atom_has_finite_xyz(cc.atoms[nb]))
        continue;
      if (adjacency[nb].size() == 1)
        terminals.push_back(nb);
      else
        anchors.push_back(nb);
    }
    if (anchors.size() != 1 || terminals.size() < 2 || terminals.size() > 3)
      continue;
    size_t anchor = anchors[0];
    std::sort(terminals.begin(), terminals.end(),
              [&](size_t a, size_t b) { return cc.atoms[a].id < cc.atoms[b].id; });
    std::vector<double> phis;
    if (terminals.size() == 2)
      phis = {0.0, pi()};
    else
      phis = {0.0, 2.0 * pi() / 3.0, 4.0 * pi() / 3.0};
    for (size_t i = 0; i != terminals.size(); ++i) {
      size_t term = terminals[i];
      double dist = 1.5;
      std::vector<Restraints::Bond>::const_iterator bond =
          cc.rt.find_bond(cc.atoms[center].id, cc.atoms[term].id);
      if (bond != cc.rt.bonds.end()) {
        double value = std::isfinite(bond->value) ? bond->value : bond->value_nucleus;
        if (std::isfinite(value))
          dist = value;
      }
      double theta = rad(default_angle_for_center(adjacency, center));
      if (const auto* angle = find_angle_restraint(cc, atom_index, anchor, center, term))
        if (std::isfinite(angle->value))
          theta = rad(angle->value);
      cc.atoms[term].xyz = position_from_axis_angle(cc.atoms[anchor].xyz, cc.atoms[center].xyz,
                                                    dist, theta, phis[i]);
    }
  }
}

double get_bond_dist_or_default(const ChemComp& cc,
                                const std::vector<std::vector<size_t>>& adjacency,
                                const std::vector<std::vector<double>>& bond_dist,
                                size_t a, size_t b) {
  for (size_t i = 0; i != adjacency[a].size(); ++i)
    if (adjacency[a][i] == b)
      return std::isfinite(bond_dist[a][i]) ? bond_dist[a][i] : 1.5;
  return 1.5;
}

const Restraints::Angle* find_angle_by_indices(const ChemComp& cc,
                                               const std::map<std::string, size_t>& atom_index,
                                               size_t a, size_t center, size_t b) {
  return find_angle_restraint(cc, atom_index, a, center, b);
}

double angle_rad_or_default(const ChemComp& cc,
                            const std::map<std::string, size_t>& atom_index,
                            const std::vector<std::vector<size_t>>& adjacency,
                            size_t a, size_t center, size_t b) {
  if (const auto* angle = find_angle_by_indices(cc, atom_index, a, center, b))
    if (std::isfinite(angle->value))
      return rad(angle->value);
  return rad(default_angle_for_center(adjacency, center));
}

const Restraints::Bond* find_bond_by_indices(const ChemComp& cc, size_t a, size_t b) {
  return cc.rt.find_bond(cc.atoms[a].id, cc.atoms[b].id) != cc.rt.bonds.end()
             ? &*cc.rt.find_bond(cc.atoms[a].id, cc.atoms[b].id)
             : nullptr;
}

const Restraints::Torsion* find_torsion_with_hydrogen(const ChemComp& cc,
                                                      const std::map<std::string, size_t>& atom_index,
                                                      size_t h_idx, size_t center_idx,
                                                      size_t heavy_idx,
                                                      size_t* tau_end_idx,
                                                      int* period) {
  for (const auto& tor : cc.rt.torsions) {
    auto it1 = atom_index.find(tor.id1.atom);
    auto it2 = atom_index.find(tor.id2.atom);
    auto it3 = atom_index.find(tor.id3.atom);
    auto it4 = atom_index.find(tor.id4.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end() ||
        it3 == atom_index.end() || it4 == atom_index.end())
      continue;
    size_t i1 = it1->second;
    size_t i2 = it2->second;
    size_t i3 = it3->second;
    size_t i4 = it4->second;
    if (i1 == h_idx && i2 == center_idx && i3 == heavy_idx && is_heavy_atom(cc.atoms[i4])) {
      if (tau_end_idx) *tau_end_idx = i4;
      if (period) *period = tor.period;
      return &tor;
    }
    if (i4 == h_idx && i3 == center_idx && i2 == heavy_idx && is_heavy_atom(cc.atoms[i1])) {
      if (tau_end_idx) *tau_end_idx = i1;
      if (period) *period = tor.period;
      return &tor;
    }
  }
  return nullptr;
}

const Restraints::Plane* find_plane_with_atoms(const ChemComp& cc,
                                               const std::map<std::string, size_t>& atom_index,
                                               size_t center_idx,
                                               size_t h_idx,
                                               size_t heavy_idx,
                                               size_t* other_idx) {
  for (const auto& plane : cc.rt.planes) {
    bool has_center = false;
    bool has_h = false;
    bool has_heavy = false;
    size_t found_other = SIZE_MAX;
    for (const auto& atom_id : plane.ids) {
      auto it = atom_index.find(atom_id.atom);
      if (it == atom_index.end())
        continue;
      size_t idx = it->second;
      if (idx == center_idx)
        has_center = true;
      else if (idx == h_idx)
        has_h = true;
      else if (idx == heavy_idx)
        has_heavy = true;
      else if (is_heavy_atom(cc.atoms[idx]) && atom_has_finite_xyz(cc.atoms[idx]))
        found_other = idx;
    }
    if (has_center && has_h && has_heavy && found_other != SIZE_MAX) {
      if (other_idx) *other_idx = found_other;
      return &plane;
    }
  }
  return nullptr;
}

const Restraints::Chirality* find_chirality_for_center(const ChemComp& cc,
                                                       const std::map<std::string, size_t>& atom_index,
                                                       size_t center_idx) {
  for (const auto& chir : cc.rt.chirs) {
    auto it = atom_index.find(chir.id_ctr.atom);
    if (it != atom_index.end() && it->second == center_idx)
      return &chir;
  }
  return nullptr;
}

double missing_angle_2h_tetrahedral_chemcomp(const ChemComp& cc,
                                             const std::map<std::string, size_t>& atom_index,
                                             size_t center_idx,
                                             const std::vector<size_t>& hs,
                                             double alpha, double theta) {
  if (hs.size() == 1)
    return 2 * pi() - alpha - theta;
  if (hs.size() == 2) {
    if (const auto* hh = find_angle_by_indices(cc, atom_index, hs[0], center_idx, hs[1])) {
      double xh = std::cos(alpha);
      double zh = std::sin(0.5 * hh->radians());
      double yh = std::sqrt(std::max(0.0, 1 - xh * xh - zh * zh));
      double st = std::sin(theta);
      double ct = std::cos(theta);
      Vec3 ah(xh, -yh * st, zh * st);
      return std::acos((ct * ah.x + st * ah.y) / ah.length());
    }
  }
  return 0.0;
}

void place_hydrogens_from_restraints(ChemComp& cc,
                                     const std::map<std::string, size_t>& atom_index,
                                     const std::vector<std::vector<size_t>>& adjacency,
                                     const std::vector<std::vector<double>>& bond_dist) {
  for (size_t center = 0; center != cc.atoms.size(); ++center) {
    if (!is_heavy_atom(cc.atoms[center]) || !atom_has_finite_xyz(cc.atoms[center]))
      continue;

    std::vector<size_t> known;
    std::vector<size_t> hs;
    for (size_t nb : adjacency[center]) {
      if (cc.atoms[nb].is_hydrogen())
        hs.push_back(nb);
      else if (atom_has_finite_xyz(cc.atoms[nb]))
        known.push_back(nb);
    }
    if (hs.empty())
      continue;

    for (size_t h : hs)
      cc.atoms[h].xyz = Position(NAN, NAN, NAN);

    auto h_dist = [&](size_t h) {
      return get_bond_dist_or_default(cc, adjacency, bond_dist, center, h);
    };

    if (known.empty()) {
      cc.atoms[hs[0]].xyz = cc.atoms[center].xyz + Position(h_dist(hs[0]), 0, 0);
      if (hs.size() > 1) {
        double theta = pi();
        if (const auto* ang = find_angle_by_indices(cc, atom_index, hs[1], center, hs[0]))
          theta = ang->radians();
        cc.atoms[hs[1]].xyz = cc.atoms[center].xyz +
                              Position(h_dist(hs[1]) * std::cos(theta),
                                       h_dist(hs[1]) * std::sin(theta), 0);
      }
      if (hs.size() > 2) {
        double theta0 = rad(109.47122);
        double theta1 = rad(109.47122);
        if (const auto* ang0 = find_angle_by_indices(cc, atom_index, hs[2], center, hs[0]))
          theta0 = ang0->radians();
        if (const auto* ang1 = find_angle_by_indices(cc, atom_index, hs[2], center, hs[1]))
          theta1 = ang1->radians();
        auto pos = position_from_two_angles_local(cc.atoms[center].xyz,
                                                  cc.atoms[hs[0]].xyz, cc.atoms[hs[1]].xyz,
                                                  h_dist(hs[2]), theta0, theta1);
        cc.atoms[hs[2]].xyz = pos.first;
        if (hs.size() == 4)
          cc.atoms[hs[3]].xyz = pos.second;
      }
      continue;
    }

    if (known.size() == 1) {
      size_t heavy = known[0];
      size_t h0 = hs[0];
      const auto* angle = find_angle_by_indices(cc, atom_index, h0, center, heavy);
      if (!angle)
        continue;
      double theta = angle->radians();
      if (hs.size() == 3) {
        double hh_sum = 0.0;
        int hh_count = 0;
        for (size_t i = 0; i != hs.size(); ++i)
          for (size_t j = i + 1; j != hs.size(); ++j)
            if (const auto* hh = find_angle_by_indices(cc, atom_index, hs[i], center, hs[j])) {
              hh_sum += hh->radians();
              ++hh_count;
            }
        if (hh_count > 0)
          theta = cone_theta_from_hh_angle(hh_sum / hh_count);
      }
      double tau = 0.0;
      int period = 0;
      size_t tau_end = SIZE_MAX;
      if (find_plane_with_atoms(cc, atom_index, center, h0, heavy, &tau_end))
        tau = 0.0;
      else if (const auto* tor = find_torsion_with_hydrogen(cc, atom_index, h0, center, heavy,
                                                            &tau_end, &period))
        tau = rad(tor->value);

      if (tau_end != SIZE_MAX)
        cc.atoms[h0].xyz = position_from_angle_and_torsion_local(
            cc.atoms[tau_end].xyz, cc.atoms[heavy].xyz, cc.atoms[center].xyz,
            h_dist(h0), theta, tau);
      else
        cc.atoms[h0].xyz = arbitrary_position_from_angle_local(
            cc.atoms[heavy].xyz, cc.atoms[center].xyz, h_dist(h0), theta);
      if (!atom_has_finite_xyz(cc.atoms[h0]))
        continue;
      if (hs.size() == 2) {
        Vec3 axis = (cc.atoms[heavy].xyz - cc.atoms[center].xyz).normalized();
        Vec3 perpendicular = get_vector_to_line(cc.atoms[h0].xyz, cc.atoms[center].xyz, axis);
        cc.atoms[hs[1]].xyz = cc.atoms[h0].xyz + Position(2 * perpendicular);
      } else if (hs.size() == 3) {
        int idx = 0;
        for (int i : {1, 2}) {
          size_t tau_tmp = SIZE_MAX;
          if (find_torsion_with_hydrogen(cc, atom_index, hs[i], center, heavy, &tau_tmp, &period)) {
            idx = i;
            cc.atoms[hs[i]].xyz = cc.atoms[h0].xyz;
          }
        }
        Vec3 axis = (cc.atoms[heavy].xyz - cc.atoms[center].xyz).normalized();
        Vec3 v1 = cc.atoms[h0].xyz - cc.atoms[center].xyz;
        Vec3 v2 = rotate_about_axis(v1, axis, rad(120));
        Vec3 v3 = rotate_about_axis(v1, axis, rad(-120));
        cc.atoms[hs[(idx + 1) % 3]].xyz = cc.atoms[center].xyz + Position(v2);
        cc.atoms[hs[(idx + 2) % 3]].xyz = cc.atoms[center].xyz + Position(v3);
      }
      continue;
    }

    if (known.size() == 2) {
      if (hs.size() >= 3)
        continue;
      const auto* ang1 = find_angle_by_indices(cc, atom_index, hs[0], center, known[0]);
      const auto* ang2 = find_angle_by_indices(cc, atom_index, hs[0], center, known[1]);
      const auto* ang3 = find_angle_by_indices(cc, atom_index, known[0], center, known[1]);
      double theta1 = ang1 ? ang1->radians() : 0.0;
      double theta2 = ang2 ? ang2->radians() : 0.0;
      double theta3 = ang3 ? ang3->radians() : 0.0;
      if (!ang3) {
        const Restraints::Bond* aptr = find_bond_by_indices(cc, center, known[1]);
        const Restraints::Bond* bptr = find_bond_by_indices(cc, known[0], known[1]);
        const Restraints::Bond* cptr = find_bond_by_indices(cc, known[0], center);
        if (!aptr || !bptr || !cptr)
          continue;
        double av = std::isfinite(aptr->value) ? aptr->value : aptr->value_nucleus;
        double bv = std::isfinite(bptr->value) ? bptr->value : bptr->value_nucleus;
        double cv = std::isfinite(cptr->value) ? cptr->value : cptr->value_nucleus;
        theta3 = angle_in_triangle(av, bv, cv);
      }
      if (theta1 == 0.0 && known[0] < cc.atoms.size() && cc.atoms[known[0]].el.is_metal())
        theta1 = missing_angle_2h_tetrahedral_chemcomp(cc, atom_index, center, hs, theta2, theta3);
      if (theta2 == 0.0 && known[1] < cc.atoms.size() && cc.atoms[known[1]].el.is_metal())
        theta2 = missing_angle_2h_tetrahedral_chemcomp(cc, atom_index, center, hs, theta1, theta3);
      if (theta1 == 0.0 || theta2 == 0.0)
        continue;

      if (theta1 + theta2 + theta3 > rad(357.0)) {
        Vec3 v12 = cc.atoms[known[0]].xyz - cc.atoms[center].xyz;
        Vec3 v13 = cc.atoms[known[1]].xyz - cc.atoms[center].xyz;
        double cur_theta3 = v12.angle(v13);
        double ratio = (2 * pi() - cur_theta3) / (theta1 + theta2);
        Vec3 axis = v13.cross(v12).normalized();
        Vec3 v14 = rotate_about_axis(v12, axis, theta1 * ratio);
        cc.atoms[hs[0]].xyz = cc.atoms[center].xyz +
                              Position(h_dist(hs[0]) / v14.length() * v14);
        continue;
      }

      double hh_half = 0.0;
      if (hs.size() == 2)
        if (const auto* hh = find_angle_by_indices(cc, atom_index, hs[0], center, hs[1]))
          hh_half = 0.5 * hh->radians();
      if (hh_half == 0.0)
        hh_half = calculate_tetrahedral_delta(theta3, theta1, theta2);

      double c0 = std::cos(theta3);
      double c1 = std::cos(theta1);
      double c2 = std::cos(theta2);
      double den = 1 / (1 - c0 * c0);
      double a = den * (c1 - c0 * c2);
      double b = den * (c2 - c0 * c1);
      Vec3 u10 = (cc.atoms[known[0]].xyz - cc.atoms[center].xyz).normalized();
      Vec3 u20 = (cc.atoms[known[1]].xyz - cc.atoms[center].xyz).normalized();
      Vec3 v = u10.cross(u20);
      if (std::isnan(v.x))
        continue;
      Vec3 d = a * u10 + b * u20;
      double dist_sin = h_dist(hs[0]) * std::sin(hh_half);
      double dist_cos = h_dist(hs[0]) * std::cos(hh_half);
      Vec3 v0s = v.changed_magnitude(dist_sin);
      Vec3 d0c = d.changed_magnitude(dist_cos);
      cc.atoms[hs[0]].xyz = cc.atoms[center].xyz + Position(d0c + v0s);
      Position other_pos = cc.atoms[center].xyz + Position(d0c - v0s);

      if (hs.size() == 1) {
        const auto* chir = find_chirality_for_center(cc, atom_index, center);
        if (chir && chir->sign != ChiralityType::Both) {
          double vol = calculate_chiral_volume(cc.atoms[center].xyz,
                                               cc.atoms[known[0]].xyz,
                                               cc.atoms[known[1]].xyz,
                                               cc.atoms[hs[0]].xyz);
          if (chir->is_wrong(vol))
            cc.atoms[hs[0]].xyz = other_pos;
        }
      } else {
        cc.atoms[hs[1]].xyz = other_pos;
      }
      continue;
    }

    if (hs.size() > 1)
      continue;
    Vec3 u10 = (cc.atoms[known[0]].xyz - cc.atoms[center].xyz).normalized();
    Vec3 u20 = (cc.atoms[known[1]].xyz - cc.atoms[center].xyz).normalized();
    Vec3 u30 = (cc.atoms[known[2]].xyz - cc.atoms[center].xyz).normalized();
    auto cos_tetrahedral = [&](size_t n) {
      const auto* angle = find_angle_by_indices(cc, atom_index, known[n], center, hs[0]);
      return angle ? std::cos(angle->radians()) : -1.0 / 3.0;
    };
    SMat33<double> m{1., 1., 1., u10.dot(u20), u10.dot(u30), u20.dot(u30)};
    Vec3 rhs(cos_tetrahedral(0), cos_tetrahedral(1), cos_tetrahedral(2));
    double det = m.determinant();
    if (std::fabs(det) < 1e-12)
      continue;
    Vec3 abc = m.inverse_(det).multiply(rhs);
    Vec3 h_dir = abc.x * u10 + abc.y * u20 + abc.z * u30;
    cc.atoms[hs[0]].xyz = cc.atoms[center].xyz + Position(h_dir.changed_magnitude(h_dist(hs[0])));
  }
}

bool place_terminal_heavy_atom(ChemComp& cc,
                               const std::map<std::string, size_t>& atom_index,
                               const std::vector<std::vector<size_t>>& adjacency,
                               const std::vector<std::vector<double>>& bond_dist,
                               size_t target) {
  if (atom_has_finite_xyz(cc.atoms[target]))
    return false;
  size_t center = SIZE_MAX;
  for (size_t nb : adjacency[target])
    if (is_heavy_atom(cc.atoms[nb]) && atom_has_finite_xyz(cc.atoms[nb])) {
      center = nb;
      break;
    }
  if (center == SIZE_MAX)
    return false;
  double dist = get_bond_dist_or_default(cc, adjacency, bond_dist, target, center);

  std::vector<size_t> anchors;
  for (size_t anchor_idx : adjacency[center])
    if (anchor_idx != target && is_heavy_atom(cc.atoms[anchor_idx]) &&
        atom_has_finite_xyz(cc.atoms[anchor_idx]))
      anchors.push_back(anchor_idx);

  if (!anchors.empty()) {
    size_t anchor = anchors[0];
    double theta = rad(default_angle_for_center(adjacency, center));
    if (const auto* angle = find_angle_by_indices(cc, atom_index, anchor, center, target))
      if (std::isfinite(angle->value))
        theta = rad(angle->value);
    for (size_t torsion_anchor : adjacency[anchor]) {
      if (torsion_anchor == center || !is_heavy_atom(cc.atoms[torsion_anchor]) ||
          !atom_has_finite_xyz(cc.atoms[torsion_anchor]))
        continue;
      if (const auto* tor = find_torsion_restraint(cc, atom_index,
                                                   torsion_anchor, anchor,
                                                   center, target)) {
        double tau = std::isfinite(tor->value) ? rad(tor->value) : pi();
        cc.atoms[target].xyz = position_from_angle_and_torsion_local(
            cc.atoms[torsion_anchor].xyz, cc.atoms[anchor].xyz,
            cc.atoms[center].xyz, dist, theta, tau);
        return atom_has_finite_xyz(cc.atoms[target]);
      }
    }
    if (anchors.size() >= 2) {
      size_t a = anchors[0];
      size_t b = anchors[1];
      auto pair = position_from_two_angles_local(cc.atoms[center].xyz,
                                                 cc.atoms[a].xyz,
                                                 cc.atoms[b].xyz,
                                                 dist,
                                                 angle_rad_or_default(cc, atom_index, adjacency, a, center, target),
                                                 angle_rad_or_default(cc, atom_index, adjacency, b, center, target));
      if (std::isfinite(pair.first.x) && std::isfinite(pair.second.x)) {
        double s1 = placement_contact_score(cc, pair.first, target, center);
        double s2 = placement_contact_score(cc, pair.second, target, center);
        cc.atoms[target].xyz = s1 <= s2 ? pair.first : pair.second;
        return true;
      }
    }
    cc.atoms[target].xyz = arbitrary_position_from_angle_local(
        cc.atoms[anchor].xyz, cc.atoms[center].xyz, dist, theta);
    return atom_has_finite_xyz(cc.atoms[target]);
  }

  cc.atoms[target].xyz = cc.atoms[center].xyz + Position(dist, 0, 0);
  return atom_has_finite_xyz(cc.atoms[target]);
}

bool place_bridge_heavy_atom(ChemComp& cc,
                             const std::map<std::string, size_t>& atom_index,
                             const std::vector<std::vector<size_t>>& adjacency,
                             const std::vector<std::vector<double>>& bond_dist,
                             size_t target) {
  if (atom_has_finite_xyz(cc.atoms[target]) || !is_heavy_atom(cc.atoms[target]))
    return false;
  std::vector<size_t> centers;
  for (size_t nb : adjacency[target])
    if (is_heavy_atom(cc.atoms[nb]) && atom_has_finite_xyz(cc.atoms[nb]))
      centers.push_back(nb);
  if (centers.size() < 2)
    return false;

  Position best_pos(NAN, NAN, NAN);
  double best_score = INFINITY;

  for (size_t ia = 0; ia != centers.size(); ++ia) {
    for (size_t ib = ia + 1; ib != centers.size(); ++ib) {
      size_t a = centers[ia];
      size_t b = centers[ib];
      double da = get_bond_dist_or_default(cc, adjacency, bond_dist, target, a);
      double db = get_bond_dist_or_default(cc, adjacency, bond_dist, target, b);
      for (size_t anchor : adjacency[a]) {
        if (anchor == target || anchor == b || !is_heavy_atom(cc.atoms[anchor]) ||
            !atom_has_finite_xyz(cc.atoms[anchor]))
          continue;
        double theta = angle_rad_or_default(cc, atom_index, adjacency, anchor, a, target);
        double dc = cc.atoms[a].xyz.dist(cc.atoms[anchor].xyz);
        double dct_sq = distance_sq_from_angle(da, dc, theta);
        if (!(dct_sq > 1e-8))
          continue;
        std::pair<Position, Position> pair = trilaterate_local(cc.atoms[a].xyz, da * da,
                                                               cc.atoms[b].xyz, db * db,
                                                               cc.atoms[anchor].xyz, dct_sq);
        Position planar = trilaterate_in_plane(cc.atoms[a].xyz, da * da,
                                               cc.atoms[b].xyz, db * db,
                                               cc.atoms[anchor].xyz, dct_sq);
        Position options[3] = {pair.first, pair.second, planar};
        for (int k = 0; k != 3; ++k) {
          const Position& pos = options[k];
          if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z))
            continue;
          double score = placement_contact_score(cc, pos, target, a) +
                         placement_contact_score(cc, pos, target, b);
          for (size_t other : centers) {
            if (other == a || other == b)
              continue;
            double d = pos.dist(cc.atoms[other].xyz);
            double want = get_bond_dist_or_default(cc, adjacency, bond_dist, target, other);
            score += 10.0 * std::fabs(d - want);
          }
          if (score < best_score) {
            best_score = score;
            best_pos = pos;
          }
        }
      }
      for (size_t anchor : adjacency[b]) {
        if (anchor == target || anchor == a || !is_heavy_atom(cc.atoms[anchor]) ||
            !atom_has_finite_xyz(cc.atoms[anchor]))
          continue;
        double theta = angle_rad_or_default(cc, atom_index, adjacency, anchor, b, target);
        double dc = cc.atoms[b].xyz.dist(cc.atoms[anchor].xyz);
        double dct_sq = distance_sq_from_angle(db, dc, theta);
        if (!(dct_sq > 1e-8))
          continue;
        std::pair<Position, Position> pair = trilaterate_local(cc.atoms[a].xyz, da * da,
                                                               cc.atoms[b].xyz, db * db,
                                                               cc.atoms[anchor].xyz, dct_sq);
        Position planar = trilaterate_in_plane(cc.atoms[a].xyz, da * da,
                                               cc.atoms[b].xyz, db * db,
                                               cc.atoms[anchor].xyz, dct_sq);
        Position options[3] = {pair.first, pair.second, planar};
        for (int k = 0; k != 3; ++k) {
          const Position& pos = options[k];
          if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z))
            continue;
          double score = placement_contact_score(cc, pos, target, a) +
                         placement_contact_score(cc, pos, target, b);
          for (size_t other : centers) {
            if (other == a || other == b)
              continue;
            double d = pos.dist(cc.atoms[other].xyz);
            double want = get_bond_dist_or_default(cc, adjacency, bond_dist, target, other);
            score += 10.0 * std::fabs(d - want);
          }
          if (score < best_score) {
            best_score = score;
            best_pos = pos;
          }
        }
      }
    }
  }
  if (!std::isfinite(best_pos.x))
    return false;
  cc.atoms[target].xyz = best_pos;
  return true;
}

bool place_trivalent_junction_atom(ChemComp& cc,
                                   const std::map<std::string, size_t>& atom_index,
                                   const std::vector<std::vector<size_t>>& adjacency,
                                   const std::vector<std::vector<double>>& bond_dist,
                                   size_t target) {
  if (atom_has_finite_xyz(cc.atoms[target]) || !is_heavy_atom(cc.atoms[target]))
    return false;
  std::vector<size_t> centers;
  for (size_t nb : adjacency[target])
    if (is_heavy_atom(cc.atoms[nb]) && atom_has_finite_xyz(cc.atoms[nb]))
      centers.push_back(nb);
  if (centers.size() != 3)
    return false;
  double d1 = get_bond_dist_or_default(cc, adjacency, bond_dist, target, centers[0]);
  double d2 = get_bond_dist_or_default(cc, adjacency, bond_dist, target, centers[1]);
  double d3 = get_bond_dist_or_default(cc, adjacency, bond_dist, target, centers[2]);
  std::pair<Position, Position> pair = trilaterate_local(cc.atoms[centers[0]].xyz, d1 * d1,
                                                         cc.atoms[centers[1]].xyz, d2 * d2,
                                                         cc.atoms[centers[2]].xyz, d3 * d3);
  Position planar = trilaterate_in_plane(cc.atoms[centers[0]].xyz, d1 * d1,
                                         cc.atoms[centers[1]].xyz, d2 * d2,
                                         cc.atoms[centers[2]].xyz, d3 * d3);
  Position best_pos(NAN, NAN, NAN);
  double best_score = INFINITY;
  Position options[3] = {pair.first, pair.second, planar};
  for (int k = 0; k != 3; ++k) {
    const Position& pos = options[k];
    if (!std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z))
      continue;
    double score = placement_contact_score(cc, pos, target, centers[0]);
    for (size_t i = 0; i != centers.size(); ++i) {
      for (size_t j = i + 1; j != centers.size(); ++j) {
        const Restraints::Angle* ang =
            find_angle_restraint(cc, atom_index, centers[i], target, centers[j]);
        if (ang && std::isfinite(ang->value)) {
          double actual = calculate_angle(cc.atoms[centers[i]].xyz, pos, cc.atoms[centers[j]].xyz);
          score += 5.0 * std::fabs(deg(actual) - ang->value);
        }
      }
    }
    for (const Restraints::Plane& plane : cc.rt.planes) {
      bool has_target = false;
      std::vector<Position> others;
      for (const Restraints::AtomId& atom_id : plane.ids) {
        std::map<std::string, size_t>::const_iterator it = atom_index.find(atom_id.atom);
        if (it == atom_index.end())
          continue;
        if (it->second == target) {
          has_target = true;
        } else if (atom_has_finite_xyz(cc.atoms[it->second])) {
          others.push_back(cc.atoms[it->second].xyz);
        }
      }
      if (!has_target || others.size() < 3)
        continue;
      Position centroid;
      Vec3 normal;
      if (!plane_from_positions(others, centroid, normal))
        continue;
      score += 20.0 * std::fabs((pos - centroid).dot(normal));
    }
    if (score < best_score) {
      best_score = score;
      best_pos = pos;
    }
  }
  if (!std::isfinite(best_pos.x))
    return false;
  cc.atoms[target].xyz = best_pos;
  return true;
}

void regrow_terminal_heavy_atoms(ChemComp& cc,
                                 const std::map<std::string, size_t>& atom_index,
                                 const std::vector<std::vector<size_t>>& adjacency,
                                 const std::vector<std::vector<double>>& bond_dist,
                                 const std::vector<bool>& ring_atom) {
  std::vector<int> heavy_leaf_count(cc.atoms.size(), 0);
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    if (!is_heavy_atom(cc.atoms[i]) || ring_atom[i])
      continue;
    int heavy_n = 0;
    size_t parent = SIZE_MAX;
    for (size_t nb : adjacency[i])
      if (is_heavy_atom(cc.atoms[nb])) {
        ++heavy_n;
        parent = nb;
      }
    if (heavy_n == 1 && parent != SIZE_MAX)
      ++heavy_leaf_count[parent];
  }
  std::vector<size_t> targets;
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    if (!is_heavy_atom(cc.atoms[i]) || ring_atom[i])
      continue;
    int heavy_n = 0;
    size_t parent = SIZE_MAX;
    for (size_t nb : adjacency[i])
      if (is_heavy_atom(cc.atoms[nb])) {
        ++heavy_n;
        parent = nb;
      }
    if (heavy_n == 1 && parent != SIZE_MAX && heavy_leaf_count[parent] == 1 &&
        (ring_atom[parent] || adjacency[parent].size() > 2))
      targets.push_back(i);
  }
  for (size_t idx : targets)
    cc.atoms[idx].xyz = Position(NAN, NAN, NAN);
  bool progress = true;
  while (progress) {
    progress = false;
    for (size_t idx : targets)
      progress = place_terminal_heavy_atom(cc, atom_index, adjacency, bond_dist, idx) || progress;
  }
}

void regrow_bridge_heavy_atoms(ChemComp& cc,
                               const std::map<std::string, size_t>& atom_index,
                               const std::vector<std::vector<size_t>>& adjacency,
                               const std::vector<std::vector<double>>& bond_dist,
                               const std::vector<bool>& ring_atom) {
  std::vector<size_t> targets;
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    if (!is_heavy_atom(cc.atoms[i]) || ring_atom[i])
      continue;
    int heavy_n = 0;
    int hydrogen_n = 0;
    bool anchored = false;
    for (size_t nb : adjacency[i]) {
      if (is_heavy_atom(cc.atoms[nb])) {
        ++heavy_n;
        for (size_t anchor : adjacency[nb]) {
          if (anchor != i && is_heavy_atom(cc.atoms[anchor]) &&
              atom_has_finite_xyz(cc.atoms[anchor])) {
            anchored = true;
            break;
          }
        }
      } else if (cc.atoms[nb].is_hydrogen()) {
        ++hydrogen_n;
      }
    }
    if (heavy_n == 2 && hydrogen_n == 1 && anchored)
      targets.push_back(i);
  }
  for (size_t idx : targets)
    cc.atoms[idx].xyz = Position(NAN, NAN, NAN);
  bool progress = true;
  while (progress) {
    progress = false;
    for (size_t idx : targets)
      progress = place_bridge_heavy_atom(cc, atom_index, adjacency, bond_dist, idx) || progress;
  }
}

void regrow_trivalent_junction_atoms(ChemComp& cc,
                                     const std::map<std::string, size_t>& atom_index,
                                     const std::vector<std::vector<size_t>>& adjacency,
                                     const std::vector<std::vector<double>>& bond_dist,
                                     const std::vector<bool>& ring_atom) {
  std::vector<size_t> targets;
  for (size_t i = 0; i != cc.atoms.size(); ++i) {
    if (!is_heavy_atom(cc.atoms[i]) || ring_atom[i] || !atom_in_plane_restraint(cc, atom_index, i))
      continue;
    int heavy_n = 0;
    int hydrogen_n = 0;
    double max_bond_error = 0.0;
    bool has_terminal_heavy_neighbor = false;
    for (size_t nb : adjacency[i]) {
      if (is_heavy_atom(cc.atoms[nb])) {
        ++heavy_n;
        if (adjacency[nb].size() == 1)
          has_terminal_heavy_neighbor = true;
        if (atom_has_finite_xyz(cc.atoms[i]) && atom_has_finite_xyz(cc.atoms[nb])) {
          double want = get_bond_dist_or_default(cc, adjacency, bond_dist, i, nb);
          max_bond_error = std::max(max_bond_error,
                                    std::fabs(cc.atoms[i].xyz.dist(cc.atoms[nb].xyz) - want));
        }
      } else if (cc.atoms[nb].is_hydrogen()) {
        ++hydrogen_n;
      }
    }
    if (heavy_n == 3 && hydrogen_n == 0 && !has_terminal_heavy_neighbor &&
        max_bond_error > 0.5)
      targets.push_back(i);
  }
  for (size_t idx : targets)
    cc.atoms[idx].xyz = Position(NAN, NAN, NAN);
  bool progress = true;
  while (progress) {
    progress = false;
    for (size_t idx : targets)
      progress = place_trivalent_junction_atom(cc, atom_index, adjacency, bond_dist, idx) || progress;
  }
}

} // namespace

int generate_chemcomp_xyz_from_restraints(ChemComp& cc) {
  size_t n = cc.atoms.size();
  if (n == 0)
    return 0;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i != n; ++i)
    atom_index[cc.atoms[i].id] = i;

  std::vector<std::vector<size_t>> adjacency(n);
  std::vector<std::vector<double>> bond_dist(n);
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    adjacency[it1->second].push_back(it2->second);
    adjacency[it2->second].push_back(it1->second);
    double dist = std::isfinite(bond.value) ? bond.value : bond.value_nucleus;
    bond_dist[it1->second].push_back(dist);
    bond_dist[it2->second].push_back(dist);
  }

  auto get_bond_dist = [&](size_t a, size_t b) {
    for (size_t i = 0; i != adjacency[a].size(); ++i)
      if (adjacency[a][i] == b)
        return std::isfinite(bond_dist[a][i]) ? bond_dist[a][i] : 1.5;
    return 1.5;
  };

  auto get_angle_rad = [&](size_t a, size_t center, size_t b) {
    if (const auto* angle = find_angle_restraint(cc, atom_index, a, center, b))
      if (std::isfinite(angle->value))
        return rad(angle->value);
    return rad(default_angle_for_center(adjacency, center));
  };

  auto count_finite = [&](bool heavy_only) {
    int count = 0;
    for (const auto& atom : cc.atoms)
      if ((!heavy_only || is_heavy_atom(atom)) && atom_has_finite_xyz(atom))
        ++count;
    return count;
  };

  for (auto& atom : cc.atoms)
    atom.xyz = Position(NAN, NAN, NAN);

  std::vector<std::vector<size_t>> small_cycles = detect_small_cycles(adjacency);
  std::vector<std::vector<size_t>> plane_fragments =
      detect_plane_fragments(cc, atom_index, adjacency);
  std::vector<int> plane_membership(n, 0);
  for (const Restraints::Plane& plane : cc.rt.planes)
    for (const Restraints::AtomId& atom_id : plane.ids) {
      std::map<std::string, size_t>::const_iterator it = atom_index.find(atom_id.atom);
      if (it != atom_index.end() && is_heavy_atom(cc.atoms[it->second]))
        ++plane_membership[it->second];
    }
  std::vector<bool> ring_atom(n, false);
  for (const std::vector<size_t>& cycle : small_cycles)
    for (size_t idx : cycle)
      ring_atom[idx] = true;

  auto place_candidate = [&](size_t target, const Position& pos) {
    cc.atoms[target].xyz = pos;
  };
  auto finite_pos = [](const Position& pos) {
    return std::isfinite(pos.x) && std::isfinite(pos.y) && std::isfinite(pos.z);
  };

  auto choose_two_angle_solution = [&](size_t target, size_t center,
                                       const std::pair<Position, Position>& pair) {
    double s1 = placement_contact_score(cc, pair.first, target, center);
    double s2 = placement_contact_score(cc, pair.second, target, center);
    place_candidate(target, s1 <= s2 ? pair.first : pair.second);
  };

  auto try_place = [&](size_t target) {
    if (!is_heavy_atom(cc.atoms[target]))
      return false;
    if (atom_has_finite_xyz(cc.atoms[target]))
      return false;
    for (size_t i = 0; i != adjacency[target].size(); ++i) {
      size_t center = adjacency[target][i];
      if (!atom_has_finite_xyz(cc.atoms[center]))
        continue;
      double dist = get_bond_dist(target, center);

      for (size_t anchor_idx : adjacency[center]) {
        if (anchor_idx == target || !atom_has_finite_xyz(cc.atoms[anchor_idx]))
          continue;
        double theta = get_angle_rad(anchor_idx, center, target);
        for (size_t torsion_anchor : adjacency[anchor_idx]) {
          if (torsion_anchor == center || !atom_has_finite_xyz(cc.atoms[torsion_anchor]))
            continue;
          if (const auto* tor = find_torsion_restraint(cc, atom_index,
                                                       torsion_anchor, anchor_idx,
                                                       center, target)) {
            double tau = std::isfinite(tor->value) ? rad(tor->value) : pi();
            place_candidate(target, position_from_angle_and_torsion_local(
                cc.atoms[torsion_anchor].xyz, cc.atoms[anchor_idx].xyz,
                cc.atoms[center].xyz, dist, theta, tau));
            if (atom_has_finite_xyz(cc.atoms[target]))
              return true;
          }
        }
      }

      std::vector<size_t> anchors;
      for (size_t anchor_idx : adjacency[center]) {
        if (anchor_idx != target && atom_has_finite_xyz(cc.atoms[anchor_idx]))
          anchors.push_back(anchor_idx);
      }
      if (anchors.size() >= 2) {
        size_t a = anchors[0];
        size_t b = anchors[1];
        auto pair = position_from_two_angles_local(cc.atoms[center].xyz,
                                                   cc.atoms[a].xyz,
                                                   cc.atoms[b].xyz,
                                                   dist,
                                                   get_angle_rad(a, center, target),
                                                   get_angle_rad(b, center, target));
        if (finite_pos(pair.first) && finite_pos(pair.second)) {
          choose_two_angle_solution(target, center, pair);
          if (atom_has_finite_xyz(cc.atoms[target]))
            return true;
        }
      }

      if (!anchors.empty()) {
        size_t anchor = anchors[0];
        place_candidate(target, arbitrary_position_from_angle_local(
            cc.atoms[anchor].xyz, cc.atoms[center].xyz, dist,
            get_angle_rad(anchor, center, target)));
        if (atom_has_finite_xyz(cc.atoms[target]))
          return true;
      }

      cc.atoms[target].xyz = cc.atoms[center].xyz + Position(dist, 0, 0);
      if (atom_has_finite_xyz(cc.atoms[target]))
        return true;
    }
    return false;
  };

  int component = 0;
  int heavy_count = 0;
  for (const auto& atom : cc.atoms)
    if (is_heavy_atom(atom))
      ++heavy_count;
  while (count_finite(true) < heavy_count) {
    size_t seed = SIZE_MAX;
    for (size_t i = 0; i != n; ++i)
      if (is_heavy_atom(cc.atoms[i]) && !atom_has_finite_xyz(cc.atoms[i])) {
        seed = i;
        break;
      }
    if (seed == SIZE_MAX)
      break;

    bool seeded_plane = false;
    int best_plane_score = -1;
    size_t best_plane = SIZE_MAX;
    for (size_t fi = 0; fi != plane_fragments.size(); ++fi) {
      const std::vector<size_t>& fragment = plane_fragments[fi];
      bool all_unplaced = true;
      int score = 0;
      for (size_t idx : fragment) {
        if (atom_has_finite_xyz(cc.atoms[idx])) {
          all_unplaced = false;
          break;
        }
        if (plane_membership[idx] > 1)
          score += plane_membership[idx] - 1;
      }
      if (all_unplaced && plane_fragments.size() <= 3 && fragment.size() <= 6 &&
          score > best_plane_score) {
        best_plane_score = score;
        best_plane = fi;
      }
    }
    if (plane_fragments.size() <= 3) {
      if (best_plane != SIZE_MAX && best_plane_score > 0 &&
          seed_plane_fragment_coordinates(cc, atom_index, adjacency, bond_dist,
                                          plane_fragments[best_plane], component * 8.0)) {
        seeded_plane = true;
      } else {
        for (const std::vector<size_t>& fragment : plane_fragments) {
          if (std::find(fragment.begin(), fragment.end(), seed) != fragment.end() &&
              seed_plane_fragment_coordinates(cc, atom_index, adjacency, bond_dist,
                                              fragment, component * 8.0)) {
            seeded_plane = true;
            break;
          }
        }
      }
    }
    if (seeded_plane) {
      bool progress = true;
      while (progress) {
        progress = false;
        for (size_t i = 0; i != n; ++i)
          progress = place_bridge_heavy_atom(cc, atom_index, adjacency, bond_dist, i) || progress;
        for (size_t i = 0; i != n; ++i)
          progress = try_place(i) || progress;
      }
      ++component;
      continue;
    }

    bool seeded_cycle = false;
    for (const auto& cycle : small_cycles) {
      if (std::find(cycle.begin(), cycle.end(), seed) != cycle.end() &&
          seed_cycle_coordinates(cc, adjacency, cycle, component * 8.0)) {
        seeded_cycle = true;
        break;
      }
    }
    if (seeded_cycle) {
      bool progress = true;
      while (progress) {
        progress = false;
        for (size_t i = 0; i != n; ++i)
          progress = place_bridge_heavy_atom(cc, atom_index, adjacency, bond_dist, i) || progress;
        for (size_t i = 0; i != n; ++i)
          progress = try_place(i) || progress;
      }
      ++component;
      continue;
    }

    cc.atoms[seed].xyz = Position(component * 8.0, 0, 0);
    bool seeded_pair = false;
    for (size_t i = 0; i != adjacency[seed].size(); ++i) {
      size_t nb = adjacency[seed][i];
      if (atom_has_finite_xyz(cc.atoms[nb]))
        continue;
      cc.atoms[nb].xyz = cc.atoms[seed].xyz + Position(get_bond_dist(seed, nb), 0, 0);
      seeded_pair = true;
      for (size_t next : adjacency[nb]) {
        if (next == seed || atom_has_finite_xyz(cc.atoms[next]))
          continue;
        cc.atoms[next].xyz = arbitrary_position_from_angle_local(
            cc.atoms[seed].xyz, cc.atoms[nb].xyz, get_bond_dist(nb, next),
            get_angle_rad(seed, nb, next));
        break;
      }
      break;
    }
    if (!seeded_pair && adjacency[seed].empty())
      cc.atoms[seed].xyz = Position(component * 8.0, 0, 0);

    bool progress = true;
    while (progress) {
      progress = false;
      for (size_t i = 0; i != n; ++i)
        progress = place_bridge_heavy_atom(cc, atom_index, adjacency, bond_dist, i) || progress;
      for (size_t i = 0; i != n; ++i)
        progress = try_place(i) || progress;
    }
    ++component;
  }

  spread_terminal_children(cc, adjacency, atom_index);
  regularize_small_cycles(cc, adjacency, small_cycles);
  enforce_plane_restraints(cc, atom_index);
  enforce_chirality_restraints(cc, atom_index, adjacency, ring_atom);
  regrow_trivalent_junction_atoms(cc, atom_index, adjacency, bond_dist, ring_atom);
  regrow_bridge_heavy_atoms(cc, atom_index, adjacency, bond_dist, ring_atom);
  regrow_terminal_heavy_atoms(cc, atom_index, adjacency, bond_dist, ring_atom);
  int heavy_placed = count_finite(true);
  if (heavy_placed > 0)
    refine_chemcomp_xyz(cc);

  place_hydrogens_from_restraints(cc, atom_index, adjacency, bond_dist);

  int placed = count_finite(false);
  cc.has_coordinates = placed > 0;

  return placed;
}

namespace {

struct ChemCompRefineTarget {
  ChemComp& cc;
  std::vector<size_t> atom_indices;
  std::vector<int> param_offset;  // -1 if atom is not refined

  struct Point {
    enum Type { Bond, Angle };
    Type type;
    size_t idx1, idx2, idx3;  // atom indices (idx3 unused for bonds)
    double target_value;       // ideal distance or angle (radians)
    double weight;             // 1/esd

    double get_y() const { return target_value; }
    double get_weight() const { return weight; }
  };

  std::vector<Point> points;

  ChemCompRefineTarget(ChemComp& cc_) : cc(cc_), param_offset(cc_.atoms.size(), -1) {
    for (size_t i = 0; i < cc.atoms.size(); ++i) {
      if (atom_has_finite_xyz(cc.atoms[i]) && is_heavy_atom(cc.atoms[i])) {
        param_offset[i] = static_cast<int>(atom_indices.size() * 3);
        atom_indices.push_back(i);
      }
    }

    std::map<std::string, size_t> atom_index;
    for (size_t i = 0; i < cc.atoms.size(); ++i)
      atom_index[cc.atoms[i].id] = i;

    for (const auto& bond : cc.rt.bonds) {
      auto it1 = atom_index.find(bond.id1.atom);
      auto it2 = atom_index.find(bond.id2.atom);
      if (it1 == atom_index.end() || it2 == atom_index.end())
        continue;
      size_t i1 = it1->second, i2 = it2->second;
      if (param_offset[i1] < 0 || param_offset[i2] < 0)
        continue;
      double value = std::isfinite(bond.value) ? bond.value : bond.value_nucleus;
      if (!std::isfinite(value) || value <= 0)
        continue;
      double esd = std::isfinite(bond.esd) ? bond.esd : 0.02;
      if (esd <= 0)
        esd = 0.02;
      points.push_back({Point::Bond, i1, i2, 0, value, 1.0 / esd});
    }

    for (const auto& angle : cc.rt.angles) {
      auto it1 = atom_index.find(angle.id1.atom);
      auto it2 = atom_index.find(angle.id2.atom);
      auto it3 = atom_index.find(angle.id3.atom);
      if (it1 == atom_index.end() || it2 == atom_index.end() || it3 == atom_index.end())
        continue;
      size_t i1 = it1->second, i2 = it2->second, i3 = it3->second;
      if (param_offset[i1] < 0 || param_offset[i2] < 0 || param_offset[i3] < 0)
        continue;
      if (!std::isfinite(angle.value))
        continue;
      double target_rad = rad(angle.value);
      double esd = std::isfinite(angle.esd) && angle.esd > 0 ? angle.esd : 1.0;
      double weight = 1.0 / rad(esd);
      points.push_back({Point::Angle, i1, i2, i3, target_rad, weight});
    }
  }

  std::vector<double> get_parameters() const {
    std::vector<double> params(atom_indices.size() * 3);
    for (size_t i = 0; i < atom_indices.size(); ++i) {
      const Position& p = cc.atoms[atom_indices[i]].xyz;
      params[i * 3 + 0] = p.x;
      params[i * 3 + 1] = p.y;
      params[i * 3 + 2] = p.z;
    }
    return params;
  }

  void set_parameters(const std::vector<double>& p) {
    for (size_t i = 0; i < atom_indices.size(); ++i) {
      Position& pos = cc.atoms[atom_indices[i]].xyz;
      pos.x = p[i * 3 + 0];
      pos.y = p[i * 3 + 1];
      pos.z = p[i * 3 + 2];
    }
  }

  double compute_value(const Point& p) const {
    if (p.type == Point::Bond) {
      return cc.atoms[p.idx1].xyz.dist(cc.atoms[p.idx2].xyz);
    } else {
      return calculate_angle(cc.atoms[p.idx1].xyz, cc.atoms[p.idx2].xyz,
                             cc.atoms[p.idx3].xyz);
    }
  }

  double compute_value_and_derivatives(const Point& p, std::vector<double>& dy_da) const {
    std::fill(dy_da.begin(), dy_da.end(), 0.0);

    if (p.type == Point::Bond) {
      const Position& p1 = cc.atoms[p.idx1].xyz;
      const Position& p2 = cc.atoms[p.idx2].xyz;
      Vec3 delta = p2 - p1;
      double d = delta.length();
      if (d < 1e-15)
        return d;
      Vec3 grad = delta * (1.0 / d);  // dd/dr2
      int off1 = param_offset[p.idx1];
      int off2 = param_offset[p.idx2];
      // dd/dr1 = -grad
      dy_da[off1 + 0] = -grad.x;
      dy_da[off1 + 1] = -grad.y;
      dy_da[off1 + 2] = -grad.z;
      // dd/dr2 = +grad
      dy_da[off2 + 0] = grad.x;
      dy_da[off2 + 1] = grad.y;
      dy_da[off2 + 2] = grad.z;
      return d;
    } else {
      // Angle p.idx1 - p.idx2 - p.idx3
      const Position& r1 = cc.atoms[p.idx1].xyz;
      const Position& r2 = cc.atoms[p.idx2].xyz;
      const Position& r3 = cc.atoms[p.idx3].xyz;
      Vec3 u = r1 - r2;
      Vec3 v = r3 - r2;
      double lu = u.length();
      double lv = v.length();
      if (lu < 1e-15 || lv < 1e-15)
        return 0.0;
      double cos_theta = u.dot(v) / (lu * lv);
      cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
      double theta = std::acos(cos_theta);
      double sin_theta = std::sin(theta);
      if (sin_theta < 1e-10)
        sin_theta = 1e-10;

      // dtheta/dr1 = (cos_theta * u/lu - v/lv) / (lu * sin_theta)
      Vec3 dtheta_dr1 = (cos_theta / lu * u - v / lv) * (1.0 / (lu * sin_theta));
      // dtheta/dr3 = (cos_theta * v/lv - u/lu) / (lv * sin_theta)
      Vec3 dtheta_dr3 = (cos_theta / lv * v - u / lu) * (1.0 / (lv * sin_theta));
      // dtheta/dr2 = -(dtheta/dr1 + dtheta/dr3)
      Vec3 dtheta_dr2 = -(dtheta_dr1 + dtheta_dr3);

      int off1 = param_offset[p.idx1];
      int off2 = param_offset[p.idx2];
      int off3 = param_offset[p.idx3];
      dy_da[off1 + 0] = dtheta_dr1.x;
      dy_da[off1 + 1] = dtheta_dr1.y;
      dy_da[off1 + 2] = dtheta_dr1.z;
      dy_da[off2 + 0] = dtheta_dr2.x;
      dy_da[off2 + 1] = dtheta_dr2.y;
      dy_da[off2 + 2] = dtheta_dr2.z;
      dy_da[off3 + 0] = dtheta_dr3.x;
      dy_da[off3 + 1] = dtheta_dr3.y;
      dy_da[off3 + 2] = dtheta_dr3.z;
      return theta;
    }
  }
};

} // namespace

double refine_chemcomp_xyz(ChemComp& cc) {
  ChemCompRefineTarget target(cc);
  if (target.atom_indices.empty() || target.points.empty())
    return 0.0;
  LevMar levmar;
  levmar.eval_limit = 200;
  return levmar.fit(target);
}

} // namespace gemmi
