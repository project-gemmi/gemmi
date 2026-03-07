// Copyright 2026 Global Phasing Ltd.

#include "gemmi/chemcomp_xyz.hpp"
#include "gemmi/calculate.hpp"
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

bool atom_has_finite_xyz(const ChemComp::Atom& atom) {
  return std::isfinite(atom.xyz.x) && std::isfinite(atom.xyz.y) && std::isfinite(atom.xyz.z);
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
                                  const std::map<std::string, size_t>& atom_index) {
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
    reflect_across_plane(cc.atoms[a3].xyz, cc.atoms[ctr].xyz,
                         cc.atoms[a1].xyz, cc.atoms[a2].xyz);
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

bool seed_cycle_coordinates(ChemComp& cc, const std::vector<std::vector<size_t>>& adjacency,
                            const std::vector<size_t>& cycle, double x_shift) {
  if (cycle.size() < 5 || cycle.size() > 6)
    return false;
  for (size_t idx : cycle)
    if (atom_has_finite_xyz(cc.atoms[idx]))
      return false;

  std::unordered_set<size_t> ring_set(cycle.begin(), cycle.end());
  std::vector<size_t> ordered;
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
      if (auto bond = cc.rt.find_bond(cc.atoms[center].id, cc.atoms[term].id);
          bond != cc.rt.bonds.end()) {
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

  auto count_finite = [&]() {
    int count = 0;
    for (const auto& atom : cc.atoms)
      if (atom_has_finite_xyz(atom))
        ++count;
    return count;
  };

  for (auto& atom : cc.atoms)
    atom.xyz = Position(NAN, NAN, NAN);

  std::vector<std::vector<size_t>> small_cycles = detect_small_cycles(adjacency);

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
  while (count_finite() < static_cast<int>(n)) {
    size_t seed = SIZE_MAX;
    for (size_t i = 0; i != n; ++i)
      if (!atom_has_finite_xyz(cc.atoms[i])) {
        seed = i;
        break;
      }
    if (seed == SIZE_MAX)
      break;

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
        progress = try_place(i) || progress;
    }
    ++component;
  }

  spread_terminal_children(cc, adjacency, atom_index);
  enforce_plane_restraints(cc, atom_index);
  enforce_chirality_restraints(cc, atom_index);

  int placed = count_finite();
  cc.has_coordinates = placed > 0;
  return placed;
}

} // namespace gemmi
