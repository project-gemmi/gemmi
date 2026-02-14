// Copyright 2026 Global Phasing Ltd.
//
// Shared graph/ring helpers for AceDRG-style ChemComp processing.

#include "gemmi/ace_graph.hpp"

#include <algorithm>
#include <set>
#include "gemmi/calculate.hpp"

namespace gemmi {

std::map<std::string, std::vector<std::string>>
make_neighbor_names(const ChemComp& cc) {
  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }
  return neighbors;
}

bool atoms_in_same_ring_by_alt_path(
    const std::string& atom1, const std::string& atom2,
    const std::map<std::string, std::vector<std::string>>& neighbors) {
  auto from_it = neighbors.find(atom1);
  if (from_it == neighbors.end())
    return false;

  // BFS from atom1 to atom2, excluding the direct bond between them.
  std::set<std::string> visited;
  std::vector<std::string> queue;
  visited.insert(atom1);
  for (const std::string& nb : from_it->second) {
    if (nb != atom2) {
      queue.push_back(nb);
      visited.insert(nb);
    }
  }

  while (!queue.empty()) {
    std::string current = queue.back();
    queue.pop_back();
    if (current == atom2)
      return true;
    auto it = neighbors.find(current);
    if (it == neighbors.end())
      continue;
    for (const std::string& nb : it->second) {
      if (visited.insert(nb).second)
        queue.push_back(nb);
    }
  }
  return false;
}

AceBondAdjacency build_bond_adjacency(
    const ChemComp& cc, const std::map<std::string, size_t>& atom_index) {
  AceBondAdjacency adj(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;
    bool aromatic = (bond.type == BondType::Aromatic ||
                     bond.type == BondType::Deloc || bond.aromatic);
    adj[idx1].push_back({idx2, bond.type, aromatic});
    adj[idx2].push_back({idx1, bond.type, aromatic});
  }
  return adj;
}

std::vector<std::vector<int>> build_neighbors(const AceBondAdjacency& adj) {
  std::vector<std::vector<int>> neighbors(adj.size());
  for (size_t i = 0; i < adj.size(); ++i) {
    neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i])
      neighbors[i].push_back(static_cast<int>(nb.idx));
  }
  return neighbors;
}

AceGraphView make_ace_graph_view(const ChemComp& cc) {
  AceGraphView out;
  out.atom_index = cc.make_atom_index();
  out.adjacency = build_bond_adjacency(cc, out.atom_index);
  out.neighbors = build_neighbors(out.adjacency);
  return out;
}

std::vector<size_t> neighbor_indices_except(
    const AceBondAdjacency& adj, size_t center_idx, size_t exclude_idx) {
  std::vector<size_t> out;
  out.reserve(adj[center_idx].size());
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx)
      out.push_back(nb.idx);
  return out;
}

std::vector<size_t> non_hydrogen_neighbors(
    const ChemComp& cc, const AceBondAdjacency& adj,
    size_t center_idx, size_t exclude_idx) {
  std::vector<size_t> out;
  out.reserve(adj[center_idx].size());
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx && !cc.atoms[nb.idx].is_hydrogen())
      out.push_back(nb.idx);
  return out;
}

std::vector<size_t> hydrogen_neighbors(
    const ChemComp& cc, const AceBondAdjacency& adj,
    size_t center_idx, size_t exclude_idx) {
  std::vector<size_t> out;
  out.reserve(adj[center_idx].size());
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx && cc.atoms[nb.idx].is_hydrogen())
      out.push_back(nb.idx);
  return out;
}

bool is_carbonyl_carbon(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx) {
  if (cc.atoms[idx].el != El::C)
    return false;
  for (const auto& nb : adj[idx]) {
    if (cc.atoms[nb.idx].el == El::O &&
        (nb.type == BondType::Double || nb.type == BondType::Deloc))
      return true;
  }
  return false;
}

bool share_ring_ids(const std::vector<int>& rings_a,
                    const std::vector<int>& rings_b) {
  for (int ring_a : rings_a)
    for (int ring_b : rings_b)
      if (ring_a == ring_b)
        return true;
  return false;
}

int shared_ring_size_from_ring_ids(const std::vector<int>& rings_a, int min_ring_a,
                                   const std::vector<int>& rings_b, int min_ring_b) {
  if (!share_ring_ids(rings_a, rings_b))
    return 0;
  int size = 999;
  if (min_ring_a > 0)
    size = std::min(size, min_ring_a);
  if (min_ring_b > 0)
    size = std::min(size, min_ring_b);
  return size == 999 ? 1 : size;
}

int expected_valence_for_nonmetal(Element el) {
  if (el == El::O) return 2;
  if (el == El::N) return 3;
  if (el == El::S) return 2;
  if (el == El::C) return 4;
  if (el == El::P) return 3;
  return 0;
}

bool has_metal_and_non_metal_heavy_neighbor(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx) {
  bool has_metal = false;
  bool has_non_metal_heavy = false;
  for (const auto& nb : adj[idx]) {
    if (cc.atoms[nb.idx].el.is_metal())
      has_metal = true;
    else if (!cc.atoms[nb.idx].is_hydrogen())
      has_non_metal_heavy = true;
  }
  return has_metal && has_non_metal_heavy;
}

float sum_non_metal_bond_order(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx) {
  float sum_bo = 0.0f;
  for (const auto& nb : adj[idx]) {
    if (!cc.atoms[nb.idx].el.is_metal())
      sum_bo += order_of_bond_type(nb.type);
  }
  return sum_bo;
}

bool has_non_hydrogen_neighbor(
    const ChemComp& cc, const std::vector<int>& neighbor_indices) {
  for (int nb : neighbor_indices)
    if (!cc.atoms[nb].is_hydrogen())
      return true;
  return false;
}

}  // namespace gemmi
