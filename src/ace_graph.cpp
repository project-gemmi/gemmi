// Copyright 2026 Global Phasing Ltd.
//
// Shared graph/ring helpers for AceDRG-style ChemComp processing.

#include "gemmi/ace_graph.hpp"
#include "gemmi/acedrg_tables.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>

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
    bool aromatic = is_aromatic_or_deloc(bond.type);
    adj[idx1].push_back({idx2, bond.type, aromatic});
    adj[idx2].push_back({idx1, bond.type, aromatic});
  }
  // AceDRG populates connAtoms in two passes: first non-H neighbors (from
  // CIF bond table), then H neighbors (added in H-generation).  Stable-sort
  // to put non-H before H while keeping within-group bond-file order.
  for (auto& nbs : adj)
    std::stable_sort(nbs.begin(), nbs.end(),
                     [&](const AceBondNeighbor& a, const AceBondNeighbor& b) {
                       bool ah = cc.atoms[a.idx].is_hydrogen();
                       bool bh = cc.atoms[b.idx].is_hydrogen();
                       return !ah && bh;  // non-H before H
                     });
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
  if (el == El::Si) return 4;
  if (el == El::Ge) return 4;
  if (el == El::As) return 3;
  if (el == El::Se) return 2;
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

bool compute_metal_neighbor_valence_charge(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx,
    int& out_charge) {
  const Element el = cc.atoms[idx].el;
  if (el.is_metal() || el == El::H)
    return false;
  if (!has_metal_and_non_metal_heavy_neighbor(cc, adj, idx))
    return false;
  int expected_valence = expected_valence_for_nonmetal(el);
  if (expected_valence == 0)
    return false;
  float sum_bo = sum_non_metal_bond_order(cc, adj, idx);
  int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
  out_charge = -rem_v;
  return true;
}

bool has_non_hydrogen_neighbor(
    const ChemComp& cc, const std::vector<int>& neighbor_indices) {
  for (int nb : neighbor_indices)
    if (!cc.atoms[nb].is_hydrogen())
      return true;
  return false;
}

void add_angles_from_bonds_if_missing(ChemComp& cc) {
  if (!cc.rt.angles.empty())
    return;

  auto atom_index = cc.make_atom_index();
  std::vector<std::vector<size_t>> neighbors(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    neighbors[it1->second].push_back(it2->second);
    neighbors[it2->second].push_back(it1->second);
  }

  std::set<std::tuple<std::string, std::string, std::string>> seen;
  for (size_t center = 0; center < neighbors.size(); ++center) {
    if (cc.atoms[center].el.is_metal())
      continue;
    auto& nbs = neighbors[center];
    if (nbs.size() < 2)
      continue;
    for (size_t i = 0; i + 1 < nbs.size(); ++i) {
      for (size_t j = i + 1; j < nbs.size(); ++j) {
        if (cc.atoms[nbs[i]].el.is_metal() && cc.atoms[nbs[j]].el.is_metal())
          continue;
        const std::string& a1 = cc.atoms[nbs[i]].id;
        const std::string& a3 = cc.atoms[nbs[j]].id;
        std::string first = a1;
        std::string third = a3;
        if (third < first)
          std::swap(first, third);
        auto key = std::make_tuple(cc.atoms[center].id, first, third);
        if (!seen.insert(key).second)
          continue;
        cc.rt.angles.push_back({{1, a1}, {1, cc.atoms[center].id},
                                {1, a3}, NAN, NAN});
      }
    }
  }
}

int rdkit_twice_bond_type(BondType bt) {
  switch (bt) {
    case BondType::Single: return 2;
    case BondType::Double: return 4;
    case BondType::Triple: return 6;
    case BondType::Aromatic:
    case BondType::Deloc: return 3;
    default: return 0;
  }
}

int rdkit_cip_bond_count(const ChemComp& cc, const AceBondAdjacency& adj,
                         const AceBondNeighbor& nb) {
  if (nb.type == BondType::Double && cc.atoms[nb.idx].el == El::P) {
    size_t deg = adj[nb.idx].size();
    if (deg == 3 || deg == 4)
      return 1;
  }
  return rdkit_twice_bond_type(nb.type);
}

struct CipSortableRef {
  std::vector<int>* entry = nullptr;
  size_t atom_idx = 0;
  unsigned curr_rank = 0;
};

bool cip_ref_less(const CipSortableRef& a, const CipSortableRef& b) {
  return *a.entry < *b.entry;
}

void find_cip_tied_segments(std::vector<CipSortableRef>& sorted,
                            std::vector<std::pair<size_t, size_t>>& segments,
                            unsigned& num_ranks) {
  segments.clear();
  num_ranks = static_cast<unsigned>(sorted.size());
  if (sorted.empty())
    return;
  CipSortableRef* current = &sorted[0];
  unsigned running_rank = 0;
  current->curr_rank = running_rank;
  bool in_equal = false;
  for (size_t i = 1; i < sorted.size(); ++i) {
    CipSortableRef& next = sorted[i];
    if (*current->entry == *next.entry) {
      next.curr_rank = running_rank;
      --num_ranks;
      if (!in_equal) {
        in_equal = true;
        segments.push_back(std::make_pair(i - 1, size_t(0)));
      }
    } else {
      ++running_rank;
      next.curr_rank = running_rank;
      current = &next;
      if (in_equal) {
        segments.back().second = i;
        in_equal = false;
      }
    }
  }
  if (in_equal)
    segments.back().second = sorted.size() - 1;
}

std::vector<unsigned> compute_rdkit_legacy_cip_ranks(
    const ChemComp& cc, const AceBondAdjacency& adj) {
  size_t n = cc.atoms.size();
  std::vector<unsigned> ranks(n, 0);
  if (n == 0)
    return ranks;

  std::vector<std::vector<int>> cip_entries(n);
  std::vector<CipSortableRef> sorted;
  sorted.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    cip_entries[i].reserve(16);
    int atomic_num = cc.atoms[i].el.atomic_number() % 128;
    int mass = 512;
    unsigned long invariant = static_cast<unsigned long>(atomic_num);
    invariant = (invariant << 10) | static_cast<unsigned long>(mass);
    invariant = (invariant << 10);
    cip_entries[i].push_back(static_cast<int>(invariant));
    sorted.push_back({&cip_entries[i], i, 0});
  }

  std::sort(sorted.begin(), sorted.end(), cip_ref_less);
  std::vector<std::pair<size_t, size_t>> needs_sorting;
  unsigned num_ranks = 0;
  find_cip_tied_segments(sorted, needs_sorting, num_ranks);
  for (size_t i = 0; i < sorted.size(); ++i)
    ranks[sorted[i].atom_idx] = sorted[i].curr_rank;

  for (size_t i = 0; i < n; ++i) {
    cip_entries[i][0] = cc.atoms[i].el.atomic_number();
    cip_entries[i].push_back(static_cast<int>(ranks[i]));
  }

  const size_t cip_rank_index = 2;
  std::vector<std::vector<std::pair<int, size_t>>> weighted_neighbors(n);
  for (size_t i = 0; i < n; ++i) {
    weighted_neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i]) {
      if (cc.atoms[i].el.is_metal() || cc.atoms[nb.idx].el.is_metal())
        continue;
      weighted_neighbors[i].push_back(
          std::make_pair(rdkit_cip_bond_count(cc, adj, nb), nb.idx));
    }
  }

  unsigned max_iterations = static_cast<unsigned>(n / 2 + 1);
  unsigned iter = 0;
  int last_num_ranks = -1;
  while (!needs_sorting.empty() && iter < max_iterations &&
         (last_num_ranks < 0 ||
          static_cast<unsigned>(last_num_ranks) < num_ranks)) {
    for (size_t i = 0; i < n; ++i) {
      std::vector<std::pair<int, size_t>>& neighbors = weighted_neighbors[i];
      if (neighbors.size() > 1) {
        std::sort(neighbors.begin(), neighbors.end(),
                  [&](const std::pair<int, size_t>& a,
                      const std::pair<int, size_t>& b) {
                    return ranks[a.second] > ranks[b.second];
                  });
      }
      std::vector<int>& cip = cip_entries[i];
      for (const std::pair<int, size_t>& nb : neighbors)
        if (nb.first > 0)
          cip.insert(cip.end(), nb.first, static_cast<int>(ranks[nb.second] + 1));
    }

    last_num_ranks = static_cast<int>(num_ranks);
    for (const std::pair<size_t, size_t>& seg : needs_sorting)
      std::sort(sorted.begin() + seg.first, sorted.begin() + seg.second + 1, cip_ref_less);
    find_cip_tied_segments(sorted, needs_sorting, num_ranks);
    for (size_t i = 0; i < sorted.size(); ++i)
      ranks[sorted[i].atom_idx] = sorted[i].curr_rank;

    if (static_cast<unsigned>(last_num_ranks) != num_ranks) {
      for (size_t i = 0; i < n; ++i) {
        cip_entries[i].resize(cip_rank_index + 1);
        cip_entries[i][cip_rank_index] = static_cast<int>(ranks[i]);
      }
    }
    ++iter;
  }
  return ranks;
}

void sort_neighbors_by_rdkit_cip_rank(std::vector<size_t>& neighbors,
                                      const std::vector<unsigned>& cip_ranks) {
  std::stable_sort(neighbors.begin(), neighbors.end(),
                   [&](size_t a, size_t b) { return cip_ranks[a] > cip_ranks[b]; });
}

bool is_oxygen_column(Element el) {
  return el == El::O || el == El::S || el == El::Se ||
         el == El::Te || el == El::Po;
}

std::pair<size_t, size_t> find_ring_sharing_pair(
    const AceBondAdjacency& adj, const std::vector<CodAtomInfo>& atom_info,
    size_t side1, size_t side2) {
  for (const auto& nb1 : adj[side1]) {
    if (nb1.idx == side2) continue;
    for (const auto& nb2 : adj[side2]) {
      if (nb2.idx == side1 || nb2.idx == nb1.idx) continue;
      if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings))
        return {nb1.idx, nb2.idx};
    }
  }
  return {SIZE_MAX, SIZE_MAX};
}

std::map<std::pair<size_t, size_t>, RingParity>
build_ring_bond_parity(const AceBondAdjacency& adj,
                       const std::vector<CodAtomInfo>& atom_info) {
  std::map<std::pair<size_t, size_t>, RingParity> bond_ring_parity;
  int max_ring_id = 0;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      max_ring_id = std::max(max_ring_id, rid + 1);
  std::vector<std::vector<size_t>> rings(max_ring_id);
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      rings[rid].push_back(i);
  for (auto& ring : rings) {
    if (ring.empty())
      continue;
    std::set<size_t> ring_set(ring.begin(), ring.end());
    size_t start = *std::min_element(ring.begin(), ring.end());
    std::vector<size_t> traversal = {start};
    size_t prev = SIZE_MAX;
    size_t cur = start;
    while (traversal.size() < ring.size()) {
      bool found = false;
      for (const auto& nb : adj[cur]) {
        if (nb.idx != prev && ring_set.count(nb.idx) &&
            std::find(traversal.begin(), traversal.end(), nb.idx) == traversal.end()) {
          traversal.push_back(nb.idx);
          prev = cur;
          cur = nb.idx;
          found = true;
          break;
        }
      }
      if (!found)
        break;
    }
    for (size_t i = 0; i + 1 < traversal.size(); ++i) {
      auto key = std::minmax(traversal[i], traversal[i + 1]);
      bond_ring_parity[key] = (i % 2 == 0) ? RingParity::Even : RingParity::Odd;
    }
    if (traversal.size() == ring.size()) {
      auto ckey = std::minmax(traversal.back(), traversal.front());
      if (bond_ring_parity.find(ckey) == bond_ring_parity.end())
        bond_ring_parity[ckey] = RingParity::NoFlip;
    }
  }
  return bond_ring_parity;
}

std::map<std::pair<size_t, size_t>, RingFlip>
build_ring_bond_flip(const ChemComp& cc,
                     const AceBondAdjacency& adj,
                     const std::vector<CodAtomInfo>& atom_info) {
  struct RingWalk {
    std::string rep;
    std::vector<size_t> seq;
  };
  std::map<int, std::vector<size_t>> ring_atoms;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      ring_atoms[rid].push_back(i);

  std::vector<RingWalk> walks;
  walks.reserve(ring_atoms.size());
  for (const auto& kv : ring_atoms) {
    const std::vector<size_t>& r_atoms = kv.second;
    if (r_atoms.size() < 3)
      continue;
    std::set<size_t> rset(r_atoms.begin(), r_atoms.end());
    std::vector<size_t> linked;
    size_t start = *std::min_element(r_atoms.begin(), r_atoms.end());
    linked.push_back(start);
    size_t cur = start;
    int guard = 1;
    while (linked.size() < r_atoms.size() && guard < (int)r_atoms.size()) {
      bool advanced = false;
      for (const auto& nb : adj[cur]) {
        if (rset.count(nb.idx) == 0)
          continue;
        if (std::find(linked.begin(), linked.end(), nb.idx) != linked.end())
          continue;
        linked.push_back(nb.idx);
        cur = nb.idx;
        advanced = true;
        break;
      }
      if (!advanced)
        break;
      ++guard;
    }
    if (linked.size() != r_atoms.size())
      continue;
    std::vector<std::string> names;
    names.reserve(linked.size());
    for (size_t idx : linked)
      names.push_back(cc.atoms[idx].id);
    std::sort(names.begin(), names.end());
    std::string rep;
    for (const auto& n : names)
      rep += n;
    walks.push_back({std::move(rep), std::move(linked)});
  }
  std::sort(walks.begin(), walks.end(),
            [](const RingWalk& a, const RingWalk& b) { return a.rep < b.rep; });

  std::map<std::pair<size_t, size_t>, RingFlip> flips;
  for (const RingWalk& rw : walks) {
    for (size_t i = 0; i + 1 < rw.seq.size(); ++i) {
      auto key = std::minmax(rw.seq[i], rw.seq[i + 1]);
      if (flips.find(key) != flips.end())
        continue;
      flips[key] = (i % 2 == 0) ? RingFlip::Even : RingFlip::Odd;
    }
  }
  return flips;
}

SugarRingInfo detect_sugar_rings(const ChemComp& cc, const AceBondAdjacency& adj,
                                 const std::vector<CodAtomInfo>& atom_info) {
  SugarRingInfo out;

  auto co_token = [&](size_t idx) -> std::string {
    int nC = 0, nO = 0;
    for (const auto& nb : adj[idx]) {
      if (cc.atoms[nb.idx].el == El::C)
        ++nC;
      else if (cc.atoms[nb.idx].el == El::O)
        ++nO;
    }
    if (nC == 1 && nO == 2) return "CCO2";
    if (nC == 2 && nO == 1) return "CC2O";
    if (nC == 1 && nO == 1) return "CCO";
    if (nC == 2) return "CC2";
    return "";
  };
  auto ring_shape = [&](const std::set<size_t>& rset, size_t o_idx, size_t start,
                        std::vector<size_t>* out_order = nullptr) -> std::string {
    std::vector<size_t> order;
    order.reserve(rset.size());
    order.push_back(o_idx);
    order.push_back(start);
    size_t prev = o_idx;
    size_t curr = start;
    while (order.size() < rset.size()) {
      size_t next = SIZE_MAX;
      for (const auto& nb : adj[curr]) {
        if (nb.idx != prev && rset.count(nb.idx)) {
          next = nb.idx;
          break;
        }
      }
      if (next == SIZE_MAX)
        return "";
      order.push_back(next);
      prev = curr;
      curr = next;
    }
    bool closes = false;
    for (const auto& nb : adj[order.back()]) {
      if (nb.idx == o_idx) {
        closes = true;
        break;
      }
    }
    if (!closes)
      return "";
    if (out_order)
      *out_order = order;
    int o_c = 0;
    for (const auto& nb : adj[o_idx])
      if (cc.atoms[nb.idx].el == El::C)
        ++o_c;
    if (o_c != 2)
      return "";
    std::string s = "(OC2)";
    for (size_t i = 1; i < order.size(); ++i) {
      std::string t = co_token(order[i]);
      if (t.empty())
        return "";
      s += "(" + t + ")";
    }
    return s;
  };
  static const char* k_ace_sugar_shapes[] = {
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CC2O)(CC2O)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)(CCO)",
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CCO2)(CC2)(CC2O)(CC2O)(CC2O)"
  };
  static const std::set<std::string> allowed_shapes(
      std::begin(k_ace_sugar_shapes), std::end(k_ace_sugar_shapes));
  std::map<int, std::vector<size_t>> ring_atoms;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      ring_atoms[rid].push_back(i);
  for (auto& kv : ring_atoms) {
    auto& atoms = kv.second;
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    int n = (int)atoms.size();
    if (n != 5 && n != 6)
      continue;
    std::set<size_t> rset(atoms.begin(), atoms.end());
    size_t oxy_idx = SIZE_MAX;
    int oxy = 0, carb = 0, sp3 = 0;
    bool ok = true;
    for (size_t idx : rset) {
      if (atom_info[idx].bonding_idx == 3)
        ++sp3;
      if (atom_info[idx].bonding_idx != 3) {
        ok = false;
        break;
      }
      if (cc.atoms[idx].el == El::O) {
        ++oxy;
        oxy_idx = idx;
      } else if (cc.atoms[idx].el == El::C) {
        ++carb;
      } else {
        ok = false;
        break;
      }
      int ring_deg = 0;
      for (const auto& nb : adj[idx])
        if (rset.count(nb.idx))
          ++ring_deg;
      if (ring_deg != 2) {
        ok = false;
        break;
      }
    }
    if (!ok || sp3 != n || oxy != 1 || carb != n - 1 || oxy_idx == SIZE_MAX)
      continue;
    std::vector<size_t> o_ring_nbs;
    for (const auto& nb : adj[oxy_idx])
      if (rset.count(nb.idx))
        o_ring_nbs.push_back(nb.idx);
    if (o_ring_nbs.size() != 2)
      continue;
    std::vector<size_t> seq1, seq2;
    std::string shape1 = ring_shape(rset, oxy_idx, o_ring_nbs[0], &seq1);
    std::string shape2 = ring_shape(rset, oxy_idx, o_ring_nbs[1], &seq2);
    bool ok1 = allowed_shapes.count(shape1) != 0;
    bool ok2 = allowed_shapes.count(shape2) != 0;
    if (!ok1 && !ok2)
      continue;
    if (ok1 && ok2) {
      std::string t1 = co_token(o_ring_nbs[0]);
      std::string t2 = co_token(o_ring_nbs[1]);
      out.ring_seq[kv.first] = (t1 == "CC2O" && t2 == "CC2O") ? seq2 : seq1;
    } else {
      out.ring_seq[kv.first] = ok1 ? seq1 : seq2;
    }
    out.ring_sets.emplace(kv.first, std::move(rset));
  }
  for (const auto& kv : out.ring_sets) {
    const auto& rset = kv.second;
    for (size_t idx : rset)
      for (const auto& nb : adj[idx])
        if (rset.count(nb.idx) && idx < nb.idx)
          out.ring_bonds.insert(std::make_pair(idx, nb.idx));
  }
  return out;
}

bool is_pyranose_ring_like_acedrg(const ChemComp& cc, const AceBondAdjacency& adj,
                                  const std::vector<CodAtomInfo>& atom_info,
                                  const std::vector<size_t>& ring_atoms,
                                  const std::set<size_t>& ring_set) {
  std::map<std::string, int> fmt;
  fmt["O"] = 0;
  fmt["C"] = 0;
  fmt["connectC"] = 0;
  fmt["connectO"] = 0;

  for (size_t idx : ring_atoms) {
    if (cc.atoms[idx].el == El::C && atom_info[idx].bonding_idx == 3) {
      ++fmt["C"];
      int iCO = 0;
      for (const auto& nb : adj[idx]) {
        if (ring_set.count(nb.idx) != 0)
          continue;
        if (cc.atoms[nb.idx].el == El::O)
          ++iCO;
        else if (cc.atoms[nb.idx].el == El::C)
          ++fmt["connectC"];
      }
      if (iCO > 1)
        break;
      if (iCO == 1)
        ++fmt["connectO"];
    } else if (cc.atoms[idx].el == El::O) {
      if (adj[idx].size() == 2) {
        int iOC = 0;
        for (const auto& nb : adj[idx])
          if (cc.atoms[nb.idx].el == El::C && ring_set.count(nb.idx) != 0)
            ++iOC;
        if (iOC == 2)
          ++fmt["O"];
        else
          break;
      }
    } else {
      break;
    }

    if (fmt["O"] == 1 && fmt["C"] == 5 && fmt["connectO"] > 1)
      return true;
  }
  return false;
}

}  // namespace gemmi
