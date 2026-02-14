// Copyright 2026 Global Phasing Ltd.
//
// Shared graph/ring helpers for AceDRG-style ChemComp processing.

#ifndef GEMMI_ACE_GRAPH_HPP_
#define GEMMI_ACE_GRAPH_HPP_

#include <map>
#include <string>
#include <vector>
#include "gemmi/chemcomp.hpp"

namespace gemmi {

struct AceBondNeighbor {
  size_t idx;
  BondType type;
  bool aromatic;
};

using AceBondAdjacency = std::vector<std::vector<AceBondNeighbor>>;

struct AceGraphView {
  std::map<std::string, size_t> atom_index;
  AceBondAdjacency adjacency;
  std::vector<std::vector<int>> neighbors;
};

std::map<std::string, std::vector<std::string>>
make_neighbor_names(const ChemComp& cc);

bool atoms_in_same_ring_by_alt_path(
    const std::string& atom1, const std::string& atom2,
    const std::map<std::string, std::vector<std::string>>& neighbors);

AceBondAdjacency build_bond_adjacency(
    const ChemComp& cc, const std::map<std::string, size_t>& atom_index);

std::vector<std::vector<int>> build_neighbors(const AceBondAdjacency& adj);

AceGraphView make_ace_graph_view(const ChemComp& cc);

std::vector<size_t> neighbor_indices_except(
    const AceBondAdjacency& adj, size_t center_idx, size_t exclude_idx);

std::vector<size_t> non_hydrogen_neighbors(
    const ChemComp& cc, const AceBondAdjacency& adj,
    size_t center_idx, size_t exclude_idx = static_cast<size_t>(-1));

std::vector<size_t> hydrogen_neighbors(
    const ChemComp& cc, const AceBondAdjacency& adj,
    size_t center_idx, size_t exclude_idx = static_cast<size_t>(-1));

bool is_carbonyl_carbon(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx);

bool share_ring_ids(const std::vector<int>& rings_a,
                    const std::vector<int>& rings_b);

int shared_ring_size_from_ring_ids(const std::vector<int>& rings_a, int min_ring_a,
                                   const std::vector<int>& rings_b, int min_ring_b);

int expected_valence_for_nonmetal(Element el);

bool has_metal_and_non_metal_heavy_neighbor(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx);

float sum_non_metal_bond_order(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx);

bool has_non_hydrogen_neighbor(
    const ChemComp& cc, const std::vector<int>& neighbor_indices);

}  // namespace gemmi

#endif  // GEMMI_ACE_GRAPH_HPP_
