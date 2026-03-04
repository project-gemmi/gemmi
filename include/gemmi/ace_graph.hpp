// Copyright 2026 Global Phasing Ltd.
//
// Shared graph/ring helpers for AceDRG-style ChemComp processing.

#ifndef GEMMI_ACE_GRAPH_HPP_
#define GEMMI_ACE_GRAPH_HPP_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "gemmi/chemcomp.hpp"

namespace gemmi {

struct CodAtomInfo;

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

inline bool is_aromatic_or_deloc(BondType type) {
  return type == BondType::Aromatic || type == BondType::Deloc;
}

std::map<std::string, std::vector<std::string>>
make_neighbor_names(const ChemComp& cc);

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

bool compute_metal_neighbor_valence_charge(
    const ChemComp& cc, const AceBondAdjacency& adj, size_t idx,
    int& out_charge);

bool has_non_hydrogen_neighbor(
    const ChemComp& cc, const std::vector<int>& neighbor_indices);

void add_angles_from_bonds_if_missing(ChemComp& cc);

std::vector<unsigned> compute_rdkit_legacy_cip_ranks(
    const ChemComp& cc, const AceBondAdjacency& adj);

void sort_neighbors_by_rdkit_cip_rank(
    std::vector<size_t>& neighbors, const std::vector<unsigned>& cip_ranks);

bool is_oxygen_column(Element el);

std::pair<size_t, size_t> find_ring_sharing_pair(
    const AceBondAdjacency& adj, const std::vector<CodAtomInfo>& atom_info,
    size_t side1, size_t side2);

struct RingInfo {
  std::vector<int> atoms;
  std::string rep;
  std::string s_rep;
  bool is_aromatic = false;
  bool is_aromatic_permissive = false;
};

void detect_rings_acedrg(const std::vector<std::vector<int>>& neighbors,
                         std::vector<CodAtomInfo>& atoms,
                         std::vector<RingInfo>& rings);

void set_ring_aromaticity_from_bonds(const AceBondAdjacency& adj,
                                     const std::vector<CodAtomInfo>& atoms,
                                     std::vector<RingInfo>& rings,
                                     int verbose = 0);

void set_atoms_ring_rep_s(std::vector<CodAtomInfo>& atoms,
                          const std::vector<RingInfo>& rings);

void append_ring_annotation(std::string& s,
                            const std::map<std::string, std::string>& ring_rep_s);

enum class RingParity { Even, Odd, NoFlip };
enum class RingFlip { Even, Odd };

std::map<std::pair<size_t, size_t>, RingParity>
build_ring_bond_parity(const AceBondAdjacency& adj,
                       const std::vector<CodAtomInfo>& atom_info);

std::map<std::pair<size_t, size_t>, RingFlip>
build_ring_bond_flip(const ChemComp& cc,
                     const AceBondAdjacency& adj,
                     const std::vector<CodAtomInfo>& atom_info);

struct SugarRingInfo {
  std::map<int, std::set<size_t>> ring_sets;
  std::set<std::pair<size_t, size_t>> ring_bonds;
  std::map<int, std::vector<size_t>> ring_seq;
};

SugarRingInfo detect_sugar_rings(const ChemComp& cc, const AceBondAdjacency& adj,
                                 const std::vector<CodAtomInfo>& atom_info);

bool is_pyranose_ring_like_acedrg(const ChemComp& cc, const AceBondAdjacency& adj,
                                  const std::vector<CodAtomInfo>& atom_info,
                                  const std::vector<size_t>& ring_atoms,
                                  const std::set<size_t>& ring_set);

}  // namespace gemmi

#endif  // GEMMI_ACE_GRAPH_HPP_
