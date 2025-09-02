// Copyright Global Phasing Ltd.
//
// DSSP (Define Secondary Structure of Proteins) implementation.

#ifndef GEMMI_DSSP_HPP_
#define GEMMI_DSSP_HPP_

#include "topo.hpp"
#include "neighbor.hpp"
#include <vector>
#include <string>

namespace gemmi {

// Secondary structure types as defined in DSSP
enum class SecondaryStructure : char {
  Loop = '~',        // Loop/coil
  Break = '=',       // Break
  Bend = 'S',        // Bend
  Turn = 'T',        // Turn
  Helix_PP = 'P',    // Polyproline helix
  Helix_5 = 'I',     // Pi helix (5-turn)
  Helix_3 = 'G',     // 3-10 helix
  Strand = 'E',      // Extended strand
  Bridge = 'B',      // Beta bridge
  Helix_4 = 'H'      // Alpha helix (4-turn)
};

// Turn types
enum class TurnType {
  Turn_3 = 3,
  Turn_4 = 4,
  Turn_5 = 5,
  Turn_PP = 6
};

// Helix positions
enum class HelixPosition {
  None = 0,
  Start,
  Middle,
  End,
  StartAndEnd
};

// Bridge types
enum class BridgeType {
  None = 0,
  Parallel,
  AntiParallel
};

struct Bridge {
    size_t partner1;
    size_t partner2;
    BridgeType type;
};

// Hydrogen bond modes
enum class HydrogenMode {
  Existing = 0,  // Use existing hydrogen atoms from structure
  Calculate      // Calculate hydrogen positions (original DSSP method)
};

// Hydrogen bond definition
enum class HBondDefinition {
  Energy = 0,   // Energy-based (original DSSP)
  Geometry      // Geometry-based (distance + angle)
};

struct GEMMI_DLL HBond {
  Topo::ResInfo *donor = nullptr, *acceptor = nullptr;
  char alt1 = '\0';
  char alt2 = '\0';
  double energy = 0;
  bool is_valid = false;
};

// Per-residue secondary structure information
struct GEMMI_DLL SecondaryStructureInfo {
  SecondaryStructure ss_type = SecondaryStructure::Loop;
  std::vector<size_t> parallel_bridges;
  std::vector<size_t> antiparallel_bridges;
  std::vector<SecondaryStructureInfo*> break_partners;
  std::array<HelixPosition, 4> helix_positions = {HelixPosition::None, HelixPosition::None,
                                                  HelixPosition::None, HelixPosition::None};
  bool has_break = false;
  bool nturn_acceptor = false;

  void set_helix_position(TurnType turn, HelixPosition pos) {
    helix_positions[static_cast<int>(turn) - 3] = pos;
  }

  HelixPosition get_helix_position(TurnType turn) const {
    return helix_positions[static_cast<int>(turn) - 3];
  }

  void add_bridge(size_t partner_idx, BridgeType type) {
    if (type == BridgeType::Parallel) {
      parallel_bridges.push_back(partner_idx);
    } else if (type == BridgeType::AntiParallel) {
      antiparallel_bridges.push_back(partner_idx);
    }
  }

  bool has_bridges(BridgeType type) const {
    return type == BridgeType::Parallel ? !parallel_bridges.empty()
                                       : !antiparallel_bridges.empty();
  }
};

// DSSP options/parameters
struct GEMMI_DLL DsspOptions {
  HydrogenMode hydrogen_mode = HydrogenMode::Calculate;
  HBondDefinition hbond_definition = HBondDefinition::Energy;
  double cutoff = 0.9;  // nm
  bool pi_helix_preference = true;
  bool search_polyproline = true;
  bool shortened_pp_stretch = false;
  double hbond_energy_cutoff = -0.5;  // kcal/mol
  double min_ca_distance = 9.0;  // Angstrom
  double bend_angle_min = 70.0;  // degrees
  double max_peptide_bond_distance = 2.5;  // Angstrom
};

// Main DSSP calculator class
struct GEMMI_DLL DsspCalculator {
  explicit DsspCalculator(const DsspOptions& opts = DsspOptions{}) : options(opts) {}

  // Calculate secondary structure for a chain
  std::string calculate_secondary_structure(NeighborSearch& ns, Topo::ChainInfo& chain_info);

  // Get detailed secondary structure information
  const std::vector<SecondaryStructureInfo>& get_detailed_info() const { return ss_info; }

  DsspOptions options;
  std::vector<SecondaryStructureInfo> ss_info;
  std::vector<Bridge> bridges_;

  // Calculate hydrogen bonds
  void calculate_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& chain_info);
  void calculate_hbond_energy(Topo::ResInfo* donor, Topo::ResInfo* acceptor);
  void calculate_hbond_geometry(Topo::ResInfo* donor, Topo::ResInfo* acceptor);

  // Pattern recognition functions
  void find_bridges_and_strands(Topo::ChainInfo& chain_info);
  void find_turns_and_helices(Topo::ChainInfo& chain_info);
  void find_bends_and_breaks(Topo::ChainInfo& chain_info);
  void find_polyproline_helices(Topo::ChainInfo& chain_info);

  // Utility functions
  bool has_hbond_between(Topo::ResInfo* donor, Topo::ResInfo* acceptor) const;
  bool no_chain_breaks_between(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;
  BridgeType calculate_bridge_type(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;

  // Generate final secondary structure string
  std::string generate_ss_string() const;
};

// Convenience function for simple use cases
GEMMI_DLL std::string calculate_dssp(NeighborSearch& ns, Topo::ChainInfo& cinfo,
                                     const DsspOptions& opts = DsspOptions{});

} // namespace gemmi
#endif
