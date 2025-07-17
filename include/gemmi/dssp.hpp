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

// Per-residue information for DSSP calculation
struct GEMMI_DLL ResidueInfo {
  Topo::ResInfo* res_info = nullptr;
  std::array<Atom*, 5> backbone_atoms = {nullptr, nullptr, nullptr, nullptr, nullptr}; // CA, C, O, N, H
  std::array<bool, 5> has_backbone_atom = {false, false, false, false, false};
  std::array<Topo::ResInfo*, 2> donors = {nullptr, nullptr};
  std::array<Topo::ResInfo*, 2> acceptors = {nullptr, nullptr};
  std::array<double, 2> donor_energies = {0.0, 0.0};
  std::array<double, 2> acceptor_energies = {0.0, 0.0};
  ResidueInfo* prev_residue = nullptr;
  ResidueInfo* next_residue = nullptr;
  bool is_proline = false;

  enum BackboneAtom { CA = 0, C = 1, O = 2, N = 3, H = 4 };

  void set_backbone_atom(BackboneAtom type, Atom* atom) {
    backbone_atoms[type] = atom;
    has_backbone_atom[type] = true;
  }

  Atom* get_backbone_atom(BackboneAtom type) const {
    return backbone_atoms[type];
  }

  bool has_atom(BackboneAtom type) const {
    return has_backbone_atom[type];
  }
};

// DSSP options/parameters
struct GEMMI_DLL DsspOptions {
  HydrogenMode hydrogen_mode = HydrogenMode::Existing;
  HBondDefinition hbond_definition = HBondDefinition::Energy;
  bool use_neighbor_search = true;
  double cutoff = 0.9;  // nm
  bool pi_helix_preference = true;
  bool search_polyproline = true;
  bool shortened_pp_stretch = false;
  bool clear_defective_residues = false;
  double hbond_energy_cutoff = -0.5;  // kcal/mol
  double min_ca_distance = 9.0;  // Angstrom
  double bend_angle_min = 70.0;  // degrees
  double max_peptide_bond_distance = 2.5;  // Angstrom
};

// Main DSSP calculator class
struct GEMMI_DLL DsspCalculator {
  explicit DsspCalculator(const DsspOptions& opts = DsspOptions{}) : options(opts) {}

  // Calculate secondary structure for a chain
  std::string calculate_secondary_structure(NeighborSearch& ns, Topo::ChainInfo& cinfo);

  // Get detailed secondary structure information
  const std::vector<SecondaryStructureInfo>& get_detailed_info() const { return ss_info; }

  DsspOptions options;
  std::vector<ResidueInfo> residue_info;
  std::vector<SecondaryStructureInfo> ss_info;

  // Setup residue information from topology
  void setup_residue_info(Topo::ChainInfo& cinfo);

  // Calculate hydrogen bonds
  void calculate_hydrogen_bonds(NeighborSearch& ns);
  void calculate_hbond_energy(ResidueInfo* donor, ResidueInfo* acceptor);
  void calculate_hbond_geometry(ResidueInfo* donor, ResidueInfo* acceptor);

  // Pattern recognition functions
  void find_bridges_and_strands();
  void find_turns_and_helices();
  void find_bends_and_breaks();
  void find_polyproline_helices();

  // Utility functions
  bool has_hbond_between(size_t donor_idx, size_t acceptor_idx) const;
  bool no_chain_breaks_between(size_t res1_idx, size_t res2_idx) const;
  BridgeType calculate_bridge_type(size_t res1_idx, size_t res2_idx) const;

  // Generate final secondary structure string
  std::string generate_ss_string() const;
};

// Convenience function for simple use cases
GEMMI_DLL std::string calculate_dssp(NeighborSearch& ns, Topo::ChainInfo& cinfo,
                                     const DsspOptions& opts = DsspOptions{});

// Legacy function for backwards compatibility
GEMMI_DLL std::vector<HBond>
dssp_determine_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& cinfo);

} // namespace gemmi
#endif
