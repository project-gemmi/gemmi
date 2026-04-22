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

/// @brief Secondary structure type codes as defined in DSSP.
enum class SecondaryStructure : char {
  Loop = '~',        /// Coil/loop (unassigned secondary structure)
  Break = '=',       /// Chain break
  Bend = 'S',        /// Bend (local geometry criterion)
  Turn = 'T',        /// Hydrogen-bonded turn
  Helix_PP = 'P',    /// Polyproline II helix
  Helix_5 = 'I',     /// π-helix (5-turn)
  Helix_3 = 'G',     /// 3₁₀-helix (3-turn)
  Strand = 'E',      /// Extended β-strand
  Bridge = 'B',      /// Isolated β-bridge
  Helix_4 = 'H'      /// α-helix (4-turn)
};

/// @brief Hydrogen bond turn types, indexed by residue separation (i to i+n).
enum class TurnType {
  Turn_3 = 3,   /// 3-residue turn (i to i+3 hydrogen bond)
  Turn_4 = 4,   /// 4-residue turn (i to i+4 hydrogen bond)
  Turn_5 = 5,   /// 5-residue turn (i to i+5 hydrogen bond)
  Turn_PP = 6   /// Polyproline II type turn (i to i+6)
};

/// @brief Position of a residue within a helix.
enum class HelixPosition {
  None = 0,       /// Not part of a helix
  Start,          /// First residue of a helix
  Middle,         /// Interior residue of a helix
  End,            /// Last residue of a helix
  StartAndEnd     /// Single-residue helix (acts as both start and end)
};

/// @brief Type of β-bridge partnership between residues.
enum class BridgeType {
  None = 0,       /// No bridge
  Parallel,       /// Parallel β-bridge
  AntiParallel    /// Antiparallel β-bridge
};

/// @brief A simple β-bridge record between two residues.
struct Bridge {
    size_t partner1;     /// Chain index of the first bridge residue
    size_t partner2;     /// Chain index of the second bridge residue
    BridgeType type;     /// Type of bridge (Parallel or AntiParallel)
};

/// @brief Source of hydrogen atoms for hydrogen bond detection.
enum class HydrogenMode {
  Existing = 0,  /// Use hydrogen atoms already present in the structure
  Calculate      /// Compute idealized hydrogen positions (original DSSP method, for structures lacking explicit H)
};

/// @brief Criterion for hydrogen bond acceptance.
enum class HBondDefinition {
  Energy = 0,   /// Use Kabsch-Sander energy criterion (E < cutoff in kcal/mol)
  Geometry      /// Use geometric criteria (distance and angle thresholds)
};

/// @brief Record of one hydrogen bond detected by DSSP.
struct GEMMI_DLL HBond {
  Topo::ResInfo *donor = nullptr;       /// Pointer to the donor residue (NH group)
  Topo::ResInfo *acceptor = nullptr;    /// Pointer to the acceptor residue (C=O group)
  char alt1 = '\0';                      /// Alternate conformer of the donor atom
  char alt2 = '\0';                      /// Alternate conformer of the acceptor atom
  double energy = 0;                     /// Hydrogen bond energy in kcal/mol (for Energy mode); NaN if not applicable
  bool is_valid = false;                 /// True if this bond meets the cutoff criterion
};

/// @brief Per-residue secondary structure annotation and structural features.
struct GEMMI_DLL SecondaryStructureInfo {
  SecondaryStructure ss_type = SecondaryStructure::Loop;  /// Assigned secondary structure code
  std::vector<size_t> parallel_bridges;                   /// Chain indices of residues forming parallel β-bridges with this residue
  std::vector<size_t> antiparallel_bridges;               /// Chain indices of residues forming antiparallel β-bridges with this residue
  std::vector<SecondaryStructureInfo*> break_partners;    /// Pointers to SecondaryStructureInfo of chain-break partners
  std::array<HelixPosition, 4> helix_positions = {HelixPosition::None, HelixPosition::None,
                                                  HelixPosition::None, HelixPosition::None};  /// Helix position for each turn type (indexed by TurnType-3)
  bool has_break = false;                                 /// True if there is a chain break before this residue
  bool nturn_acceptor = false;                            /// True if this residue is the acceptor of an n-turn hydrogen bond

  /// @brief Set the helix position for a given turn type.
  /// @param turn TurnType value (Turn_3, Turn_4, Turn_5, or Turn_PP)
  /// @param pos HelixPosition value
  void set_helix_position(TurnType turn, HelixPosition pos) {
    helix_positions[static_cast<int>(turn) - 3] = pos;
  }

  /// @brief Retrieve the helix position for a given turn type.
  /// @param turn TurnType value
  /// @return HelixPosition at the specified turn
  HelixPosition get_helix_position(TurnType turn) const {
    return helix_positions[static_cast<int>(turn) - 3];
  }

  /// @brief Record a β-bridge partnership with another residue.
  /// @param partner_idx Chain index of the bridge partner
  /// @param type BridgeType (Parallel or AntiParallel)
  void add_bridge(size_t partner_idx, BridgeType type) {
    if (type == BridgeType::Parallel) {
      parallel_bridges.push_back(partner_idx);
    } else if (type == BridgeType::AntiParallel) {
      antiparallel_bridges.push_back(partner_idx);
    }
  }

  /// @brief Check if this residue has any bridges of the given type.
  /// @param type BridgeType to check
  /// @return True if this residue has at least one bridge of the given type
  bool has_bridges(BridgeType type) const {
    return type == BridgeType::Parallel ? !parallel_bridges.empty()
                                       : !antiparallel_bridges.empty();
  }
};

/// @brief Configuration parameters for the DSSP secondary structure algorithm.
struct GEMMI_DLL DsspOptions {
  HydrogenMode hydrogen_mode = HydrogenMode::Calculate;  /// How to source hydrogen atoms (Existing or Calculate)
  HBondDefinition hbond_definition = HBondDefinition::Energy;  /// Criterion for hydrogen bond acceptance (Energy or Geometry)
  double cutoff = 0.9;  /// Distance cutoff in nm for hydrogen bond search (default 0.9)
  bool pi_helix_preference = true;  /// If true, prefer π-helix over α-helix at ambiguous positions
  bool search_polyproline = true;  /// If true, detect polyproline II helices
  bool shortened_pp_stretch = false;  /// If true, allow polyproline stretches shorter than canonical
  double hbond_energy_cutoff = -0.5;  /// Energy threshold in kcal/mol for hydrogen bond acceptance (default -0.5)
  double min_ca_distance = 9.0;  /// Minimum Cα–Cα distance in Angstrom to consider residues potentially bonded (default 9.0)
  double bend_angle_min = 70.0;  /// Minimum bend angle in degrees to assign Bend secondary structure (default 70.0)
  double max_peptide_bond_distance = 2.5;  /// Maximum C–N distance in Angstrom considered a valid peptide bond (default 2.5)
};

/// @brief Main DSSP engine that performs secondary structure assignment for a protein chain.
/// Detects and annotates secondary structure elements (helices, strands, turns, etc.)
/// based on hydrogen bond patterns and geometry.
struct GEMMI_DLL DsspCalculator {
  /// @brief Initialize the calculator with options.
  /// @param opts DsspOptions configuration (default: standard DSSP parameters)
  explicit DsspCalculator(const DsspOptions& opts = DsspOptions{}) : options(opts) {}

  /// @brief Run the full DSSP secondary structure assignment for a chain.
  /// Requires a populated NeighborSearch for spatial lookups.
  /// @param ns NeighborSearch object with the protein atoms indexed
  /// @param chain_info Topo::ChainInfo containing residue topology
  /// @return Secondary structure string, one character per residue
  std::string calculate_secondary_structure(NeighborSearch& ns, Topo::ChainInfo& chain_info);

  /// @brief Access the per-residue secondary structure information.
  /// @return Const reference to the vector of SecondaryStructureInfo
  const std::vector<SecondaryStructureInfo>& get_detailed_info() const { return ss_info; }

  DsspOptions options;                              /// Configuration options for this calculation
  std::vector<SecondaryStructureInfo> ss_info;      /// Per-residue secondary structure and structural features (populated after calculate_secondary_structure)
  std::vector<Bridge> bridges_;                      /// List of all detected β-bridges

  /// @brief Detect and record all hydrogen bonds in the chain.
  /// @param ns NeighborSearch object with atoms indexed
  /// @param chain_info Topo::ChainInfo containing residue topology
  void calculate_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& chain_info);

  /// @brief Compute and store Kabsch-Sander hydrogen bond energy between donor and acceptor.
  /// @param donor Topo::ResInfo pointer for the NH (donor) residue
  /// @param acceptor Topo::ResInfo pointer for the C=O (acceptor) residue
  void calculate_hbond_energy(Topo::ResInfo* donor, Topo::ResInfo* acceptor);

  /// @brief Evaluate geometric hydrogen bond criterion (distance and angle) between donor and acceptor.
  /// @param donor Topo::ResInfo pointer for the NH (donor) residue
  /// @param acceptor Topo::ResInfo pointer for the C=O (acceptor) residue
  void calculate_hbond_geometry(Topo::ResInfo* donor, Topo::ResInfo* acceptor);

  /*
  // Pattern recognition functions
  void find_bridges_and_strands(Topo::ChainInfo& chain_info);
  void find_turns_and_helices(Topo::ChainInfo& chain_info);
  void find_bends_and_breaks(Topo::ChainInfo& chain_info);
  void find_polyproline_helices(Topo::ChainInfo& chain_info);
  */

  /// @brief Check if a valid hydrogen bond exists from donor to acceptor.
  /// @param donor Topo::ResInfo pointer for the donor
  /// @param acceptor Topo::ResInfo pointer for the acceptor
  /// @return True if a hydrogen bond exists between the residues
  bool has_hbond_between(Topo::ResInfo* donor, Topo::ResInfo* acceptor) const;

  /// @brief Check for chain breaks between two residue indices.
  /// @param chain_info Topo::ChainInfo containing residue topology
  /// @param res1_idx Chain index of the first residue
  /// @param res2_idx Chain index of the second residue
  /// @return True if no chain breaks exist between res1_idx and res2_idx
  bool no_chain_breaks_between(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;

  /// @brief Determine the type of β-bridge between two residues.
  /// @param chain_info Topo::ChainInfo containing residue topology
  /// @param res1_idx Chain index of the first residue
  /// @param res2_idx Chain index of the second residue
  /// @return BridgeType (Parallel, AntiParallel, or None)
  BridgeType calculate_bridge_type(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;

  // Generate final secondary structure string
  //std::string generate_ss_string() const;
};

/// @brief Convenience function to calculate secondary structure with default or custom options.
/// Creates a DsspCalculator with the given options and runs the full secondary structure assignment.
/// @param ns NeighborSearch object with atoms indexed for spatial lookups
/// @param cinfo Topo::ChainInfo containing the chain topology
/// @param opts DsspOptions configuration (default: standard DSSP parameters)
/// @return Secondary structure string, one character per residue
GEMMI_DLL std::string calculate_dssp(NeighborSearch& ns, Topo::ChainInfo& cinfo,
                                     const DsspOptions& opts = DsspOptions{});

} // namespace gemmi
#endif
