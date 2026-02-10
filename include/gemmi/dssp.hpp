//! @file
//! @brief DSSP (Define Secondary Structure of Proteins) implementation.

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

//! @brief Secondary structure types as defined in DSSP.
enum class SecondaryStructure : char {
  Loop = '~',        //!< Loop/coil
  Break = '=',       //!< Break
  Bend = 'S',        //!< Bend
  Turn = 'T',        //!< Turn
  Helix_PP = 'P',    //!< Polyproline helix
  Helix_5 = 'I',     //!< Pi helix (5-turn)
  Helix_3 = 'G',     //!< 3-10 helix
  Strand = 'E',      //!< Extended strand
  Bridge = 'B',      //!< Beta bridge
  Helix_4 = 'H'      //!< Alpha helix (4-turn)
};

//! @brief Turn types.
enum class TurnType {
  Turn_3 = 3,   //!< 3-turn
  Turn_4 = 4,   //!< 4-turn (alpha helix)
  Turn_5 = 5,   //!< 5-turn (pi helix)
  Turn_PP = 6   //!< Polyproline turn
};

//! @brief Helix positions.
enum class HelixPosition {
  None = 0,      //!< Not in helix
  Start,         //!< Helix start
  Middle,        //!< Helix middle
  End,           //!< Helix end
  StartAndEnd    //!< Single-turn helix
};

//! @brief Bridge types.
enum class BridgeType {
  None = 0,       //!< No bridge
  Parallel,       //!< Parallel beta bridge
  AntiParallel    //!< Antiparallel beta bridge
};

//! @brief Beta bridge between residues.
struct Bridge {
    size_t partner1;    //!< First partner residue index
    size_t partner2;    //!< Second partner residue index
    BridgeType type;    //!< Bridge type
};

//! @brief Hydrogen bond modes.
enum class HydrogenMode {
  Existing = 0,  //!< Use existing hydrogen atoms from structure
  Calculate      //!< Calculate hydrogen positions (original DSSP method)
};

//! @brief Hydrogen bond definition.
enum class HBondDefinition {
  Energy = 0,   //!< Energy-based (original DSSP)
  Geometry      //!< Geometry-based (distance + angle)
};

//! @brief Hydrogen bond between residues.
struct GEMMI_DLL HBond {
  Topo::ResInfo *donor = nullptr, *acceptor = nullptr;  //!< Donor and acceptor residues
  char alt1 = '\0';     //!< Altloc of donor
  char alt2 = '\0';     //!< Altloc of acceptor
  double energy = 0;    //!< H-bond energy (kcal/mol)
  bool is_valid = false;  //!< Whether H-bond is valid
};

//! @brief Per-residue secondary structure information.
struct GEMMI_DLL SecondaryStructureInfo {
  SecondaryStructure ss_type = SecondaryStructure::Loop;  //!< Secondary structure type
  std::vector<size_t> parallel_bridges;              //!< Parallel bridge partners
  std::vector<size_t> antiparallel_bridges;          //!< Antiparallel bridge partners
  std::vector<SecondaryStructureInfo*> break_partners;  //!< Chain break partners
  std::array<HelixPosition, 4> helix_positions = {HelixPosition::None, HelixPosition::None,
                                                  HelixPosition::None, HelixPosition::None};  //!< Helix positions for 3,4,5,PP turns
  bool has_break = false;      //!< Chain break after this residue
  bool nturn_acceptor = false;  //!< N-turn acceptor

  //! @brief Set helix position for turn type.
  //! @param turn Turn type
  //! @param pos Helix position
  void set_helix_position(TurnType turn, HelixPosition pos) {
    helix_positions[static_cast<int>(turn) - 3] = pos;
  }

  //! @brief Get helix position for turn type.
  //! @param turn Turn type
  //! @return Helix position
  HelixPosition get_helix_position(TurnType turn) const {
    return helix_positions[static_cast<int>(turn) - 3];
  }

  //! @brief Add bridge partner.
  //! @param partner_idx Partner residue index
  //! @param type Bridge type
  void add_bridge(size_t partner_idx, BridgeType type) {
    if (type == BridgeType::Parallel) {
      parallel_bridges.push_back(partner_idx);
    } else if (type == BridgeType::AntiParallel) {
      antiparallel_bridges.push_back(partner_idx);
    }
  }

  //! @brief Check if has bridges of type.
  //! @param type Bridge type
  //! @return True if has bridges
  bool has_bridges(BridgeType type) const {
    return type == BridgeType::Parallel ? !parallel_bridges.empty()
                                       : !antiparallel_bridges.empty();
  }
};

//! @brief DSSP options/parameters.
struct GEMMI_DLL DsspOptions {
  HydrogenMode hydrogen_mode = HydrogenMode::Calculate;  //!< Hydrogen calculation mode
  HBondDefinition hbond_definition = HBondDefinition::Energy;  //!< H-bond definition
  double cutoff = 0.9;                  //!< Distance cutoff (nm)
  bool pi_helix_preference = true;      //!< Prefer pi helix over alpha helix
  bool search_polyproline = true;       //!< Search for polyproline helices
  bool shortened_pp_stretch = false;    //!< Use shortened PP stretch criteria
  double hbond_energy_cutoff = -0.5;    //!< H-bond energy cutoff (kcal/mol)
  double min_ca_distance = 9.0;         //!< Minimum CA distance (Angstrom)
  double bend_angle_min = 70.0;         //!< Bend angle minimum (degrees)
  double max_peptide_bond_distance = 2.5;  //!< Maximum peptide bond distance (Angstrom)
};

//! @brief Main DSSP calculator class.
struct GEMMI_DLL DsspCalculator {
  //! @brief Construct DSSP calculator.
  //! @param opts DSSP options
  explicit DsspCalculator(const DsspOptions& opts = DsspOptions{}) : options(opts) {}

  //! @brief Calculate secondary structure for a chain.
  //! @param ns Neighbor search structure
  //! @param chain_info Chain information
  //! @return Secondary structure string
  std::string calculate_secondary_structure(NeighborSearch& ns, Topo::ChainInfo& chain_info);

  //! @brief Get detailed secondary structure information.
  //! @return Vector of per-residue SS info
  const std::vector<SecondaryStructureInfo>& get_detailed_info() const { return ss_info; }

  DsspOptions options;                          //!< DSSP options
  std::vector<SecondaryStructureInfo> ss_info;  //!< Per-residue SS info
  std::vector<Bridge> bridges_;                 //!< Beta bridges

  //! @brief Calculate hydrogen bonds.
  //! @param ns Neighbor search structure
  //! @param chain_info Chain information
  void calculate_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& chain_info);

  //! @brief Calculate H-bond energy.
  //! @param donor Donor residue
  //! @param acceptor Acceptor residue
  void calculate_hbond_energy(Topo::ResInfo* donor, Topo::ResInfo* acceptor);

  //! @brief Calculate H-bond geometry.
  //! @param donor Donor residue
  //! @param acceptor Acceptor residue
  void calculate_hbond_geometry(Topo::ResInfo* donor, Topo::ResInfo* acceptor);

  /*
  // Pattern recognition functions
  void find_bridges_and_strands(Topo::ChainInfo& chain_info);
  void find_turns_and_helices(Topo::ChainInfo& chain_info);
  void find_bends_and_breaks(Topo::ChainInfo& chain_info);
  void find_polyproline_helices(Topo::ChainInfo& chain_info);
  */

  //! @brief Check if H-bond exists between residues.
  //! @param donor Donor residue
  //! @param acceptor Acceptor residue
  //! @return True if H-bond exists
  bool has_hbond_between(Topo::ResInfo* donor, Topo::ResInfo* acceptor) const;

  //! @brief Check if no chain breaks between residues.
  //! @param chain_info Chain information
  //! @param res1_idx First residue index
  //! @param res2_idx Second residue index
  //! @return True if no chain breaks
  bool no_chain_breaks_between(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;

  //! @brief Calculate bridge type between residues.
  //! @param chain_info Chain information
  //! @param res1_idx First residue index
  //! @param res2_idx Second residue index
  //! @return Bridge type
  BridgeType calculate_bridge_type(Topo::ChainInfo& chain_info, size_t res1_idx, size_t res2_idx) const;

  // Generate final secondary structure string
  //std::string generate_ss_string() const;
};

//! @brief Convenience function for simple use cases.
//! @param ns Neighbor search structure
//! @param cinfo Chain information
//! @param opts DSSP options
//! @return Secondary structure string
GEMMI_DLL std::string calculate_dssp(NeighborSearch& ns, Topo::ChainInfo& cinfo,
                                     const DsspOptions& opts = DsspOptions{});

} // namespace gemmi
#endif
