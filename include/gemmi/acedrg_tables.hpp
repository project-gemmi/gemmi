// Copyright 2025 Global Phasing Ltd.
//
// AcedrgTables - COD/CSD-based atom classification and restraint value lookup
// Port of AceDRG codClassify system to gemmi.

#ifndef GEMMI_ACEDRG_TABLES_HPP_
#define GEMMI_ACEDRG_TABLES_HPP_

#include <string>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <algorithm>
#include "chemcomp.hpp"
#include "elem.hpp"
#include "fail.hpp"
#include "util.hpp"

namespace gemmi {

// Hybridization states used in atom classification
enum class Hybridization {
  SP1,    // sp hybridization (linear)
  SP2,    // sp2 hybridization (trigonal planar)
  SP3,    // sp3 hybridization (tetrahedral)
  SPD5,   // d-orbital involvement (5-coordinate)
  SPD6,   // d-orbital involvement (6-coordinate)
  SPD7,   // d-orbital involvement (7-coordinate)
  SPD8,   // d-orbital involvement (8-coordinate)
  SP_NON  // non-standard/unknown hybridization
};

GEMMI_DLL const char* hybridization_to_string(Hybridization h);
GEMMI_DLL Hybridization hybridization_from_string(const std::string& s);

// Metal coordination geometry types
enum class CoordGeometry {
  LINEAR,           // CN=2: 180°
  TRIGONAL_PLANAR,  // CN=3: 120°
  T_SHAPED,         // CN=3: 90°, 180°
  TETRAHEDRAL,      // CN=4: 109.47°
  SQUARE_PLANAR,    // CN=4: 90°, 180°
  TRIGONAL_BIPYRAMIDAL, // CN=5: 90°, 120°, 180°
  SQUARE_PYRAMIDAL, // CN=5: 90°, 180°
  OCTAHEDRAL,       // CN=6: 90°
  TRIGONAL_PRISM,   // CN=6: alternative
  PENTAGONAL_BIPYRAMIDAL, // CN=7
  CAPPED_OCTAHEDRAL,      // CN=7
  SQUARE_ANTIPRISM,       // CN=8
  UNKNOWN
};

// Classification information for a single atom
struct CodAtomInfo {
  int index;                // Index in ChemComp.atoms
  std::string id;           // Atom id (name)
  int hashing_value;        // 0-1000+ hash code
  Element el;               // Element
  Hybridization hybrid;     // Hybridization state
  std::string cod_class;    // Full COD class (e.g., "C[6a](C[6a]C[6a])(C[6a])(H)")
  std::string cod_class_no_charge;  // COD class computed without formal charges (for COD table lookup)
  std::string cod_main;     // COD main type (codAtmMain)
  std::string cod_root;     // COD root type (codAtmRoot)
  std::string nb_symb;      // codNBSymb
  std::string nb2_symb;     // codNB2Symb
  std::string nb3_symb;     // codNB3Symb
  std::string nb1nb2_sp;    // codNB1NB2_SP
  std::vector<int> conn_atoms_no_metal; // Non-metal neighbors (index list)
  int connectivity;         // Number of bonded atoms
  int metal_connectivity;   // Number of metal neighbors
  int min_ring_size;        // Minimum ring size (0 = not in ring)
  bool is_aromatic;         // In aromatic ring
  bool is_metal;            // Is a metal atom
  int excess_electrons;     // Formal charge/lone pair info
  float charge;             // Formal/partial charge
  float par_charge;         // Partial charge (AceDRG parCharge)
  int bonding_idx;          // AceDRG bonding index (1=sp1,2=sp2,3=sp3,...)

  // Ring bookkeeping (AceDRG-style).
  std::map<std::string, int> ring_rep;
  std::map<std::string, std::string> ring_rep_s;
  std::vector<int> in_rings;

  CodAtomInfo()
    : index(-1), hashing_value(0), el(El::X), hybrid(Hybridization::SP_NON),
      connectivity(0), metal_connectivity(0), min_ring_size(0),
      is_aromatic(false), is_metal(false), excess_electrons(0), charge(0.0f),
      par_charge(0.0f), bonding_idx(0) {}
};

// Statistical value with count
struct ValueStats {
  double value = NAN;
  double sigma = NAN;
  int count = 0;
  int level = 0;  // match specificity: 0=none, 1-4=aggregated, 10=full

  ValueStats() = default;
  ValueStats(double v, double s, int c, int lvl = 0) : value(v), sigma(s), count(c), level(lvl) {}
};

// Protonated hydrogen distances (both electron cloud and nucleus)
struct ProtHydrDist {
  double electron_val = NAN;
  double electron_sigma = NAN;
  double nucleus_val = NAN;
  double nucleus_sigma = NAN;
};

// Metal bond entry
struct MetalBondEntry {
  Element metal = El::X;
  Element ligand = El::X;
  int metal_coord = 0;
  int ligand_coord = 0;
  std::string ligand_class;
  double pre_value = NAN;
  double pre_sigma = NAN;
  int pre_count = 0;
  double value = NAN;
  double sigma = NAN;
  int count = 0;
};

// Metal coordination angle entry
struct MetalAngleEntry {
  Element metal = El::X;
  int coord_number = 0;
  CoordGeometry geometry = CoordGeometry::UNKNOWN;
  double angle = NAN;
  double sigma = NAN;
};

struct TorsionEntry {
  double value = 0.0;
  int period = 0;
  int priority = 0;
  std::string id;
};

// Periodic table row and group information
inline int element_row(Element el) {
  int n = el.ordinal();
  if (n == 0) return 0;
  if (n <= 2) return 1;
  if (n <= 10) return 2;
  if (n <= 18) return 3;
  if (n <= 36) return 4;
  if (n <= 54) return 5;
  if (n <= 86) return 6;
  return 7;
}


// ============================================================================
// Main AcedrgTables class
// ============================================================================

struct GEMMI_DLL AcedrgTables {
public:
  AcedrgTables() = default;

  // Load all tables from directory
  void load_tables(const std::string& tables_dir);

  // Process a ChemComp - fill all missing restraint values
  void fill_restraints(ChemComp& cc) const;

  // Assign CCP4 atom energy types (type_energy) following AceDRG rules
  void assign_ccp4_types(ChemComp& cc) const;
  // Adjust charges for atoms bonded to metals using AceDRG valence rules
  void apply_metal_charge_corrections(ChemComp& cc) const;

  bool lookup_pep_tors(const std::string& a1, const std::string& a2,
                                 const std::string& a3, const std::string& a4,
                                 TorsionEntry& out) const;

  // Individual lookups - returns match level (10=full, 4+=neighbor matched, 0-3=aggregated)
  int fill_bond(const ChemComp& cc,
                const std::vector<CodAtomInfo>& atom_info,
                Restraints::Bond& bond) const;
  void fill_angle(const ChemComp& cc,
                            const std::vector<CodAtomInfo>& atom_info,
                            Restraints::Angle& angle) const;

  // Atom classification - returns info for all atoms
  std::vector<CodAtomInfo> classify_atoms(const ChemComp& cc) const;

  // Compute acedrg_type string (like acedrg --typeOut)
  // Format: CentralElement(Neighbor1_desc)(Neighbor2_desc)...
  // where each neighbor description = neighbor element + sorted neighbor's other neighbors
  std::string compute_acedrg_type(const CodAtomInfo& atom,
                                  const std::vector<CodAtomInfo>& atoms,
                                  const std::vector<std::vector<int>>& neighbors) const;
  std::vector<std::string> compute_acedrg_types(const ChemComp& cc) const;

  // Configuration
  double upper_bond_sigma = 0.02;
  double lower_bond_sigma = 0.01;
  double upper_angle_sigma = 3.0;
  double lower_angle_sigma = 1.5;
  int min_observations_angle = 3;  // AceDRG default for angles
  int min_observations_angle_fallback = 3;
  int min_observations_bond = 3;   // AceDRG default for bonds (aNumTh=3)
  int metal_class_min_count = 5; // AceDRG uses >5 for metal class selection
  int verbose = 0;  // Debug output level (0=off, 1=basic, 2=detailed)

  static constexpr int HASH_SIZE = 1000;
  // Table directory
  std::string tables_dir_;
  bool tables_loaded_ = false;

  // Hash code tables
  std::map<int, std::string> digit_keys_;  // hash -> footprint
  std::map<int, int> linked_hash_;         // hash -> linked hash

  // HRS (High-Resolution Summary) bond tables
  // Key: hash1, hash2, hybrid_pair, in_ring
  struct BondHRSKey {
    int hash1, hash2;
    std::string hybrid_pair;
    std::string in_ring;
    bool operator<(const BondHRSKey& o) const {
      if (hash1 != o.hash1) return hash1 < o.hash1;
      if (hash2 != o.hash2) return hash2 < o.hash2;
      if (hybrid_pair != o.hybrid_pair) return hybrid_pair < o.hybrid_pair;
      return in_ring < o.in_ring;
    }
  };
  std::map<BondHRSKey, ValueStats> bond_hrs_;

  // HRS angle tables
  // Key: hash1, hash2, hash3, value_key (ring:hybr_tuple)
  struct AngleHRSKey {
    int hash1, hash2, hash3;
    std::string value_key;
    bool operator<(const AngleHRSKey& o) const {
      if (hash1 != o.hash1) return hash1 < o.hash1;
      if (hash2 != o.hash2) return hash2 < o.hash2;
      if (hash3 != o.hash3) return hash3 < o.hash3;
      return value_key < o.value_key;
    }
  };
  std::map<AngleHRSKey, ValueStats> angle_hrs_;

  struct Ccp4AtomInfo {
    Element el = El::X;
    std::string chem_type;
    std::string ccp4_type;
    int bonding_idx = 0;
    std::map<std::string, int> ring_rep;
    std::vector<int> conn_atoms;
    std::vector<int> conn_atoms_no_metal;
    std::vector<int> conn_h_atoms;
    float par_charge = 0.0f;
    int formal_charge = 0;
  };

  struct Ccp4BondEntry {
    double length = NAN;
    double sigma = NAN;
  };

  std::map<std::string, std::map<std::string, std::map<std::string, Ccp4BondEntry>>> ccp4_bonds_;

  static int ccp4_material_type(Element el);
  static void set_one_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx);
  static void set_hydro_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx);
  static void set_org_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx);
  static std::string bond_order_key(BondType type);
  void load_ccp4_bonds(const std::string& path);
  std::vector<std::string> compute_ccp4_types(const ChemComp& cc,
                                              const std::vector<CodAtomInfo>& atom_info,
                                              const std::vector<std::vector<int>>& neighbors) const;
  bool search_ccp4_bond(const std::string& type1,
                        const std::string& type2,
                        const std::string& order,
                        ValueStats& out) const;

  // Detailed indexed bond tables from allOrgBondTables/*.table
  // Level 0: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2, a1NB, a2NB, a1TypeM, a2TypeM
  using BondIdx1D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>>>>>;
  BondIdx1D bond_idx_1d_;

  // Exact match with full COD class (a1TypeF/a2TypeF) and main types
  using BondIdxFull = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    ValueStats>>>>>>>>>>>>;
  BondIdxFull bond_idx_full_;

  // Levels 3-6: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2, a1NB, a2NB (no atom types)
  using BondIdx2D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>>>;
  BondIdx2D bond_idx_2d_;

  // Level Nb2D: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2 (no nb1nb2_sp)
  using BondNb2D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>;
  BondNb2D bond_nb2d_;

  // Level Nb2DType: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2, type1, type2 (with atom types)
  using BondNb2DType = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>>>;
  BondNb2DType bond_nb2d_type_;

  // Levels 9-11: Hash+Sp fallback structures
  // Level 9: ha1, ha2, hybrComb, inRing
  using BondHaSp2D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>;
  BondHaSp2D bond_hasp_2d_;

  // Level 10: ha1, ha2, hybrComb only
  using BondHaSp1D = std::map<int, std::map<int,
    std::map<std::string, std::vector<ValueStats>>>>;
  BondHaSp1D bond_hasp_1d_;

  // Level 11: ha1, ha2 only
  using BondHaSp0D = std::map<int, std::map<int, std::vector<ValueStats>>>;
  BondHaSp0D bond_hasp_0d_;

  // Bond file index: maps (ha1, ha2) -> table file number
  std::map<int, std::map<int, int>> bond_file_index_;

  // Atom type code mapping: coded -> full type string
  std::map<std::string, std::string> atom_type_codes_;

  // Detailed indexed angle tables from allOrgAngleTables/*.table
  // Angles have 3 hashes (center, flank1, flank2)
  using AngleTypes = std::map<std::string,
    std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>;
  using AngleNB = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleTypes>>>;
  using AngleNB2 = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleNB>>>;
  using AngleRoots = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleNB2>>>;

  // Level 1D: Full detail with atom types
  using AngleIdx1D = std::map<int, std::map<int, std::map<int,
    std::map<std::string, AngleRoots>>>>;
  AngleIdx1D angle_idx_1d_;

  using AngleNBNoTypes = std::map<std::string,
    std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>;
  using AngleNB2NoTypes = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleNBNoTypes>>>;
  using AngleRootsNoTypes = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleNB2NoTypes>>>;

  // Level 2D: No atom types
  using AngleIdx2D = std::map<int, std::map<int, std::map<int,
    std::map<std::string, AngleRootsNoTypes>>>>;
  AngleIdx2D angle_idx_2d_;

  using AngleNB2Only = std::map<std::string,
    std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>;
  using AngleRootsNB2 = std::map<std::string,
    std::map<std::string, std::map<std::string, AngleNB2Only>>>;

  // Level 3D: Hash + valueKey + roots + NB2 only
  using AngleIdx3D = std::map<int, std::map<int, std::map<int,
    std::map<std::string, AngleRootsNB2>>>>;
  AngleIdx3D angle_idx_3d_;

  using AngleRootsOnly = std::map<std::string,
    std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>;

  // Level 4D: Hash + valueKey + roots
  using AngleIdx4D = std::map<int, std::map<int, std::map<int,
    std::map<std::string, AngleRootsOnly>>>>;
  AngleIdx4D angle_idx_4d_;

  // Level 5D: Hash + valueKey only
  using AngleIdx5D = std::map<int, std::map<int, std::map<int,
    std::map<std::string, std::vector<ValueStats>>>>>; // 5 closes
  AngleIdx5D angle_idx_5d_;

  // Level 6D: Hash only (4 levels)
  using AngleIdx6D = std::map<int, std::map<int, std::map<int,
    std::vector<ValueStats>>>>; // 4 closes
  AngleIdx6D angle_idx_6d_;

  // Angle file index: maps (ha1, ha2, ha3) -> table file number
  std::map<int, std::map<int, std::map<int, int>>> angle_file_index_;

  // Element + hybridization based fallback bonds
  using ENBonds = std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>;
  ENBonds en_bonds_;

  // Metal bond tables
  std::vector<MetalBondEntry> metal_bonds_;
  std::map<std::string, double> covalent_radii_;
  std::vector<MetalAngleEntry> metal_angles_;
  std::map<Element, std::map<int, CoordGeometry>> metal_coord_geo_;
  std::map<std::string, TorsionEntry> pep_tors_;

  // Protonated hydrogen distances: maps type (e.g., "H_sp3_C") -> ProtHydrDist
  std::map<std::string, ProtHydrDist> prot_hydr_dists_;

  // Internal helper functions

  // Loading functions
  void load_hash_codes(const std::string& path);
  void load_bond_hrs(const std::string& path);
  void load_angle_hrs(const std::string& path);
  void load_metal_tables(const std::string& dir);
  void load_covalent_radii(const std::string& path);
  void load_en_bonds(const std::string& path);
  void load_atom_type_codes(const std::string& path);
  void load_bond_index(const std::string& path);
  void load_bond_tables(const std::string& dir);
  void load_pep_tors(const std::string& path);
  void load_prot_hydr_dists(const std::string& path);
  void load_angle_index(const std::string& path);
  void load_angle_tables(const std::string& dir);

  // Atom classification helpers
  struct BondInfo {
    int neighbor_idx;
    BondType type;
    bool aromatic = false;
  };
  struct RingInfo {
    std::vector<int> atoms;
    std::string rep;
    std::string s_rep;
    bool is_aromatic = false;
    bool is_aromatic_permissive = false;
  };
  struct SortMap {
    std::string key;
    int val = 0;
  };
  struct SortMap2 {
    std::string key;
    int val = 0;
    int nNB = 0;
  };
  struct NB1stFam {
    std::string name;
    std::vector<std::string> NB2ndList;
    int repN = 1;
  };

  void set_atom_cod_class_name_new2(CodAtomInfo& atom,
                                    const CodAtomInfo& ori_atom,
                                    int lev,
                                    const std::vector<CodAtomInfo>& atoms,
                                    const std::vector<std::vector<int>>& neighbors) const;
  void set_special_3nb_symb2(CodAtomInfo& atom,
                             const std::vector<CodAtomInfo>& atoms,
                             const std::vector<std::vector<int>>& neighbors) const;
  void cod_class_to_atom2(const std::string& cod_class,
                          CodAtomInfo& atom) const;
  void set_atoms_nb1nb2_sp(std::vector<CodAtomInfo>& atoms,
                           const std::vector<std::vector<int>>& neighbors) const;
  void set_atoms_nb_symb_from_neighbors(std::vector<CodAtomInfo>& atoms,
                                        const std::vector<std::vector<int>>& neighbors) const;
  void set_atoms_bonding_and_chiral_center(std::vector<CodAtomInfo>& atoms,
                                           const std::vector<std::vector<int>>& neighbors) const;
 private:
  std::vector<std::vector<BondInfo>> build_adjacency(const ChemComp& cc) const;
  std::vector<std::vector<int>> build_neighbors(const std::vector<std::vector<BondInfo>>& adj) const;
  void set_ring_aromaticity_from_bonds(const std::vector<std::vector<BondInfo>>& adj,
                                       const std::vector<CodAtomInfo>& atoms,
                                       std::vector<RingInfo>& rings) const;
  void detect_rings_acedrg(const std::vector<std::vector<int>>& neighbors,
                           std::vector<CodAtomInfo>& atoms,
                           std::vector<RingInfo>& rings) const;
  void check_one_path_acedrg(const std::vector<std::vector<int>>& neighbors,
                             std::vector<CodAtomInfo>& atoms,
                             std::vector<RingInfo>& rings,
                             std::map<std::string, int>& ring_index,
                             int ori_idx,
                             int cur_idx,
                             int prev_idx,
                             int cur_lev,
                             int max_ring,
                             std::map<int, std::string>& seen_atom_ids,
                             std::map<int, std::string>& atom_ids_in_path) const;
  void set_atoms_ring_rep_s(std::vector<CodAtomInfo>& atoms,
                            std::vector<RingInfo>& rings) const;
  int get_num_oxy_connect(const std::vector<CodAtomInfo>& atoms,
                          const CodAtomInfo& atom,
                          const std::vector<std::vector<int>>& neighbors) const;
  int get_min_ring2_from_cod_class(const std::string& cod_class) const;
  bool cod_class_is_aromatic(const std::string& cod_class) const;
  void get_small_family(const std::string& in_str, NB1stFam& fam) const;
  bool are_in_same_ring(const CodAtomInfo& a1, const CodAtomInfo& a2) const;
  int angle_ring_size(const CodAtomInfo& center,
                      const CodAtomInfo& a1,
                      const CodAtomInfo& a3) const;
  void order_bond_atoms(const CodAtomInfo& a1, const CodAtomInfo& a2,
                        const CodAtomInfo*& first,
                        const CodAtomInfo*& second) const;
  void order_angle_flanks(const CodAtomInfo& a1, const CodAtomInfo& a3,
                          const CodAtomInfo*& flank1,
                          const CodAtomInfo*& flank3) const;
  Hybridization hybrid_from_bonding_idx(int bonding_idx, bool is_metal,
                                        int connectivity) const;
  void compute_hash(CodAtomInfo& atom) const;
  static int str_to_int(const std::string& s);
  static bool compare_no_case(const std::string& first, const std::string& second);
  static bool compare_no_case2(const std::string& first, const std::string& second);
  static bool desc_sort_map_key(const SortMap& a, const SortMap& b);
  static bool desc_sort_map_key2(const SortMap2& a, const SortMap2& b);

  // Bond search helpers
  ValueStats search_bond_multilevel(const CodAtomInfo& a1,
                                    const CodAtomInfo& a2) const;
  ValueStats search_bond_hrs(const CodAtomInfo& a1, const CodAtomInfo& a2,
                             bool in_ring) const;
  ValueStats search_bond_en(const CodAtomInfo& a1, const CodAtomInfo& a2) const;
  ProtHydrDist search_prot_hydr_dist(const CodAtomInfo& h_atom,
                                     const CodAtomInfo& heavy_atom) const;
  ValueStats search_metal_bond(const CodAtomInfo& metal,
                               const CodAtomInfo& ligand,
                               const std::vector<CodAtomInfo>& atoms) const;
  // Angle search helpers
  ValueStats search_angle_multilevel(const CodAtomInfo& a1,
                                     const CodAtomInfo& center,
                                     const CodAtomInfo& a3) const;
  ValueStats search_angle_hrs(const CodAtomInfo& a1, const CodAtomInfo& center,
                              const CodAtomInfo& a3, int ring_size) const;
  std::vector<double> get_metal_angles(Element metal, int coord_number) const;

  // Statistical aggregation
  ValueStats aggregate_stats(const std::vector<ValueStats>& values) const;

  // Utility: clamp sigma to reasonable range
  double clamp_bond_sigma(double sigma) const {
    return std::max(lower_bond_sigma, std::min(upper_bond_sigma, sigma));
  }
  double clamp_angle_sigma(double sigma) const {
    return std::max(lower_angle_sigma, std::min(upper_angle_sigma, sigma));
  }
};

/// Run the full restraint-generation pipeline on a ChemComp:
/// chemical-group adjustments, protonation, fill_restraints,
/// torsion/chirality/plane generation, and CCP4 type assignment.
/// \param atom_stereo  maps atom names to pdbx_stereo_config strings
///                     (needed for chirality generation).
GEMMI_DLL void prepare_chemcomp(ChemComp& cc, const AcedrgTables& tables,
                                const std::map<std::string, std::string>& atom_stereo = {});

} // namespace gemmi

#endif
