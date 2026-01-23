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
#include <fstream>
#include <sstream>
#include "chemcomp.hpp"
#include "elem.hpp"
#include "fail.hpp"
#include "util.hpp"
#include "read_cif.hpp"

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

inline const char* hybridization_to_string(Hybridization h) {
  switch (h) {
    case Hybridization::SP1: return "SP1";
    case Hybridization::SP2: return "SP2";
    case Hybridization::SP3: return "SP3";
    case Hybridization::SPD5: return "SPD5";
    case Hybridization::SPD6: return "SPD6";
    case Hybridization::SPD7: return "SPD7";
    case Hybridization::SPD8: return "SPD8";
    case Hybridization::SP_NON: default: return "SP-NON";
  }
}

inline Hybridization hybridization_from_string(const std::string& s) {
  if (s == "SP1") return Hybridization::SP1;
  if (s == "SP2") return Hybridization::SP2;
  if (s == "SP3") return Hybridization::SP3;
  if (s == "SPD5") return Hybridization::SPD5;
  if (s == "SPD6") return Hybridization::SPD6;
  if (s == "SPD7") return Hybridization::SPD7;
  if (s == "SPD8") return Hybridization::SPD8;
  return Hybridization::SP_NON;
}

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

  ValueStats() = default;
  ValueStats(double v, double s, int c) : value(v), sigma(s), count(c) {}
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

inline int element_group(Element el) {
  // Lookup table for periodic table groups (1-18) by element ordinal (0-118)
  // Lanthanides (57-71) and actinides (89-103) are assigned to group 3
  static constexpr int groups[119] = {
    // 0: unknown
    0,
    // 1-2: H, He
    1, 18,
    // 3-10: Li-Ne
    1, 2, 13, 14, 15, 16, 17, 18,
    // 11-18: Na-Ar
    1, 2, 13, 14, 15, 16, 17, 18,
    // 19-36: K-Kr
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    // 37-54: Rb-Xe
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    // 55-56: Cs, Ba
    1, 2,
    // 57-71: La-Lu (lanthanides)
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    // 72-86: Hf-Rn
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    // 87-88: Fr, Ra
    1, 2,
    // 89-103: Ac-Lr (actinides)
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    // 104-118: Rf-Og
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
  };
  int n = el.ordinal();
  return (n >= 0 && n <= 118) ? groups[n] : 0;
}

// ============================================================================
// Main AcedrgTables class
// ============================================================================

class AcedrgTables {
public:
  AcedrgTables();

  // Load all tables from directory
  void load_tables(const std::string& tables_dir);

  // Check if tables are loaded
  bool tables_loaded() const { return tables_loaded_; }

  // Get the tables directory
  const std::string& tables_dir() const { return tables_dir_; }

  // Process a ChemComp - fill all missing restraint values
  void fill_restraints(ChemComp& cc) const;

  // Assign CCP4 atom energy types (type_energy) following AceDRG rules
  void assign_ccp4_types(ChemComp& cc) const;

  bool lookup_pep_tors(const std::string& a1, const std::string& a2,
                       const std::string& a3, const std::string& a4,
                       TorsionEntry& out) const;

  // Individual lookups
  void fill_bond(const ChemComp& cc,
                 const std::vector<CodAtomInfo>& atom_info,
                 Restraints::Bond& bond) const;
  void fill_angle(const ChemComp& cc,
                  const std::vector<CodAtomInfo>& atom_info,
                  Restraints::Angle& angle) const;

  // Atom classification - returns info for all atoms
  std::vector<CodAtomInfo> classify_atoms(const ChemComp& cc) const;

  // Configuration
  double upper_bond_sigma = 0.2;
  double lower_bond_sigma = 0.02;
  double upper_angle_sigma = 3.0;
  double lower_angle_sigma = 1.5;
  int min_observations_angle = 3;  // AceDRG default for angles
  int min_observations_angle_fallback = 3;
  int min_observations_bond = 4;   // AceDRG default for bonds
  int metal_class_min_count = 5; // AceDRG uses >5 for metal class selection
  int verbose = 0;  // Debug output level (0=off, 1=basic, 2=detailed)

private:
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
  std::vector<MetalAngleEntry> metal_angles_;
  std::map<Element, std::map<int, CoordGeometry>> metal_coord_geo_;
  std::map<std::string, TorsionEntry> pep_tors_;

  // Internal helper functions

  // Loading functions
  void load_hash_codes(const std::string& path);
  void load_bond_hrs(const std::string& path);
  void load_angle_hrs(const std::string& path);
  void load_metal_tables(const std::string& dir);
  void load_en_bonds(const std::string& path);
  void load_atom_type_codes(const std::string& path);
  void load_bond_index(const std::string& path);
  void load_bond_tables(const std::string& dir);
  void load_pep_tors(const std::string& path);
  void load_angle_index(const std::string& path);
  void load_angle_tables(const std::string& dir);

  // Atom classification helpers
  struct BondInfo {
    int neighbor_idx;
    BondType type;
  };
  struct RingInfo {
    std::vector<int> atoms;
    std::string rep;
    std::string s_rep;
    bool is_aromatic = false;
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

  std::vector<std::vector<BondInfo>> build_adjacency(const ChemComp& cc) const;
  std::vector<std::vector<int>> build_neighbors(const std::vector<std::vector<BondInfo>>& adj) const;
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
  void set_atoms_bonding_and_chiral_center(std::vector<CodAtomInfo>& atoms,
                                           const std::vector<std::vector<int>>& neighbors) const;
  int get_num_oxy_connect(const std::vector<CodAtomInfo>& atoms,
                          const CodAtomInfo& atom,
                          const std::vector<std::vector<int>>& neighbors) const;
  int get_min_ring2_from_cod_class(const std::string& cod_class) const;
  bool cod_class_is_aromatic(const std::string& cod_class) const;
  void get_small_family(const std::string& in_str, NB1stFam& fam) const;
  static std::string trim_spaces(const std::string& s);
  static std::vector<std::string> split(const std::string& s, char delim);
  static int str_to_int(const std::string& s);
  static bool compare_no_case(const std::string& first, const std::string& second);
  static bool compare_no_case2(const std::string& first, const std::string& second);
  static bool des_sort_map_key(const SortMap& a, const SortMap& b);
  static bool des_sort_map_key2(const SortMap2& a, const SortMap2& b);
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

  // Bond search helpers
  ValueStats search_bond_multilevel(const CodAtomInfo& a1,
                                    const CodAtomInfo& a2) const;
  ValueStats search_bond_hrs(const CodAtomInfo& a1, const CodAtomInfo& a2,
                             bool in_ring) const;
  ValueStats search_bond_en(const CodAtomInfo& a1, const CodAtomInfo& a2) const;
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

// ============================================================================
// Implementation - Loading functions
// ============================================================================

inline AcedrgTables::AcedrgTables() {}

inline void AcedrgTables::load_tables(const std::string& tables_dir) {
  tables_dir_ = tables_dir;

  // Load hash code mapping
  load_hash_codes(tables_dir + "/allOrgLinkedHashCode.table");

  // Load HRS (summary) tables
  load_bond_hrs(tables_dir + "/allOrgBondsHRS.table");
  load_angle_hrs(tables_dir + "/allOrgAnglesHRS.table");

  // Load element+hybridization fallback
  load_en_bonds(tables_dir + "/allOrgBondEN.table");

  // Load metal tables
  load_metal_tables(tables_dir);

  // Load CCP4 energetic library bonds (for AceDRG fallback)
  load_ccp4_bonds(tables_dir + "/ener_lib.cif");

  // Load detailed indexed tables
  load_atom_type_codes(tables_dir + "/allAtomTypesFromMolsCoded.list");
  load_bond_index(tables_dir + "/allOrgBondTables/bond_idx.table");
  load_bond_tables(tables_dir + "/allOrgBondTables");
  load_angle_index(tables_dir + "/allOrgAngleTables/angle_idx.table");
  load_angle_tables(tables_dir + "/allOrgAngleTables");
  load_pep_tors(tables_dir + "/pep_tors.table");

  tables_loaded_ = true;
}

inline void AcedrgTables::load_hash_codes(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash_code, linked_hash;
    std::string footprint;

    if (iss >> hash_code >> footprint >> linked_hash) {
      digit_keys_[hash_code] = footprint;
      linked_hash_[hash_code] = linked_hash;
    }
  }
}

inline void AcedrgTables::load_bond_hrs(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    fail("Cannot open bond HRS table: ", path);

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash1, hash2;
    std::string hybrid_pair, in_ring;
    double value, sigma;
    int count;

    if (iss >> hash1 >> hash2 >> hybrid_pair >> in_ring >> value >> sigma >> count) {
      BondHRSKey key;
      key.hash1 = std::min(hash1, hash2);
      key.hash2 = std::max(hash1, hash2);
      key.hybrid_pair = hybrid_pair;
      key.in_ring = (in_ring == "Y" || in_ring == "y") ? "Y" : "N";
      bond_hrs_[key] = ValueStats(value, sigma, count);
    }
  }
}

inline void AcedrgTables::load_angle_hrs(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    fail("Cannot open angle HRS table: ", path);

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash1, hash2, hash3;
    std::string value_key;
    std::string a1_cod, a2_cod, a3_cod;
    double value1, sigma1;
    int count1;
    double value2, sigma2;
    int count2;

    if (iss >> hash1 >> hash2 >> hash3 >> value_key
        >> a1_cod >> a2_cod >> a3_cod
        >> value1 >> sigma1 >> count1
        >> value2 >> sigma2 >> count2) {
      AngleHRSKey key;
      // For angles, center is always hash2, but we canonicalize flanking atoms
      bool swap_flanks = hash1 > hash3;
      key.hash1 = swap_flanks ? hash3 : hash1;
      key.hash2 = hash2;
      key.hash3 = swap_flanks ? hash1 : hash3;
      if (swap_flanks) {
        size_t colon = value_key.find(':');
        if (colon != std::string::npos) {
          std::string ring_part = value_key.substr(0, colon);
          std::string hybrid_part = value_key.substr(colon + 1);
          std::vector<std::string> parts = split(hybrid_part, '_');
          if (parts.size() == 3) {
            std::swap(parts[0], parts[2]);
            hybrid_part = parts[0] + "_" + parts[1] + "_" + parts[2];
          }
          value_key = ring_part + ":" + hybrid_part;
        }
      }
      key.value_key = value_key;
      angle_hrs_[key] = ValueStats(value1, sigma1, count1);
    }
  }
}

inline void AcedrgTables::load_en_bonds(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string elem1, sp1, elem2, sp2;
    double value, sigma;
    int count;

    if (iss >> elem1 >> sp1 >> elem2 >> sp2 >> value >> sigma >> count) {
      // Canonicalize order
      if (elem1 > elem2 || (elem1 == elem2 && sp1 > sp2)) {
        std::swap(elem1, elem2);
        std::swap(sp1, sp2);
      }
      en_bonds_[elem1][sp1][elem2][sp2].emplace_back(value, sigma, count);
    }
  }
}

inline void AcedrgTables::load_metal_tables(const std::string& dir) {
  // Load allMetalBonds.table
  std::ifstream f(dir + "/allMetalBonds.table");
  if (!f)
    return;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token)
      tokens.push_back(token);

    if (tokens.size() == 10 || tokens.size() == 11) {
      MetalBondEntry entry;
      entry.metal = Element(tokens[0]);
      entry.metal_coord = str_to_int(tokens[1]);
      entry.ligand = Element(tokens[2]);
      entry.ligand_coord = str_to_int(tokens[3]);
      entry.pre_value = std::stod(tokens[4]);
      entry.pre_sigma = std::stod(tokens[5]);
      entry.pre_count = str_to_int(tokens[6]);
      if (tokens.size() == 11) {
        entry.ligand_class = tokens[7];
        entry.value = std::stod(tokens[8]);
        entry.sigma = std::stod(tokens[9]);
        entry.count = str_to_int(tokens[10]);
      } else {
        entry.ligand_class = "NONE";
        entry.value = std::stod(tokens[7]);
        entry.sigma = std::stod(tokens[8]);
        entry.count = str_to_int(tokens[9]);
      }
      metal_bonds_.push_back(entry);
    }
  }

  // Load metal coordination geometry
  std::ifstream f2(dir + "/allMetalDefCoordGeos.table");
  if (f2) {
    while (std::getline(f2, line)) {
      if (line.empty() || line[0] == '#')
        continue;

      std::istringstream iss(line);
      std::string metal_str, geo_str;
      int coord;

      if (iss >> metal_str >> coord >> geo_str) {
        Element metal(metal_str);
        CoordGeometry geo = CoordGeometry::UNKNOWN;
        if (geo_str == "LINEAR") geo = CoordGeometry::LINEAR;
        else if (geo_str == "TRIGONAL_PLANAR") geo = CoordGeometry::TRIGONAL_PLANAR;
        else if (geo_str == "T_SHAPED") geo = CoordGeometry::T_SHAPED;
        else if (geo_str == "TETRAHEDRAL") geo = CoordGeometry::TETRAHEDRAL;
        else if (geo_str == "SQUARE_PLANAR") geo = CoordGeometry::SQUARE_PLANAR;
        else if (geo_str == "TRIGONAL_BIPYRAMIDAL") geo = CoordGeometry::TRIGONAL_BIPYRAMIDAL;
        else if (geo_str == "SQUARE_PYRAMIDAL") geo = CoordGeometry::SQUARE_PYRAMIDAL;
        else if (geo_str == "OCTAHEDRAL") geo = CoordGeometry::OCTAHEDRAL;
        else if (geo_str == "TRIGONAL_PRISM") geo = CoordGeometry::TRIGONAL_PRISM;
        else if (geo_str == "PENTAGONAL_BIPYRAMIDAL") geo = CoordGeometry::PENTAGONAL_BIPYRAMIDAL;
        else if (geo_str == "CAPPED_OCTAHEDRAL") geo = CoordGeometry::CAPPED_OCTAHEDRAL;
        else if (geo_str == "SQUARE_ANTIPRISM") geo = CoordGeometry::SQUARE_ANTIPRISM;
        metal_coord_geo_[metal][coord] = geo;
      }
    }
  }

  // Load metal coordination angles
  std::ifstream f3(dir + "/allMetalCoordGeoAngles.table");
  if (f3) {
    while (std::getline(f3, line)) {
      if (line.empty() || line[0] == '#')
        continue;

      std::istringstream iss(line);
      std::string metal_str, geo_str;
      int coord;
      double angle, sigma;

      if (iss >> metal_str >> coord >> geo_str >> angle >> sigma) {
        MetalAngleEntry entry;
        entry.metal = Element(metal_str);
        entry.coord_number = coord;
        entry.angle = angle;
        entry.sigma = sigma;
        metal_angles_.push_back(entry);
      }
    }
  }
}

inline void AcedrgTables::load_atom_type_codes(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: code<whitespace>full_type
    std::istringstream iss(line);
    std::string code, full_type;
    if (iss >> code >> full_type) {
      atom_type_codes_[code] = full_type;
    }
  }
}

inline void AcedrgTables::load_bond_index(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: ha1 ha2 fileNum
    std::istringstream iss(line);
    int ha1, ha2, file_num;
    if (iss >> ha1 >> ha2 >> file_num) {
      bond_file_index_[ha1][ha2] = file_num;
    }
  }
}

inline void AcedrgTables::load_bond_tables(const std::string& dir) {
  // Load each bond table file referenced in the index
  std::set<int> loaded_files;

  for (const auto& ha1_pair : bond_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      int file_num = ha2_pair.second;
      if (loaded_files.count(file_num))
        continue;
      loaded_files.insert(file_num);

      std::string path = dir + "/" + std::to_string(file_num) + ".table";
      std::ifstream f(path);
      if (!f)
        continue;

      std::string line;
      while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#')
          continue;

        // Format: ha1 ha2 hybrComb inRing a1NB2 a2NB2 a1NB a2NB atomCode1 atomCode2
        //         value sigma count val2 sig2 cnt2
        std::istringstream iss(line);
        int ha1, ha2;
        std::string hybr_comb, in_ring, a1_nb2, a2_nb2, a1_nb, a2_nb;
        std::string atom_code1, atom_code2;
        double value, sigma;
        int count;
        double value2, sigma2;
        int count2;

        if (!(iss >> ha1 >> ha2 >> hybr_comb >> in_ring
                  >> a1_nb2 >> a2_nb2 >> a1_nb >> a2_nb
                  >> atom_code1 >> atom_code2
                  >> value >> sigma >> count
                  >> value2 >> sigma2 >> count2))
          continue;

        // Get main atom types from codes
        std::string a1_type_m, a2_type_m;
        std::string a1_type_f, a2_type_f;
        auto it1 = atom_type_codes_.find(atom_code1);
        auto it2 = atom_type_codes_.find(atom_code2);
        if (it1 != atom_type_codes_.end()) {
          a1_type_f = it1->second;
          // Extract main type (before '{' if present)
          a1_type_m = a1_type_f;
          size_t brace = a1_type_m.find('{');
          if (brace != std::string::npos)
            a1_type_m = a1_type_m.substr(0, brace);
        }
        if (it2 != atom_type_codes_.end()) {
          a2_type_f = it2->second;
          a2_type_m = a2_type_f;
          size_t brace = a2_type_m.find('{');
          if (brace != std::string::npos)
            a2_type_m = a2_type_m.substr(0, brace);
        }

        // Extract root types (element + ring annotation, e.g., "C[5a]")
        std::string a1_root, a2_root;
        if (!a1_type_m.empty()) {
          size_t paren = a1_type_m.find('(');
          a1_root = (paren != std::string::npos) ? a1_type_m.substr(0, paren) : a1_type_m;
        }
        if (!a2_type_m.empty()) {
          size_t paren = a2_type_m.find('(');
          a2_root = (paren != std::string::npos) ? a2_type_m.substr(0, paren) : a2_type_m;
        }

        ValueStats vs(value, sigma, count);
        ValueStats vs1d(value2, sigma2, count2);

        // Populate 1D structure (full detail)
        bond_idx_1d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb][a1_type_m][a2_type_m].push_back(vs1d);

        // Store full COD-class stats for exact matches (AceDRG uses full codClass first).
        if (!a1_type_f.empty() && !a2_type_f.empty()) {
          bond_idx_full_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                         [a1_nb][a2_nb][a1_type_m][a2_type_m]
                         [a1_type_f][a2_type_f] = vs;
        }

        // Populate 2D structure (no atom types)
        bond_idx_2d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb].push_back(vs);

        // Populate Nb2D structure (nb2 only, no nb1nb2_sp)
        bond_nb2d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2].push_back(vs);

        // Populate Nb2DType structure (nb2 + root types, no nb1nb2_sp)
        if (!a1_root.empty() && !a2_root.empty())
          bond_nb2d_type_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                         [a1_root][a2_root].push_back(vs1d);

        // Populate HaSp structures for fallback
        bond_hasp_2d_[ha1][ha2][hybr_comb][in_ring].push_back(vs);
        bond_hasp_1d_[ha1][ha2][hybr_comb].push_back(vs);
        bond_hasp_0d_[ha1][ha2].push_back(vs);
      }
    }
  }
}

inline void AcedrgTables::load_angle_index(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: ha1 ha2 ha3 fileNum
    std::istringstream iss(line);
    int ha1, ha2, ha3, file_num;
    if (iss >> ha1 >> ha2 >> ha3 >> file_num) {
      angle_file_index_[ha1][ha2][ha3] = file_num;
    }
  }
}

inline void AcedrgTables::load_angle_tables(const std::string& dir) {
  // Load each angle table file referenced in the index
  std::set<int> loaded_files;

  for (const auto& ha1_pair : angle_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      for (const auto& ha3_pair : ha2_pair.second) {
        int file_num = ha3_pair.second;
        if (loaded_files.count(file_num))
          continue;
        loaded_files.insert(file_num);

        std::string path = dir + "/" + std::to_string(file_num) + ".table";
        std::ifstream f(path);
        if (!f)
          continue;

        std::string line;
        while (std::getline(f, line)) {
          if (line.empty() || line[0] == '#')
            continue;

          // 34-column format:
          // 1-3: ha1 ha2 ha3
          // 4: valueKey (ring:hybr_tuple, e.g. "0:SP2_SP2_SP3")
          // 5-7: a1_root a2_root a3_root
          // 8-10: a1_nb2 a2_nb2 a3_nb2
          // 11-13: a1_nb a2_nb a3_nb
          // 14-16: a1_code a2_code a3_code
          // 17-19: AxC value, sigma, count (unused)
          // 20-22: AxM value, sigma, count (1D)
          // 23-25: A_NB value, sigma, count (2D)
          // 26-28: A_NB2 value, sigma, count (3D)
          // 29-31: a1R/a2R/a3R value, sigma, count (4D)
          // 32-34: R3A value, sigma, count (5D)
          std::istringstream iss(line);
          int ha1, ha2, ha3;
          std::string value_key;
          std::string a1_root, a2_root, a3_root;
          std::string a1_nb2, a2_nb2, a3_nb2;
          std::string a1_nb, a2_nb, a3_nb;
          std::string a1_code, a2_code, a3_code;

          if (!(iss >> ha1 >> ha2 >> ha3 >> value_key
                    >> a1_root >> a2_root >> a3_root
                    >> a1_nb2 >> a2_nb2 >> a3_nb2
                    >> a1_nb >> a2_nb >> a3_nb
                    >> a1_code >> a2_code >> a3_code))
            continue;

          // Read 6 sets of value/sigma/count
          double values[6], sigmas[6];
          int counts[6];
          bool ok = true;
          for (int lvl = 0; lvl < 6 && ok; ++lvl) {
            if (!(iss >> values[lvl] >> sigmas[lvl] >> counts[lvl]))
              ok = false;
          }
          if (!ok)
            continue;

          // Get main atom types from codes
          std::string a1_type, a2_type, a3_type;
          auto it1 = atom_type_codes_.find(a1_code);
          auto it2 = atom_type_codes_.find(a2_code);
          auto it3 = atom_type_codes_.find(a3_code);
          if (it1 != atom_type_codes_.end()) {
            a1_type = it1->second;
            size_t brace = a1_type.find('{');
            if (brace != std::string::npos)
              a1_type = a1_type.substr(0, brace);
          }
          if (it2 != atom_type_codes_.end()) {
            a2_type = it2->second;
            size_t brace = a2_type.find('{');
            if (brace != std::string::npos)
              a2_type = a2_type.substr(0, brace);
          }
          if (it3 != atom_type_codes_.end()) {
            a3_type = it3->second;
            size_t brace = a3_type.find('{');
            if (brace != std::string::npos)
              a3_type = a3_type.substr(0, brace);
          }

          // Populate structures at each level with corresponding pre-computed values.
          // AceDRG keeps only the first entry for each key (no aggregation).
          // Level 1D: full detail with atom types
          ValueStats vs1(values[1], sigmas[1], counts[1]);
          auto& angle_1d_vec = angle_idx_1d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb]
                               [a1_type][a2_type][a3_type];
          if (angle_1d_vec.empty())
            angle_1d_vec.push_back(vs1);

          // Level 2D: no atom types
          ValueStats vs2(values[2], sigmas[2], counts[2]);
          auto& angle_2d_vec = angle_idx_2d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb];
          if (angle_2d_vec.empty())
            angle_2d_vec.push_back(vs2);

          // Level 3D: hash + valueKey + NB2 only
          ValueStats vs3(values[3], sigmas[3], counts[3]);
          auto& angle_3d_vec = angle_idx_3d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2];
          if (angle_3d_vec.empty())
            angle_3d_vec.push_back(vs3);

          // Level 4D: hash + valueKey + roots
          ValueStats vs4(values[4], sigmas[4], counts[4]);
          auto& angle_4d_vec = angle_idx_4d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root];
          if (angle_4d_vec.empty())
            angle_4d_vec.push_back(vs4);

          // Level 5D: hash + valueKey only
          ValueStats vs5(values[5], sigmas[5], counts[5]);
          auto& angle_5d_vec = angle_idx_5d_[ha1][ha2][ha3][value_key];
          if (angle_5d_vec.empty())
            angle_5d_vec.push_back(vs5);

          // Level 6D: hash only (leave empty for 34-column data)
        }
      }
    }
  }
}

inline void AcedrgTables::load_pep_tors(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    std::string tors_id, label, a1, a2, a3, a4;
    int period = 0;
    int idx = 0;
    double value = 0.0;
    if (!(iss >> tors_id >> label >> a1 >> a2 >> a3 >> a4 >> period >> idx >> value))
      continue;
    pep_tors_[a1 + "_" + a2 + "_" + a3 + "_" + a4] = {value, period};
  }
}

// ============================================================================
// Implementation - Atom classification
// ============================================================================

inline std::vector<CodAtomInfo> AcedrgTables::classify_atoms(const ChemComp& cc) const {
  std::vector<CodAtomInfo> atoms(cc.atoms.size());

  // Initialize basic info
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    atoms[i].index = static_cast<int>(i);
    atoms[i].id = cc.atoms[i].id;
    atoms[i].el = cc.atoms[i].el;
    atoms[i].is_metal = is_metal(cc.atoms[i].el);
    atoms[i].charge = cc.atoms[i].charge;
    atoms[i].par_charge = cc.atoms[i].charge;
  }

  // Build adjacency and neighbor lists
  std::vector<std::vector<BondInfo>> adj = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adj);

  // Connectivity counts
  for (size_t i = 0; i < atoms.size(); ++i) {
    atoms[i].connectivity = static_cast<int>(neighbors[i].size());
    atoms[i].conn_atoms_no_metal.clear();
    int metal_conn = 0;
    for (int nb : neighbors[i]) {
      if (atoms[nb].is_metal) {
        ++metal_conn;
      } else {
        atoms[i].conn_atoms_no_metal.push_back(nb);
      }
    }
    atoms[i].metal_connectivity = metal_conn;
  }

  // Seed aromaticity from input bonds
  for (const auto& bond : cc.rt.bonds) {
    if (!bond.aromatic)
      continue;
    auto it1 = cc.find_atom(bond.id1.atom);
    auto it2 = cc.find_atom(bond.id2.atom);
    if (it1 != cc.atoms.end()) {
      int idx1 = static_cast<int>(it1 - cc.atoms.begin());
      atoms[idx1].is_aromatic = true;
    }
    if (it2 != cc.atoms.end()) {
      int idx2 = static_cast<int>(it2 - cc.atoms.begin());
      atoms[idx2].is_aromatic = true;
    }
  }

  // Detect rings and populate ring representations
  std::vector<RingInfo> rings;
  detect_rings_acedrg(neighbors, atoms, rings);
  set_atoms_ring_rep_s(atoms, rings);

  // Build COD class names (AceDRG style)
  for (size_t i = 0; i < atoms.size(); ++i)
    set_atom_cod_class_name_new2(atoms[i], atoms[i], 2, atoms, neighbors);

  for (size_t i = 0; i < atoms.size(); ++i)
    set_special_3nb_symb2(atoms[i], atoms, neighbors);

  for (size_t i = 0; i < atoms.size(); ++i)
    cod_class_to_atom2(atoms[i].cod_class, atoms[i]);

  // Hybridization and NB1/NB2_SP
  set_atoms_bonding_and_chiral_center(atoms, neighbors);
  set_atoms_nb1nb2_sp(atoms, neighbors);

  // nb2_symb (codNB2Symb) stays derived from cod_class (AceDRG behavior).

  // Ring props from codClass and hash codes
  for (auto& atom : atoms) {
    atom.min_ring_size = get_min_ring2_from_cod_class(atom.cod_class);
    atom.is_aromatic = cod_class_is_aromatic(atom.cod_class);
    atom.hybrid = hybrid_from_bonding_idx(atom.bonding_idx, atom.is_metal,
                                          atom.connectivity);
    compute_hash(atom);
  }

  return atoms;
}

// AceDrg hash used in *HRS.table files
inline void AcedrgTables::compute_hash(CodAtomInfo& atom) const {
  static const int primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229
  };

  // d1: aromaticity (0 or 1)
  int d1 = atom.is_aromatic ? 1 : 0;

  // d2: min ring size mapped to 2-7
  int d2;
  switch (atom.min_ring_size) {
    case 0: d2 = 2; break;
    case 3: d2 = 3; break;
    case 4: d2 = 4; break;
    case 5: d2 = 5; break;
    case 6: d2 = 6; break;
    default: d2 = 7; break;
  }

  // d3: connectivity + 8
  int d3 = 8 + atom.connectivity;

  // d4: periodic row + 16
  int d4 = 16 + element_row(atom.el);

  // d5: periodic group + 24
  int d5 = 24 + element_group(atom.el);

  // Compute hash as product of primes mod HASH_SIZE
  int64_t prime_product = static_cast<int64_t>(primes[d1]) *
                          primes[d2] * primes[d3] * primes[d4] * primes[d5];
  atom.hashing_value = static_cast<int>(prime_product % HASH_SIZE);

  // If we have hash tables loaded, resolve collisions
  if (!digit_keys_.empty()) {
    std::string footprint = std::to_string(d1) + "_" + std::to_string(d2) +
                            "_" + std::to_string(d3) + "_" + std::to_string(d4) +
                            "_" + std::to_string(d5);

    int pseudo_hash = atom.hashing_value;
    auto it = digit_keys_.find(pseudo_hash);

    if (it != digit_keys_.end()) {
      if (it->second == footprint) {
        atom.hashing_value = pseudo_hash;
      } else {
        // Follow linked hash chain
        bool found = false;
        while (!found) {
          auto link_it = linked_hash_.find(pseudo_hash);
          if (link_it == linked_hash_.end() || link_it->second == -1) {
            break;
          }
          pseudo_hash = link_it->second;
          auto key_it = digit_keys_.find(pseudo_hash);
          if (key_it != digit_keys_.end() && key_it->second == footprint) {
            atom.hashing_value = pseudo_hash;
            found = true;
          }
        }
      }
    }
  }
}

inline std::vector<std::vector<AcedrgTables::BondInfo>>
AcedrgTables::build_adjacency(const ChemComp& cc) const {
  std::vector<std::vector<BondInfo>> adj(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = cc.find_atom(bond.id1.atom);
    auto it2 = cc.find_atom(bond.id2.atom);
    if (it1 != cc.atoms.end() && it2 != cc.atoms.end()) {
      int idx1 = static_cast<int>(it1 - cc.atoms.begin());
      int idx2 = static_cast<int>(it2 - cc.atoms.begin());
      adj[idx1].push_back({idx2, bond.type});
      adj[idx2].push_back({idx1, bond.type});
    }
  }
  return adj;
}

inline std::vector<std::vector<int>> AcedrgTables::build_neighbors(
    const std::vector<std::vector<BondInfo>>& adj) const {
  std::vector<std::vector<int>> neighbors(adj.size());
  for (size_t i = 0; i < adj.size(); ++i) {
    neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i])
      neighbors[i].push_back(nb.neighbor_idx);
  }
  return neighbors;
}

inline void AcedrgTables::detect_rings_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) const {
  rings.clear();
  std::map<std::string, int> ring_index;

  for (size_t i = 0; i < atoms.size(); ++i) {
    if (atoms[i].is_metal || atoms[i].el == El::H)
      continue;
    std::map<int, std::string> seen;
    std::map<int, std::string> path;
    check_one_path_acedrg(neighbors, atoms, rings, ring_index,
                          static_cast<int>(i), static_cast<int>(i), -999,
                          1, 7, seen, path);
  }
}

inline void AcedrgTables::check_one_path_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings,
    std::map<std::string, int>& ring_index,
    int ori_idx, int cur_idx, int prev_idx, int cur_lev, int max_ring,
    std::map<int, std::string>& seen_atom_ids,
    std::map<int, std::string>& atom_ids_in_path) const {
  if (cur_lev >= max_ring)
    return;

  bool path_collision = false;
  for (int nb : neighbors[cur_idx]) {
    if (nb != ori_idx && nb != prev_idx &&
        seen_atom_ids.find(nb) != seen_atom_ids.end()) {
      path_collision = true;
      break;
    }
  }

  if (!path_collision) {
    for (int nb : neighbors[cur_idx]) {
      if (nb == ori_idx && nb != prev_idx && cur_lev > 2 && atoms[nb].el != El::H) {
        atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
        std::list<std::string> all_ids;
        std::list<std::string> all_seris;
        std::vector<int> ring_atoms;
        for (const auto& it : atom_ids_in_path) {
          all_seris.push_back(std::to_string(it.first));
          all_ids.push_back(it.second);
          ring_atoms.push_back(it.first);
        }
        all_seris.sort(compare_no_case);
        all_ids.sort(compare_no_case);

        std::string rep;
        for (const auto& id : all_ids)
          rep += id;

        std::string s_rep;
        int nrs = 0;
        for (const auto& seri : all_seris) {
          if (nrs == 0)
            s_rep += seri;
          else
            s_rep += "_" + seri;
          ++nrs;
        }

        atoms[ori_idx].ring_rep[rep] = static_cast<int>(atom_ids_in_path.size());

        if (ring_index.find(s_rep) == ring_index.end()) {
          RingInfo ring;
          ring.atoms = ring_atoms;
          ring.rep = rep;
          ring.s_rep = s_rep;
          ring.is_aromatic = true;
          for (int idx : ring_atoms)
            if (!atoms[idx].is_aromatic)
              ring.is_aromatic = false;

          int idx = static_cast<int>(rings.size());
          ring_index[s_rep] = idx;
          rings.push_back(ring);
          for (int atom_idx : ring_atoms)
            atoms[atom_idx].in_rings.push_back(idx);
        }

        atom_ids_in_path.erase(cur_idx);
        path_collision = true;
        break;
      }
    }
  }

  if (!path_collision) {
    int next_lev = cur_lev + 1;
    seen_atom_ids[cur_idx] = atoms[cur_idx].id;
    if (next_lev < max_ring) {
      if (cur_lev == 1) {
        seen_atom_ids.clear();
        atom_ids_in_path.clear();
        seen_atom_ids[cur_idx] = atoms[cur_idx].id;
        atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
      }
      atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
      for (int nb : neighbors[cur_idx]) {
        if (nb != prev_idx && !atoms[nb].is_metal) {
          check_one_path_acedrg(neighbors, atoms, rings, ring_index,
                                ori_idx, nb, cur_idx, next_lev, max_ring,
                                seen_atom_ids, atom_ids_in_path);
        }
      }
      atom_ids_in_path.erase(cur_idx);
      seen_atom_ids.erase(cur_idx);
    }
    atom_ids_in_path.erase(cur_idx);
    seen_atom_ids.erase(cur_idx);
  }
}

inline void AcedrgTables::set_atoms_ring_rep_s(
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) const {
  for (const auto& ring : rings) {
    std::string size = std::to_string(ring.atoms.size());
    std::string rep_id;
    std::list<std::string> all_seris;
    for (int idx : ring.atoms)
      all_seris.push_back(std::to_string(idx));
    all_seris.sort(compare_no_case);
    int nrs = 0;
    for (const auto& seri : all_seris) {
      if (nrs == 0)
        rep_id += seri;
      else
        rep_id += "_" + seri;
      ++nrs;
    }

    for (int idx : ring.atoms) {
      if (ring.is_aromatic)
        atoms[idx].ring_rep_s[rep_id] = size + "a";
      else
        atoms[idx].ring_rep_s[rep_id] = size;
    }
  }
}

// Build ring size_map for cod_class annotation.
// For fused aromatic rings, collect all sizes (e.g., C[5a,6a] for indole fusion).
inline void build_ring_size_map(const std::map<std::string, std::string>& ring_rep_s,
                                std::map<std::string, int>& size_map) {
  for (const auto& it : ring_rep_s) {
    const std::string& ring_str = it.second;
    size_map[ring_str] += 1;
  }
}

inline void AcedrgTables::set_atom_cod_class_name_new2(
    CodAtomInfo& atom, const CodAtomInfo& ori_atom, int lev,
    const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  if (lev == 1) {
    atom.cod_class.clear();
    atom.cod_class.append(atom.el.name());

    if (!atom.ring_rep_s.empty()) {
      std::map<std::string, int> size_map;
      build_ring_size_map(atom.ring_rep_s, size_map);

      atom.cod_class.append("[");
      int i = 0;
      int j = static_cast<int>(size_map.size());
      for (const auto& it : size_map) {
        std::string size = it.first;
        std::string num = std::to_string(it.second);
        if (it.second >= 3)
          atom.cod_class.append(num + "x" + size);
        else if (it.second == 2)
          atom.cod_class.append(size + "," + size);
        else
          atom.cod_class.append(size);
        if (i != j - 1)
          atom.cod_class.append(",");
        else
          atom.cod_class.append("]");
        ++i;
      }
    }

    std::string t_str;
    std::list<std::string> t_str_list;
    std::map<std::string, int> comps;
    for (int nb : neighbors[atom.index]) {
      if (nb == ori_atom.index)
        continue;
      std::string nb_type = atoms[nb].el.name();
      if (!atoms[nb].ring_rep_s.empty()) {
        std::map<std::string, int> size_map;
        build_ring_size_map(atoms[nb].ring_rep_s, size_map);
        nb_type.append("[");
        int i = 0;
        int j = static_cast<int>(size_map.size());
        for (const auto& it : size_map) {
          std::string size = it.first;
          std::string num = std::to_string(it.second);
          if (it.second >= 3)
            nb_type.append(num + "x" + size);
          else if (it.second == 2)
            nb_type.append(size + "," + size);
          else
            nb_type.append(size);
          if (i != j - 1)
            nb_type.append(",");
          else
            nb_type.append("]");
          ++i;
        }
      }
      comps[nb_type] += 1;
    }

    std::vector<SortMap> sorted;
    for (const auto& it : comps) {
      SortMap sm;
      sm.key = it.first;
      sm.val = it.second;
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), des_sort_map_key);
    for (const auto& sm : sorted) {
      std::string s1 = sm.key + std::to_string(sm.val);
      std::string s2;
      for (int i = 0; i < sm.val; ++i)
        s2.append(sm.key);
      if (s1.size() < s2.size())
        t_str_list.push_back(s1);
      else
        t_str_list.push_back(s2);
    }
    for (const auto& s : t_str_list)
      t_str.append(s);

    atom.cod_class.append(t_str);
  } else if (lev == 2) {
    atom.cod_class.clear();
    atom.cod_class.append(atom.el.name());

    if (!atom.ring_rep_s.empty()) {
      std::map<std::string, int> size_map;
      build_ring_size_map(atom.ring_rep_s, size_map);

      atom.cod_class.append("[");
      int i = 0;
      int j = static_cast<int>(size_map.size());
      for (const auto& it : size_map) {
        std::string size = it.first;
        std::string num = std::to_string(it.second);
        if (it.second >= 3)
          atom.cod_class.append(num + "x" + size);
        else if (it.second == 2)
          atom.cod_class.append(size + "," + size);
        else
          atom.cod_class.append(size);
        if (i != j - 1)
          atom.cod_class.append(",");
        else
          atom.cod_class.append("]");
        ++i;
      }
    }

    int low_lev = lev - 1;
    std::map<std::string, std::vector<int>> id_map;
    for (int nb : neighbors[atom.index]) {
      CodAtomInfo nb_atom = atoms[nb];
      set_atom_cod_class_name_new2(nb_atom, ori_atom, low_lev, atoms, neighbors);
      auto& entry = id_map[nb_atom.cod_class];
      if (entry.empty()) {
        entry.push_back(1);
        entry.push_back(static_cast<int>(neighbors[nb].size()));
      } else {
        entry[0] += 1;
      }
    }

    std::vector<SortMap2> sorted;
    for (const auto& it : id_map) {
      SortMap2 sm;
      sm.key = it.first;
      sm.val = it.second[0];
      sm.nNB = it.second[1];
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), des_sort_map_key2);
    for (const auto& sm : sorted) {
      if (sm.val == 1)
        atom.cod_class.append("(" + sm.key + ")");
      else
        atom.cod_class.append("(" + sm.key + ")" + std::to_string(sm.val));
    }
  }
}

inline void AcedrgTables::set_special_3nb_symb2(
    CodAtomInfo& atom, const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  if (atom.ring_rep.empty())
    return;

  std::vector<int> ser_num_nb123;
  std::map<std::string, int> nb3_props;

  for (int nb1 : neighbors[atom.index]) {
    if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb1) == ser_num_nb123.end())
      ser_num_nb123.push_back(nb1);
    for (int nb2 : neighbors[nb1]) {
      if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb2) == ser_num_nb123.end() &&
          nb2 != atom.index) {
        ser_num_nb123.push_back(nb2);
      }
    }
  }

  for (int nb1 : neighbors[atom.index]) {
    if (atoms[nb1].ring_rep.empty())
      continue;
    for (int nb2 : neighbors[nb1]) {
      if (atoms[nb2].ring_rep.empty())
        continue;
      for (int nb3 : neighbors[nb2]) {
        if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb3) == ser_num_nb123.end() &&
            nb3 != atom.index) {
          std::string prop = atoms[nb3].el.name();
          prop.append("<" + std::to_string(neighbors[nb3].size()) + ">");
          nb3_props[prop] += 1;
          ser_num_nb123.push_back(nb3);
        }
      }
    }
  }

  std::list<std::string> comps;
  for (const auto& it : nb3_props) {
    std::string id = std::to_string(it.second) + "|" + it.first;
    comps.push_back(id);
  }
  comps.sort(compare_no_case2);

  if (!comps.empty()) {
    std::string all3 = "{";
    unsigned i = 0;
    unsigned n = comps.size();
    for (const auto& id : comps) {
      if (i < n - 1)
        all3.append(id + ",");
      else
        all3.append(id);
      ++i;
    }
    all3.append("}");
    atom.cod_class.append(all3);
  }
}

inline void AcedrgTables::cod_class_to_atom2(const std::string& cod_class,
                                             CodAtomInfo& atom) const {
  std::string t_cod = trim_spaces(cod_class);
  atom.cod_class = t_cod;
  atom.nb_symb.clear();
  atom.nb2_symb.clear();
  atom.nb3_symb.clear();

  std::vector<std::string> two_parts;
  if (t_cod.find('{') != std::string::npos) {
    two_parts = split(t_cod, '{');
    if (two_parts.size() == 2) {
      atom.cod_main = two_parts[0];
      std::vector<std::string> nb3 = split(two_parts[1], '}');
      if (!nb3.empty())
        atom.nb3_symb = nb3[0];
    } else {
      atom.cod_main = t_cod;
    }
  } else {
    atom.cod_main = t_cod;
  }

  std::vector<std::string> atm_strs = split(atom.cod_main, '(');
  if (!atm_strs.empty()) {
    atom.cod_root = trim_spaces(atm_strs[0]);
  }

  std::vector<NB1stFam> all_nbs;
  for (size_t i = 1; i < atm_strs.size(); ++i) {
    std::string tS = trim_spaces(atm_strs[i]);
    std::vector<std::string> nb1 = split(tS, ')');
    NB1stFam fam;
    if (nb1.size() > 1) {
      fam.repN = str_to_int(nb1[1]);
      if (fam.repN == 0)
        fam.repN = 1;
    } else {
      fam.repN = 1;
    }

    std::string tS1 = trim_spaces(nb1[0]);
    get_small_family(tS1, fam);
    all_nbs.push_back(fam);
  }

  for (const auto& fam : all_nbs) {
    for (int j = 0; j < fam.repN; ++j) {
      std::string sN = std::to_string(static_cast<int>(fam.NB2ndList.size()) + 1);
      atom.nb_symb += fam.name + "-" + sN + ":";
      atom.nb2_symb += sN + ":";
    }
  }
}

inline void AcedrgTables::set_atoms_nb1nb2_sp(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  for (auto& atom : atoms) {
    std::vector<std::string> nb1_nb2_sp_set;
    for (int nb1 : neighbors[atom.index]) {
      std::string nb1_main = atoms[nb1].cod_root;
      std::vector<int> nb2_sp_set;
      for (int nb2 : neighbors[nb1])
        nb2_sp_set.push_back(atoms[nb2].bonding_idx);
      std::sort(nb2_sp_set.begin(), nb2_sp_set.end(), std::greater<int>());
      std::string nb2_sp_str;
      for (size_t i = 0; i < nb2_sp_set.size(); ++i) {
        nb2_sp_str.append(std::to_string(nb2_sp_set[i]));
        if (i != nb2_sp_set.size() - 1)
          nb2_sp_str.append("_");
      }
      nb1_nb2_sp_set.emplace_back(nb1_main + "-" + nb2_sp_str);
    }
    // Sort alphabetically by the string (same order as AceDRG tables)
    std::sort(nb1_nb2_sp_set.begin(), nb1_nb2_sp_set.end(),
              [](const auto& a, const auto& b) {
                return compare_no_case(a, b);
              });
    // Build nb1nb2_sp in alphabetical order
    atom.nb1nb2_sp.clear();
    for (size_t i = 0; i < nb1_nb2_sp_set.size(); ++i) {
      atom.nb1nb2_sp.append(nb1_nb2_sp_set[i]);
      if (i != nb1_nb2_sp_set.size() - 1)
        atom.nb1nb2_sp.append(":");
    }
  }
}

inline void AcedrgTables::set_atoms_bonding_and_chiral_center(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  std::map<int, std::vector<int>> num_conn_map;

  for (auto& atom : atoms) {
    int t_len = 0;
    int t_m_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
      else
        t_m_len++;
    }
    if (atom.metal_connectivity > 0)
      t_m_len = atom.metal_connectivity;

    if (atom.el == El::C) {
      if (t_len == 3 && t_m_len == 2)
        t_len = 3;
    }

    num_conn_map[atom.index].push_back(t_len);
    num_conn_map[atom.index].push_back(t_m_len);

    if (t_len > 4) {
      atom.bonding_idx = t_len;
    } else if (atom.el == El::C || atom.el == El::Si || atom.el == El::Ge) {
      if (t_len == 4) {
        atom.bonding_idx = 3;
      } else if (t_len == 3) {
        atom.bonding_idx = 2;
      } else if (t_len == 2) {
        if (get_num_oxy_connect(atoms, atom, neighbors) == 1)
          atom.bonding_idx = 1;
        else if (t_m_len == 1 || atom.charge == -1.0f)
          atom.bonding_idx = 2;
        else if (t_m_len == 1 || atom.charge == -2.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 1;
      }
    } else if (atom.el == El::N || atom.el == El::As) {
      if (t_len == 4 || t_len == 3) {
        atom.bonding_idx = 3;
      } else if (t_len == 2) {
        if (atom.charge == 1.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 2;
      } else if (t_len == 1) {
        atom.bonding_idx = 1;
      }
    } else if (atom.el == El::B) {
      if (t_len == 4) {
        atom.bonding_idx = 3;
      } else if (t_len == 3) {
        atom.bonding_idx = 2;
      } else if (t_len == 2) {
        if (atom.charge == 1.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 2;
      } else if (t_len == 1) {
        atom.bonding_idx = 1;
      }
    } else if (atom.el == El::O) {
      if (neighbors[atom.index].size() == 2) {
        atom.bonding_idx = 3;
      } else if (neighbors[atom.index].size() == 1) {
        if (neighbors[neighbors[atom.index][0]].size() != 1)
          atom.bonding_idx = 2;
        else
          atom.bonding_idx = 3;
      } else {
        atom.bonding_idx = 3;
      }
    } else if (atom.el == El::P) {
      if (t_len == 4 || t_len == 3 || t_len == 2 || t_len == 5)
        atom.bonding_idx = 3;
    } else if (atom.el == El::S) {
      if (t_len == 2 || t_len == 3 || t_len == 4)
        atom.bonding_idx = 3;
      else if (t_len == 6)
        atom.bonding_idx = 5;
      else if (t_len == 1)
        atom.bonding_idx = 3;
    } else if (atom.el == El::Se) {
      if (t_len == 4 || t_len == 3 || t_len == 2) {
        atom.bonding_idx = (t_len == 3) ? 2 : 3;
      } else if (t_len == 6) {
        atom.bonding_idx = 5;
      } else if (t_len == 1) {
        atom.bonding_idx = 3;
      }
    } else if (atom.el == El::Br) {
      if (t_len == 3)
        atom.bonding_idx = 3;
    }
  }

  for (auto& atom : atoms) {
    int t_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
    }

    if (atom.el == El::O) {
      if (t_len == 2 && atom.par_charge == 0.0f) {
        bool l_sp2 = false;
        for (int nb : neighbors[atom.index]) {
          if (atoms[nb].bonding_idx == 2) {
            l_sp2 = true;
            break;
          }
        }
        if (l_sp2)
          atom.bonding_idx = 2;
      }
    }
  }

  std::map<int, int> pre_bonding;
  for (const auto& atom : atoms)
    pre_bonding[atom.index] = atom.bonding_idx;

  for (auto& atom : atoms) {
    int t_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
    }
    if (atom.el == El::N || atom.el == El::As) {
      if (t_len == 3) {
        if (atom.charge == 0.0f) {
          bool l_sp2 = false;
          for (int nb : neighbors[atom.index]) {
            if (pre_bonding[nb] == 2 && atoms[nb].el != El::O) {
              l_sp2 = true;
              break;
            }
          }
          if (l_sp2) {
            if (num_conn_map[atom.index][1] != 0)
              atom.bonding_idx = 3;
            else
              atom.bonding_idx = 2;
          } else {
            atom.bonding_idx = 3;
          }
        } else if (atom.charge == 1.0f) {
          atom.bonding_idx = 2;
        }
      }
    }
    if (atom.el == El::S) {
      if (t_len == 2 && atom.charge == 0.0f) {
        bool l_sp2 = false;
        for (int nb : neighbors[atom.index]) {
          if (pre_bonding[nb] == 2 && atoms[nb].el != El::O) {
            l_sp2 = true;
            break;
          }
        }
        if (l_sp2)
          atom.bonding_idx = 2;
      }
    }
    if (atom.el == El::C) {
      if (t_len == 3 && atom.charge == -1.0f) {
        std::vector<int> sp2_set;
        for (int nb : neighbors[atom.index]) {
          if (atoms[nb].bonding_idx == 2)
            sp2_set.push_back(nb);
        }
        if (sp2_set.size() == 2)
          atom.bonding_idx = 2;
        else
          atom.bonding_idx = 3;
      }
    }
  }
}

inline int AcedrgTables::get_num_oxy_connect(const std::vector<CodAtomInfo>& atoms,
                                             const CodAtomInfo& atom,
                                             const std::vector<std::vector<int>>& neighbors) const {
  int nO = 0;
  for (int nb : neighbors[atom.index])
    if (atoms[nb].el == El::O)
      nO++;
  return nO;
}

inline int AcedrgTables::get_min_ring2_from_cod_class(const std::string& cod_class) const {
  int r_size = 0;
  if (!cod_class.empty()) {
    std::vector<std::string> tmp1 = split(cod_class, '(');
    if (!tmp1.empty()) {
      if (tmp1[0].find('[') != std::string::npos) {
        std::vector<std::string> tmp2 = split(tmp1[0], '[');
        if (tmp2.size() > 1) {
          if (tmp2[1].find(',') != std::string::npos) {
            std::vector<std::string> tmp3 = split(tmp2[1], ',');
            if (!tmp3.empty()) {
              if (tmp3[0].find('x') != std::string::npos) {
                std::vector<std::string> tmp4 = split(tmp3[0], 'x');
                if (tmp4.size() > 1)
                  r_size = str_to_int(tmp4[1]);
              } else {
                r_size = str_to_int(tmp3[0]);
              }
            }
          } else {
            std::vector<std::string> tmp3 = split(tmp2[1], ']');
            if (!tmp3.empty()) {
              if (tmp3[0].find('x') != std::string::npos) {
                std::vector<std::string> tmp4 = split(tmp3[0], 'x');
                if (tmp4.size() > 1)
                  r_size = str_to_int(tmp4[1]);
              } else {
                r_size = str_to_int(tmp3[0]);
              }
            }
          }
        }
      }
    }
  }
  return r_size;
}

inline bool AcedrgTables::cod_class_is_aromatic(const std::string& cod_class) const {
  if (cod_class.empty())
    return false;
  std::vector<std::string> secs = split(cod_class, '(');
  if (secs.empty())
    return false;
  if (secs[0].find('[') != std::string::npos) {
    std::vector<std::string> rs = split(secs[0], '[');
    if (rs.size() > 1)
      return rs[1].find('a') != std::string::npos;
  }
  return false;
}

inline void AcedrgTables::get_small_family(const std::string& in_str, NB1stFam& fam) const {
  fam.name.clear();
  fam.NB2ndList.clear();
  std::vector<std::string> ch_list = {",", "x"};
  std::string name_str;
  bool l_r = false;
  int n_rep = 1;
  for (size_t i = 0; i < in_str.size(); ++i) {
    char c = in_str[i];
    if (std::isalpha(static_cast<unsigned char>(c))) {
      if (std::toupper(c) == c) {
        if (!name_str.empty()) {
          if (fam.name.empty()) {
            fam.name = name_str;
            n_rep = 1;
          } else {
            for (int l = 0; l < n_rep; ++l)
              fam.NB2ndList.push_back(name_str);
            n_rep = 1;
          }
        }
        name_str = c;
      } else {
        name_str += c;
      }
    } else if (c == '[') {
      name_str += c;
      l_r = true;
    } else if (c == ']') {
      name_str += c;
      l_r = false;
    } else if (std::find(ch_list.begin(), ch_list.end(), std::string(1, c)) != ch_list.end()) {
      name_str += c;
    } else if (std::isdigit(static_cast<unsigned char>(c))) {
      if (l_r) {
        name_str += c;
      } else {
        n_rep = std::stoi(in_str.substr(i, 1));
      }
    }
  }

  if (!name_str.empty()) {
    if (fam.name.empty()) {
      fam.name = name_str;
    } else {
      for (int l = 0; l < n_rep; ++l)
        fam.NB2ndList.push_back(name_str);
    }
  }
}

inline std::string AcedrgTables::trim_spaces(const std::string& s) {
  size_t start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos)
    return "";
  size_t end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

inline std::vector<std::string> AcedrgTables::split(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == delim) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

inline int AcedrgTables::str_to_int(const std::string& s) {
  std::istringstream iss(s);
  int value = 0;
  iss >> value;
  return value;
}

inline bool AcedrgTables::compare_no_case(const std::string& first,
                                          const std::string& second) {
  size_t i = 0;
  while (i < first.length() && i < second.length()) {
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(first[i])));
    char b = static_cast<char>(std::toupper(static_cast<unsigned char>(second[i])));
    if (a < b)
      return true;
    if (a > b)
      return false;
    ++i;
  }
  return first.length() > second.length();
}

inline bool AcedrgTables::compare_no_case2(const std::string& first,
                                           const std::string& second) {
  if (first.length() > second.length())
    return true;
  if (first.length() < second.length())
    return false;
  for (size_t i = 0; i < first.length() && i < second.length(); ++i) {
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(first[i])));
    char b = static_cast<char>(std::toupper(static_cast<unsigned char>(second[i])));
    if (a < b)
      return true;
    if (a > b)
      return false;
  }
  return true;
}

inline bool AcedrgTables::des_sort_map_key(const SortMap& a, const SortMap& b) {
  return a.key.length() > b.key.length();
}

inline bool AcedrgTables::des_sort_map_key2(const SortMap2& a, const SortMap2& b) {
  if (a.key.length() > b.key.length())
    return true;
  if (a.key.length() == b.key.length())
    return a.nNB > b.nNB;
  return false;
}

inline bool AcedrgTables::are_in_same_ring(const CodAtomInfo& a1,
                                           const CodAtomInfo& a2) const {
  for (int ring_idx : a1.in_rings)
    if (std::find(a2.in_rings.begin(), a2.in_rings.end(), ring_idx) != a2.in_rings.end())
      return true;
  return false;
}

inline int AcedrgTables::angle_ring_size(const CodAtomInfo& center,
                                         const CodAtomInfo& a1,
                                         const CodAtomInfo& a3) const {
  for (const auto& it : center.ring_rep) {
    if (a1.ring_rep.find(it.first) != a1.ring_rep.end() &&
        a3.ring_rep.find(it.first) != a3.ring_rep.end())
      return it.second;
  }
  return 0;
}

inline void AcedrgTables::order_bond_atoms(const CodAtomInfo& a1, const CodAtomInfo& a2,
                                           const CodAtomInfo*& first,
                                           const CodAtomInfo*& second) const {
  if (a1.hashing_value < a2.hashing_value) {
    first = &a1;
    second = &a2;
    return;
  }
  if (a1.hashing_value > a2.hashing_value) {
    first = &a2;
    second = &a1;
    return;
  }
  if (compare_no_case2(a1.cod_main, a2.cod_main)) {
    first = &a1;
    second = &a2;
  } else {
    first = &a2;
    second = &a1;
  }
}

inline void AcedrgTables::order_angle_flanks(const CodAtomInfo& a1, const CodAtomInfo& a3,
                                             const CodAtomInfo*& flank1,
                                             const CodAtomInfo*& flank3) const {
  if (a1.hashing_value < a3.hashing_value) {
    flank1 = &a1;
    flank3 = &a3;
    return;
  }
  if (a1.hashing_value > a3.hashing_value) {
    flank1 = &a3;
    flank3 = &a1;
    return;
  }
  if (a1.cod_main.size() > a3.cod_main.size()) {
    flank1 = &a1;
    flank3 = &a3;
    return;
  }
  if (a1.cod_main.size() < a3.cod_main.size()) {
    flank1 = &a3;
    flank3 = &a1;
    return;
  }
  if (compare_no_case2(a1.cod_main, a3.cod_main)) {
    flank1 = &a1;
    flank3 = &a3;
  } else {
    flank1 = &a3;
    flank3 = &a1;
  }
}

inline Hybridization AcedrgTables::hybrid_from_bonding_idx(int bonding_idx,
                                                           bool is_metal_atom,
                                                           int connectivity) const {
  if (is_metal_atom) {
    if (connectivity <= 4) return Hybridization::SPD5;
    if (connectivity == 5) return Hybridization::SPD5;
    if (connectivity == 6) return Hybridization::SPD6;
    if (connectivity == 7) return Hybridization::SPD7;
    return Hybridization::SPD8;
  }
  if (bonding_idx <= 0)
    return Hybridization::SP_NON;
  if (bonding_idx == 1)
    return Hybridization::SP1;
  if (bonding_idx == 2)
    return Hybridization::SP2;
  if (bonding_idx == 3)
    return Hybridization::SP3;
  if (bonding_idx == 5)
    return Hybridization::SPD5;
  if (bonding_idx == 6)
    return Hybridization::SPD6;
  if (bonding_idx == 7)
    return Hybridization::SPD7;
  if (bonding_idx >= 8)
    return Hybridization::SPD8;
  return Hybridization::SP_NON;
}

// ============================================================================
// CCP4 atom type assignment (AceDRG)
// ============================================================================

inline int AcedrgTables::ccp4_material_type(Element el) {
  switch (el.elem) {
    case El::H: case El::D:
      return 1;
    case El::C: case El::N: case El::O: case El::P: case El::S: case El::Se:
      return 2;
    case El::Li: case El::Na: case El::K: case El::Rb: case El::Cs: case El::Fr:
      return 3;
    case El::Be: case El::Mg: case El::Ca: case El::Sr: case El::Ba: case El::Ra:
      return 4;
    case El::Sc: case El::Y: case El::Ti: case El::Zr: case El::Hf: case El::Rf:
    case El::V: case El::Nb: case El::Ta: case El::Db:
    case El::Cr: case El::Mo: case El::W: case El::Sg:
    case El::Mn: case El::Tc: case El::Re: case El::Bh:
    case El::Fe: case El::Ru: case El::Os: case El::Hs:
    case El::Co: case El::Rh: case El::Ir: case El::Mt:
    case El::Ni: case El::Pd: case El::Pt: case El::Ds:
    case El::Cu: case El::Ag: case El::Au: case El::Rg:
    case El::Zn: case El::Cd: case El::Hg: case El::Cn:
      return 5;
    case El::Al: case El::Ga: case El::In: case El::Tl: case El::Sn:
    case El::Pb: case El::Bi:
      return 6;
    case El::B: case El::Si: case El::Ge: case El::As: case El::Sb:
    case El::Te: case El::Po:
      return 7;
    case El::F: case El::Cl: case El::Br: case El::I: case El::At:
      return 8;
    case El::La: case El::Ce: case El::Pr: case El::Nd: case El::Pm:
    case El::Sm: case El::Eu: case El::Gd: case El::Tb: case El::Dy:
    case El::Ho: case El::Er: case El::Tm: case El::Yb: case El::Lu:
    case El::Ac: case El::Th: case El::Pa: case El::U: case El::Np:
    case El::Pu: case El::Am: case El::Cm: case El::Bk: case El::Cf:
    case El::Es: case El::Fm: case El::Md: case El::No: case El::Lr:
      return 9;
    case El::He: case El::Ne: case El::Ar: case El::Kr: case El::Xe: case El::Rn:
      return 10;
    default:
      return 0;
  }
}

inline std::string AcedrgTables::bond_order_key(BondType type) {
  std::string s = bond_type_to_string(type);
  if (s.empty())
    s = "single";
  s = to_upper(s);
  if (s.size() > 4)
    s.resize(4);
  return s;
}

inline void AcedrgTables::load_ccp4_bonds(const std::string& path) {
  try {
    cif::Document doc = read_cif_gz(path);
    if (doc.blocks.empty())
      return;
    for (const auto& row : doc.blocks[0].find("_lib_bond.",
                    {"atom_type_1", "atom_type_2", "type",
                     "length", "value_esd"})) {
      std::string type1 = row.str(0);
      std::string type2 = row.str(1);
      std::string order = row.str(2);
      double length = cif::as_number(row[3], NAN);
      double sigma = cif::as_number(row[4], NAN);
      if (type1.empty() || order.empty() || std::isnan(length))
        continue;
      std::string order_key = to_upper(order);
      if (order_key.size() > 4)
        order_key.resize(4);
      ccp4_bonds_[type1][type2][order_key] = {length, sigma};
    }
  } catch (std::exception&) {
    return;
  }
}

inline bool AcedrgTables::search_ccp4_bond(const std::string& type1,
                                           const std::string& type2,
                                           const std::string& order,
                                           ValueStats& out) const {
  auto it1 = ccp4_bonds_.find(type1);
  if (it1 == ccp4_bonds_.end())
    return false;
  auto it2 = it1->second.find(type2);
  if (it2 == it1->second.end())
    return false;
  auto it3 = it2->second.find(order);
  if (it3 == it2->second.end())
    return false;
  out.value = it3->second.length;
  out.sigma = it3->second.sigma;
  out.count = 1;
  return true;
}

inline std::vector<std::string> AcedrgTables::compute_ccp4_types(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int>>& neighbors) const {
  std::vector<Ccp4AtomInfo> atoms;
  atoms.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    Ccp4AtomInfo info;
    info.el = cc.atoms[i].el;
    info.chem_type = cc.atoms[i].el.name();
    info.ccp4_type = info.chem_type;
    info.bonding_idx = atom_info[i].bonding_idx;
    info.ring_rep = atom_info[i].ring_rep;
    info.conn_atoms = neighbors[i];
    info.conn_atoms_no_metal.clear();
    for (int nb : neighbors[i]) {
      if (cc.atoms[nb].is_hydrogen())
        info.conn_h_atoms.push_back(nb);
      if (!cc.atoms[nb].el.is_metal())
        info.conn_atoms_no_metal.push_back(nb);
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(cc.atoms[i].charge));
    atoms.emplace_back(std::move(info));
  }

  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);

  std::vector<std::string> out;
  out.reserve(atoms.size());
  for (const auto& atom : atoms)
    out.push_back(atom.ccp4_type);
  return out;
}

inline void AcedrgTables::set_one_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                            size_t idx) {
  int ntype = ccp4_material_type(atoms[idx].el);
  switch (ntype) {
    case 1:
      set_hydro_ccp4_type(atoms, idx);
      break;
    case 2:
      set_org_ccp4_type(atoms, idx);
      break;
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
      atoms[idx].ccp4_type = atoms[idx].chem_type;
      break;
    default:
      atoms[idx].ccp4_type = atoms[idx].chem_type;
      break;
  }
  atoms[idx].ccp4_type = to_upper(atoms[idx].ccp4_type);
}

inline void AcedrgTables::set_hydro_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                              size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  atom.ccp4_type = "H";
  if (atom.conn_atoms.size() == 1) {
    int nb = atom.conn_atoms[0];
    if (atoms[nb].chem_type == "S")
      atom.ccp4_type = "HSH1";
  }
}

inline void AcedrgTables::set_org_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                            size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  int r5 = 0;
  int r6 = 0;
  for (const auto& item : atom.ring_rep) {
    if (item.second == 5)
      r5 += 1;
    if (item.second == 6)
      r6 += 1;
  }

  const size_t nconn = atom.conn_atoms_no_metal.size();
  const size_t nh = atom.conn_h_atoms.size();

  if (atom.chem_type == "C") {
    if (atom.bonding_idx == 2) {
      if (r5 && r6) {
        atom.ccp4_type = "CR56";
      } else if (r5 == 2) {
        atom.ccp4_type = "CR55";
      } else if (r6 == 2) {
        atom.ccp4_type = "CR66";
      } else if (r5 == 1) {
        if (nh == 1) {
          atom.ccp4_type = "CR15";
        } else if (nh == 0) {
          atom.ccp4_type = "CR5";
        }
      } else if (r6 == 1) {
        if (nh == 1) {
          atom.ccp4_type = "CR16";
        } else if (nh == 0) {
          atom.ccp4_type = "CR6";
        }
      } else {
        if (nh == 1) {
          atom.ccp4_type = "C1";
        } else if (nh == 2) {
          atom.ccp4_type = "C2";
        } else if (nh == 0) {
          atom.ccp4_type = "C";
        }
      }
    } else if (atom.bonding_idx == 3) {
      if (nh == 0) {
        atom.ccp4_type = "CT";
      } else if (nh == 1) {
        atom.ccp4_type = "CH1";
      } else if (nh == 2) {
        atom.ccp4_type = "CH2";
      } else if (nh == 3) {
        atom.ccp4_type = "CH3";
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "CSP";
    }
  } else if (atom.chem_type == "N") {
    if (atom.bonding_idx == 2) {
      if (nconn == 3) {
        if (nh == 1) {
          atom.ccp4_type = "NH1";
        } else if (nh == 2) {
          atom.ccp4_type = "NH2";
        } else if (nh == 0) {
          atom.ccp4_type = "NH0";
        } else {
          atom.ccp4_type = "N";
        }
      } else if (nconn == 2) {
        if (nh == 1) {
          atom.ccp4_type = "N21";
        } else if (nh == 0) {
          atom.ccp4_type = "N20";
        } else {
          atom.ccp4_type = "N";
        }
      }
    } else if (atom.bonding_idx == 3) {
      if (nconn == 4) {
        if (nh == 1) {
          atom.ccp4_type = "NT1";
        } else if (nh == 2) {
          atom.ccp4_type = "NT2";
        } else if (nh == 3) {
          atom.ccp4_type = "NT3";
        } else if (nh == 4) {
          atom.ccp4_type = "NT4";
        } else if (nh == 0) {
          atom.ccp4_type = "NT";
        } else {
          atom.ccp4_type = "N";
        }
      } else if (nconn == 3) {
        if (nh == 1) {
          atom.ccp4_type = "N31";
        } else if (nh == 2) {
          atom.ccp4_type = "N32";
        } else if (nh == 3) {
          atom.ccp4_type = "N33";
        } else if (nh == 0) {
          atom.ccp4_type = "N30";
        } else {
          atom.ccp4_type = "N3";
        }
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "NSP";
    }
  } else if (atom.chem_type == "P") {
    atom.ccp4_type = (nconn == 4 ? "P" : "P1");
  } else if (atom.chem_type == "O") {
    bool lP = false;
    bool lS = false;
    bool lB = false;
    for (int nb : atom.conn_atoms) {
      if (atoms[nb].chem_type == "P") lP = true;
      if (atoms[nb].chem_type == "S") lS = true;
      if (atoms[nb].chem_type == "B") lB = true;
    }

    bool has_par_charge = std::fabs(atom.par_charge) > 1e-6f;
    bool has_formal_charge = atom.formal_charge != 0;

    if (atom.bonding_idx == 2) {
      if (has_par_charge) {
        if (lP) {
          atom.ccp4_type = "OP";
        } else if (lS) {
          atom.ccp4_type = "OS";
        } else if (lB) {
          atom.ccp4_type = "OB";
        } else {
          atom.ccp4_type = "OC";
        }
      } else if (nconn == 2) {
        if (nh == 1) {
          atom.ccp4_type = "OH1";
        } else if (nh == 2) {
          atom.ccp4_type = "OH2";
        } else {
          atom.ccp4_type = "O";
        }
      } else {
        if (has_formal_charge) {
          if (lP) {
            atom.ccp4_type = "OP";
          } else if (lS) {
            atom.ccp4_type = "OS";
          } else if (lB) {
            atom.ccp4_type = "OB";
          } else {
            atom.ccp4_type = "OC";
          }
        } else {
          atom.ccp4_type = "O";
        }
      }
    } else if (atom.bonding_idx == 3) {
      bool lC = false;
      for (int nb : atom.conn_atoms)
        if (atoms[nb].chem_type == "C")
          lC = true;
      if (lC && nh == 1 && nconn == 2) {
        atom.ccp4_type = "OH1";
      } else if (nh == 2) {
        atom.ccp4_type = "OH2";
      } else if (nconn == 2) {
        if (has_par_charge) {
          atom.ccp4_type = "OC2";
        } else if (nh == 1) {
          atom.ccp4_type = "OH1";
        } else {
          atom.ccp4_type = "O2";
        }
      }
    } else if (nconn == 1) {
      if (has_formal_charge) {
        if (lP) {
          atom.ccp4_type = "OP";
        } else if (lS) {
          atom.ccp4_type = "OS";
        } else if (lB) {
          atom.ccp4_type = "OB";
        } else {
          atom.ccp4_type = "OC";
        }
      } else {
        atom.ccp4_type = "O";
      }
    } else {
      atom.ccp4_type = "O";
    }
  } else if (atom.chem_type == "S") {
    if (nconn == 3 || nconn == 4) {
      if (nh == 0) {
        atom.ccp4_type = "S3";
      } else if (nh == 1) {
        atom.ccp4_type = "SH1";
      }
    } else if (nconn == 2) {
      if (nh == 0) {
        atom.ccp4_type = "S2";
      } else {
        atom.ccp4_type = "SH1";
      }
    } else if (nconn == 1) {
      atom.ccp4_type = "S1";
    } else if (nh == 0) {
      atom.ccp4_type = "S";
    } else if (nh == 1) {
      atom.ccp4_type = "SH1";
    } else {
      atom.ccp4_type = "S";
    }
  } else if (atom.chem_type == "Se") {
    // Selenium - similar to sulfur
    if (nconn == 3 || nconn == 4) {
      atom.ccp4_type = "SE";
    } else if (nconn == 2) {
      atom.ccp4_type = "SE";
    } else if (nconn == 1) {
      atom.ccp4_type = "SE";
    } else {
      atom.ccp4_type = "SE";
    }
  } else {
    atom.ccp4_type = atom.chem_type;
  }
}

inline void AcedrgTables::assign_ccp4_types(ChemComp& cc) const {
  if (cc.atoms.empty())
    return;

  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);
  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);

  std::vector<Ccp4AtomInfo> atoms;
  atoms.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    Ccp4AtomInfo info;
    info.el = cc.atoms[i].el;
    info.chem_type = cc.atoms[i].el.name();
    info.ccp4_type = info.chem_type;
    info.bonding_idx = atom_info[i].bonding_idx;
    info.ring_rep = atom_info[i].ring_rep;
    info.conn_atoms = neighbors[i];
    info.conn_atoms_no_metal.clear();
    for (int nb : neighbors[i]) {
      if (cc.atoms[nb].is_hydrogen())
        info.conn_h_atoms.push_back(nb);
      if (!cc.atoms[nb].el.is_metal())
        info.conn_atoms_no_metal.push_back(nb);
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(cc.atoms[i].charge));
    atoms.emplace_back(std::move(info));
  }

  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);

  for (size_t i = 0; i < atoms.size(); ++i)
    cc.atoms[i].chem_type = atoms[i].ccp4_type;
}

// ============================================================================
// Implementation - Bond search
// ============================================================================

inline void AcedrgTables::fill_restraints(ChemComp& cc) const {
  if (!tables_loaded_)
    return;

  // Classify atoms
  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);
  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);
  std::vector<std::string> ccp4_types;
  if (!ccp4_bonds_.empty())
    ccp4_types = compute_ccp4_types(cc, atom_info, neighbors);

  // Print atom classification if verbose (level 1 = bonds only, level 2+ = atoms too)
  if (verbose >= 1) {
    std::fprintf(stderr, "  Atom classification:\n");
    for (size_t i = 0; i < atom_info.size(); ++i) {
      const auto& a = atom_info[i];
      std::fprintf(stderr, "    %s: el=%s conn=%d ring=%d arom=%d hybr=%s "
                   "bonding_idx=%d hash=%d cod_class=%s\n",
                   a.id.c_str(), element_name(a.el), a.connectivity,
                   a.min_ring_size, a.is_aromatic ? 1 : 0,
                   hybridization_to_string(a.hybrid), a.bonding_idx,
                   a.hashing_value, a.cod_class.c_str());
      if (verbose >= 2)
        std::fprintf(stderr, "        cod_main=%s nb1nb2_sp=%s nb=%s nb2=%s\n",
                     a.cod_main.c_str(), a.nb1nb2_sp.c_str(),
                     a.nb_symb.c_str(), a.nb2_symb.c_str());
    }
  }

  // Fill bonds
  for (auto& bond : cc.rt.bonds) {
    if (std::isnan(bond.value)) {
      fill_bond(cc, atom_info, bond);
      // CCP4 energetic library fallback
      if (std::isnan(bond.value) && !ccp4_types.empty()) {
        auto it1 = cc.find_atom(bond.id1.atom);
        auto it2 = cc.find_atom(bond.id2.atom);
        if (it1 != cc.atoms.end() && it2 != cc.atoms.end()) {
          int idx1 = static_cast<int>(it1 - cc.atoms.begin());
          int idx2 = static_cast<int>(it2 - cc.atoms.begin());
          std::string order = bond_order_key(bond.type);
          ValueStats vs;
          if (!search_ccp4_bond(ccp4_types[idx1], ccp4_types[idx2], order, vs) &&
              !search_ccp4_bond(ccp4_types[idx2], ccp4_types[idx1], order, vs)) {
            ValueStats v1, v2;
            bool has1 = search_ccp4_bond(ccp4_types[idx1], ".", order, v1);
            bool has2 = search_ccp4_bond(ccp4_types[idx2], ".", order, v2);
            if (has1 && has2) {
              bond.value = 0.5 * (v1.value + v2.value);
              bond.esd = 0.02;
            }
          } else {
            bond.value = vs.value;
            bond.esd = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
          }
        }
      }
    }
  }

  // Fill angles
  for (auto& angle : cc.rt.angles) {
    if (std::isnan(angle.value)) {
      fill_angle(cc, atom_info, angle);
    }
  }

  // AceDRG adjustment: enforce 360-degree sum for sp2 centers with 3 angles.
  std::map<std::string, std::vector<size_t>> center_angles;
  for (size_t i = 0; i < cc.rt.angles.size(); ++i)
    center_angles[cc.rt.angles[i].id2.atom].push_back(i);
  for (const auto& entry : center_angles) {
    if (entry.second.size() != 3)
      continue;
    auto it = cc.find_atom(entry.first);
    if (it == cc.atoms.end())
      continue;
    int idx = static_cast<int>(it - cc.atoms.begin());
    if (atom_info[idx].hybrid != Hybridization::SP2)
      continue;
    double sum = 0.0;
    for (size_t idx_ang : entry.second)
      sum += cc.rt.angles[idx_ang].value;
    double diff = (360.0 - sum) / 3.0;
    if (std::fabs(diff) > 0.01) {
      double new_sum = 0.0;
      for (size_t idx_ang : entry.second) {
        cc.rt.angles[idx_ang].value += diff;
        new_sum += cc.rt.angles[idx_ang].value;
      }
      cc.rt.angles[entry.second[0]].value += (360.0 - new_sum);
    } else {
      cc.rt.angles[entry.second[0]].value += diff;
    }
  }
}

inline void AcedrgTables::fill_bond(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Bond& bond) const {

  auto it1 = cc.find_atom(bond.id1.atom);
  auto it2 = cc.find_atom(bond.id2.atom);
  if (it1 == cc.atoms.end() || it2 == cc.atoms.end())
    return;

  int idx1 = static_cast<int>(it1 - cc.atoms.begin());
  int idx2 = static_cast<int>(it2 - cc.atoms.begin());

  const CodAtomInfo& a1 = atom_info[idx1];
  const CodAtomInfo& a2 = atom_info[idx2];

  const char* source = "none";

  // Check for metal bond
  if (a1.is_metal || a2.is_metal) {
    const CodAtomInfo& metal = a1.is_metal ? a1 : a2;
    const CodAtomInfo& ligand = a1.is_metal ? a2 : a1;

    ValueStats vs = search_metal_bond(metal, ligand, atom_info);
    if (vs.count > 0) {
      bond.value = vs.value;
      bond.esd = clamp_bond_sigma(vs.sigma);
      source = "metal";
      if (verbose)
        std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                     bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                     a1.hashing_value, a2.hashing_value,
                     hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                     source, bond.value, bond.esd, vs.count);
      return;
    }
  }

  // Try detailed multilevel search first (if tables loaded)
  ValueStats vs = search_bond_multilevel(a1, a2);
  if (vs.count >= min_observations_bond) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    source = "multilevel";
    if (verbose)
      std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                   bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                   a1.hashing_value, a2.hashing_value,
                   hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                   source, bond.value, bond.esd, vs.count);
    return;
  }

  // Try HRS (summary) table
  vs = search_bond_hrs(a1, a2, are_in_same_ring(a1, a2));
  if (vs.count >= min_observations_bond) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    source = "HRS";
    if (verbose)
      std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                   bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                   a1.hashing_value, a2.hashing_value,
                   hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                   source, bond.value, bond.esd, vs.count);
    return;
  }

  // Try element+hybridization fallback
  vs = search_bond_en(a1, a2);
  if (vs.count > 0) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    source = "EN";
    if (verbose)
      std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                   bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                   a1.hashing_value, a2.hashing_value,
                   hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                   source, bond.value, bond.esd, vs.count);
    return;
  }

  // Ultimate fallback: sum of covalent radii
  bond.value = a1.el.covalent_r() + a2.el.covalent_r();
  bond.esd = upper_bond_sigma;
  source = "covalent_radii";
  if (verbose)
    std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f)\n",
                 bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                 a1.hashing_value, a2.hashing_value,
                 hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                 source, bond.value, bond.esd);
}

inline ValueStats AcedrgTables::search_bond_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& a2) const {

  if (verbose >= 2)
    std::fprintf(stderr, "      search_bond_multilevel: input a1=%s(hash=%d) a2=%s(hash=%d)\n",
                 a1.id.c_str(), a1.hashing_value, a2.id.c_str(), a2.hashing_value);

  const CodAtomInfo* left = nullptr;
  const CodAtomInfo* right = nullptr;
  order_bond_atoms(a1, a2, left, right);

  if (verbose >= 2)
    std::fprintf(stderr, "      after order: left=%s(hash=%d) right=%s(hash=%d)\n",
                 left->id.c_str(), left->hashing_value, right->id.c_str(), right->hashing_value);

  // Build lookup keys
  int ha1 = left->hashing_value;
  int ha2 = right->hashing_value;

  std::string h1 = hybridization_to_string(left->hybrid);
  std::string h2 = hybridization_to_string(right->hybrid);
  if (h1 > h2) std::swap(h1, h2);
  std::string hybr_comb = h1 + "_" + h2;

  std::string in_ring = are_in_same_ring(a1, a2) ? "Y" : "N";

  // Use neighbor symbols as additional keys
  const std::string& a1_nb2 = left->nb2_symb;
  const std::string& a2_nb2 = right->nb2_symb;
  const std::string& a1_nb = left->nb1nb2_sp;
  const std::string& a2_nb = right->nb1nb2_sp;

  // Use COD main type as atom type (for 1D lookup)
  const std::string& a1_type = left->cod_main;
  const std::string& a2_type = right->cod_main;

  if (verbose >= 2)
    std::fprintf(stderr, "      lookup: hash=%d/%d hybr=%s ring=%s nb1nb2_sp=%s/%s nb2=%s/%s type=%s/%s\n",
                 ha1, ha2, hybr_comb.c_str(), in_ring.c_str(),
                 a1_nb.c_str(), a2_nb.c_str(), a1_nb2.c_str(), a2_nb2.c_str(),
                 a1_type.c_str(), a2_type.c_str());

  // Exact match with full COD class (a1C/a2C) before any aggregation.
  {
    auto it1 = bond_idx_full_.find(ha1);
    if (it1 != bond_idx_full_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(hybr_comb);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(in_ring);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_nb2);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_nb2);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a1_nb);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a2_nb);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a1_type);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a2_type);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(left->cod_class);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(right->cod_class);
                          if (it12 != it11->second.end()) {
                            const ValueStats& vs = it12->second;
                            if (vs.count >= min_observations_bond) {
                              if (verbose >= 2)
                                std::fprintf(stderr, "      matched: level=full\n");
                              return vs;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  const auto* map_1d = [&]() -> const std::map<std::string, std::map<std::string, std::vector<ValueStats>>>* {
    auto it1 = bond_idx_1d_.find(ha1);
    if (it1 == bond_idx_1d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    auto it5 = it4->second.find(a1_nb2);
    if (it5 == it4->second.end()) return nullptr;
    auto it6 = it5->second.find(a2_nb2);
    if (it6 == it5->second.end()) return nullptr;
    auto it7 = it6->second.find(a1_nb);
    if (it7 == it6->second.end()) return nullptr;
    auto it8 = it7->second.find(a2_nb);
    if (it8 == it7->second.end()) return nullptr;
    return &it8->second;
  }();

  const auto* map_2d = [&]() -> const std::map<std::string, std::map<std::string, std::vector<ValueStats>>>* {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 == bond_idx_2d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    auto it5 = it4->second.find(a1_nb2);
    if (it5 == it4->second.end()) return nullptr;
    auto it6 = it5->second.find(a2_nb2);
    if (it6 == it5->second.end()) return nullptr;
    return &it6->second;
  }();

  const auto* map_nb2 = [&]() -> const std::map<std::string, std::map<std::string,
                                      std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>>* {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 == bond_idx_2d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    return &it4->second;
  }();

  // AceDRG-like multilevel fallback ordering (levels 0..11).
  for (int level = 0; level < 12; ++level) {
    ValueStats vs;
    if (level == 0) {
      if (map_1d) {
        auto it1 = map_1d->find(a1_type);
        if (it1 != map_1d->end()) {
          auto it2 = it1->second.find(a2_type);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            vs = it2->second.front();
          }
        }
      }
    } else if (level == 1 || level == 2) {
      if (map_1d) {
        std::vector<ValueStats> values;
        for (const auto& it1 : *map_1d) {
          if (level == 1 && it1.first == a1_type) {
            for (const auto& it2 : it1.second) {
              if (!it2.second.empty())
                values.push_back(it2.second.front());
            }
          } else if (it1.first != a1_type) {
            auto it2 = it1.second.find(a2_type);
            if (it2 != it1.second.end() && !it2->second.empty())
              values.push_back(it2->second.front());
          }
        }
        if (static_cast<int>(values.size()) >= min_observations_bond)
          vs = aggregate_stats(values);
      }
    } else if (level == 3) {
      if (map_2d) {
        auto it1 = map_2d->find(a1_nb);
        if (it1 != map_2d->end()) {
          auto it2 = it1->second.find(a2_nb);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            if (static_cast<int>(it2->second.size()) >= min_observations_bond)
              vs = aggregate_stats(it2->second);
          }
        }
      }
    } else if (level == 4 || level == 5) {
      if (map_2d) {
        std::vector<ValueStats> values;
        for (const auto& it1 : *map_2d) {
          if (level == 4 && it1.first == a1_nb) {
            for (const auto& it2 : it1.second) {
              for (const auto& v : it2.second)
                values.push_back(v);
            }
          } else if (it1.first != a1_nb) {
            auto it2 = it1.second.find(a2_nb);
            if (it2 != it1.second.end()) {
              for (const auto& v : it2->second)
                values.push_back(v);
            }
          }
        }
        if (static_cast<int>(values.size()) >= min_observations_bond)
          vs = aggregate_stats(values);
      }
    } else if (level == 6) {
      if (map_2d) {
        std::vector<ValueStats> values;
        for (const auto& it1 : *map_2d) {
          for (const auto& it2 : it1.second) {
            for (const auto& v : it2.second)
              values.push_back(v);
          }
        }
        if (static_cast<int>(values.size()) >= min_observations_bond)
          vs = aggregate_stats(values);
      }
    } else if (level == 7 || level == 8) {
      if (map_nb2) {
        std::vector<ValueStats> values;
        auto it1 = map_nb2->find(a1_nb2);
        if (level == 7 && it1 != map_nb2->end()) {
          for (const auto& it2 : it1->second) {
            for (const auto& it3 : it2.second) {
              for (const auto& it4 : it3.second)
                for (const auto& v : it4.second)
                  values.push_back(v);
            }
          }
        }
        for (const auto& it2 : *map_nb2) {
          if (it2.first == a1_nb2)
            continue;
          auto it3 = it2.second.find(a2_nb2);
          if (it3 != it2.second.end()) {
            for (const auto& it4 : it3->second) {
              for (const auto& it5 : it4.second)
                for (const auto& v : it5.second)
                  values.push_back(v);
            }
          }
        }
        if (static_cast<int>(values.size()) >= min_observations_bond)
          vs = aggregate_stats(values);
      }
    } else if (level == 9) {
      auto it1 = bond_hasp_2d_.find(ha1);
      if (it1 != bond_hasp_2d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end()) {
          auto it3 = it2->second.find(hybr_comb);
          if (it3 != it2->second.end()) {
            auto it4 = it3->second.find(in_ring);
            if (it4 != it3->second.end() && !it4->second.empty())
              vs = aggregate_stats(it4->second);
          }
        }
      }
    } else if (level == 10) {
      auto it1 = bond_hasp_1d_.find(ha1);
      if (it1 != bond_hasp_1d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end()) {
          auto it3 = it2->second.find(hybr_comb);
          if (it3 != it2->second.end() && !it3->second.empty())
            vs = aggregate_stats(it3->second);
        }
      }
    } else if (level == 11) {
      auto it1 = bond_hasp_0d_.find(ha1);
      if (it1 != bond_hasp_0d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end() && !it2->second.empty())
          vs = aggregate_stats(it2->second);
      }
    }

    if (vs.count >= min_observations_bond) {
      if (verbose >= 2)
        std::fprintf(stderr, "      matched: level=%d\n", level);
      return vs;
    }
  }

  // No match found in detailed tables
  if (verbose >= 2)
    std::fprintf(stderr, "      matched: level=none (no multilevel match)\n");
  return ValueStats();
}

inline ValueStats AcedrgTables::search_bond_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& a2, bool in_ring) const {

  BondHRSKey key;
  key.hash1 = std::min(a1.hashing_value, a2.hashing_value);
  key.hash2 = std::max(a1.hashing_value, a2.hashing_value);

  // Build hybrid pair string
  std::string h1 = hybridization_to_string(a1.hybrid);
  std::string h2 = hybridization_to_string(a2.hybrid);
  if (h1 > h2) std::swap(h1, h2);
  key.hybrid_pair = h1 + "_" + h2;
  key.in_ring = in_ring ? "Y" : "N";

  auto it = bond_hrs_.find(key);
  if (it != bond_hrs_.end()) {
    return it->second;
  }

  return ValueStats();
}

inline ValueStats AcedrgTables::search_bond_en(const CodAtomInfo& a1,
    const CodAtomInfo& a2) const {

  std::string elem1 = a1.el.name();
  std::string elem2 = a2.el.name();
  std::string sp1 = hybridization_to_string(a1.hybrid);
  std::string sp2 = hybridization_to_string(a2.hybrid);

  // Canonicalize order
  if (elem1 > elem2 || (elem1 == elem2 && sp1 > sp2)) {
    std::swap(elem1, elem2);
    std::swap(sp1, sp2);
  }

  auto it1 = en_bonds_.find(elem1);
  if (it1 == en_bonds_.end()) return ValueStats();

  auto it2 = it1->second.find(sp1);
  if (it2 == it1->second.end()) return ValueStats();

  auto it3 = it2->second.find(elem2);
  if (it3 == it2->second.end()) return ValueStats();

  auto it4 = it3->second.find(sp2);
  if (it4 == it3->second.end()) return ValueStats();

  if (it4->second.empty()) return ValueStats();

  return aggregate_stats(it4->second);
}

inline ValueStats AcedrgTables::search_metal_bond(const CodAtomInfo& metal,
    const CodAtomInfo& ligand, const std::vector<CodAtomInfo>& atoms) const {

  int ha1 = metal.connectivity;
  int ha2 = static_cast<int>(ligand.conn_atoms_no_metal.size());

  std::vector<std::string> nb_ids;
  nb_ids.reserve(ligand.conn_atoms_no_metal.size());
  for (int nb : ligand.conn_atoms_no_metal) {
    if (nb >= 0 && nb < static_cast<int>(atoms.size()))
      nb_ids.push_back(atoms[nb].el.name());
  }
  std::sort(nb_ids.begin(), nb_ids.end(), compare_no_case);

  std::string ligand_class;
  for (const auto& id : nb_ids)
    ligand_class.append("_" + id);

  const MetalBondEntry* class_entry = nullptr;
  const MetalBondEntry* pre_entry = nullptr;

  for (const auto& entry : metal_bonds_) {
    if (entry.metal != metal.el || entry.ligand != ligand.el)
      continue;
    if (entry.metal_coord != ha1 || entry.ligand_coord != ha2)
      continue;
    if (!pre_entry)
      pre_entry = &entry;
    if (!ligand_class.empty() && entry.ligand_class == ligand_class)
      class_entry = &entry;
  }

  if (class_entry) {
    if (class_entry->count > metal_class_min_count)
      return ValueStats(class_entry->value, class_entry->sigma,
                        class_entry->count);
    if (pre_entry)
      return ValueStats(pre_entry->pre_value, pre_entry->pre_sigma,
                        pre_entry->pre_count);
  }

  if (pre_entry)
    return ValueStats(pre_entry->pre_value, pre_entry->pre_sigma,
                      pre_entry->pre_count);

  return ValueStats();
}

inline bool AcedrgTables::lookup_pep_tors(const std::string& a1,
    const std::string& a2, const std::string& a3, const std::string& a4,
    TorsionEntry& out) const {
  auto it = pep_tors_.find(a1 + "_" + a2 + "_" + a3 + "_" + a4);
  if (it == pep_tors_.end())
    return false;
  out = it->second;
  return true;
}

// ============================================================================
// Implementation - Angle search
// ============================================================================

inline ValueStats AcedrgTables::search_angle_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3) const {

  // Build lookup keys - canonicalize flanking atoms
  int ha1, ha3;
  const CodAtomInfo *flank1, *flank3;
  order_angle_flanks(a1, a3, flank1, flank3);
  ha1 = center.hashing_value;
  int ha2 = flank1->hashing_value;
  ha3 = flank3->hashing_value;

  // Build hybridization tuple
  std::string h1 = hybridization_to_string(center.hybrid);
  std::string h2 = hybridization_to_string(flank3->hybrid);
  std::string h3 = hybridization_to_string(flank1->hybrid);
  std::string hybr_tuple = h1 + "_" + h2 + "_" + h3;

  // Build valueKey (ring:hybr_tuple)
  int ring_val = angle_ring_size(center, *flank1, *flank3);
  std::string value_key = std::to_string(ring_val) + ":" + hybr_tuple;

  // Get neighbor symbols
  const std::string& a1_nb2 = center.nb2_symb;
  const std::string& a2_nb2 = flank1->nb2_symb;
  const std::string& a3_nb2 = flank3->nb2_symb;
  const std::string& a1_root = center.cod_root;
  const std::string& a2_root = flank1->cod_root;
  const std::string& a3_root = flank3->cod_root;
  const std::string& a1_nb = center.nb_symb;
  const std::string& a2_nb = flank1->nb_symb;
  const std::string& a3_nb = flank3->nb_symb;
  const std::string& a1_type = center.cod_main;
  const std::string& a2_type = flank1->cod_main;
  const std::string& a3_type = flank3->cod_main;

  if (verbose >= 2) {
    std::fprintf(stderr,
                 "      angle key %s-%s-%s: hash=%d/%d/%d value=%s "
                 "nb2=%s/%s/%s nb=%s/%s/%s type=%s/%s/%s\n",
                 a1.id.c_str(), center.id.c_str(), a3.id.c_str(),
                 ha1, ha2, ha3, value_key.c_str(),
                 a1_nb2.c_str(), a2_nb2.c_str(), a3_nb2.c_str(),
                 a1_nb.c_str(), a2_nb.c_str(), a3_nb.c_str(),
                 a1_type.c_str(), a2_type.c_str(), a3_type.c_str());
  }

  // Level 1D: Try exact match with atom types
  {
    auto it1 = angle_idx_1d_.find(ha1);
    if (it1 != angle_idx_1d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(a1_nb);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(a2_nb);
                          if (it12 != it11->second.end()) {
                            auto it13 = it12->second.find(a3_nb);
                            if (it13 != it12->second.end()) {
                              auto it14 = it13->second.find(a1_type);
                              if (it14 != it13->second.end()) {
                                auto it15 = it14->second.find(a2_type);
                                if (it15 != it14->second.end()) {
                                  auto it16 = it15->second.find(a3_type);
                                  if (it16 != it15->second.end() && !it16->second.empty()) {
                                    ValueStats vs = it16->second.front();
                                    if (vs.count >= min_observations_angle) {
                                      if (verbose >= 2)
                                        std::fprintf(stderr,
                                                     "      matched angle %s-%s-%s: level=1D\n",
                                                     a1.id.c_str(), center.id.c_str(),
                                                     a3.id.c_str());
                                      return vs;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 2D: No atom types
  {
    auto it1 = angle_idx_2d_.find(ha1);
    if (it1 != angle_idx_2d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(a1_nb);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(a2_nb);
                          if (it12 != it11->second.end()) {
                            auto it13 = it12->second.find(a3_nb);
                            if (it13 != it12->second.end() && !it13->second.empty()) {
                              ValueStats vs = it13->second.front();
                              if (vs.count >= min_observations_angle_fallback) {
                                if (verbose >= 2)
                                  std::fprintf(stderr,
                                               "      matched angle %s-%s-%s: level=2D\n",
                                               a1.id.c_str(), center.id.c_str(),
                                               a3.id.c_str());
                                return vs;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 3D: Hash + valueKey + NB2 only
  {
    auto it1 = angle_idx_3d_.find(ha1);
    if (it1 != angle_idx_3d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end() && !it10->second.empty()) {
                        ValueStats vs = it10->second.front();
                        if (vs.count >= min_observations_angle_fallback) {
                          if (verbose >= 2)
                            std::fprintf(stderr,
                                         "      matched angle %s-%s-%s: level=3D\n",
                                         a1.id.c_str(), center.id.c_str(),
                                         a3.id.c_str());
                          return vs;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 4D: Hash + valueKey only
  {
    auto it1 = angle_idx_4d_.find(ha1);
    if (it1 != angle_idx_4d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end() && !it7->second.empty()) {
                  ValueStats vs = it7->second.front();
                  if (vs.count >= min_observations_angle_fallback) {
                    if (verbose >= 2)
                      std::fprintf(stderr,
                                   "      matched angle %s-%s-%s: level=4D\n",
                                   a1.id.c_str(), center.id.c_str(),
                                   a3.id.c_str());
                    return vs;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 5D: Hash + valueKey only
  {
    auto it1 = angle_idx_5d_.find(ha1);
    if (it1 != angle_idx_5d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end() && !it4->second.empty()) {
            ValueStats vs = it4->second.front();
            if (vs.count >= min_observations_angle_fallback) {
              if (verbose >= 2)
                std::fprintf(stderr,
                             "      matched angle %s-%s-%s: level=5D\n",
                             a1.id.c_str(), center.id.c_str(),
                             a3.id.c_str());
              return vs;
            }
          }
        }
      }
    }
  }

  // Level 6D: Hash only
  {
    auto it1 = angle_idx_6d_.find(ha1);
    if (it1 != angle_idx_6d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end() && !it3->second.empty()) {
          if (verbose >= 2)
            std::fprintf(stderr, "      matched angle %s-%s-%s: level=6D\n",
                         a1.id.c_str(), center.id.c_str(), a3.id.c_str());
          return it3->second.front();
        }
      }
    }
  }

  // No match found in detailed tables
  return ValueStats();
}

inline void AcedrgTables::fill_angle(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Angle& angle) const {

  auto it1 = cc.find_atom(angle.id1.atom);
  auto it2 = cc.find_atom(angle.id2.atom);  // center
  auto it3 = cc.find_atom(angle.id3.atom);

  if (it1 == cc.atoms.end() || it2 == cc.atoms.end() || it3 == cc.atoms.end())
    return;

  int idx1 = static_cast<int>(it1 - cc.atoms.begin());
  int idx2 = static_cast<int>(it2 - cc.atoms.begin());
  int idx3 = static_cast<int>(it3 - cc.atoms.begin());

  const CodAtomInfo& a1 = atom_info[idx1];
  const CodAtomInfo& center = atom_info[idx2];
  const CodAtomInfo& a3 = atom_info[idx3];

  // Check for metal center
  if (center.is_metal) {
    std::vector<double> ideal_angles = get_metal_angles(center.el,
                                                        center.connectivity);
    if (!ideal_angles.empty()) {
      // Use the most common angle for this geometry
      angle.value = ideal_angles[0];
      angle.esd = clamp_angle_sigma(3.0);
      return;
    }
  }

  int ring_size = angle_ring_size(center, a1, a3);

  // Try detailed multilevel search first
  ValueStats vs = search_angle_multilevel(a1, center, a3);
  if (vs.count >= min_observations_angle) {
    angle.value = vs.value;
    angle.esd = clamp_angle_sigma(vs.sigma);
    return;
  }

  // Try HRS table
  vs = search_angle_hrs(a1, center, a3, ring_size);
  if (vs.count >= min_observations_angle) {
    angle.value = vs.value;
    angle.esd = clamp_angle_sigma(vs.sigma);
    return;
  }

  // Fallback based on hybridization
  switch (center.hybrid) {
    case Hybridization::SP1:
      angle.value = 180.0;
      break;
    case Hybridization::SP2:
      angle.value = 120.0;
      break;
    case Hybridization::SP3:
      angle.value = 109.47;  // tetrahedral
      break;
    case Hybridization::SPD5:
    case Hybridization::SPD6:
      angle.value = 90.0;    // octahedral
      break;
    default:
      angle.value = 109.47;  // default tetrahedral
      break;
  }
  angle.esd = upper_angle_sigma;
}

inline ValueStats AcedrgTables::search_angle_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3, int ring_size) const {

  AngleHRSKey key;
  key.hash1 = center.hashing_value;
  key.hash2 = std::min(a1.hashing_value, a3.hashing_value);
  key.hash3 = std::max(a1.hashing_value, a3.hashing_value);

  // Build hybrid tuple
  std::string h1 = hybridization_to_string(center.hybrid);
  std::string h2;
  std::string h3;
  if (a1.hashing_value > a3.hashing_value) {
    h2 = hybridization_to_string(a1.hybrid);
    h3 = hybridization_to_string(a3.hybrid);
  } else {
    h2 = hybridization_to_string(a3.hybrid);
    h3 = hybridization_to_string(a1.hybrid);
  }
  std::string hybr_tuple = h1 + "_" + h2 + "_" + h3;
  key.value_key = std::to_string(ring_size) + ":" + hybr_tuple;

  auto it = angle_hrs_.find(key);
  if (it != angle_hrs_.end()) {
    return it->second;
  }

  return ValueStats();
}

inline std::vector<double> AcedrgTables::get_metal_angles(Element metal,
    int coord_number) const {

  std::vector<double> angles;

  // First check if we have specific angles in the table
  for (const auto& entry : metal_angles_) {
    if (entry.metal == metal && entry.coord_number == coord_number) {
      angles.push_back(entry.angle);
    }
  }

  if (!angles.empty()) {
    return angles;
  }

  // Fall back to geometry defaults
  auto geo_it = metal_coord_geo_.find(metal);
  if (geo_it != metal_coord_geo_.end()) {
    auto cn_it = geo_it->second.find(coord_number);
    if (cn_it != geo_it->second.end()) {
      switch (cn_it->second) {
        case CoordGeometry::LINEAR:
          angles.push_back(180.0);
          break;
        case CoordGeometry::TRIGONAL_PLANAR:
          angles.push_back(120.0);
          break;
        case CoordGeometry::T_SHAPED:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::TETRAHEDRAL:
          angles.push_back(109.47);
          break;
        case CoordGeometry::SQUARE_PLANAR:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::TRIGONAL_BIPYRAMIDAL:
          angles.push_back(90.0);
          angles.push_back(120.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::SQUARE_PYRAMIDAL:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::OCTAHEDRAL:
          angles.push_back(90.0);
          break;
        default:
          break;
      }
    }
  }

  // Ultimate fallback based on coordination number
  if (angles.empty()) {
    switch (coord_number) {
      case 2: angles.push_back(180.0); break;
      case 3: angles.push_back(120.0); break;
      case 4: angles.push_back(109.47); break;
      case 5: angles.push_back(90.0); angles.push_back(120.0); break;
      case 6: angles.push_back(90.0); break;
      default: angles.push_back(90.0); break;
    }
  }

  return angles;
}

// ============================================================================
// Implementation - Statistical helpers
// ============================================================================

inline ValueStats AcedrgTables::aggregate_stats(
    const std::vector<ValueStats>& values) const {

  if (values.empty())
    return ValueStats();

  if (values.size() == 1)
    return values[0];

  // Weighted average by count
  double sum_val = 0;
  double sum_weight = 0;
  int total_count = 0;

  for (const auto& v : values) {
    if (v.count > 0) {
      sum_val += v.value * v.count;
      sum_weight += v.count;
      total_count += v.count;
    }
  }

  if (sum_weight == 0)
    return ValueStats();

  double mean = sum_val / sum_weight;

  // Compute pooled standard deviation
  double sum_sq = 0;
  for (const auto& v : values) {
    if (v.count > 0) {
      double diff = v.value - mean;
      sum_sq += v.count * (v.sigma * v.sigma + diff * diff);
    }
  }

  double sigma = std::sqrt(sum_sq / sum_weight);

  return ValueStats(mean, sigma, total_count);
}

} // namespace gemmi

#endif
