// Copyright 2025 Global Phasing Ltd.
//
// AcedrgTables - COD/CSD-based atom classification and restraint value lookup
// Port of AceDRG codClassify system to gemmi.

#ifndef GEMMI_ACEDRG_TABLES_HPP_
#define GEMMI_ACEDRG_TABLES_HPP_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <sstream>
#include "chemcomp.hpp"
#include "elem.hpp"
#include "fail.hpp"

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
    case Hybridization::SP_NON: return "SP-NON";
  }
  return "SP-NON";
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
  int hashing_value;        // 0-1000+ hash code
  Element el;               // Element
  Hybridization hybrid;     // Hybridization state
  std::string cod_class;    // Full COD class (e.g., "C[6a](C[6a]C[6a])(C[6a])(H)")
  std::string cod_main;     // Main atom type without neighbors
  std::string nb_symb;      // 1st neighbor symbol string
  std::string nb2_symb;     // 2nd neighbor symbol string
  int connectivity;         // Number of bonded atoms
  int min_ring_size;        // Minimum ring size (0 = not in ring)
  bool is_aromatic;         // In aromatic ring
  bool is_metal;            // Is a metal atom
  int excess_electrons;     // Formal charge/lone pair info
  std::string nb1nb2_sp;    // Hybridization info of neighbors
  std::string nb1nb2_ex_elec; // Excess electron info of neighbors

  CodAtomInfo()
    : index(-1), hashing_value(0), el(El::X), hybrid(Hybridization::SP_NON),
      connectivity(0), min_ring_size(0), is_aromatic(false), is_metal(false),
      excess_electrons(0) {}
};

// Statistical value with count
struct ValueStats {
  double value;
  double sigma;
  int count;

  ValueStats() : value(NAN), sigma(NAN), count(0) {}
  ValueStats(double v, double s, int c) : value(v), sigma(s), count(c) {}
};

// Metal bond entry
struct MetalBondEntry {
  Element metal = El::X;
  Element ligand = El::X;
  int coord_number = 0;
  std::string ligand_class;
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
  int n = el.ordinal();
  if (n == 0) return 0;
  // Period 1
  if (n == 1) return 1;  // H
  if (n == 2) return 18; // He
  // Period 2
  if (n >= 3 && n <= 4) return n - 2;     // Li, Be
  if (n >= 5 && n <= 10) return n + 8;    // B-Ne
  // Period 3
  if (n >= 11 && n <= 12) return n - 10;  // Na, Mg
  if (n >= 13 && n <= 18) return n;       // Al-Ar
  // Period 4
  if (n >= 19 && n <= 20) return n - 18;  // K, Ca
  if (n >= 21 && n <= 30) return n - 18;  // Sc-Zn
  if (n >= 31 && n <= 36) return n - 18;  // Ga-Kr
  // Period 5
  if (n >= 37 && n <= 38) return n - 36;  // Rb, Sr
  if (n >= 39 && n <= 48) return n - 36;  // Y-Cd
  if (n >= 49 && n <= 54) return n - 36;  // In-Xe
  // Period 6
  if (n >= 55 && n <= 56) return n - 54;  // Cs, Ba
  if (n >= 57 && n <= 71) return 3;       // La-Lu (lanthanides, group 3)
  if (n >= 72 && n <= 80) return n - 68;  // Hf-Hg
  if (n >= 81 && n <= 86) return n - 68;  // Tl-Rn
  // Period 7
  if (n >= 87 && n <= 88) return n - 86;  // Fr, Ra
  if (n >= 89 && n <= 103) return 3;      // Ac-Lr (actinides, group 3)
  if (n >= 104 && n <= 112) return n - 100; // Rf-Cn
  if (n >= 113 && n <= 118) return n - 100; // Nh-Og
  return 0;
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
  int min_observations = 5;   // Minimum count for exact match
  int min_observations_fallback = 4; // Minimum for fallback levels

private:
  static constexpr int HASH_SIZE = 1000;

  // Table directory
  std::string tables_dir_;
  bool tables_loaded_ = false;

  // Hash code tables
  std::map<int, std::string> digit_keys_;  // hash -> footprint
  std::map<int, int> linked_hash_;         // hash -> linked hash

  // HRS (High-Resolution Summary) bond tables
  // Key: hash1, hash2, hybrid_pair, aromatic
  struct BondHRSKey {
    int hash1, hash2;
    std::string hybrid_pair;
    bool aromatic;
    bool operator<(const BondHRSKey& o) const {
      if (hash1 != o.hash1) return hash1 < o.hash1;
      if (hash2 != o.hash2) return hash2 < o.hash2;
      if (hybrid_pair != o.hybrid_pair) return hybrid_pair < o.hybrid_pair;
      return aromatic < o.aromatic;
    }
  };
  std::map<BondHRSKey, ValueStats> bond_hrs_;

  // HRS angle tables
  // Key: hash1, hash2, hash3, hybrid_tuple, aromatic
  struct AngleHRSKey {
    int hash1, hash2, hash3;
    std::string hybrid_tuple;
    bool aromatic;
    bool operator<(const AngleHRSKey& o) const {
      if (hash1 != o.hash1) return hash1 < o.hash1;
      if (hash2 != o.hash2) return hash2 < o.hash2;
      if (hash3 != o.hash3) return hash3 < o.hash3;
      if (hybrid_tuple != o.hybrid_tuple) return hybrid_tuple < o.hybrid_tuple;
      return aromatic < o.aromatic;
    }
  };
  std::map<AngleHRSKey, ValueStats> angle_hrs_;

  // Detailed indexed bond tables from allOrgBondTables/*.table
  // Level 0: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2, a1NB, a2NB, a1TypeM, a2TypeM
  using BondIdx1D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>>>>>;
  BondIdx1D bond_idx_1d_;

  // Levels 3-6: ha1, ha2, hybrComb, inRing, a1NB2, a2NB2, a1NB, a2NB (no atom types)
  using BondIdx2D = std::map<int, std::map<int,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>>>>>;
  BondIdx2D bond_idx_2d_;

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
  // Similar multi-level structure for angles (to be implemented)
  // For now, keep existing HRS-based angle lookup

  // Element + hybridization based fallback bonds
  using ENBonds = std::map<std::string, std::map<std::string,
    std::map<std::string, std::map<std::string,
    std::vector<ValueStats>>>>>;
  ENBonds en_bonds_;

  // Metal bond tables
  std::vector<MetalBondEntry> metal_bonds_;
  std::vector<MetalAngleEntry> metal_angles_;
  std::map<Element, std::map<int, CoordGeometry>> metal_coord_geo_;

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

  // Atom classification helpers
  void detect_rings(const ChemComp& cc,
                    std::vector<std::vector<int>>& rings) const;
  struct BondInfo {
    int neighbor_idx;
    BondType type;
  };
  std::vector<std::vector<BondInfo>> build_adjacency(const ChemComp& cc) const;
  int compute_min_ring(int atom_idx,
                       const std::vector<std::vector<int>>& rings) const;
  Hybridization determine_hybridization(const ChemComp& cc,
                                        int atom_idx,
                                        int connectivity,
                                        bool in_ring,
                                        bool is_aromatic) const;
  void compute_hash(CodAtomInfo& atom) const;
  void compute_nb_symbols(const std::vector<std::vector<BondInfo>>& adj,
                          std::vector<CodAtomInfo>& atoms) const;
  void compute_cod_class(const ChemComp& cc,
                         const std::vector<std::vector<BondInfo>>& adj,
                         std::vector<CodAtomInfo>& atoms) const;
  void refine_hybridization(const std::vector<std::vector<BondInfo>>& adj,
                            std::vector<CodAtomInfo>& atoms) const;
  int count_non_metal_neighbors(const std::vector<BondInfo>& neighbors,
                                const std::vector<CodAtomInfo>& atoms) const;

  // Bond search helpers
  ValueStats search_bond_multilevel(const CodAtomInfo& a1, const CodAtomInfo& a2,
                                    bool aromatic) const;
  ValueStats search_bond_hrs(const CodAtomInfo& a1, const CodAtomInfo& a2,
                             bool aromatic) const;
  ValueStats search_bond_en(const CodAtomInfo& a1, const CodAtomInfo& a2) const;
  ValueStats search_metal_bond(const CodAtomInfo& metal,
                               const CodAtomInfo& ligand,
                               int coord_number) const;

  // Angle search helpers
  ValueStats search_angle_hrs(const CodAtomInfo& a1, const CodAtomInfo& center,
                              const CodAtomInfo& a3, bool aromatic) const;
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

  // Load detailed indexed tables
  load_atom_type_codes(tables_dir + "/allAtomTypesFromMolsCoded.list");
  load_bond_index(tables_dir + "/allOrgBondTables/bond_idx.table");
  load_bond_tables(tables_dir + "/allOrgBondTables");

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
    std::string hybrid_pair, arom_str;
    double value, sigma;
    int count;

    if (iss >> hash1 >> hash2 >> hybrid_pair >> arom_str >> value >> sigma >> count) {
      BondHRSKey key;
      key.hash1 = std::min(hash1, hash2);
      key.hash2 = std::max(hash1, hash2);
      key.hybrid_pair = hybrid_pair;
      key.aromatic = (arom_str == "Y" || arom_str == "y");
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
    std::string hybrid_tuple, arom_str;
    double value, sigma;
    int count;

    if (iss >> hash1 >> hash2 >> hash3 >> hybrid_tuple >> arom_str
        >> value >> sigma >> count) {
      AngleHRSKey key;
      // For angles, center is always hash2, but we canonicalize flanking atoms
      key.hash1 = std::min(hash1, hash3);
      key.hash2 = hash2;
      key.hash3 = std::max(hash1, hash3);
      key.hybrid_tuple = hybrid_tuple;
      key.aromatic = (arom_str == "Y" || arom_str == "y");
      angle_hrs_[key] = ValueStats(value, sigma, count);
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
    std::string metal_str, ligand_str, ligand_class;
    int coord;
    double value, sigma;
    int count;

    if (iss >> metal_str >> ligand_str >> coord >> ligand_class >> value >> sigma >> count) {
      MetalBondEntry entry;
      entry.metal = Element(metal_str);
      entry.ligand = Element(ligand_str);
      entry.coord_number = coord;
      entry.ligand_class = ligand_class;
      entry.value = value;
      entry.sigma = sigma;
      entry.count = count;
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

        if (!(iss >> ha1 >> ha2 >> hybr_comb >> in_ring
                  >> a1_nb2 >> a2_nb2 >> a1_nb >> a2_nb
                  >> atom_code1 >> atom_code2
                  >> value >> sigma >> count))
          continue;

        // Get main atom types from codes
        std::string a1_type_m, a2_type_m;
        auto it1 = atom_type_codes_.find(atom_code1);
        auto it2 = atom_type_codes_.find(atom_code2);
        if (it1 != atom_type_codes_.end()) {
          // Extract main type (before '{' if present)
          a1_type_m = it1->second;
          size_t brace = a1_type_m.find('{');
          if (brace != std::string::npos)
            a1_type_m = a1_type_m.substr(0, brace);
        }
        if (it2 != atom_type_codes_.end()) {
          a2_type_m = it2->second;
          size_t brace = a2_type_m.find('{');
          if (brace != std::string::npos)
            a2_type_m = a2_type_m.substr(0, brace);
        }

        ValueStats vs(value, sigma, count);

        // Populate 1D structure (full detail)
        bond_idx_1d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb][a1_type_m][a2_type_m].push_back(vs);

        // Populate 2D structure (no atom types)
        bond_idx_2d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb].push_back(vs);

        // Populate HaSp structures for fallback
        bond_hasp_2d_[ha1][ha2][hybr_comb][in_ring].push_back(vs);
        bond_hasp_1d_[ha1][ha2][hybr_comb].push_back(vs);
        bond_hasp_0d_[ha1][ha2].push_back(vs);
      }
    }
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
    atoms[i].el = cc.atoms[i].el;
    atoms[i].is_metal = is_metal(cc.atoms[i].el);
  }

  // Build connectivity
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = cc.find_atom(bond.id1.atom);
    auto it2 = cc.find_atom(bond.id2.atom);
    if (it1 != cc.atoms.end() && it2 != cc.atoms.end()) {
      int idx1 = static_cast<int>(it1 - cc.atoms.begin());
      int idx2 = static_cast<int>(it2 - cc.atoms.begin());
      atoms[idx1].connectivity++;
      atoms[idx2].connectivity++;
    }
  }

  // Detect rings and set ring-related properties
  std::vector<std::vector<int>> rings;
  detect_rings(cc, rings);

  // Set min ring size
  for (size_t i = 0; i < atoms.size(); ++i) {
    atoms[i].min_ring_size = compute_min_ring(static_cast<int>(i), rings);
  }

  // Determine aromaticity from bonds
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

  // Build adjacency list (needed for hybridization refinement and later steps)
  std::vector<std::vector<BondInfo>> adj = build_adjacency(cc);

  // Determine hybridization - Phase 1: element-based rules
  for (size_t i = 0; i < atoms.size(); ++i) {
    atoms[i].hybrid = determine_hybridization(cc, static_cast<int>(i),
                                              atoms[i].connectivity,
                                              atoms[i].min_ring_size > 0,
                                              atoms[i].is_aromatic);
  }

  // Refine hybridization - Phase 2: neighbor-based refinement
  refine_hybridization(adj, atoms);

  // Compute hash codes
  for (auto& atom : atoms) {
    compute_hash(atom);
  }

  // Compute neighbor symbols
  compute_nb_symbols(adj, atoms);

  // Compute full COD class
  compute_cod_class(cc, adj, atoms);

  return atoms;
}

inline void AcedrgTables::detect_rings(const ChemComp& cc,
                                      std::vector<std::vector<int>>& rings) const {
  for (const auto& bond : cc.rt.bonds) {
    std::vector<Restraints::AtomId> path =
        cc.rt.find_shortest_path(bond.id1, bond.id2, {}, 2);
    if (path.size() < 3)
      continue;

    std::vector<int> ring;
    ring.reserve(path.size());
    for (const auto& id : path) {
      auto it = cc.find_atom(id.atom);
      if (it != cc.atoms.end())
        ring.push_back(static_cast<int>(it - cc.atoms.begin()));
    }
    if (ring.size() < 3)
      continue;
    std::sort(ring.begin(), ring.end());
    ring.erase(std::unique(ring.begin(), ring.end()), ring.end());
    if (ring.size() >= 3)
      rings.push_back(std::move(ring));
  }

  // Remove duplicate rings (same atoms in different order)
  std::vector<std::vector<int>> unique_rings;
  for (auto& ring : rings) {
    std::sort(ring.begin(), ring.end());
    bool found = false;
    for (const auto& ur : unique_rings) {
      if (ur == ring) {
        found = true;
        break;
      }
    }
    if (!found) {
      unique_rings.push_back(ring);
    }
  }
  rings = std::move(unique_rings);
}

inline int AcedrgTables::compute_min_ring(int atom_idx,
    const std::vector<std::vector<int>>& rings) const {
  int min_size = 0;
  for (const auto& ring : rings) {
    if (std::find(ring.begin(), ring.end(), atom_idx) != ring.end()) {
      int size = static_cast<int>(ring.size());
      if (min_size == 0 || size < min_size) {
        min_size = size;
      }
    }
  }
  return min_size;
}

inline Hybridization AcedrgTables::determine_hybridization(
    const ChemComp& cc, int atom_idx, int connectivity,
    bool in_ring, bool is_aromatic) const {

  Element el = cc.atoms[atom_idx].el;

  // Metals have special hybridization
  if (is_metal(el)) {
    if (connectivity <= 4) return Hybridization::SPD5;
    if (connectivity == 5) return Hybridization::SPD5;
    if (connectivity == 6) return Hybridization::SPD6;
    if (connectivity == 7) return Hybridization::SPD7;
    return Hybridization::SPD8;
  }

  // Count double and triple bonds to this atom
  int double_count = 0;
  int triple_count = 0;

  for (const auto& bond : cc.rt.bonds) {
    auto it1 = cc.find_atom(bond.id1.atom);
    auto it2 = cc.find_atom(bond.id2.atom);
    if (it1 == cc.atoms.end() || it2 == cc.atoms.end())
      continue;

    int idx1 = static_cast<int>(it1 - cc.atoms.begin());
    int idx2 = static_cast<int>(it2 - cc.atoms.begin());

    if (idx1 != atom_idx && idx2 != atom_idx)
      continue;

    if (bond.type == BondType::Double)
      double_count++;
    else if (bond.type == BondType::Triple)
      triple_count++;
  }

  // Determine hybridization based on element and bonding
  if (triple_count > 0) {
    return Hybridization::SP1;
  }

  if (is_aromatic) {
    return Hybridization::SP2;
  }

  if (double_count > 0) {
    return Hybridization::SP2;
  }

  // Default based on connectivity
  if (connectivity == 2) {
    // Could be SP1 (linear) or SP2/SP3 (bent)
    // Check if it's a linear arrangement based on element
    if (el == El::C || el == El::N) {
      // Need more context - default to SP3 for now
      return in_ring ? Hybridization::SP2 : Hybridization::SP3;
    }
  }

  if (connectivity == 3) {
    return Hybridization::SP2;
  }

  return Hybridization::SP3;
}

inline int AcedrgTables::count_non_metal_neighbors(
    const std::vector<BondInfo>& neighbors,
    const std::vector<CodAtomInfo>& atoms) const {
  int count = 0;
  for (const auto& nb : neighbors)
    if (!atoms[nb.neighbor_idx].is_metal)
      ++count;
  return count;
}

// Phase 2 hybridization refinement based on neighbor states.
// This implements the second pass from AceDRG's setAtomsBondingAndChiralCenter().
inline void AcedrgTables::refine_hybridization(
    const std::vector<std::vector<BondInfo>>& adj,
    std::vector<CodAtomInfo>& atoms) const {

  // Phase 2a: Oxygen atoms with 2 non-metal connections.
  // If any neighbor is SP2, the oxygen is also SP2 (e.g., ester oxygen).
  for (size_t i = 0; i < atoms.size(); ++i) {
    if (atoms[i].el != El::O || atoms[i].is_metal)
      continue;
    if (count_non_metal_neighbors(adj[i], atoms) != 2)
      continue;
    for (const auto& nb : adj[i]) {
      if (!atoms[nb.neighbor_idx].is_metal &&
          atoms[nb.neighbor_idx].hybrid == Hybridization::SP2) {
        atoms[i].hybrid = Hybridization::SP2;
        break;
      }
    }
  }

  // Save hybridization state before Phase 2b to avoid order-dependency
  std::vector<Hybridization> pre_hybrid(atoms.size());
  for (size_t i = 0; i < atoms.size(); ++i)
    pre_hybrid[i] = atoms[i].hybrid;

  // Phase 2b: Nitrogen/Arsenic with 3 non-metal connections.
  // If any non-oxygen neighbor was SP2, set to SP2 (e.g., amide nitrogen).
  for (size_t i = 0; i < atoms.size(); ++i) {
    if (atoms[i].el != El::N && atoms[i].el != El::As)
      continue;
    if (atoms[i].is_metal)
      continue;
    if (count_non_metal_neighbors(adj[i], atoms) != 3)
      continue;
    for (const auto& nb : adj[i]) {
      if (atoms[nb.neighbor_idx].is_metal)
        continue;
      if (atoms[nb.neighbor_idx].el == El::O)
        continue;
      if (pre_hybrid[nb.neighbor_idx] == Hybridization::SP2) {
        atoms[i].hybrid = Hybridization::SP2;
        break;
      }
    }
  }
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
  long prime_product = static_cast<long>(primes[d1]) *
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

inline void AcedrgTables::compute_nb_symbols(
    const std::vector<std::vector<BondInfo>>& adj,
    std::vector<CodAtomInfo>& atoms) const {
  // For each atom, build neighbor symbols
  for (size_t i = 0; i < atoms.size(); ++i) {
    std::vector<std::string> nb1_parts;
    std::vector<std::string> nb2_parts;

    for (const auto& nb : adj[i]) {
      const CodAtomInfo& nb_atom = atoms[nb.neighbor_idx];
      // First neighbor: element + connectivity
      std::string nb1 = std::string(nb_atom.el.name()) + "-" +
                        std::to_string(nb_atom.connectivity) + ":";
      nb1_parts.push_back(nb1);

      // Second neighbor: just connectivity count
      nb2_parts.push_back(std::to_string(nb_atom.connectivity) + ":");
    }

    // Sort for canonical form
    std::sort(nb1_parts.begin(), nb1_parts.end());
    std::sort(nb2_parts.begin(), nb2_parts.end());

    // Concatenate
    for (const auto& p : nb1_parts) atoms[i].nb_symb += p;
    for (const auto& p : nb2_parts) atoms[i].nb2_symb += p;
  }
}

inline void AcedrgTables::compute_cod_class(const ChemComp& cc,
    const std::vector<std::vector<BondInfo>>& adj,
    std::vector<CodAtomInfo>& atoms) const {
  // Build COD class for each atom
  // Format: Element[ring_size aromatic](neighbor1)(neighbor2)...
  for (size_t i = 0; i < atoms.size(); ++i) {
    std::string& cod = atoms[i].cod_class;
    cod = cc.atoms[i].el.name();

    // Ring annotation
    if (atoms[i].min_ring_size > 0) {
      cod += "[";
      cod += std::to_string(atoms[i].min_ring_size);
      if (atoms[i].is_aromatic) cod += "a";
      cod += "]";
    }

    // Main type (without neighbors)
    atoms[i].cod_main = cod;

    // Neighbor groups
    std::vector<std::string> nb_groups;
    for (const auto& nb : adj[i]) {
      const CodAtomInfo& nb_atom = atoms[nb.neighbor_idx];
      std::string nb_str = cc.atoms[nb.neighbor_idx].el.name();
      if (nb_atom.min_ring_size > 0) {
        nb_str += "[";
        nb_str += std::to_string(nb_atom.min_ring_size);
        if (nb_atom.is_aromatic) nb_str += "a";
        nb_str += "]";
      }
      nb_groups.push_back(nb_str);
    }

    // Sort neighbor groups for canonical form
    std::sort(nb_groups.begin(), nb_groups.end());

    // Add to class
    for (const auto& ng : nb_groups) {
      cod += "(" + ng + ")";
    }
  }
}

// ============================================================================
// Implementation - Bond search
// ============================================================================

inline void AcedrgTables::fill_restraints(ChemComp& cc) const {
  if (!tables_loaded_)
    return;

  // Classify atoms
  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);

  // Fill bonds
  for (auto& bond : cc.rt.bonds) {
    if (std::isnan(bond.value)) {
      fill_bond(cc, atom_info, bond);
    }
  }

  // Fill angles
  for (auto& angle : cc.rt.angles) {
    if (std::isnan(angle.value)) {
      fill_angle(cc, atom_info, angle);
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

  // Check for metal bond
  if (a1.is_metal || a2.is_metal) {
    const CodAtomInfo& metal = a1.is_metal ? a1 : a2;
    const CodAtomInfo& ligand = a1.is_metal ? a2 : a1;
    int coord = metal.connectivity;

    ValueStats vs = search_metal_bond(metal, ligand, coord);
    if (vs.count >= min_observations) {
      bond.value = vs.value;
      bond.esd = clamp_bond_sigma(vs.sigma);
      return;
    }
  }

  // Try detailed multilevel search first (if tables loaded)
  ValueStats vs = search_bond_multilevel(a1, a2, bond.aromatic);
  if (vs.count >= min_observations) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    return;
  }

  // Try HRS (summary) table
  vs = search_bond_hrs(a1, a2, bond.aromatic);
  if (vs.count >= min_observations) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    return;
  }

  // Try element+hybridization fallback
  vs = search_bond_en(a1, a2);
  if (vs.count > 0) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    return;
  }

  // Ultimate fallback: sum of covalent radii
  bond.value = a1.el.covalent_r() + a2.el.covalent_r();
  bond.esd = upper_bond_sigma;
}

inline ValueStats AcedrgTables::search_bond_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& a2, bool /*aromatic*/) const {

  // Build lookup keys
  int ha1 = std::min(a1.hashing_value, a2.hashing_value);
  int ha2 = std::max(a1.hashing_value, a2.hashing_value);

  std::string h1 = hybridization_to_string(a1.hybrid);
  std::string h2 = hybridization_to_string(a2.hybrid);
  if (h1 > h2) std::swap(h1, h2);
  std::string hybr_comb = h1 + "_" + h2;

  std::string in_ring = (a1.min_ring_size > 0 || a2.min_ring_size > 0) ? "Y" : "N";

  // Use neighbor symbols as additional keys
  const std::string& a1_nb2 = a1.nb2_symb;
  const std::string& a2_nb2 = a2.nb2_symb;
  const std::string& a1_nb = a1.nb_symb;
  const std::string& a2_nb = a2.nb_symb;

  // Use COD main type as atom type (simplified)
  const std::string& a1_type = a1.cod_main;
  const std::string& a2_type = a2.cod_main;

  // Level 0: Try exact match in 1D structure
  {
    auto it1 = bond_idx_1d_.find(ha1);
    if (it1 != bond_idx_1d_.end()) {
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
                      if (it10 != it9->second.end() && !it10->second.empty()) {
                        ValueStats vs = aggregate_stats(it10->second);
                        if (vs.count >= min_observations)
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

  // Level 3: Try 2D structure (no atom types)
  {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 != bond_idx_2d_.end()) {
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
                  if (it8 != it7->second.end() && !it8->second.empty()) {
                    ValueStats vs = aggregate_stats(it8->second);
                    if (vs.count >= min_observations_fallback)
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

  // Level 9: Try HaSp2D (hash + hybridization + ring)
  {
    auto it1 = bond_hasp_2d_.find(ha1);
    if (it1 != bond_hasp_2d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(hybr_comb);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(in_ring);
          if (it4 != it3->second.end() && !it4->second.empty()) {
            ValueStats vs = aggregate_stats(it4->second);
            if (vs.count >= min_observations_fallback)
              return vs;
          }
        }
      }
    }
  }

  // Level 10: Try HaSp1D (hash + hybridization only)
  {
    auto it1 = bond_hasp_1d_.find(ha1);
    if (it1 != bond_hasp_1d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(hybr_comb);
        if (it3 != it2->second.end() && !it3->second.empty()) {
          ValueStats vs = aggregate_stats(it3->second);
          if (vs.count >= min_observations_fallback)
            return vs;
        }
      }
    }
  }

  // Level 11: Try HaSp0D (hash codes only)
  {
    auto it1 = bond_hasp_0d_.find(ha1);
    if (it1 != bond_hasp_0d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end() && !it2->second.empty()) {
        return aggregate_stats(it2->second);
      }
    }
  }

  // No match found in detailed tables
  return ValueStats();
}

inline ValueStats AcedrgTables::search_bond_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& a2, bool aromatic) const {

  BondHRSKey key;
  key.hash1 = std::min(a1.hashing_value, a2.hashing_value);
  key.hash2 = std::max(a1.hashing_value, a2.hashing_value);

  // Build hybrid pair string
  std::string h1 = hybridization_to_string(a1.hybrid);
  std::string h2 = hybridization_to_string(a2.hybrid);
  if (h1 > h2) std::swap(h1, h2);
  key.hybrid_pair = h1 + "_" + h2;
  key.aromatic = aromatic;

  auto it = bond_hrs_.find(key);
  if (it != bond_hrs_.end()) {
    return it->second;
  }

  // Try without aromatic flag
  key.aromatic = !aromatic;
  it = bond_hrs_.find(key);
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
    const CodAtomInfo& ligand, int coord_number) const {

  std::vector<ValueStats> matches;

  for (const auto& entry : metal_bonds_) {
    if (entry.metal == metal.el && entry.ligand == ligand.el) {
      if (entry.coord_number == coord_number ||
          entry.coord_number == 0) { // 0 means any coordination
        matches.emplace_back(entry.value, entry.sigma, entry.count);
      }
    }
  }

  if (!matches.empty()) {
    return aggregate_stats(matches);
  }

  return ValueStats();
}

// ============================================================================
// Implementation - Angle search
// ============================================================================

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

  // Determine if any atom is aromatic
  bool aromatic = a1.is_aromatic || center.is_aromatic || a3.is_aromatic;

  // Try HRS table
  ValueStats vs = search_angle_hrs(a1, center, a3, aromatic);
  if (vs.count >= min_observations) {
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
    const CodAtomInfo& center, const CodAtomInfo& a3, bool aromatic) const {

  AngleHRSKey key;
  key.hash1 = std::min(a1.hashing_value, a3.hashing_value);
  key.hash2 = center.hashing_value;
  key.hash3 = std::max(a1.hashing_value, a3.hashing_value);

  // Build hybrid tuple
  std::string h1 = hybridization_to_string(a1.hybrid);
  std::string h2 = hybridization_to_string(center.hybrid);
  std::string h3 = hybridization_to_string(a3.hybrid);
  if (h1 > h3) std::swap(h1, h3);
  key.hybrid_tuple = h1 + "_" + h2 + "_" + h3;
  key.aromatic = aromatic;

  auto it = angle_hrs_.find(key);
  if (it != angle_hrs_.end()) {
    return it->second;
  }

  // Try without aromatic flag
  key.aromatic = !aromatic;
  it = angle_hrs_.find(key);
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
