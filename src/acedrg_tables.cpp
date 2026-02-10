// Copyright 2025 Global Phasing Ltd.
//
// AcedrgTables - COD/CSD-based atom classification and restraint value lookup
// Port of AceDRG codClassify system to gemmi.

#include "gemmi/acedrg_tables.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <set>
#include "gemmi/atof.hpp"
#include "gemmi/atox.hpp"
#include "gemmi/fileutil.hpp"
#include "gemmi/elem.hpp"
#include "gemmi/fail.hpp"
#include "gemmi/util.hpp"
#include "gemmi/read_cif.hpp"

namespace gemmi {

namespace {

// Look up a key in a map, returning pointer to value or nullptr.
template<typename Map, typename Key>
const typename Map::mapped_type* find_val(const Map& m, const Key& k) {
  auto it = m.find(k);
  return it != m.end() ? &it->second : nullptr;
}

// Returns the prefix of s before the first occurrence of ch, or s itself.
std::string prefix_before(const std::string& s, char ch) {
  size_t pos = s.find(ch);
  return pos != std::string::npos ? s.substr(0, pos) : s;
}

// Expected valence for common non-metal elements (used in metal charge adjustment).
int get_expected_valence(Element el) {
  if (el == El::O) return 2;
  if (el == El::N) return 3;
  if (el == El::S) return 2;
  if (el == El::C) return 4;
  if (el == El::P) return 3;
  return 0;
}

// Skip blank lines, comments, and empty lines when reading table files.
bool is_skip_line(const char* line) {
  return line[0] == '\n' || line[0] == '#' || line[0] == '\0';
}

bool compare_no_case(const std::string& first,
                     const std::string& second) {
  size_t i = 0;
  while (i < first.length() && i < second.length()) {
    char a = alpha_up(first[i]);
    char b = alpha_up(second[i]);
    if (a < b)
      return true;
    if (a > b)
      return false;
    ++i;
  }
  return first.length() > second.length();
}

bool compare_no_case2(const std::string& first,
                      const std::string& second) {
  if (first.length() > second.length())
    return true;
  if (first.length() < second.length())
    return false;
  for (size_t i = 0; i < first.length() && i < second.length(); ++i) {
    char a = alpha_up(first[i]);
    char b = alpha_up(second[i]);
    if (a < b)
      return true;
    if (a > b)
      return false;
  }
  return false;
}

struct SortMap {
  std::string key;
  int val = 0;
  int nNB = 0;
};

bool desc_sort_map_key(const SortMap& a, const SortMap& b) {
  if (a.key.length() != b.key.length())
    return a.key.length() > b.key.length();
  return a.nNB > b.nNB;
}

}  // namespace

const char* hybridization_to_string(Hybridization h) {
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

Hybridization hybridization_from_string(const std::string& s) {
  if (s == "SP1") return Hybridization::SP1;
  if (s == "SP2") return Hybridization::SP2;
  if (s == "SP3") return Hybridization::SP3;
  if (s == "SPD5") return Hybridization::SPD5;
  if (s == "SPD6") return Hybridization::SPD6;
  if (s == "SPD7") return Hybridization::SPD7;
  if (s == "SPD8") return Hybridization::SPD8;
  return Hybridization::SP_NON;
}


void AcedrgTables::load_tables(const std::string& tables_dir, bool skip_angles) {
  using Clock = std::chrono::steady_clock;
  auto t0 = Clock::now();
  auto lap = [&](const char* label) {
    if (verbose >= 1) {
      auto t1 = Clock::now();
      std::fprintf(stderr, "  %-44s %6.3f s\n", label,
                   std::chrono::duration<double>(t1 - t0).count());
      t0 = t1;
    }
  };

  tables_dir_ = tables_dir;

  load_hash_codes(tables_dir + "/allOrgLinkedHashCode.table");
  lap("load_hash_codes");
  load_bond_hrs(tables_dir + "/allOrgBondsHRS.table");
  lap("load_bond_hrs");
  if (!skip_angles) {
    load_angle_hrs(tables_dir + "/allOrgAnglesHRS.table");
    lap("load_angle_hrs");
  }
  load_en_bonds(tables_dir + "/allOrgBondEN.table");
  lap("load_en_bonds");
  load_prot_hydr_dists(tables_dir + "/prot_hydr_dists.table");
  lap("load_prot_hydr_dists");
  load_metal_tables(tables_dir);
  lap("load_metal_tables");
  load_covalent_radii(tables_dir + "/radii.table");
  lap("load_covalent_radii");
  load_ccp4_bonds(tables_dir + "/ener_lib.cif");
  lap("load_ccp4_bonds");
  load_atom_type_codes(tables_dir + "/allAtomTypesFromMolsCoded.list");
  lap("load_atom_type_codes");
  load_bond_index(tables_dir + "/allOrgBondTables/bond_idx.table");
  lap("load_bond_index");
  load_bond_tables(tables_dir + "/allOrgBondTables");
  lap("load_bond_tables");
  if (!skip_angles) {
    load_angle_index(tables_dir + "/allOrgAngleTables/angle_idx.table");
    lap("load_angle_index");
    load_angle_tables(tables_dir + "/allOrgAngleTables");
    lap("load_angle_tables");
  }
  load_pep_tors(tables_dir + "/pep_tors.table");
  lap("load_pep_tors");

  tables_loaded_ = true;
}

void AcedrgTables::load_hash_codes(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    int hash_code, linked_hash;
    char footprint[64];

    if (std::sscanf(line, "%d %63s %d", &hash_code, footprint, &linked_hash) == 3) {
      digit_keys_[hash_code] = footprint;
      linked_hash_[hash_code] = linked_hash;
    }
  }
}

void AcedrgTables::load_bond_hrs(const std::string& path) {
  fileptr_t f = file_open(path.c_str(), "r");

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    int hash1, hash2, count;
    char hybrid_pair[64], in_ring[8];
    double value, sigma;

    if (std::sscanf(line, "%d %d %63s %7s %lf %lf %d",
                    &hash1, &hash2, hybrid_pair, in_ring, &value, &sigma, &count) == 7) {
      BondHRSKey key;
      key.hash1 = std::min(hash1, hash2);
      key.hash2 = std::max(hash1, hash2);
      key.hybrid_pair = hybrid_pair;
      key.in_ring = (in_ring[0] == 'Y' || in_ring[0] == 'y') ? "Y" : "N";
      CodStats vs(value, sigma, count);
      bond_hrs_[key] = vs;
      // AceDRG levels 9-11 use the HRS table hierarchy.
      bond_hasp_2d_[key.hash1][key.hash2][key.hybrid_pair][key.in_ring].push_back(vs);
      bond_hasp_1d_[key.hash1][key.hash2][key.hybrid_pair].push_back(vs);
      bond_hasp_0d_[key.hash1][key.hash2].push_back(vs);
    }
  }
}

void AcedrgTables::load_angle_hrs(const std::string& path) {
  fileptr_t f = file_open(path.c_str(), "r");

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    int hash1, hash2, hash3, count1, count2;
    char value_key[128], a1_cod[64], a2_cod[64], a3_cod[64];
    double value1, sigma1, value2, sigma2;

    if (std::sscanf(line, "%d %d %d %127s %63s %63s %63s %lf %lf %d %lf %lf %d",
                    &hash1, &hash2, &hash3, value_key, a1_cod, a2_cod, a3_cod,
                    &value1, &sigma1, &count1, &value2, &sigma2, &count2) == 13) {
      AngleHRSKey key;
      // File format: hash1=center, hash2/hash3=flanks. Canonicalize flanks.
      bool swap_flanks = hash2 > hash3;
      key.hash1 = hash1;  // center stays
      key.hash2 = swap_flanks ? hash3 : hash2;
      key.hash3 = swap_flanks ? hash2 : hash3;
      if (swap_flanks) {
        const char* colon = std::strchr(value_key, ':');
        if (colon) {
          std::string ring_part(value_key, colon - value_key);
          std::string hp(colon + 1);
          // Tuple format: center_flank1_flank2. Swap flanks: center_flank2_flank1
          size_t u1 = hp.find('_');
          size_t u2 = u1 != std::string::npos ? hp.find('_', u1 + 1) : std::string::npos;
          if (u2 != std::string::npos)
            hp = cat(hp.substr(0, u1), '_', hp.substr(u2 + 1), '_', hp.substr(u1 + 1, u2 - u1 - 1));
          key.value_key = cat(ring_part, ':', hp);
        } else {
          key.value_key = value_key;
        }
      } else {
        key.value_key = value_key;
      }
      auto it = angle_hrs_.find(key);
      if (it == angle_hrs_.end() || count1 > it->second.count)
        angle_hrs_[key] = CodStats(value1, sigma1, count1);
    }
  }
}

void AcedrgTables::load_en_bonds(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    char elem1[16], sp1[16], elem2[16], sp2[16];
    double value, sigma;
    int count;

    if (std::sscanf(line, "%15s %15s %15s %15s %lf %lf %d",
                    elem1, sp1, elem2, sp2, &value, &sigma, &count) == 7) {
      // Canonicalize order
      std::string e1(elem1), s1(sp1), e2(elem2), s2(sp2);
      if (e1 > e2 || (e1 == e2 && s1 > s2)) {
        std::swap(e1, e2);
        std::swap(s1, s2);
      }
      en_bonds_[e1][s1][e2][s2].emplace_back(value, sigma, count);
    }
  }
}

void AcedrgTables::load_prot_hydr_dists(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    // Remove stray commas that may appear in numeric fields
    for (char* p = line; *p; ++p)
      if (*p == ',') *p = ' ';

    char h_elem[16], heavy_elem[16], type_key[64];
    double nucleus_val, nucleus_sigma, v1, s1, electron_val, electron_sigma;

    // Format: H  C  H_sp3_C  1.092  0.010  1.093107  0.003875  0.988486  0.005251  ...
    // Column 4-5: nucleus distance (ener_lib-like)
    // Column 6-7: refined nucleus distance
    // Column 8-9: electron distance (X-ray) - this is what acedrg uses for value_dist
    if (std::sscanf(line, "%15s %15s %63s %lf %lf %lf %lf %lf %lf",
                    h_elem, heavy_elem, type_key,
                    &nucleus_val, &nucleus_sigma, &v1, &s1,
                    &electron_val, &electron_sigma) == 9) {
      ProtHydrDist& phd = prot_hydr_dists_[type_key];
      phd.electron_val = electron_val;
      phd.electron_sigma = electron_sigma;
      phd.nucleus_val = nucleus_val;
      phd.nucleus_sigma = nucleus_sigma;
    }
  }
}

void AcedrgTables::load_metal_tables(const std::string& dir) {
  // Load allMetalBonds.table
  {
    std::string path = dir + "/allMetalBonds.table";
    fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
    if (!f)
      return;

    char line[512];
    while (std::fgets(line, sizeof(line), f.get())) {
      if (is_skip_line(line))
        continue;

      char metal[8], ligand[8], ligand_class[64];
      int metal_coord, ligand_coord, pre_count, cnt;
      double pre_val, pre_sig, val, sig;

      // Try 11-field format first (with ligand_class)
      int n = std::sscanf(line, "%7s %d %7s %d %lf %lf %d %63s %lf %lf %d",
                          metal, &metal_coord, ligand, &ligand_coord,
                          &pre_val, &pre_sig, &pre_count,
                          ligand_class, &val, &sig, &cnt);
      if (n == 11 || n == 10) {
        if (n == 10) {
          // 10-field: the 8th field was parsed into ligand_class but is actually val
          std::sscanf(line, "%*s %*d %*s %*d %*f %*f %*d %lf %lf %d", &val, &sig, &cnt);
        }
        MetalBondEntry entry;
        entry.metal = Element(metal);
        entry.metal_coord = metal_coord;
        entry.ligand = Element(ligand);
        entry.ligand_coord = ligand_coord;
        entry.pre_value = pre_val;
        entry.pre_sigma = pre_sig;
        entry.pre_count = pre_count;
        entry.ligand_class = (n == 11) ? ligand_class : "NONE";
        entry.value = val;
        entry.sigma = sig;
        entry.count = cnt;
        metal_bonds_.push_back(entry);
      }
    }
  }

  // Load metal coordination geometry
  {
    std::string path = dir + "/allMetalDefCoordGeos.table";
    fileptr_t f2(std::fopen(path.c_str(), "r"), needs_fclose{true});
    if (f2) {
      char line[512];
      while (std::fgets(line, sizeof(line), f2.get())) {
        if (is_skip_line(line))
          continue;

        char metal_str[8], geo_str[64];
        int coord;

        if (std::sscanf(line, "%7s %d %63s", metal_str, &coord, geo_str) == 3) {
          Element metal(metal_str);
          static const struct { const char* name; CoordGeometry geo; } geo_names[] = {
            {"LINEAR", CoordGeometry::LINEAR},
            {"TRIGONAL_PLANAR", CoordGeometry::TRIGONAL_PLANAR},
            {"T_SHAPED", CoordGeometry::T_SHAPED},
            {"TETRAHEDRAL", CoordGeometry::TETRAHEDRAL},
            {"SQUARE_PLANAR", CoordGeometry::SQUARE_PLANAR},
            {"TRIGONAL_BIPYRAMIDAL", CoordGeometry::TRIGONAL_BIPYRAMIDAL},
            {"SQUARE_PYRAMIDAL", CoordGeometry::SQUARE_PYRAMIDAL},
            {"OCTAHEDRAL", CoordGeometry::OCTAHEDRAL},
            {"TRIGONAL_PRISM", CoordGeometry::TRIGONAL_PRISM},
            {"PENTAGONAL_BIPYRAMIDAL", CoordGeometry::PENTAGONAL_BIPYRAMIDAL},
            {"CAPPED_OCTAHEDRAL", CoordGeometry::CAPPED_OCTAHEDRAL},
            {"SQUARE_ANTIPRISM", CoordGeometry::SQUARE_ANTIPRISM},
          };
          CoordGeometry geo = CoordGeometry::UNKNOWN;
          for (const auto& g : geo_names)
            if (std::strcmp(geo_str, g.name) == 0) { geo = g.geo; break; }
          metal_coord_geo_[metal][coord] = geo;
        }
      }
    }
  }

  // Load metal coordination angles
  {
    std::string path = dir + "/allMetalCoordGeoAngles.table";
    fileptr_t f3(std::fopen(path.c_str(), "r"), needs_fclose{true});
    if (f3) {
      char line[512];
      while (std::fgets(line, sizeof(line), f3.get())) {
        if (is_skip_line(line))
          continue;

        char metal_str[8], geo_str[64];
        int coord;
        double angle, sigma;

        if (std::sscanf(line, "%7s %d %63s %lf %lf",
                        metal_str, &coord, geo_str, &angle, &sigma) == 5) {
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
}

void AcedrgTables::load_covalent_radii(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return;

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;
    char elem[16], kind[16];
    double value = NAN;
    if (std::sscanf(line, "%15s %15s %lf", elem, kind, &value) != 3)
      continue;
    if (std::strcmp(kind, "cova") != 0)
      continue;
    covalent_radii_[to_upper(elem)] = value;
  }
}

void AcedrgTables::load_atom_type_codes(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    // Format: code<whitespace>full_type
    char code[64], full_type[256];
    if (std::sscanf(line, "%63s %255s", code, full_type) == 2) {
      atom_type_codes_[code] = full_type;
    }
  }
}

void AcedrgTables::load_bond_index(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    // Format: ha1 ha2 fileNum
    int ha1, ha2, file_num;
    if (std::sscanf(line, "%d %d %d", &ha1, &ha2, &file_num) == 3) {
      bond_file_index_[ha1][ha2] = file_num;
    }
  }
}

void AcedrgTables::load_bond_tables(const std::string& dir) {
  // Load each bond table file referenced in the index
  std::set<int> loaded_files;
  int n_files = 0, n_lines = 0;

  for (const auto& ha1_pair : bond_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      int file_num = ha2_pair.second;
      if (!loaded_files.insert(file_num).second)
        continue;

      std::string path = cat(dir, '/', file_num, ".table");
      fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
      if (!f)
        continue;
      ++n_files;

      char line[512];
      while (std::fgets(line, sizeof(line), f.get())) {
        if (is_skip_line(line))
          continue;

        // Format: ha1 ha2 hybrComb inRing a1NB2 a2NB2 a1NB a2NB atomCode1 atomCode2
        //         value sigma count val2 sig2 cnt2
        const char* p = line;
        int ha1 = simple_atoi(p, &p);
        int ha2 = simple_atoi(p, &p);
        const char* s;
        s = skip_blank(p); p = skip_word(s); std::string hybr_comb(s, p);
        s = skip_blank(p); p = skip_word(s); std::string in_ring(s, p);
        s = skip_blank(p); p = skip_word(s); std::string a1_nb2(s, p);
        s = skip_blank(p); p = skip_word(s); std::string a2_nb2(s, p);
        s = skip_blank(p); p = skip_word(s); std::string a1_nb(s, p);
        s = skip_blank(p); p = skip_word(s); std::string a2_nb(s, p);
        s = skip_blank(p); p = skip_word(s);
        const char* code1_s = s; const char* code1_e = p;
        s = skip_blank(p); p = skip_word(s);
        if (s == p)  // not enough fields
          continue;
        const char* code2_s = s; const char* code2_e = p;
        double value = fast_atof(p, &p);
        double sigma = fast_atof(p, &p);
        int count = simple_atoi(p, &p);
        double value2 = fast_atof(p, &p);
        double sigma2 = fast_atof(p, &p);
        int count2 = simple_atoi(p, &p);

        ++n_lines;
        // Get atom types from codes: full, main (before '{'), root (before '(')
        auto* p1 = find_val(atom_type_codes_, std::string(code1_s, code1_e));
        auto* p2 = find_val(atom_type_codes_, std::string(code2_s, code2_e));
        std::string a1_type_f = p1 ? *p1 : std::string();
        std::string a2_type_f = p2 ? *p2 : std::string();
        std::string a1_type_m = prefix_before(a1_type_f, '{');
        std::string a2_type_m = prefix_before(a2_type_f, '{');
        std::string a1_root = prefix_before(a1_type_m, '(');
        std::string a2_root = prefix_before(a2_type_m, '(');

        CodStats vs(value, sigma, count);
        CodStats vs1d(value2, sigma2, count2);

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

        // Levels 9-11 are populated from allOrgBondsHRS.table in load_bond_hrs().
      }
    }
  }
  if (verbose >= 1)
    std::fprintf(stderr, "    bond tables: %d files, %d lines\n", n_files, n_lines);
}

void AcedrgTables::load_angle_index(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return; // Optional file

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;

    // Format: ha1 ha2 ha3 fileNum
    int ha1, ha2, ha3, file_num;
    if (std::sscanf(line, "%d %d %d %d", &ha1, &ha2, &ha3, &file_num) == 4) {
      angle_file_index_[ha1][ha2][ha3] = file_num;
    }
  }
}

void AcedrgTables::load_angle_tables(const std::string& dir) {
  // Load each angle table file referenced in the index
  std::set<int> loaded_files;
  int n_files = 0, n_lines = 0;

  for (const auto& ha1_pair : angle_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      for (const auto& ha3_pair : ha2_pair.second) {
        int file_num = ha3_pair.second;
        if (!loaded_files.insert(file_num).second)
          continue;

        std::string path = cat(dir, '/', file_num, ".table");
        fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
        if (!f)
          continue;
        ++n_files;

        char line[1024];
        while (std::fgets(line, sizeof(line), f.get())) {
          if (is_skip_line(line))
            continue;

          // 34-column format:
          // 1-3: ha1 ha2 ha3
          // 4: valueKey (ring:hybr_tuple, e.g. "0:SP2_SP2_SP3")
          // 5-7: a1_root a2_root a3_root
          // 8-10: a1_nb2 a2_nb2 a3_nb2
          // 11-13: a1_nb a2_nb a3_nb
          // 14-16: a1_code a2_code a3_code
          // 17-34: 6 sets of value sigma count
          const char* p = line;
          int ha1 = simple_atoi(p, &p);
          int ha2 = simple_atoi(p, &p);
          int ha3 = simple_atoi(p, &p);
          const char* s;
          s = skip_blank(p); p = skip_word(s); std::string value_key(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a1_root(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a2_root(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a3_root(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a1_nb2(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a2_nb2(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a3_nb2(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a1_nb(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a2_nb(s, p);
          s = skip_blank(p); p = skip_word(s); std::string a3_nb(s, p);
          s = skip_blank(p); p = skip_word(s);
          const char* code1_s = s; const char* code1_e = p;
          s = skip_blank(p); p = skip_word(s);
          const char* code2_s = s; const char* code2_e = p;
          s = skip_blank(p); p = skip_word(s);
          if (s == p)  // not enough fields
            continue;
          const char* code3_s = s; const char* code3_e = p;

          // Read 6 sets of value/sigma/count
          double values[6], sigmas[6];
          int counts[6];
          for (int lvl = 0; lvl < 6; ++lvl) {
            values[lvl] = fast_atof(p, &p);
            sigmas[lvl] = fast_atof(p, &p);
            counts[lvl] = simple_atoi(p, &p);
          }

          ++n_lines;
          // Get main atom types (before '{') from codes
          auto* q1 = find_val(atom_type_codes_, std::string(code1_s, code1_e));
          auto* q2 = find_val(atom_type_codes_, std::string(code2_s, code2_e));
          auto* q3 = find_val(atom_type_codes_, std::string(code3_s, code3_e));
          std::string a1_type = q1 ? prefix_before(*q1, '{') : std::string();
          std::string a2_type = q2 ? prefix_before(*q2, '{') : std::string();
          std::string a3_type = q3 ? prefix_before(*q3, '{') : std::string();

          // Populate structures at each level with corresponding pre-computed values.
          // AceDRG keeps only the first entry for each key (no aggregation).
          // Level 1D: full detail with atom types
          CodStats vs1(values[1], sigmas[1], counts[1]);
          auto& angle_1d_vec = angle_idx_1d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb]
                               [a1_type][a2_type][a3_type];
          if (angle_1d_vec.empty())
            angle_1d_vec.push_back(vs1);

          // Level 2D: no atom types
          CodStats vs2(values[2], sigmas[2], counts[2]);
          auto& angle_2d_vec = angle_idx_2d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb];
          if (angle_2d_vec.empty())
            angle_2d_vec.push_back(vs2);

          // Level 3D: hash + valueKey + NB2 only
          CodStats vs3(values[3], sigmas[3], counts[3]);
          auto& angle_3d_vec = angle_idx_3d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2];
          if (angle_3d_vec.empty())
            angle_3d_vec.push_back(vs3);

          // Level 4D: hash + valueKey + roots
          CodStats vs4(values[4], sigmas[4], counts[4]);
          auto& angle_4d_vec = angle_idx_4d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root];
          if (angle_4d_vec.empty())
            angle_4d_vec.push_back(vs4);

          // Level 5D: hash + valueKey only
          CodStats vs5(values[5], sigmas[5], counts[5]);
          auto& angle_5d_vec = angle_idx_5d_[ha1][ha2][ha3][value_key];
          if (angle_5d_vec.empty())
            angle_5d_vec.push_back(vs5);

          // Level 6D: hash only (leave empty for 34-column data)
        }
      }
    }
  }
  if (verbose >= 1)
    std::fprintf(stderr, "    angle tables: %d files, %d lines\n", n_files, n_lines);
}

void AcedrgTables::load_pep_tors(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "r"), needs_fclose{true});
  if (!f)
    return;

  char line[512];
  while (std::fgets(line, sizeof(line), f.get())) {
    if (is_skip_line(line))
      continue;
    char tors_id[64], label[64], a1[16], a2[16], a3[16], a4[16];
    int period = 0, idx = 0;
    double value = 0.0;
    if (std::sscanf(line, "%63s %63s %15s %15s %15s %15s %d %d %lf",
                    tors_id, label, a1, a2, a3, a4, &period, &idx, &value) != 9)
      continue;
    TorsionEntry entry;
    entry.value = value;
    entry.period = period;
    entry.priority = idx;
    entry.id = label;
    pep_tors_[cat(a1, '_', a2, '_', a3, '_', a4)] = std::move(entry);
  }
}

namespace {

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

struct NB1stFam {
  std::string name;
  std::vector<std::string> NB2ndList;
  int repN = 1;
};

std::vector<std::vector<BondInfo>>
build_adjacency(const ChemComp& cc) {
  std::vector<std::vector<BondInfo>> adj(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom);
    int idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 >= 0 && idx2 >= 0) {
      adj[idx1].push_back({idx2, bond.type, bond.aromatic});
      adj[idx2].push_back({idx1, bond.type, bond.aromatic});
    }
  }
  return adj;
}

std::vector<std::vector<int>> build_neighbors(
    const std::vector<std::vector<BondInfo>>& adj) {
  std::vector<std::vector<int>> neighbors(adj.size());
  for (size_t i = 0; i < adj.size(); ++i) {
    neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i])
      neighbors[i].push_back(nb.neighbor_idx);
  }
  return neighbors;
}

void set_ring_aromaticity_from_bonds(
    const std::vector<std::vector<BondInfo>>& adj,
    const std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings,
    int verbose) {
  // AceDRG has two-phase aromaticity:
  // - Strict (mode 0): only NoMetal pi count → isAromatic (used for COD table lookup)
  // - Permissive (mode 1): NoMetal+All pi counts → isAromaticP (used for output CIF)
  // Ring must be planar (bondingIdx==2, N allowed).
  auto count_non_mc = [&](int idx) -> int {
    int non_mc = 0;
    for (const auto& nb : adj[idx]) {
      const auto& at = atoms[nb.neighbor_idx];
      if (!at.is_metal)
        ++non_mc;
    }
    return non_mc;
  };

  auto is_atom_planar = [&](int idx) -> bool {
    const auto& atom = atoms[idx];
    if (atom.el == El::N)
      return true;
    return atom.bonding_idx == 2;
  };

  auto is_ring_planar = [&](const RingInfo& ring) -> bool {
    for (int idx : ring.atoms)
      if (!is_atom_planar(idx))
        return false;
    return true;
  };

  // Count pi electrons for an atom in a ring.
  // mode 0 (strict): C(-1)/non_mc=2 and N(+1)/non_mc=3 give 2 pi.
  // mode 1 (permissive): those same cases give 1 pi.
  // include_c_minus2: if false, skip the C charge=-2 case (used for "all" count).
  auto count_atom_pi = [&](int idx, int mode, bool include_c_minus2) -> int {
    const auto& atom = atoms[idx];
    int non_mc = count_non_mc(idx);
    int aN = 0;

    if (atom.bonding_idx == 2) {
      if (atom.charge == 0.0f) {
        if (atom.el == El::C) {
          if (non_mc == 3) {
            bool has_exo = false;
            for (int nb : atom.conn_atoms_no_metal) {
              if (atoms[nb].el == El::O &&
                  atoms[nb].conn_atoms_no_metal.size() == 1 &&
                  atoms[nb].charge == 0.0f) {
                has_exo = true;
              } else if (atoms[nb].el == El::C &&
                         atoms[nb].conn_atoms_no_metal.size() == 3) {
                int h_count = 0;
                for (int nb2 : atoms[nb].conn_atoms_no_metal)
                  if (atoms[nb2].el == El::H)
                    ++h_count;
                if (h_count >= 2)
                  has_exo = true;
              }
            }
            if (!has_exo)
              aN = 1;
          }
        } else if (atom.el == El::N) {
          if (non_mc == 2)
            aN = 1;
          else if (non_mc == 3)
            aN = 2;
        } else if (atom.el == El::B) {
          if (non_mc == 2)
            aN = 1;
        } else if (atom.el == El::O || atom.el == El::S) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.el == El::P) {
          if (non_mc == 3)
            aN = 2;
        }
      } else {
        if (atom.el == El::C) {
          if (atom.charge == -1.0f) {
            if (non_mc == 3)
              aN = 2;
            else if (non_mc == 2)
              aN = (mode == 1) ? 1 : 2;
          } else if (include_c_minus2 && atom.charge == -2.0f) {
            if (non_mc == 2)
              aN = 2;
          }
        } else if (atom.el == El::N) {
          if (atom.charge == -1.0f) {
            if (non_mc == 2)
              aN = 2;
          } else if (atom.charge == 1.0f) {
            if (non_mc == 3)
              aN = (mode == 1) ? 1 : 2;
          }
        } else if (atom.el == El::O) {
          if (atom.charge == 1.0f && non_mc == 2)
            aN = 1;
        } else if (atom.el == El::B) {
          if (atom.charge == -1.0f && non_mc == 3)
            aN = 1;
        }
      }
    } else if (atom.bonding_idx == 3 &&
               (atom.el == El::N || atom.el == El::B)) {
      if (atom.el == El::N) {
        if (atom.charge == -1.0f) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.charge == 1.0f) {
          if (non_mc == 3)
            aN = 1;
        } else {
          aN = 2;
        }
      } else if (atom.el == El::B) {
        aN = 0;
      }
    }

    return aN;
  };

  for (size_t i = 0; i < rings.size(); ++i) {
    RingInfo& ring = rings[i];
    ring.is_aromatic = false;

    if (!is_ring_planar(ring))
      continue;

    // AceDRG uses strict aromaticity (mode 0) for the COD table lookup.
    // Mode 0 differs from mode 1 for charged N+ with 3 connections:
    // mode 0 gives 2 pi electrons, mode 1 gives 1.
    // AceDRG only checks the NoMetal pi count in strict mode.
    int pi1 = 0;
    for (int idx : ring.atoms)
      pi1 += count_atom_pi(idx, 0, true);
    if (pi1 > 0 && pi1 % 4 == 2)
      ring.is_aromatic = true;
    if (verbose >= 2) {
      std::fprintf(stderr, "    ring %zu (size=%zu): pi1=%d aromatic=%d atoms:",
                   i, ring.atoms.size(), pi1, ring.is_aromatic ? 1 : 0);
      for (int idx : ring.atoms)
        std::fprintf(stderr, " %s", atoms[idx].id.c_str());
      std::fprintf(stderr, "\n");
    }
  }

  // AceDRG pyrole rule: if there are exactly 4 five-member rings with 4C+1N
  // and they are planar, mark them aromatic even if pi-count failed.
  std::vector<size_t> pyrole_rings;
  for (size_t i = 0; i < rings.size(); ++i) {
    const RingInfo& ring = rings[i];
    if (ring.atoms.size() != 5)
      continue;
    int num_c = 0;
    int num_n = 0;
    for (int idx : ring.atoms) {
      if (atoms[idx].el == El::C)
        ++num_c;
      else if (atoms[idx].el == El::N)
        ++num_n;
    }
    if (num_c == 4 && num_n == 1)
      pyrole_rings.push_back(i);
  }
  if (pyrole_rings.size() == 4) {
    for (size_t i : pyrole_rings) {
      if (is_ring_planar(rings[i]))
        rings[i].is_aromatic = true;
    }
  }

  for (auto& ring : rings) {
    ring.is_aromatic_permissive = ring.is_aromatic;
    if (ring.is_aromatic)
      continue;
    if (!is_ring_planar(ring))
      continue;
    int pi1 = 0;
    int pi2 = 0;
    for (int idx : ring.atoms) {
      pi1 += count_atom_pi(idx, 1, true);
      pi2 += count_atom_pi(idx, 1, false);
    }
    if ((pi1 > 0 && pi1 % 4 == 2) || (pi2 > 0 && pi2 % 4 == 2))
      ring.is_aromatic_permissive = true;
  }
}

void check_one_path_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings,
    std::map<std::string, int>& ring_index,
    int ori_idx, int cur_idx, int prev_idx, int cur_lev, int max_ring,
    std::map<int, std::string>& seen_atom_ids,
    std::map<int, std::string>& atom_ids_in_path);

void detect_rings_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) {
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

void check_one_path_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings,
    std::map<std::string, int>& ring_index,
    int ori_idx, int cur_idx, int prev_idx, int cur_lev, int max_ring,
    std::map<int, std::string>& seen_atom_ids,
    std::map<int, std::string>& atom_ids_in_path) {
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
        std::vector<std::string> all_ids;
        std::vector<std::string> all_seris;
        std::vector<int> ring_atoms;
        for (const auto& it : atom_ids_in_path) {
          all_seris.push_back(std::to_string(it.first));
          all_ids.push_back(it.second);
          ring_atoms.push_back(it.first);
        }
        std::sort(all_seris.begin(), all_seris.end(), compare_no_case);
        std::sort(all_ids.begin(), all_ids.end(), compare_no_case);

        std::string rep;
        for (const auto& id : all_ids)
          rep += id;

        std::string s_rep = join_str(all_seris, '_');

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

void set_atoms_ring_rep_s(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<RingInfo>& rings) {
  for (const auto& ring : rings) {
    std::string size = std::to_string(ring.atoms.size());
    std::vector<std::string> all_seris;
    for (int idx : ring.atoms)
      all_seris.push_back(std::to_string(idx));
    std::sort(all_seris.begin(), all_seris.end(), compare_no_case);
    std::string rep_id = join_str(all_seris, '_');

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
void build_ring_size_map(const std::map<std::string, std::string>& ring_rep_s,
                         std::map<std::string, int>& size_map) {
  for (const auto& it : ring_rep_s) {
    const std::string& ring_str = it.second;
    size_map[ring_str] += 1;
  }
}

// Append ring size annotation like "[5a,6a]" to string s.
void append_ring_annotation(std::string& s,
                            const std::map<std::string, std::string>& ring_rep_s) {
  std::map<std::string, int> size_map;
  build_ring_size_map(ring_rep_s, size_map);
  s += '[';
  int i = 0;
  int j = static_cast<int>(size_map.size());
  for (const auto& it : size_map) {
    const std::string& size = it.first;
    if (it.second >= 3)
      cat_to(s, it.second, 'x', size);
    else if (it.second == 2)
      cat_to(s, size, ',', size);
    else
      s.append(size);
    s += (i != j - 1) ? ',' : ']';
    ++i;
  }
}

int get_num_oxy_connect(const std::vector<CodAtomInfo>& atoms,
                        const CodAtomInfo& atom,
                        const std::vector<std::vector<int>>& neighbors) {
  int nO = 0;
  for (int nb : neighbors[atom.index])
    if (atoms[nb].el == El::O)
      nO++;
  return nO;
}

// Parse the first ring size from a COD class like "C[5a,6a](...)".
// Extracts from the bracket content: first token before ',' or ']',
// handling "NxSIZE" multiplier format.
int get_min_ring2_from_cod_class(const std::string& cod_class) {
  // Find bracket content before '('
  size_t bracket = cod_class.find('[');
  if (bracket == std::string::npos)
    return 0;
  size_t paren = cod_class.find('(');
  if (paren != std::string::npos && paren < bracket)
    return 0;
  // Extract first token: from '[' to first ',' or ']'
  size_t end = cod_class.find_first_of(",]", bracket + 1);
  if (end == std::string::npos)
    return 0;
  std::string token = cod_class.substr(bracket + 1, end - bracket - 1);
  // Handle "NxSIZE" format — take part after 'x'
  size_t x = token.find('x');
  if (x != std::string::npos)
    token = token.substr(x + 1);
  return string_to_int(token, false);
}

bool cod_class_is_aromatic(const std::string& cod_class) {
  // Check if the ring annotation (between '[' and ']', before any '(') contains 'a'
  size_t paren = cod_class.find('(');
  size_t bracket = cod_class.find('[');
  if (bracket == std::string::npos || (paren != std::string::npos && bracket > paren))
    return false;
  size_t end = cod_class.find(']', bracket + 1);
  if (end == std::string::npos)
    end = cod_class.size();
  return cod_class.find('a', bracket + 1) < end;
}

void get_small_family(const std::string& in_str, NB1stFam& fam) {
  fam.name.clear();
  fam.NB2ndList.clear();
  std::vector<std::string> ch_list = {",", "x"};
  std::string name_str;
  bool l_r = false;
  int n_rep = 1;
  for (size_t i = 0; i < in_str.size(); ++i) {
    char c = in_str[i];
    if (std::isalpha(static_cast<unsigned char>(c))) {
      if (alpha_up(c) == c) {
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
    } else if (in_vector(std::string(1, c), ch_list)) {
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

bool are_in_same_ring(const CodAtomInfo& a1, const CodAtomInfo& a2) {
  for (int ring_idx : a1.in_rings)
    if (in_vector(ring_idx, a2.in_rings))
      return true;
  return false;
}

int angle_ring_size(const CodAtomInfo& center,
                    const CodAtomInfo& a1,
                    const CodAtomInfo& a3) {
  for (const auto& it : center.ring_rep) {
    if (a1.ring_rep.find(it.first) != a1.ring_rep.end() &&
        a3.ring_rep.find(it.first) != a3.ring_rep.end())
      return it.second;
  }
  return 0;
}

// Order two atoms by: hashing_value, then cod_main (longer first, then
// case-insensitive), then atom id as tiebreaker.
void order_two_atoms(const CodAtomInfo& a1, const CodAtomInfo& a2,
                     const CodAtomInfo*& first,
                     const CodAtomInfo*& second) {
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
  } else if (compare_no_case2(a2.cod_main, a1.cod_main)) {
    first = &a2;
    second = &a1;
  } else if (!compare_no_case2(a1.id, a2.id)) {
    first = &a1;
    second = &a2;
  } else {
    first = &a2;
    second = &a1;
  }
}

Hybridization hybrid_from_bonding_idx(int bonding_idx,
                                      bool is_metal_atom,
                                      int connectivity) {
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


void set_atom_cod_class_name_new2(
    CodAtomInfo& atom, const CodAtomInfo& ori_atom, int lev,
    const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) {
  if (lev == 1) {
    atom.cod_class.clear();
    atom.cod_class.append(atom.el.name());

    if (!atom.ring_rep_s.empty())
      append_ring_annotation(atom.cod_class, atom.ring_rep_s);

    std::string t_str;
    std::vector<std::string> t_str_list;
    std::map<std::string, int> comps;
    for (int nb : neighbors[atom.index]) {
      if (nb == ori_atom.index || atoms[nb].is_metal)
        continue;
      std::string nb_type = atoms[nb].el.name();
      if (!atoms[nb].ring_rep_s.empty())
        append_ring_annotation(nb_type, atoms[nb].ring_rep_s);
      comps[nb_type] += 1;
    }

    std::vector<SortMap> sorted;
    for (const auto& it : comps) {
      SortMap sm;
      sm.key = it.first;
      sm.val = it.second;
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), desc_sort_map_key);
    for (const auto& sm : sorted) {
      std::string s1 = cat(sm.key, sm.val);
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

    if (!atom.ring_rep_s.empty())
      append_ring_annotation(atom.cod_class, atom.ring_rep_s);

    int low_lev = lev - 1;
    std::map<std::string, std::vector<int>> id_map;
    for (int nb : neighbors[atom.index]) {
      if (atoms[nb].is_metal)
        continue;
      CodAtomInfo nb_atom = atoms[nb];
      set_atom_cod_class_name_new2(nb_atom, ori_atom, low_lev, atoms, neighbors);
      auto& entry = id_map[nb_atom.cod_class];
      if (entry.empty()) {
        entry.push_back(1);
        entry.push_back((int)std::count_if(neighbors[nb].begin(), neighbors[nb].end(),
                                           [&](int nb2) { return !atoms[nb2].is_metal; }));
      } else {
        entry[0] += 1;
      }
    }

    std::vector<SortMap> sorted;
    for (const auto& it : id_map) {
      SortMap sm;
      sm.key = it.first;
      sm.val = it.second[0];
      sm.nNB = it.second[1];
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), desc_sort_map_key);
    for (const auto& sm : sorted) {
      if (sm.val == 1)
        cat_to(atom.cod_class, '(', sm.key, ')');
      else
        cat_to(atom.cod_class, '(', sm.key, ')', sm.val);
    }
  }
}

void set_special_3nb_symb2(
    CodAtomInfo& atom, const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) {
  if (atom.ring_rep.empty())
    return;

  std::vector<int> ser_num_nb123;
  std::map<std::string, int> nb3_props;

  for (int nb1 : neighbors[atom.index]) {
    if (atoms[nb1].is_metal)
      continue;
    if (!in_vector(nb1, ser_num_nb123))
      ser_num_nb123.push_back(nb1);
    for (int nb2 : neighbors[nb1]) {
      if (atoms[nb2].is_metal)
        continue;
      if (!in_vector(nb2, ser_num_nb123) &&
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
        if (atoms[nb3].is_metal)
          continue;
        if (!in_vector(nb3, ser_num_nb123) &&
            nb3 != atom.index) {
          std::string prop = atoms[nb3].el.name();
          int deg = 0;
          for (int nbx : neighbors[nb3])
            if (!atoms[nbx].is_metal)
              ++deg;
          cat_to(prop, '<', deg, '>');
          nb3_props[prop] += 1;
          ser_num_nb123.push_back(nb3);
        }
      }
    }
  }

  std::vector<std::string> comps;
  for (const auto& it : nb3_props) {
    std::string id = cat(it.second, '|', it.first);
    comps.push_back(id);
  }
  std::sort(comps.begin(), comps.end(), compare_no_case2);

  if (!comps.empty()) {
    atom.cod_class += '{';
    for (size_t i = 0; i < comps.size(); ++i) {
      if (i > 0)
        atom.cod_class += ',';
      atom.cod_class.append(comps[i]);
    }
    atom.cod_class += '}';
  }
}

void cod_class_to_atom2(const std::string& cod_class, CodAtomInfo& atom) {
  std::string t_cod = trim_str(cod_class);
  atom.cod_class = t_cod;
  atom.nb_symb.clear();
  atom.nb2_symb.clear();
  atom.nb3_symb.clear();

  size_t brace = t_cod.find('{');
  if (brace != std::string::npos) {
    atom.cod_main = t_cod.substr(0, brace);
    size_t brace_end = t_cod.find('}', brace + 1);
    if (brace_end != std::string::npos)
      atom.nb3_symb = t_cod.substr(brace + 1, brace_end - brace - 1);
  } else {
    atom.cod_main = t_cod;
  }

  std::vector<std::string> atm_strs = split_str(atom.cod_main, '(');
  if (!atm_strs.empty()) {
    atom.cod_root = trim_str(atm_strs[0]);
  }

  std::vector<NB1stFam> all_nbs;
  for (size_t i = 1; i < atm_strs.size(); ++i) {
    std::string tS = trim_str(atm_strs[i]);
    std::vector<std::string> nb1 = split_str(tS, ')');
    NB1stFam fam;
    if (nb1.size() > 1) {
      fam.repN = string_to_int(nb1[1], false);
      if (fam.repN == 0)
        fam.repN = 1;
    } else {
      fam.repN = 1;
    }

    std::string tS1 = trim_str(nb1[0]);
    get_small_family(tS1, fam);
    all_nbs.push_back(fam);
  }

  for (const auto& fam : all_nbs) {
    for (int j = 0; j < fam.repN; ++j) {
      int sN = static_cast<int>(fam.NB2ndList.size()) + 1;
      cat_to(atom.nb_symb, fam.name, '-', sN, ":");
      cat_to(atom.nb2_symb, sN, ":");
    }
  }
}

void set_atoms_nb1nb2_sp(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) {
  for (auto& atom : atoms) {
    std::vector<std::string> nb1_nb2_sp_set;
    for (int nb1 : neighbors[atom.index]) {
      if (atoms[nb1].is_metal)
        continue;
      std::string nb1_main = atoms[nb1].cod_root;
      std::vector<int> nb2_sp_set;
      for (int nb2 : neighbors[nb1]) {
        if (atoms[nb2].is_metal)
          continue;
        nb2_sp_set.push_back(atoms[nb2].bonding_idx);
      }
      std::sort(nb2_sp_set.begin(), nb2_sp_set.end(), std::greater<int>());
      std::string nb2_sp_str = join_str(nb2_sp_set, '_',
          [](int v) { return std::to_string(v); });
      nb1_nb2_sp_set.emplace_back(cat(nb1_main, '-', nb2_sp_str));
    }
    // Sort alphabetically by the string (same order as AceDRG tables)
    std::sort(nb1_nb2_sp_set.begin(), nb1_nb2_sp_set.end(),
              [](const auto& a, const auto& b) {
                return compare_no_case(a, b);
              });
    // Build nb1nb2_sp in alphabetical order
    atom.nb1nb2_sp.clear();
    for (size_t i = 0; i < nb1_nb2_sp_set.size(); ++i) {
      if (i > 0)
        atom.nb1nb2_sp += ':';
      atom.nb1nb2_sp += nb1_nb2_sp_set[i];
    }
  }
}

void set_atoms_bonding_and_chiral_center(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) {
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
      // Use t_len (non-metal connection count) like AceDRG, which
      // excludes metal bonds from connAtoms before hybridization assignment.
      if (t_len == 2) {
        atom.bonding_idx = 3;
      } else if (t_len == 1) {
        // Find the first non-metal neighbor and check its non-metal connections
        int first_nb = -1;
        for (int nb : neighbors[atom.index])
          if (!atoms[nb].is_metal) { first_nb = nb; break; }
        int nb_non_metal = 0;
        if (first_nb >= 0)
          for (int nnb : neighbors[first_nb])
            if (!atoms[nnb].is_metal)
              nb_non_metal++;
        if (nb_non_metal != 1)
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
    } else if (atom.el == El::H || atom.el == El::D) {
      atom.bonding_idx = 0;  // H has SP0 hybridization (no p orbitals)
    }
  }

  auto non_metal_count = [&](int idx) {
    return (int)std::count_if(neighbors[idx].begin(), neighbors[idx].end(),
                              [&](int nb) { return !atoms[nb].is_metal; });
  };

  for (auto& atom : atoms) {
    if (atom.el == El::O) {
      if (non_metal_count(atom.index) == 2 && atom.par_charge == 0.0f) {
        if (std::any_of(neighbors[atom.index].begin(), neighbors[atom.index].end(),
                        [&](int nb) { return atoms[nb].bonding_idx == 2; }))
          atom.bonding_idx = 2;
      }
    }
  }

  std::map<int, int> pre_bonding;
  for (const auto& atom : atoms)
    pre_bonding[atom.index] = atom.bonding_idx;

  auto has_sp2_nb_not_O = [&](int idx) {
    return std::any_of(neighbors[idx].begin(), neighbors[idx].end(),
                       [&](int nb) { return pre_bonding[nb] == 2 && atoms[nb].el != El::O; });
  };

  for (auto& atom : atoms) {
    int t_len = non_metal_count(atom.index);
    if (atom.el == El::N || atom.el == El::As) {
      if (t_len == 3) {
        if (atom.charge == 0.0f) {
          if (has_sp2_nb_not_O(atom.index)) {
            atom.bonding_idx = num_conn_map[atom.index][1] != 0 ? 3 : 2;
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
        if (has_sp2_nb_not_O(atom.index))
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

}  // namespace

// ============================================================================
// Implementation - Atom classification
// ============================================================================

std::vector<CodAtomInfo> AcedrgTables::classify_atoms(const ChemComp& cc) const {
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

  // AceDRG adjusts charges for atoms bonded to metals via valence bookkeeping.
  // Use non-metal bond orders to recompute charge for such atoms so that
  // bonding_idx and aromaticity match AceDRG behavior (e.g., AIV N2).
  for (size_t i = 0; i < atoms.size(); ++i) {
    CodAtomInfo& atom = atoms[i];
    if (atom.is_metal || atom.el == El::H)
      continue;
    bool has_metal_neighbor = atom.connectivity > static_cast<int>(atom.conn_atoms_no_metal.size());
    if (!has_metal_neighbor)
      continue;
    if (std::none_of(atom.conn_atoms_no_metal.begin(), atom.conn_atoms_no_metal.end(),
                      [&](int nb) { return atoms[nb].el != El::H; }))
      continue;
    int expected_valence = get_expected_valence(atom.el);
    if (expected_valence == 0)
      continue;
    float sum_bo = 0.0f;
    for (const BondInfo& bi : adj[i]) {
      if (!atoms[bi.neighbor_idx].is_metal)
        sum_bo += order_of_bond_type(bi.type);
    }
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    atom.charge = static_cast<float>(-rem_v);
  }

  // Bonding/planarity info is needed for AceDRG ring aromaticity rules.
  set_atoms_bonding_and_chiral_center(atoms, neighbors);

  // Detect rings and populate ring representations
  std::vector<RingInfo> rings;
  detect_rings_acedrg(neighbors, atoms, rings);
  set_ring_aromaticity_from_bonds(adj, atoms, rings, verbose);
  set_atoms_ring_rep_s(atoms, rings);

  // Build COD class names (AceDRG style)
  for (size_t i = 0; i < atoms.size(); ++i) {
    set_atom_cod_class_name_new2(atoms[i], atoms[i], 2, atoms, neighbors);
    set_special_3nb_symb2(atoms[i], atoms, neighbors);
    cod_class_to_atom2(atoms[i].cod_class, atoms[i]);
  }

  // Hybridization and NB1/NB2_SP
  set_atoms_nb1nb2_sp(atoms, neighbors);

  // Ring props from codClass and hash codes
  for (auto& atom : atoms) {
    atom.min_ring_size = get_min_ring2_from_cod_class(atom.cod_class);
    atom.is_aromatic = cod_class_is_aromatic(atom.cod_class);
    atom.hybrid = hybrid_from_bonding_idx(atom.bonding_idx, atom.is_metal,
                                          atom.connectivity);
    compute_hash(atom);
    // Store no-charge variant for COD table lookup consistency.
    // cod_class itself doesn't encode charges (it's based on ring aromaticity and
    // neighbor topology). The charge-sensitive fields are bonding_idx/hybrid/nb1nb2_sp,
    // but those are used separately in lookups. For now, cod_class_no_charge equals
    // cod_class since the COD class string format doesn't include charge info.
    atom.cod_class_no_charge = atom.cod_class;
  }

  // Phase 2: Rebuild codClass with permissive aromaticity for output.
  // AceDRG's reDoAtomCodClassNames() sets isAromatic = isAromaticP after lookups.
  // cod_class_no_charge (used for table lookups) retains the strict version.
  bool has_permissive_diff = false;
  for (const auto& ring : rings)
    if (ring.is_aromatic_permissive != ring.is_aromatic)
      has_permissive_diff = true;
  if (has_permissive_diff) {
    for (auto& ring : rings)
      ring.is_aromatic = ring.is_aromatic_permissive;
    for (auto& atom : atoms)
      atom.ring_rep_s.clear();
    set_atoms_ring_rep_s(atoms, rings);
    for (size_t i = 0; i < atoms.size(); ++i) {
      set_atom_cod_class_name_new2(atoms[i], atoms[i], 2, atoms, neighbors);
      set_special_3nb_symb2(atoms[i], atoms, neighbors);
    }
    // Don't call cod_class_to_atom2 here — it would overwrite cod_main, cod_root,
    // nb_symb, nb2_symb, nb3_symb with permissive-aromaticity values, but these
    // fields must retain strict-aromaticity values for COD table lookups.
  }

  return atoms;
}

// AceDrg hash used in *HRS.table files
void AcedrgTables::compute_hash(CodAtomInfo& atom) const {
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

  // d3: non-metal connectivity + 8
  // AceDRG excludes metal neighbors from connectivity when computing hash.
  int d3 = 8 + static_cast<int>(atom.conn_atoms_no_metal.size());

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
    std::string footprint = cat(d1, '_', d2, '_', d3, '_', d4, '_', d5);

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



// ============================================================================
// CCP4 atom type assignment (AceDRG)
// ============================================================================

namespace {

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

int ccp4_material_type(Element el) {
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

std::string bond_order_key(BondType type) {
  std::string s = bond_type_to_string(type);
  if (s.empty())
    s = "single";
  s = to_upper(s);
  if (s.size() > 4)
    s.resize(4);
  return s;
}

void set_hydro_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  atom.ccp4_type = "H";
  if (atom.conn_atoms.size() == 1) {
    int nb = atom.conn_atoms[0];
    if (atoms[nb].chem_type == "S")
      atom.ccp4_type = "HSH1";
  }
}

void set_org_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
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
      if (r5 && r6) atom.ccp4_type = "CR56";
      else if (r5 == 2) atom.ccp4_type = "CR55";
      else if (r6 == 2) atom.ccp4_type = "CR66";
      else if (r5 == 1) atom.ccp4_type = (nh == 1 ? "CR15" : nh == 0 ? "CR5" : "C");
      else if (r6 == 1) atom.ccp4_type = (nh == 1 ? "CR16" : nh == 0 ? "CR6" : "C");
      else if (nh == 1) atom.ccp4_type = "C1";
      else if (nh == 2) atom.ccp4_type = "C2";
      else if (nh == 0) atom.ccp4_type = "C";
    } else if (atom.bonding_idx == 3) {
      if (nh == 0) atom.ccp4_type = "CT";
      else if (nh == 1) atom.ccp4_type = "CH1";
      else if (nh == 2) atom.ccp4_type = "CH2";
      else if (nh == 3) atom.ccp4_type = "CH3";
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "CSP";
    }
  } else if (atom.chem_type == "N") {
    if (atom.bonding_idx == 2) {
      if (nconn == 3) {
        if (nh == 1) atom.ccp4_type = "NH1";
        else if (nh == 2) atom.ccp4_type = "NH2";
        else if (nh == 0) atom.ccp4_type = "NH0";
        else atom.ccp4_type = "N";
      } else if (nconn == 2) {
        if (nh == 1) atom.ccp4_type = "N21";
        else if (nh == 0) atom.ccp4_type = "N20";
        else atom.ccp4_type = "N";
      }
    } else if (atom.bonding_idx == 3) {
      if (nconn == 4) {
        if (nh == 1) atom.ccp4_type = "NT1";
        else if (nh == 2) atom.ccp4_type = "NT2";
        else if (nh == 3) atom.ccp4_type = "NT3";
        else if (nh == 4) atom.ccp4_type = "NT4";
        else if (nh == 0) atom.ccp4_type = "NT";
        else atom.ccp4_type = "N";
      } else if (nconn == 3) {
        if (nh == 1) atom.ccp4_type = "N31";
        else if (nh == 2) atom.ccp4_type = "N32";
        else if (nh == 3) atom.ccp4_type = "N33";
        else if (nh == 0) atom.ccp4_type = "N30";
        else atom.ccp4_type = "N3";
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "NSP";
    }
  } else if (atom.chem_type == "P") {
    atom.ccp4_type = (nconn == 4 ? "P" : "P1");
  } else if (atom.chem_type == "O") {
    auto has_nb_type = [&](const char* t) {
      return std::any_of(atom.conn_atoms.begin(), atom.conn_atoms.end(),
                         [&](int nb) { return atoms[nb].chem_type == t; });
    };
    bool lP = has_nb_type("P"), lS = has_nb_type("S"), lB = has_nb_type("B");
    bool has_par_charge = std::fabs(atom.par_charge) > 1e-6f;
    bool has_negative_charge = atom.formal_charge < 0;
    auto oc_type = [&]() -> const char* {
      return lP ? "OP" : lS ? "OS" : lB ? "OB" : "OC";
    };
    if (atom.bonding_idx == 2) {
      if (has_par_charge && atom.par_charge < 0) {
        atom.ccp4_type = oc_type();
      } else if (nconn == 2) {
        if (nh == 1) atom.ccp4_type = "OH1";
        else if (nh == 2) atom.ccp4_type = "OH2";
        else atom.ccp4_type = "O";
      } else if (has_negative_charge) {
        atom.ccp4_type = oc_type();
      } else {
        atom.ccp4_type = "O";
      }
    } else if (atom.bonding_idx == 3) {
      bool lC = has_nb_type("C");
      if (lC && nh == 1 && nconn == 2) atom.ccp4_type = "OH1";
      else if (nh == 2) atom.ccp4_type = "OH2";
      else if (nconn == 2) {
        if (has_par_charge) atom.ccp4_type = "OC2";
        else if (nh == 1) atom.ccp4_type = "OH1";
        else atom.ccp4_type = "O2";
      } else if (nconn == 1 && has_negative_charge) {
        atom.ccp4_type = oc_type();
      }
    } else if (nconn == 1) {
      atom.ccp4_type = has_negative_charge ? oc_type() : "O";
    } else {
      atom.ccp4_type = "O";
    }
  } else if (atom.chem_type == "S") {
    if (nconn == 3 || nconn == 4) atom.ccp4_type = (nh == 0 ? "S3" : "SH1");
    else if (nconn == 2) atom.ccp4_type = (nh == 0 ? "S2" : "SH1");
    else if (nconn == 1) atom.ccp4_type = "S1";
    else atom.ccp4_type = (nh == 1 ? "SH1" : "S");
  } else if (atom.chem_type == "Se") {
    atom.ccp4_type = "SE";
  } else {
    atom.ccp4_type = atom.chem_type;
  }
}

void set_one_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
  int ntype = ccp4_material_type(atoms[idx].el);
  switch (ntype) {
    case 1: set_hydro_ccp4_type(atoms, idx); break;
    case 2: set_org_ccp4_type(atoms, idx); break;
    default: atoms[idx].ccp4_type = atoms[idx].chem_type; break;
  }
  atoms[idx].ccp4_type = to_upper(atoms[idx].ccp4_type);
}

}  // namespace

void AcedrgTables::load_ccp4_bonds(const std::string& path) {
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

bool AcedrgTables::search_ccp4_bond(const std::string& type1,
                                           const std::string& type2,
                                           const std::string& order,
                                           CodStats& out) const {
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

std::vector<Ccp4AtomInfo> build_ccp4_atoms(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int>>& neighbors) {
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
    // AceDRG treats metal-coordinated carbons with 4+ total neighbors as sp3 in CCP4 types.
    if (info.el == El::C && info.bonding_idx == 2) {
      size_t total_conn = info.conn_atoms.size();
      if (total_conn >= 4 && total_conn > info.conn_atoms_no_metal.size())
        info.bonding_idx = 3;
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(cc.atoms[i].charge));
    atoms.emplace_back(std::move(info));
  }
  return atoms;
}

void assign_all_ccp4_types(std::vector<Ccp4AtomInfo>& atoms) {
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);
}

std::vector<std::string> AcedrgTables::compute_ccp4_types(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int>>& neighbors) const {
  std::vector<Ccp4AtomInfo> atoms = build_ccp4_atoms(cc, atom_info, neighbors);
  assign_all_ccp4_types(atoms);

  std::vector<std::string> out;
  out.reserve(atoms.size());
  for (const auto& atom : atoms)
    out.push_back(atom.ccp4_type);
  return out;
}

void AcedrgTables::assign_ccp4_types(ChemComp& cc) const {
  if (cc.atoms.empty())
    return;

  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);
  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);

  std::vector<Ccp4AtomInfo> atoms = build_ccp4_atoms(cc, atom_info, neighbors);

  // Valence-based charge calculation for atoms bonded to metals.
  // Acedrg computes ligand charges via valence bookkeeping: for non-metal atoms,
  // sum only non-metal bond orders, then set charge = -(expectedValence - sumBo).
  // This ensures metal-bonded oxygens (like O bonded to both C and Hg) get
  // formal_charge=-1 and thus type "OC".
  // Note: atoms bonded ONLY to metals (no non-metal neighbors except H) are excluded
  // from this adjustment - their type stays "O" even with formal charge.
  for (size_t i = 0; i < atoms.size(); ++i) {
    const Ccp4AtomInfo& info = atoms[i];
    // Skip metals, hydrogens, and atoms already charged
    if (info.el.is_metal() || info.el == El::H || info.formal_charge != 0)
      continue;
    // Check if this atom has any metal neighbors
    bool has_metal_neighbor = info.conn_atoms.size() > info.conn_atoms_no_metal.size();
    if (!has_metal_neighbor)
      continue;
    if (std::none_of(info.conn_atoms_no_metal.begin(), info.conn_atoms_no_metal.end(),
                      [&](int nb) { return !cc.atoms[nb].is_hydrogen(); }))
      continue;
    int expected_valence = get_expected_valence(info.el);
    if (expected_valence == 0)
      continue;
    // Sum bond orders to non-metal neighbors only
    float sum_bo = 0.0f;
    for (const BondInfo& bi : adjacency[i]) {
      if (!cc.atoms[bi.neighbor_idx].el.is_metal())
        sum_bo += order_of_bond_type(bi.type);
    }
    // Calculate remaining valence and set formal charge
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    if (rem_v != 0)
      atoms[i].formal_charge = -rem_v;
  }

  assign_all_ccp4_types(atoms);

  for (size_t i = 0; i < atoms.size(); ++i)
    cc.atoms[i].chem_type = atoms[i].ccp4_type;
}

// ============================================================================
// Implementation - Bond search
// ============================================================================

void AcedrgTables::fill_restraints(ChemComp& cc) const {
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
      int match_level = fill_bond(cc, atom_info, bond);
      int idx1 = cc.find_atom_index(bond.id1.atom);
      int idx2 = cc.find_atom_index(bond.id2.atom);
      auto apply_ccp4 = [&](bool override_existing) {
        if (ccp4_types.empty() || idx1 < 0 || idx2 < 0)
          return false;
        if (!override_existing && !std::isnan(bond.value))
          return false;
        std::string order = bond_order_key(bond.type);
        CodStats vs;
        if (search_ccp4_bond(ccp4_types[idx1], ccp4_types[idx2], order, vs) ||
            search_ccp4_bond(ccp4_types[idx2], ccp4_types[idx1], order, vs)) {
          bond.value = vs.value;
          bond.esd = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
          return true;
        }
        CodStats v1, v2;
        bool has1 = search_ccp4_bond(ccp4_types[idx1], ".", order, v1);
        bool has2 = search_ccp4_bond(ccp4_types[idx2], ".", order, v2);
        if (has1 && has2) {
          bond.value = 0.5 * (v1.value + v2.value);
          bond.esd = 0.02;
          return true;
        }
        return false;
      };
      // CCP4 energetic library fallback
      if (std::isnan(bond.value))
        apply_ccp4(false);
      // CCP4 DELO override for explicit delocalized bonds only.
      // Applying this to regular single/double C-O bonds causes mismatches with AceDRG.
      if (!std::isnan(bond.value) && !ccp4_types.empty() &&
          bond.type == BondType::Deloc && match_level < 4) {
        if (idx1 >= 0 && idx2 >= 0) {
          const std::string& t1 = ccp4_types[idx1];
          const std::string& t2 = ccp4_types[idx2];
          // Check for carboxylate C-O bond (C bonded to O or OC)
          bool is_c_o_bond = (t1 == "C" && (t2 == "O" || t2 == "OC")) ||
                             (t2 == "C" && (t1 == "O" || t1 == "OC"));
          if (is_c_o_bond) {
            // Check if this C has two terminal O neighbors (carboxylate pattern)
            int c_idx = (t1 == "C") ? idx1 : idx2;
            int terminal_o_count = 0;
            for (const auto& b : cc.rt.bonds) {
              if (b.id1.atom == cc.atoms[c_idx].id || b.id2.atom == cc.atoms[c_idx].id) {
                const std::string& other_name = (b.id1.atom == cc.atoms[c_idx].id) ? b.id2.atom : b.id1.atom;
                auto it_other = cc.find_atom(other_name);
                if (it_other != cc.atoms.end() && it_other->el == El::O) {
                  // Check if this oxygen is terminal (only bonded to C and possibly H)
                  int heavy_neighbors = 0;
                  for (const auto& b2 : cc.rt.bonds) {
                    if (b2.id1.atom == other_name || b2.id2.atom == other_name) {
                      const std::string& nb = (b2.id1.atom == other_name) ? b2.id2.atom : b2.id1.atom;
                      auto it_nb = cc.find_atom(nb);
                      if (it_nb != cc.atoms.end() && it_nb->el != El::H)
                        ++heavy_neighbors;
                    }
                  }
                  if (heavy_neighbors == 1)  // Only bonded to C (the carboxyl carbon)
                    ++terminal_o_count;
                }
              }
            }
            if (terminal_o_count >= 2) {
              // True carboxylate pattern - use DELO value
              CodStats vs;
              if (search_ccp4_bond("OC", "C", "DELO", vs) ||
                  search_ccp4_bond("C", "OC", "DELO", vs)) {
                bond.value = vs.value;
                bond.esd = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
              }
            }
          }
        }
      }
    }
  }

  // Populate value_nucleus/esd_nucleus for X-H bonds from prot_hydr_dists
  for (auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom);
    int idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0)
      continue;
    const CodAtomInfo& a1 = atom_info[idx1];
    const CodAtomInfo& a2 = atom_info[idx2];
    if (a1.el == El::H || a2.el == El::H) {
      const CodAtomInfo& h_atom = (a1.el == El::H) ? a1 : a2;
      const CodAtomInfo& heavy_atom = (a1.el == El::H) ? a2 : a1;
      ProtHydrDist phd = search_prot_hydr_dist(h_atom, heavy_atom);
      if (!std::isnan(phd.nucleus_val)) {
        bond.value_nucleus = phd.nucleus_val;
        bond.esd_nucleus = phd.nucleus_sigma;
      }
    }
  }

  // Fill angles (skipped when angle tables are not loaded)
  if (!angle_hrs_.empty()) {
  for (auto& angle : cc.rt.angles) {
    if (std::isnan(angle.value)) {
      fill_angle(cc, atom_info, angle);
    }
  }

  // AceDRG adjustment: enforce planar ring angle sum ((n-2)*180/n) for SP2 rings.
  auto atom_index = cc.make_atom_index();
  std::map<int, std::vector<size_t>> rings;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int ring_id : atom_info[i].in_rings)
      rings[ring_id].push_back(i);
  for (const auto& ring : rings) {
    const std::vector<size_t>& ring_atoms = ring.second;
    if (ring_atoms.size() < 3)
      continue;
    bool planar = true;
    for (size_t idx : ring_atoms) {
      if (atom_info[idx].bonding_idx == 3) {
        planar = false;
        break;
      }
    }
    if (!planar)
      continue;
    std::set<std::string> ring_ids;
    for (size_t idx : ring_atoms)
      ring_ids.insert(cc.atoms[idx].id);
    std::vector<size_t> ring_angles;
    ring_angles.reserve(ring_atoms.size());
    for (size_t i = 0; i < cc.rt.angles.size(); ++i) {
      const auto& ang = cc.rt.angles[i];
      if (ring_ids.count(ang.id2.atom) &&
          ring_ids.count(ang.id1.atom) &&
          ring_ids.count(ang.id3.atom))
        ring_angles.push_back(i);
    }
    if (ring_angles.size() != ring_atoms.size())
      continue;
    double sum = 0.0;
    for (size_t idx : ring_angles)
      sum += cc.rt.angles[idx].value;
    double mean = sum / static_cast<double>(ring_angles.size());
    double target = (static_cast<double>(ring_atoms.size()) - 2.0) *
                    180.0 / static_cast<double>(ring_atoms.size());
    double shift = target - mean;
    for (size_t idx : ring_angles)
      cc.rt.angles[idx].value += shift;
  }

  // AceDRG adjustment: enforce 360-degree sum for sp2 centers with 3 angles,
  // keeping ring angles fixed (checkRingAngleConstraints behavior).
  std::map<std::string, std::vector<size_t>> center_angles;
  for (size_t i = 0; i < cc.rt.angles.size(); ++i)
    center_angles[cc.rt.angles[i].id2.atom].push_back(i);
  for (const auto& entry : center_angles) {
    if (entry.second.size() != 3)
      continue;
    auto it = atom_index.find(entry.first);
    if (it == atom_index.end())
      continue;
    size_t center_idx = it->second;
    if (atom_info[center_idx].hybrid != Hybridization::SP2)
      continue;
    std::vector<size_t> fixed, free;
    for (size_t idx_ang : entry.second) {
      const auto& ang = cc.rt.angles[idx_ang];
      auto it1 = atom_index.find(ang.id1.atom);
      auto it3 = atom_index.find(ang.id3.atom);
      if (it1 == atom_index.end() || it3 == atom_index.end())
        continue;
      const CodAtomInfo& a1 = atom_info[it1->second];
      const CodAtomInfo& a3 = atom_info[it3->second];
      if (angle_ring_size(atom_info[center_idx], a1, a3) > 0)
        fixed.push_back(idx_ang);
      else
        free.push_back(idx_ang);
    }
    if (free.empty())
      continue;
    double fixed_sum = 0.0;
    for (size_t idx_ang : fixed)
      fixed_sum += cc.rt.angles[idx_ang].value;
    double free_sum = 0.0;
    for (size_t idx_ang : free)
      free_sum += cc.rt.angles[idx_ang].value;
    double diff = (360.0 - fixed_sum - free_sum) / static_cast<double>(free.size());
    if (std::fabs(diff) > 0.01) {
      double new_sum = 0.0;
      for (size_t idx_ang : free) {
        cc.rt.angles[idx_ang].value += diff;
        new_sum += cc.rt.angles[idx_ang].value;
      }
      cc.rt.angles[free[0]].value += (360.0 - fixed_sum - new_sum);
    } else {
      cc.rt.angles[free[0]].value += diff;
    }
  }
  } // !angle_hrs_.empty()
}

namespace {

CodStats aggregate_stats(const std::vector<CodStats>& values) {
  if (values.empty())
    return CodStats();
  if (values.size() == 1)
    return values[0];
  // Match AceDRG setValueSet() aggregation (weighted mean + pooled sigma).
  double sum_val = 0.0;
  double sum1 = 0.0;
  int total_count = 0;
  for (const auto& v : values) {
    if (v.count > 0) {
      sum_val += v.value * v.count;
      sum1 += (v.count - 1) * v.sigma * v.sigma + v.count * v.value * v.value;
      total_count += v.count;
    }
  }
  if (total_count == 0)
    return CodStats();
  double mean = sum_val / total_count;
  double sum2 = mean * sum_val;
  double sigma = 0.0;
  if (total_count > 1)
    sigma = std::sqrt(std::fabs(sum1 - sum2) / (total_count - 1));
  else
    sigma = std::sqrt(std::fabs(sum1 - sum2) / total_count);
  return CodStats(mean, sigma, total_count);
}

}  // namespace

int AcedrgTables::fill_bond(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Bond& bond) const {

  int idx1 = cc.find_atom_index(bond.id1.atom);
  int idx2 = cc.find_atom_index(bond.id2.atom);
  if (idx1 < 0 || idx2 < 0)
    return 0;

  const CodAtomInfo& a1 = atom_info[idx1];
  const CodAtomInfo& a2 = atom_info[idx2];

  const char* source = "no_match";

  auto log_bond = [&](const char* src, int count = -1) {
    if (!verbose) return;
    std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s",
                 bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                 a1.hashing_value, a2.hashing_value,
                 hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid));
    if (!std::isnan(bond.value)) {
      std::fprintf(stderr, " → %s (%.3f, %.3f", src, bond.value, bond.esd);
      if (count >= 0) std::fprintf(stderr, ", n=%d", count);
      std::fprintf(stderr, ")\n");
    } else {
      std::fprintf(stderr, " -> %s\n", src);
    }
  };

  // Check for metal bond
  if (a1.is_metal || a2.is_metal) {
    const CodAtomInfo& metal = a1.is_metal ? a1 : a2;
    const CodAtomInfo& ligand = a1.is_metal ? a2 : a1;

    std::string mkey = to_upper(metal.el.name());
    std::string lkey = to_upper(ligand.el.name());
    auto it_m = covalent_radii_.find(mkey);
    auto it_l = covalent_radii_.find(lkey);
    if (it_m != covalent_radii_.end() && it_l != covalent_radii_.end()) {
      bond.value = it_m->second + it_l->second;
      bond.esd = 0.04;
      source = "metal_cova";
      log_bond(source);
      return 10;
    }

    CodStats vs = search_metal_bond(metal, ligand, atom_info);
    if (vs.count > 0) {
      bond.value = vs.value;
      double sigma = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
      bond.esd = std::max(0.02, clamp_bond_sigma(sigma));
      source = "metal";
      log_bond(source, vs.count);
      return 10;  // metal bond - treat as high-specificity match
    }
  }

  // Try both multilevel and HRS
  bool same_ring = are_in_same_ring(a1, a2);
  CodStats vs_ml = search_bond_multilevel(a1, a2);
  CodStats vs_hrs = search_bond_hrs(a1, a2, same_ring);

  // Acedrg's logic: use multilevel when it matches with sufficient threshold.
  // Acedrg iterates from start_level upward and uses the first level that meets threshold.
  // The threshold check is done inside search_bond_multilevel (entry count for levels 1-8,
  // and no threshold for HRS-like levels 9-11). When it returns with level >= 0, the threshold was met.
  // There's no special preference for HRS over type-based (levels 0-2) matches.
  CodStats vs;
  vs.level = -1;  // no match by default; prevents false acceptance below
  if (vs_ml.level >= 0) {
    // Multilevel matched with threshold met - use it
    vs = vs_ml;
    source = "multilevel";
  } else if (vs_hrs.count > 0) {
    // Multilevel didn't match - try HRS
    vs = vs_hrs;
    source = "HRS";
  }

  bool is_hrs = (source && std::strcmp(source, "HRS") == 0);
  bool accept = false;
  if (is_hrs || vs.level >= 9) {
    accept = vs.count > 0;
  } else {
    // For multilevel levels 0-8, thresholding is handled inside search_bond_multilevel().
    accept = vs.level >= 0;
  }
  if (accept) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    log_bond(source, vs.count);
    return vs.level;  // return match level for multilevel, or 0 for HRS
  }

  // Try element+hybridization fallback
  vs = search_bond_en(a1, a2);
  if (vs.count > 0) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    source = "EN";
    log_bond(source, vs.count);
    return 0;  // fallback - no good type-specific match
  }

  // No match found - leave bond.value as NaN so CCP4 fallback in fill_restraints
  // can be applied. AceDRG's search order is: multilevel -> HRS -> EN -> CCP4.
  log_bond(source);
  return 0;  // no type-specific match
}

CodStats AcedrgTables::search_bond_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& a2) const {

  if (verbose >= 2)
    std::fprintf(stderr, "      search_bond_multilevel: input a1=%s(hash=%d) a2=%s(hash=%d)\n",
                 a1.id.c_str(), a1.hashing_value, a2.id.c_str(), a2.hashing_value);

  const CodAtomInfo* left = nullptr;
  const CodAtomInfo* right = nullptr;
  order_two_atoms(a1, a2, left, right);

  if (verbose >= 2)
    std::fprintf(stderr, "      after order: left=%s(hash=%d) right=%s(hash=%d)\n",
                 left->id.c_str(), left->hashing_value, right->id.c_str(), right->hashing_value);

  // Build lookup keys
  int ha1 = left->hashing_value;
  int ha2 = right->hashing_value;

  std::string h1 = hybridization_to_string(left->hybrid);
  std::string h2 = hybridization_to_string(right->hybrid);
  if (h1 > h2) std::swap(h1, h2);
  std::string hybr_comb = cat(h1, '_', h2);

  std::string in_ring = are_in_same_ring(a1, a2) ? "Y" : "N";
  // AceDRG: if requested ring key is missing, fall back to the other (Y/N).
  auto has_ring_key = [&](const std::string& key) -> bool {
    auto it1 = bond_idx_full_.find(ha1);
    if (it1 == bond_idx_full_.end()) return false;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return false;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return false;
    return it3->second.find(key) != it3->second.end();
  };
  if (!has_ring_key(in_ring)) {
    std::string alt = (in_ring == "Y") ? "N" : "Y";
    if (has_ring_key(alt)) {
      if (verbose >= 2)
        std::fprintf(stderr, "      ring key '%s' missing, using '%s'\n",
                     in_ring.c_str(), alt.c_str());
      in_ring = alt;
    }
  }

  // Use neighbor symbols as additional keys.
  std::string a1_nb2 = left->nb2_symb;
  std::string a2_nb2 = right->nb2_symb;
  std::string a1_nb = left->nb1nb2_sp;
  std::string a2_nb = right->nb1nb2_sp;

  // Use COD main type as atom type (for 1D lookup).
  std::string a1_type = left->cod_main;
  std::string a2_type = right->cod_main;

  // Full COD class keys (as generated, AceDRG uses codClass directly).
  std::string a1_class = left->cod_class_no_charge;
  std::string a2_class = right->cod_class_no_charge;

  int num_th = min_observations_bond;
  if (a1.el == El::As || a2.el == El::As || a1.el == El::Ge || a2.el == El::Ge)
    num_th = 1;

  if (verbose >= 2)
    std::fprintf(stderr, "      lookup: hash=%d/%d hybr=%s ring=%s nb1nb2_sp=%s/%s nb2=%s/%s type=%s/%s\n",
                 ha1, ha2, hybr_comb.c_str(), in_ring.c_str(),
                 a1_nb.c_str(), a2_nb.c_str(), a1_nb2.c_str(), a2_nb2.c_str(),
                 a1_type.c_str(), a2_type.c_str());

  // AceDRG first tries the exact full codClass match (approx level 0) before
  // entering inter-level fallback.
  bool has_a1_class_only = false;  // a1 class exists, a2 class missing
  if (auto* p = find_val(bond_idx_full_, ha1))
  if (auto* p2 = find_val(*p, ha2))
  if (auto* p3 = find_val(*p2, hybr_comb))
  if (auto* p4 = find_val(*p3, in_ring))
  if (auto* p5 = find_val(*p4, a1_nb2))
  if (auto* p6 = find_val(*p5, a2_nb2))
  if (auto* p7 = find_val(*p6, a1_nb))
  if (auto* p8 = find_val(*p7, a2_nb))
  if (auto* p9 = find_val(*p8, a1_type))
  if (auto* p10 = find_val(*p9, a2_type))
  if (auto* p11 = find_val(*p10, a1_class)) {
    if (auto* vs = find_val(*p11, a2_class)) {
      if (vs->count >= num_th) {
        CodStats exact = *vs;
        exact.level = 0;
        if (verbose >= 2)
          std::fprintf(stderr,
                       "      matched exact codClass level=0: value=%.3f sigma=%.3f count=%d\n",
                       exact.value, exact.sigma, exact.count);
        return exact;
      }
    } else {
      has_a1_class_only = true;
    }
  }


  // Pre-resolve nested map pointers for bond lookup levels.
  auto* bond_ha = find_val(bond_idx_1d_, ha1);
  auto* bond_ha2 = bond_ha ? find_val(*bond_ha, ha2) : nullptr;
  auto* bond_hybr = bond_ha2 ? find_val(*bond_ha2, hybr_comb) : nullptr;
  auto* bond_ring = bond_hybr ? find_val(*bond_hybr, in_ring) : nullptr;
  auto* bond_nb2a = bond_ring ? find_val(*bond_ring, a1_nb2) : nullptr;
  auto* bond_nb2b = bond_nb2a ? find_val(*bond_nb2a, a2_nb2) : nullptr;
  auto* bond_nba = bond_nb2b ? find_val(*bond_nb2b, a1_nb) : nullptr;
  const auto* map_1d = bond_nba ? find_val(*bond_nba, a2_nb) : nullptr;

  // Same prefix but from bond_idx_2d_ (no atom type levels).
  auto* bond2_ha = find_val(bond_idx_2d_, ha1);
  auto* bond2_ha2 = bond2_ha ? find_val(*bond2_ha, ha2) : nullptr;
  auto* bond2_hybr = bond2_ha2 ? find_val(*bond2_ha2, hybr_comb) : nullptr;
  auto* bond2_ring = bond2_hybr ? find_val(*bond2_hybr, in_ring) : nullptr;
  auto* bond2_nb2a = bond2_ring ? find_val(*bond2_ring, a1_nb2) : nullptr;
  const auto* map_2d = bond2_nb2a ? find_val(*bond2_nb2a, a2_nb2) : nullptr;

  // Keep the signal only for level-gating compatibility with AceDRG's
  // class-missing branches. Do not short-circuit to direct 2D stats here:
  // AceDRG still proceeds through start-level logic in these cases.
  (void) has_a1_class_only;

  const auto* map_nb2 = bond2_ring;  // same as bond_idx_2d_[ha1][ha2][hybr][ring]

  // Determine key availability for gating (mirror AceDRG's presence checks).
  bool has_hybr = bond2_hybr != nullptr;
  bool has_in_ring = bond2_ring != nullptr;
  bool has_a1_nb2 = bond2_nb2a != nullptr;
  bool has_a2_nb2 = map_2d != nullptr;
  auto* map_2d_a1nb = map_2d ? find_val(*map_2d, a1_nb) : nullptr;
  bool has_a1_nb = map_2d_a1nb != nullptr;
  bool has_a2_nb = has_a1_nb && find_val(*map_2d_a1nb, a2_nb) != nullptr;
  auto* map_1d_a1type = map_1d ? find_val(*map_1d, a1_type) : nullptr;
  bool has_a1_type = map_1d_a1type != nullptr;
  bool has_a2_type = has_a1_type && find_val(*map_1d_a1type, a2_type) != nullptr;

  // Determine start level based on available keys (matching acedrg's dynamic logic from codClassify.cpp)
  // acedrg iterates from start_level upward (0→1→2→...) until threshold is met
  int start_level = 0;
  if (!has_a2_type) start_level = std::max(start_level, 1);   // a2M missing
  if (!has_a1_type) start_level = std::max(start_level, 2);   // a1M missing
  if (!has_a2_nb) start_level = std::max(start_level, 4);     // a2NB missing
  if (!has_a1_nb) start_level = std::max(start_level, 5);     // a1NB missing
  if (!has_a2_nb2) start_level = std::max(start_level, 7);    // a2NB2 missing
  if (!has_a1_nb2) start_level = std::max(start_level, 8);    // a1NB2 missing
  if (!has_in_ring) start_level = std::max(start_level, 10);
  if (!has_hybr) start_level = std::max(start_level, 11);

  if (verbose >= 2)
    std::fprintf(stderr, "      has: hybr=%d ring=%d a1_nb2=%d a2_nb2=%d a1_nb=%d a2_nb=%d a1_type=%d a2_type=%d start_level=%d\n",
                 has_hybr, has_in_ring, has_a1_nb2, has_a2_nb2, has_a1_nb, has_a2_nb, has_a1_type, has_a2_type, start_level);

  // AceDRG-like multilevel fallback: iterate from start_level upward until threshold is met.
  // For the special "a2 class missing" branch, AceDRG continues with post-search overrides,
  // so we must not return early here.
  CodStats matched_vs;
  bool have_matched_vs = false;
  for (int level = start_level; level < 12; level++) {
    // Skip levels that require keys we don't have
    if (level <= 2 && !map_1d) continue;  // type levels need map_1d
    if (level >= 3 && level <= 6 && !map_2d) continue;  // nb levels need map_2d
    if ((level == 7 || level == 8) && !map_nb2) continue;  // nb2 levels need map_nb2
    if (level == 9 && !has_in_ring) continue;
    if (level == 10 && !has_hybr) continue;
    CodStats vs;
    int values_size = 0;

    if (level == 0) {
      if (map_1d) {
        auto it1 = map_1d->find(a1_type);
        if (it1 != map_1d->end()) {
          auto it2 = it1->second.find(a2_type);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            vs = it2->second.front();
            values_size = 1;
          }
        }
      }
    } else if (level == 1) {
      // Level 1: a1_type matches exactly, aggregate over all a2_types
      if (map_1d) {
        std::vector<CodStats> values;
        auto it1 = map_1d->find(a1_type);
        if (it1 != map_1d->end()) {
          for (const auto& it2 : it1->second) {
            if (!it2.second.empty())
              values.push_back(it2.second.front());
          }
        }
        // Add entries where a2_type matches but a1_type differs (AceDRG tLev==1)
        for (const auto& it1b : *map_1d) {
          if (it1b.first == a1_type)
            continue;
          auto it2 = it1b.second.find(a2_type);
          if (it2 != it1b.second.end() && !it2->second.empty())
            values.push_back(it2->second.front());
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 2) {
      // Level 2: a2_type matches exactly, aggregate over different a1_types
      if (map_1d) {
        std::vector<CodStats> values;
        for (const auto& it1 : *map_1d) {
          if (it1.first != a1_type) {
            auto it2 = it1.second.find(a2_type);
            if (it2 != it1.second.end() && !it2->second.empty())
              values.push_back(it2->second.front());
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 3) {
      if (map_2d) {
        auto it1 = map_2d->find(a1_nb);
        if (it1 != map_2d->end()) {
          auto it2 = it1->second.find(a2_nb);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            values_size = static_cast<int>(it2->second.size());
            vs = aggregate_stats(it2->second);
            if (verbose >= 2)
              std::fprintf(stderr, "      level 3: found %d entries, value=%.3f count=%d (need %d)\n",
                           values_size, vs.value, vs.count, num_th);
          }
        }
        if (verbose >= 2 && values_size == 0)
          std::fprintf(stderr, "      level 3 miss: a1_nb='%s' a2_nb='%s'\n", a1_nb.c_str(), a2_nb.c_str());
      }
    } else if (level == 4 || level == 5) {
      if (map_2d) {
        std::vector<CodStats> values;
        // AceDRG level 4: combine entries with a1NB plus entries with a2NB but not a1NB.
        // AceDRG level 5: only entries with a2NB but not a1NB.
        if (level == 4) {
          auto it1 = map_2d->find(a1_nb);
          if (it1 != map_2d->end()) {
            for (const auto& it2 : it1->second)
              for (const auto& v : it2.second)
                values.push_back(v);
          }
        }
        for (const auto& it1 : *map_2d) {
          if (it1.first == a1_nb)
            continue;
          auto it2 = it1.second.find(a2_nb);
          if (it2 != it1.second.end()) {
            for (const auto& v : it2->second)
              values.push_back(v);
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 6) {
      if (map_2d) {
        std::vector<CodStats> values;
        for (const auto& it1 : *map_2d)
          for (const auto& it2 : it1.second)
            for (const auto& v : it2.second)
              values.push_back(v);
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 7 || level == 8) {
      if (map_nb2) {
        std::vector<CodStats> values;
        if (level == 7) {
          auto it1 = map_nb2->find(a1_nb2);
          if (it1 != map_nb2->end()) {
            for (const auto& it2 : it1->second)
              for (const auto& it3 : it2.second)
                for (const auto& it4 : it3.second)
                  for (const auto& v : it4.second)
                    values.push_back(v);
          }
        }
        for (const auto& it2 : *map_nb2) {
          if (it2.first == a1_nb2)
            continue;
          auto it3 = it2.second.find(a2_nb2);
          if (it3 != it2.second.end()) {
            for (const auto& it4 : it3->second)
              for (const auto& it5 : it4.second)
                for (const auto& v : it5.second)
                  values.push_back(v);
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 9) {
      if (auto* p = find_val(bond_hasp_2d_, ha1))
      if (auto* p2 = find_val(*p, ha2))
      if (auto* p3 = find_val(*p2, hybr_comb))
      if (auto* p4 = find_val(*p3, in_ring))
      if (!p4->empty()) {
        values_size = static_cast<int>(p4->size());
        vs = aggregate_stats(*p4);
      }
    } else if (level == 10) {
      if (auto* p = find_val(bond_hasp_1d_, ha1))
      if (auto* p2 = find_val(*p, ha2))
      if (auto* p3 = find_val(*p2, hybr_comb))
      if (!p3->empty()) {
        values_size = static_cast<int>(p3->size());
        vs = aggregate_stats(*p3);
      }
    } else if (level == 11) {
      if (auto* p = find_val(bond_hasp_0d_, ha1))
      if (auto* p2 = find_val(*p, ha2))
      if (!p2->empty()) {
        values_size = static_cast<int>(p2->size());
        vs = aggregate_stats(*p2);
      }
    }

    // Check threshold (matching acedrg's logic)
    // Levels 1-8 use entry count threshold, level 0 and others use observation count
    bool threshold_met = false;
    if (level >= 1 && level <= 8) {
      threshold_met = values_size >= num_th;  // entry count
    } else if (level >= 9) {
      threshold_met = values_size > 0;        // HRS entries only require presence
    } else {
      threshold_met = vs.count >= num_th;     // observation count
    }
    if (values_size > 0 && threshold_met) {
      if (verbose >= 2)
        std::fprintf(stderr, "      matched: level=%d\n", level);
      vs.level = level;
      if (!has_a1_class_only)
        return vs;
      matched_vs = vs;
      have_matched_vs = true;
      break;
    }
  }

  // AceDRG branch (searchCodOrgBonds2_2, a2C missing):
  // after inter-level search, prefer exact (a1M,a2M)[0] if enough observations,
  // otherwise use exact (a1NB1NB2,a2NB1NB2) 2D stats regardless threshold.
  if (has_a1_class_only) {
    if (map_1d) {
      auto it1 = map_1d->find(a1_type);
      if (it1 != map_1d->end()) {
        auto it2 = it1->second.find(a2_type);
        if (it2 != it1->second.end() && !it2->second.empty()) {
          CodStats exact_1d = it2->second.front();
          if (exact_1d.count >= num_th) {
            exact_1d.level = 0;
            return exact_1d;
          }
        }
      }
    }
    if (map_2d) {
      auto it1 = map_2d->find(a1_nb);
      if (it1 != map_2d->end()) {
        auto it2 = it1->second.find(a2_nb);
        if (it2 != it1->second.end() && !it2->second.empty()) {
          CodStats exact_2d = aggregate_stats(it2->second);
          exact_2d.level = 3;
          return exact_2d;
        }
      }
    }
    if (have_matched_vs)
      return matched_vs;
  }

  if (have_matched_vs)
    return matched_vs;

  // No match found in detailed tables
  if (verbose >= 2)
    std::fprintf(stderr, "      matched: level=none (no multilevel match)\n");
  CodStats no_match;
  no_match.level = -1;  // sentinel for no match
  return no_match;
}

CodStats AcedrgTables::search_bond_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& a2, bool in_ring) const {

  const CodAtomInfo* left = nullptr;
  const CodAtomInfo* right = nullptr;
  order_two_atoms(a1, a2, left, right);

  BondHRSKey key;
  key.hash1 = left->hashing_value;
  key.hash2 = right->hashing_value;

  std::string h1 = hybridization_to_string(left->hybrid);
  std::string h2 = hybridization_to_string(right->hybrid);
  if (h1 > h2) std::swap(h1, h2);
  key.hybrid_pair = cat(h1, '_', h2);
  key.in_ring = in_ring ? "Y" : "N";

  auto it = bond_hrs_.find(key);
  if (it != bond_hrs_.end()) {
    return it->second;
  }

  return CodStats();
}

CodStats AcedrgTables::search_bond_en(const CodAtomInfo& a1,
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

  if (auto* p1 = find_val(en_bonds_, elem1))
    if (auto* p2 = find_val(*p1, sp1))
      if (auto* p3 = find_val(*p2, elem2))
        if (auto* p4 = find_val(*p3, sp2))
          if (!p4->empty())
            return aggregate_stats(*p4);
  return CodStats();
}

ProtHydrDist AcedrgTables::search_prot_hydr_dist(const CodAtomInfo& /* h_atom */,
    const CodAtomInfo& heavy_atom) const {
  // Build type key: H_sp[1|2|3][_arom]_ELEM
  // Examples: H_sp3_C, H_sp2_arom_C, H_sp2_N
  std::string hybr;
  switch (heavy_atom.hybrid) {
    case Hybridization::SP1: hybr = "sp1"; break;
    case Hybridization::SP2: hybr = heavy_atom.is_aromatic ? "sp2_arom" : "sp2"; break;
    case Hybridization::SP3: hybr = "sp3"; break;
    default: return ProtHydrDist();
  }

  std::string type_key = cat("H_", hybr, '_', heavy_atom.el.name());

  auto it = prot_hydr_dists_.find(type_key);
  if (it != prot_hydr_dists_.end()) {
    if (verbose >= 2)
      std::fprintf(stderr, "      prot_hydr_dists: found %s -> electron=%.4f nucleus=%.4f\n",
                   type_key.c_str(), it->second.electron_val, it->second.nucleus_val);
    return it->second;
  }

  // Try without aromatic qualifier if sp2_arom not found
  if (heavy_atom.is_aromatic && heavy_atom.hybrid == Hybridization::SP2) {
    type_key = std::string("H_sp2_") + heavy_atom.el.name();
    it = prot_hydr_dists_.find(type_key);
    if (it != prot_hydr_dists_.end()) {
      if (verbose >= 2)
        std::fprintf(stderr, "      prot_hydr_dists: found %s (fallback from arom) -> electron=%.4f nucleus=%.4f\n",
                     type_key.c_str(), it->second.electron_val, it->second.nucleus_val);
      return it->second;
    }
  }

  return ProtHydrDist();
}

CodStats AcedrgTables::search_metal_bond(const CodAtomInfo& metal,
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
    cat_to(ligand_class, '_', id);

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

  if (class_entry && class_entry->count > metal_class_min_count)
    return CodStats(class_entry->value, class_entry->sigma,
                      class_entry->count);
  if (pre_entry)
    return CodStats(pre_entry->pre_value, pre_entry->pre_sigma,
                      pre_entry->pre_count);
  return CodStats();
}

bool AcedrgTables::lookup_pep_tors(const std::string& a1,
    const std::string& a2, const std::string& a3, const std::string& a4,
    TorsionEntry& out) const {
  auto it = pep_tors_.find(cat(a1, '_', a2, '_', a3, '_', a4));
  if (it == pep_tors_.end())
    return false;
  out = it->second;
  return true;
}

// ============================================================================
// Implementation - Angle search
// ============================================================================

CodStats AcedrgTables::search_angle_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3) const {

  // Build lookup keys - canonicalize flanking atoms
  // Table format: ha1=center, ha2=flank_min, ha3=flank_max
  const CodAtomInfo *flank_min, *flank_max;
  order_two_atoms(a1, a3, flank_min, flank_max);
  int ha1 = center.hashing_value;
  int ha2 = flank_min->hashing_value;
  int ha3 = flank_max->hashing_value;

  // Build hybridization tuple - table format: center_flankMax_flankMin
  std::string hc = hybridization_to_string(center.hybrid);
  std::string hmin = hybridization_to_string(flank_min->hybrid);
  std::string hmax = hybridization_to_string(flank_max->hybrid);
  std::string hybr_tuple = cat(hc, '_', hmax, '_', hmin);

  // Build valueKey (ring:hybr_tuple)
  int ring_val = angle_ring_size(center, *flank_min, *flank_max);
  std::string value_key = cat(ring_val, ':', hybr_tuple);

  // Get neighbor symbols - table format: a1=center, a2=flank_min, a3=flank_max
  std::string a1_nb2 = center.nb2_symb;
  std::string a2_nb2 = flank_min->nb2_symb;
  std::string a3_nb2 = flank_max->nb2_symb;
  std::string a1_root = center.cod_root;
  std::string a2_root = flank_min->cod_root;
  std::string a3_root = flank_max->cod_root;
  std::string a1_nb = center.nb_symb;
  std::string a2_nb = flank_min->nb_symb;
  std::string a3_nb = flank_max->nb_symb;
  std::string a1_type = center.cod_main;
  std::string a2_type = flank_min->cod_main;
  std::string a3_type = flank_max->cod_main;

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

  // Helper lambda: check if vs from angle level meets threshold and return it.
  auto try_angle = [&](const std::vector<CodStats>* vec, int min_obs, const char* level_name) -> CodStats {
    if (vec && !vec->empty()) {
      CodStats vs = vec->front();
      if (vs.count >= min_obs) {
        if (verbose >= 2)
          std::fprintf(stderr, "      matched angle %s-%s-%s: level=%s\n",
                       a1.id.c_str(), center.id.c_str(), a3.id.c_str(), level_name);
        return vs;
      }
    }
    return CodStats();
  };

  // When center and a flank share the same hash, the table may store their
  // property columns (root, nb2, nb, type) in either order. Try both.
  int max_attempts = (ha1 == ha2 || ha1 == ha3) ? 2 : 1;
  for (int attempt = 0; attempt < max_attempts; ++attempt) {
    if (attempt == 1) {
      if (ha1 == ha2) {
        std::swap(a1_root, a2_root);
        std::swap(a1_nb2, a2_nb2);
        std::swap(a1_nb, a2_nb);
        std::swap(a1_type, a2_type);
      } else {  // ha1 == ha3
        std::swap(a1_root, a3_root);
        std::swap(a1_nb2, a3_nb2);
        std::swap(a1_nb, a3_nb);
        std::swap(a1_type, a3_type);
      }
    }

    // Level 1D: exact match with atom types (16 keys)
    if (auto* p = find_val(angle_idx_1d_, ha1))
    if (auto* p2 = find_val(*p, ha2))
    if (auto* p3 = find_val(*p2, ha3))
    if (auto* p4 = find_val(*p3, value_key))
    if (auto* p5 = find_val(*p4, a1_root))
    if (auto* p6 = find_val(*p5, a2_root))
    if (auto* p7 = find_val(*p6, a3_root))
    if (auto* p8 = find_val(*p7, a1_nb2))
    if (auto* p9 = find_val(*p8, a2_nb2))
    if (auto* p10 = find_val(*p9, a3_nb2))
    if (auto* p11 = find_val(*p10, a1_nb))
    if (auto* p12 = find_val(*p11, a2_nb))
    if (auto* p13 = find_val(*p12, a3_nb))
    if (auto* p14 = find_val(*p13, a1_type))
    if (auto* p15 = find_val(*p14, a2_type)) {
      CodStats vs = try_angle(find_val(*p15, a3_type), min_observations_angle, "1D");
      if (!std::isnan(vs.value))
        return vs;
    }

    // Level 2D: no atom types (13 keys)
    if (auto* p = find_val(angle_idx_2d_, ha1))
    if (auto* p2 = find_val(*p, ha2))
    if (auto* p3 = find_val(*p2, ha3))
    if (auto* p4 = find_val(*p3, value_key))
    if (auto* p5 = find_val(*p4, a1_root))
    if (auto* p6 = find_val(*p5, a2_root))
    if (auto* p7 = find_val(*p6, a3_root))
    if (auto* p8 = find_val(*p7, a1_nb2))
    if (auto* p9 = find_val(*p8, a2_nb2))
    if (auto* p10 = find_val(*p9, a3_nb2))
    if (auto* p11 = find_val(*p10, a1_nb))
    if (auto* p12 = find_val(*p11, a2_nb)) {
      CodStats vs = try_angle(find_val(*p12, a3_nb), min_observations_angle_fallback, "2D");
      if (!std::isnan(vs.value))
        return vs;
    }

    // Level 3D: hash + valueKey + roots + NB2 (10 keys)
    if (auto* p = find_val(angle_idx_3d_, ha1))
    if (auto* p2 = find_val(*p, ha2))
    if (auto* p3 = find_val(*p2, ha3))
    if (auto* p4 = find_val(*p3, value_key))
    if (auto* p5 = find_val(*p4, a1_root))
    if (auto* p6 = find_val(*p5, a2_root))
    if (auto* p7 = find_val(*p6, a3_root))
    if (auto* p8 = find_val(*p7, a1_nb2))
    if (auto* p9 = find_val(*p8, a2_nb2)) {
      CodStats vs = try_angle(find_val(*p9, a3_nb2), min_observations_angle_fallback, "3D");
      if (!std::isnan(vs.value))
        return vs;
    }
  }  // end same-hash retry loop

  // Level 4D: hash + valueKey + roots (7 keys)
  if (auto* p = find_val(angle_idx_4d_, ha1))
  if (auto* p2 = find_val(*p, ha2))
  if (auto* p3 = find_val(*p2, ha3))
  if (auto* p4 = find_val(*p3, value_key))
  if (auto* p5 = find_val(*p4, a1_root))
  if (auto* p6 = find_val(*p5, a2_root)) {
    CodStats vs = try_angle(find_val(*p6, a3_root), min_observations_angle_fallback, "4D");
    if (!std::isnan(vs.value))
      return vs;
  }

  // Level 5D: hash + valueKey (4 keys)
  if (auto* p = find_val(angle_idx_5d_, ha1))
  if (auto* p2 = find_val(*p, ha2))
  if (auto* p3 = find_val(*p2, ha3)) {
    CodStats vs = try_angle(find_val(*p3, value_key), min_observations_angle_fallback, "5D");
    if (!std::isnan(vs.value))
      return vs;
  }

  // Level 6D: hash only (3 keys)
  if (auto* p = find_val(angle_idx_6d_, ha1))
  if (auto* p2 = find_val(*p, ha2))
  if (auto* p3 = find_val(*p2, ha3))
  if (!p3->empty()) {
    if (verbose >= 2)
      std::fprintf(stderr, "      matched angle %s-%s-%s: level=6D\n",
                   a1.id.c_str(), center.id.c_str(), a3.id.c_str());
    return p3->front();
  }

  // No match found in detailed tables
  return CodStats();
}

void AcedrgTables::fill_angle(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Angle& angle) const {

  int idx1 = cc.find_atom_index(angle.id1.atom);
  int idx2 = cc.find_atom_index(angle.id2.atom);  // center
  int idx3 = cc.find_atom_index(angle.id3.atom);
  if (idx1 < 0 || idx2 < 0 || idx3 < 0)
    return;

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
  CodStats vs = search_angle_multilevel(a1, center, a3);
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

CodStats AcedrgTables::search_angle_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3, int ring_size) const {

  AngleHRSKey key;
  // Table format: hash1=center, hash2=min(flank), hash3=max(flank)
  key.hash1 = center.hashing_value;
  key.hash2 = std::min(a1.hashing_value, a3.hashing_value);
  key.hash3 = std::max(a1.hashing_value, a3.hashing_value);

  // Build hybrid tuple - table format: center_flankMax_flankMin
  std::string hc = hybridization_to_string(center.hybrid);
  std::string hmin, hmax;
  if (a1.hashing_value <= a3.hashing_value) {
    hmin = hybridization_to_string(a1.hybrid);
    hmax = hybridization_to_string(a3.hybrid);
  } else {
    hmin = hybridization_to_string(a3.hybrid);
    hmax = hybridization_to_string(a1.hybrid);
  }
  std::string hybr_tuple = cat(hc, '_', hmax, '_', hmin);
  key.value_key = cat(ring_size, ':', hybr_tuple);

  auto it = angle_hrs_.find(key);
  if (it != angle_hrs_.end()) {
    return it->second;
  }

  return CodStats();
}

std::vector<double> AcedrgTables::get_metal_angles(Element metal,
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

// Helper to compact element list: ["H","H","H"] -> "H3"
// Only compacts runs of 3+ identical elements (acedrg convention)
static std::string compact_element_list(const std::vector<std::string>& elems) {
  std::string result;
  for (size_t i = 0; i < elems.size(); ) {
    size_t count = 1;
    while (i + count < elems.size() && elems[i + count] == elems[i])
      ++count;
    if (count >= 3) {
      // Compact runs of 3+
      result += elems[i];
      cat_to(result, count);
    } else {
      // Don't compact runs of 1-2, just repeat the element
      for (size_t j = 0; j < count; ++j)
        result += elems[i];
    }
    i += count;
  }
  return result;
}

std::string AcedrgTables::compute_acedrg_type(
    const CodAtomInfo& atom,
    const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  std::string result = element_name(atom.el);

  // Build neighbor descriptions: neighbor_element + neighbor's_other_neighbors
  // Note: AceDRG excludes metals from the atom typing system entirely
  std::map<std::string, int> neighbor_groups;
  for (int nb_idx : neighbors[atom.index]) {
    if (atoms[nb_idx].is_metal)  // Skip metals in neighbor description
      continue;
    std::string desc = element_name(atoms[nb_idx].el);
    // Collect second neighbors (excluding the central atom and metals)
    std::vector<std::string> second_nbs;
    for (int nb2_idx : neighbors[nb_idx]) {
      if (nb2_idx != atom.index && !atoms[nb2_idx].is_metal)
        second_nbs.push_back(element_name(atoms[nb2_idx].el));
    }
    // Sort alphabetically
    std::sort(second_nbs.begin(), second_nbs.end());
    // Compact runs: H,H,H -> H3
    desc += compact_element_list(second_nbs);
    neighbor_groups[desc]++;
  }

  // Sort groups: by length desc, then alphabetically
  std::vector<std::pair<std::string, int>> sorted_groups(
      neighbor_groups.begin(), neighbor_groups.end());
  std::sort(sorted_groups.begin(), sorted_groups.end(),
    [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
      if (a.first.size() != b.first.size())
        return a.first.size() > b.first.size();  // longer first
      return a.first < b.first;  // then alphabetically
    });

  for (const auto& group : sorted_groups) {
    result += "(" + group.first + ")";
    if (group.second > 1)
      cat_to(result, group.second);
  }
  return result;
}

std::vector<std::string> AcedrgTables::compute_acedrg_types(const ChemComp& cc) const {
  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);

  std::vector<std::string> types;
  types.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    types.push_back(atom_info[i].cod_class);
  }
  return types;
}

} // namespace gemmi
