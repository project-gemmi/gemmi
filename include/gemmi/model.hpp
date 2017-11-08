// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if
#include <cmath>      // for NAN
#include <cstring>    // for size_t
#include <map>        // for map
#include <memory>     // for unique_ptr
#include <string>
#include <vector>

#include <linalg.h>
#include "elem.hpp"
#include "unitcell.hpp"

namespace gemmi {

namespace impl {

template<typename T>
T* find_or_null(std::vector<T>& vec, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(), [&name](const T& m) {
      return m.name == name;
  });
  return it != vec.end() ? &*it : nullptr;
}

template<typename T>
T* find_or_add(std::vector<T>& vec, const std::string& name) {
  T* ret = find_or_null(vec, name);
  if (ret)
    return ret;
  vec.emplace_back(name);
  return &vec.back();
}

} // namespace impl


struct Structure;
struct Model;
struct Chain;
struct Residue;

// File format with macromolecular model.
enum class CoorFormat { Unknown, Pdb, Cif, Json };

enum class EntityType : unsigned char {
  Unknown,
  Polymer,
  NonPolymer,
  // _entity.type macrolide is in PDBx/mmCIF, but no PDB entry uses it
  //Macrolide,
  Water
};

// number of different _entity_poly.type values in the PDB in mid-2017:
//   168923 polypeptide(L)
//   9905   polydeoxyribonucleotide
//   4559   polyribonucleotide
//   156    polydeoxyribonucleotide/polyribonucleotide hybrid
//   57     polypeptide(D)
//   18     polysaccharide(D)
//   4      other
//   2      peptide nucleic acid
//   1      cyclic-pseudo-peptide
//   0      polysaccharide(L)  (never used but present in mmcif_pdbx_v50.dic)
enum class PolymerType : unsigned char {
  PeptideL,
  PeptideD,
  Dna,
  Rna,
  DnaRnaHybrid,
  SaccharideD,
  SaccharideL,
  Pna, // artificial thing
  CyclicPseudoPeptide,
  Other,
  NA // not applicable or not determined
};


struct SequenceItem {
  int num;
  std::string mon;
  explicit SequenceItem(std::string m) noexcept : num(-1), mon(m) {}
  SequenceItem(int n, std::string m) noexcept : num(n), mon(m) {}
};

using Sequence = std::vector<SequenceItem>;

typedef linalg::mat<double,4,4> Mat4x4;

struct Entity {
  std::string id;  // it does not need to be number according to mmCIF spec
  EntityType type;
  PolymerType polymer_type;
  Sequence sequence;

  std::string type_as_string() {
    switch (type) {
      case EntityType::Polymer: return "polymer";
      case EntityType::NonPolymer: return "non-polymer";
      case EntityType::Water: return "water";
      default /*EntityType::Unknown*/: return "?";
    }
  }
  std::string polymer_type_as_string() {
    switch (polymer_type) {
      case PolymerType::PeptideL: return "polypeptide(L)";
      case PolymerType::PeptideD: return "polypeptide(D)";
      case PolymerType::Dna: return "polydeoxyribonucleotide";
      case PolymerType::Rna: return "polyribonucleotide";
      case PolymerType::DnaRnaHybrid:
        return "polydeoxyribonucleotide/polyribonucleotide hybrid";
      case PolymerType::SaccharideD: return "polysaccharide(D)";
      case PolymerType::SaccharideL: return "polysaccharide(L)";
      case PolymerType::Other: return "other";
      case PolymerType::Pna: return "peptide nucleic acid";
      case PolymerType::CyclicPseudoPeptide: return "cyclic-pseudo-peptide";
      default /*PolymerType::NA*/: return "?";
    }
  }
};

inline EntityType entity_type_from_string(const std::string& t) {
  if (t == "polymer")     return EntityType::Polymer;
  if (t == "non-polymer") return EntityType::NonPolymer;
  if (t == "water")       return EntityType::Water;
  return EntityType::Unknown;
}

inline PolymerType polymer_type_from_string(const std::string& t) {
  if (t == "polypeptide(L)")          return PolymerType::PeptideL;
  if (t == "polydeoxyribonucleotide") return PolymerType::Dna;
  if (t == "polyribonucleotide")      return PolymerType::Rna;
  if (t == "polydeoxyribonucleotide/polyribonucleotide hybrid")
                                      return PolymerType::DnaRnaHybrid;
  if (t == "polypeptide(D)")          return PolymerType::PeptideD;
  if (t == "polysaccharide(D)")       return PolymerType::SaccharideD;
  if (t == "other")                   return PolymerType::Other;
  if (t == "peptide nucleic acid")    return PolymerType::Pna;
  if (t == "cyclic-pseudo-peptide")   return PolymerType::CyclicPseudoPeptide;
  if (t == "polysaccharide(L)")       return PolymerType::SaccharideL;
  return PolymerType::NA;
}

struct Atom {
  std::string name;
  char group;
  char altloc; // 0 if not set
  signed char charge;  // [-8, +8]
  gemmi::Element element = gemmi::El::X;
  Position pos;
  float occ;
  float b_iso;
  float u11=0, u22=0, u33=0, u12=0, u13=0, u23=0;
  Residue* parent = nullptr;
};

inline double calculate_distance_sq(const Atom* a, const Atom* b) {
  if (!a || !b)
    return NAN;
  double dx = a->pos.x - b->pos.x;
  double dy = a->pos.y - b->pos.y;
  double dz = a->pos.z - b->pos.z;
  return dx*dx + dy*dy + dz*dz;
}


struct ResidueInfo {
  // simple classification, ELS is something else
  enum Kind { UNKNOWN=0, AA=1, NA=2, RNA=(2|4), DNA=(2|8), HOH=16, ELS=32 };
  Kind kind;
  // PDB format has non-standard residues (modified AA) marked as HETATM.
  bool pdb_standard;
  // rough count of hydrogens used to estimate mass with implicit hydrogens
  short hydrogen_count;

  bool is_water() const { return kind == HOH; }
  bool is_nucleic() const { return kind & NA; }
  bool is_amino() const { return kind == AA; }
};

// hyderogen_count needs to be verified
inline const ResidueInfo find_tabulated_residue(const std::string& name) {
  if (name.size() == 3) {
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(name.c_str())) {
      case ID("ALA"): return { ResidueInfo::AA,  true,  5 };
      case ID("ARG"): return { ResidueInfo::AA,  true, 13 };
      case ID("ASN"): return { ResidueInfo::AA,  true,  6 };
      case ID("ABA"): return { ResidueInfo::AA,  false, 7 };
      case ID("ASP"): return { ResidueInfo::AA,  true,  4 };
      case ID("ASX"): return { ResidueInfo::AA,  true,  4 };
      case ID("CYS"): return { ResidueInfo::AA,  true,  4 };
      case ID("CSH"): return { ResidueInfo::AA,  false, 5 };
      case ID("GLN"): return { ResidueInfo::AA,  true,  8 };
      case ID("GLU"): return { ResidueInfo::AA,  true,  6 };
      case ID("GLX"): return { ResidueInfo::AA,  true,  8 };
      case ID("GLY"): return { ResidueInfo::AA,  true,  3 };
      case ID("HIS"): return { ResidueInfo::AA,  true,  8 };
      case ID("ILE"): return { ResidueInfo::AA,  true, 11 };
      case ID("LEU"): return { ResidueInfo::AA,  true, 11 };
      case ID("LYS"): return { ResidueInfo::AA,  true, 13 };
      case ID("MET"): return { ResidueInfo::AA,  true,  9 };
      case ID("MSE"): return { ResidueInfo::AA,  false, 9 };
      case ID("ORN"): return { ResidueInfo::AA,  false,10 };
      case ID("PHE"): return { ResidueInfo::AA,  true,  9 };
      case ID("PRO"): return { ResidueInfo::AA,  true,  7 };
      case ID("SER"): return { ResidueInfo::AA,  true,  5 };
      case ID("THR"): return { ResidueInfo::AA,  true,  7 };
      case ID("TRP"): return { ResidueInfo::AA,  true, 10 };
      case ID("TYR"): return { ResidueInfo::AA,  true,  9 };
      case ID("UNK"): return { ResidueInfo::AA,  true,  0 };
      case ID("VAL"): return { ResidueInfo::AA,  true,  9 };
      case ID("SEC"): return { ResidueInfo::AA,  true,  6 };
      case ID("PYL"): return { ResidueInfo::AA,  true, 19 };
      case ID("HOH"): return { ResidueInfo::HOH, false, 2 };
      case ID("WAT"): return { ResidueInfo::HOH, false, 2 };
      case ID("H20"): return { ResidueInfo::HOH, false, 2 };
      case ID("HEM"): return { ResidueInfo::ELS, false,30 };
      case ID("SO4"): return { ResidueInfo::ELS, false, 0 };
      case ID("SUL"): return { ResidueInfo::ELS, false, 0 };
    }
#undef ID
  } else if (name.size() == 1) {
    switch (name[0]) {
      case 'A': return { ResidueInfo::RNA, true, 13 };
      case 'C': return { ResidueInfo::RNA, true, 13 };
      case 'G': return { ResidueInfo::RNA, true, 13 };
      case 'I': return { ResidueInfo::RNA, true, 12 };
      case 'T': return { ResidueInfo::RNA, true, 14 };
      case 'U': return { ResidueInfo::RNA, true, 12 };
    }
  } else if (name.size() == 2) {
    if (name[0] == 'D' || name[0] == '+')
      switch (name[1]) {
        case 'A': return { ResidueInfo::DNA, true, 13 };
        case 'C': return { ResidueInfo::DNA, true, 13 };
        case 'G': return { ResidueInfo::DNA, true, 13 };
        case 'I': return { ResidueInfo::DNA, true, 12 };
        case 'T': return { ResidueInfo::DNA, true, 14 };
        case 'U': return { ResidueInfo::DNA, true, 12 };
      }
  }
  return { ResidueInfo::UNKNOWN, false, 0 };
}


struct ResidueId {
  // traditional residue sequence numbers are coupled with insertion codes
  struct SNIC {
    int seq_num;
    char ins_code;
    bool operator==(const SNIC& o) const {
      return seq_num == o.seq_num && ins_code == o.ins_code;
    }
    bool operator!=(const SNIC& o) const { return !operator==(o); }
  };

  ResidueId(int id, SNIC auth_id, std::string rname) noexcept
    : seq_id(id), snic(auth_id), name(rname) {}
  ResidueId(SNIC id, std::string rname) noexcept : snic(id), name(rname) {}
  ResidueId(int id, std::string rname) noexcept : seq_id(id), name(rname) {}

  enum { UnknownId=-1000 };
  int seq_id = UnknownId;
  SNIC snic = {UnknownId, '\0'};
  //bool in_main_conformer/is_point_mut
  //uint32_t segment_id; // number or 4 characters
  std::string segment; // normally up to 4 characters in the PDB file
  std::string name;
};

struct Residue : public ResidueId {
  bool is_cis = false;  // bond to the next residue marked as cis
  std::vector<Atom> atoms;
  // Maps Connection::id (from Model::connections) to Atom::name.
  std::map<std::string, std::string> conn;
  Chain* parent = nullptr;

  explicit Residue(const ResidueId& rid) noexcept : ResidueId(rid) {}
  int seq_id_for_pdb() const {
    return snic.seq_num != UnknownId ? snic.seq_num : seq_id;
  }
  ResidueInfo get_info() const { return find_tabulated_residue(name); }
  std::vector<Atom>& children() { return atoms; }
  const std::vector<Atom>& children() const { return atoms; }
  bool matches(const ResidueId& rid) const;
  const Atom* find_by_element(El el) const {
    for (const Atom& a : atoms)
      if (a.element == el)
        return &a;
    return nullptr;
  }
  const Atom* find_by_name(const std::string& name) const {
    for (const Atom& a : atoms)
      if (a.name == name)
        return &a;
    return nullptr;
  }
  const Atom* find_by_name_and_elem(const std::string& name, El elem) const {
    for (const Atom& a : atoms)
      if (a.name == name && a.element == elem &&
          (a.altloc == 0 || a.altloc == 'A'))
        return &a;
    return nullptr;
  }
  const Atom* get_ca() const {
    static const std::string CA("CA");
    return find_by_name_and_elem(CA, El::C);
  }
  const Residue* prev_bonded_aa() const;
  const Residue* next_bonded_aa() const;
  double calculate_omega(const Residue& next) const;
};

struct Chain {
  std::string name;
  std::string auth_name;
  std::vector<Residue> residues;
  Entity *entity = nullptr;
  Model* parent = nullptr;
  int force_pdb_serial = 0;

  explicit Chain(std::string cname) noexcept : name(cname) {}
  Residue* find_residue(const ResidueId& rid);
  Residue* find_or_add_residue(const ResidueId& rid);
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }
  const std::string& name_for_pdb() const {
    return auth_name.empty() ? name : auth_name;
  }
};

// A connection. Corresponds to _struct_conn.
// Symmetry operators are not trusted and not stored.
// We assume that the nearest symmetry mate is connected.
struct Connection {
  enum Type { Covale, CoveleBase, CovalePhosphate, CovaleSugar, Disulf,
              Hydrog, MetalC, Mismat, ModRes, SaltBr };
  std::string id;  // the id is refered by Residue::conn;
  Type type;
  // The pointers get invalidated by some changes to the model.
  Residue* res1;
  Residue* res2;
};

struct Model {
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  std::vector<Connection> connections;
  Structure* parent = nullptr;
  explicit Model(std::string mname) noexcept : name(mname) {}

  Chain* find_chain(const std::string& chain_name) {
    return impl::find_or_null(chains, chain_name);
  }
  Chain* find_or_add_chain(const std::string& chain_name) {
    return impl::find_or_add(chains, chain_name);
  }
  std::vector<Chain>& children() { return chains; }
  const std::vector<Chain>& children() const { return chains; }
};

struct NcsOp {
  std::string id;
  bool given;
  Mat4x4 transform;
};

struct Structure {
  std::string name;
  gemmi::UnitCell cell;
  std::string sg_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::vector<std::unique_ptr<Entity>> entities;

  // Store ORIGXn / _database_PDB_matrix.origx*
  Mat4x4 origx = linalg::identity;

  // Minimal metadata with keys being mmcif tags: _entry.id, _exptl.method, ...
  std::map<std::string, std::string> info;

  const char* get_info(const std::string& tag, const char* def=nullptr) const {
    auto it = info.find(tag);
    return it != info.end() ? it->second.c_str() : def;
  }
  Model* find_model(const std::string& model_name) {
    return impl::find_or_null(models, model_name);
  }
  Model* find_or_add_model(const std::string& model_name) {
    return impl::find_or_add(models, model_name);
  }

  Entity* find_entity(const std::string& ent_id) {
    for (auto& ent : entities)
      if (ent->id == ent_id)
        return ent.get();
    return nullptr;
  }
  Entity* find_or_add_entity(const std::string& ent_id) {
    Entity* ent = find_entity(ent_id);
    if (!ent) {
      ent = new Entity{ent_id, EntityType::Unknown, PolymerType::NA, {}};
      entities.emplace_back(ent);
    }
    return ent;
  }
  const std::vector<Chain>& get_chains() const {
    // We don't handle yet a corner case (ever happening?)
    // in which the first model is lacking a chain.
    return models.at(0).chains;
  }

  double get_ncs_multiplier() const {
    int given = std::count_if(ncs.begin(), ncs.end(),
                              [](const NcsOp& o) { return o.given; });
    return (ncs.size() + 1.0) / (given + 1.0);  // +1 b/c identity not included
  }

  std::vector<Model>& children() { return models; }
  const std::vector<Model>& children() const { return models; }
  void finish();
};


inline bool Residue::matches(const ResidueId& rid) const {
  return seq_id == rid.seq_id &&
         (seq_id != Residue::UnknownId || snic == rid.snic) &&
         segment == rid.segment &&
         name == rid.name;
}

// TODO: handle alternative conformations (point mutations)
inline const Residue* Residue::prev_bonded_aa() const {
  if (parent && this != parent->residues.data() &&
      calculate_distance_sq(get_ca(), (this - 1)->get_ca()) < 6.0 * 6.0)
    return this - 1;
  return nullptr;
}

inline const Residue* Residue::next_bonded_aa() const {
  const Residue* next = this + 1;
  if (parent && next != parent->residues.data() + parent->residues.size() &&
      calculate_distance_sq(get_ca(), next->get_ca()) < 6.0 * 6.0)
    return next;
  return nullptr;
}

inline double Residue::calculate_omega(const Residue& next) const {
  const Atom* CA = get_ca();
  const Atom* C = find_by_name_and_elem("C", El::C);
  const Atom* nextN = next.find_by_name_and_elem("N", El::N);
  const Atom* nextCA = next.get_ca();
  if (CA && C && nextN && nextCA)
    return calculate_dihedral(CA->pos, C->pos, nextN->pos, nextCA->pos);
  return NAN;
}

inline Residue* Chain::find_residue(const ResidueId& rid) {
  auto it = std::find_if(residues.begin(), residues.end(),
                         [&](const Residue& r) { return r.matches(rid); });
  return it != residues.end() ? &*it : nullptr;
}

inline Residue* Chain::find_or_add_residue(const ResidueId& rid) {
  Residue* r = find_residue(rid);
  if (r)
    return r;
  residues.emplace_back(rid);
  return &residues.back();
}

template<class T> void add_backlinks(T& obj) {
  for (auto& child : obj.children()) {
    child.parent = &obj;
    add_backlinks(child);
  }
}
template<> inline void add_backlinks(Atom&) {}


template<class T> size_t count_atom_sites(const T& obj) {
  size_t sum = 0;
  for (const auto& child : obj.children())
    sum += count_atom_sites(child);
  return sum;
}
template<> inline size_t count_atom_sites(const Residue& res) {
  return res.atoms.size();
}


template<class T> double count_occupancies(const T& obj) {
  double sum = 0;
  for (const auto& child : obj.children())
    sum += count_occupancies(child);
  return sum;
}
template<> inline double count_occupancies(const Atom& atom) {
  return atom.occ;
}

inline void Structure::finish() {
  add_backlinks(*this);
  // if "entities" were not specifed, deduce them based on sequence
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
