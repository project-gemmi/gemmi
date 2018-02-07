// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if
#include <cmath>      // for NAN
#include <cstring>    // for size_t
#include <iterator>   // for back_inserter
#include <map>        // for map
#include <string>
#include <vector>

#include "elem.hpp"
#include "unitcell.hpp"
#include "symmetry.hpp"

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
T& find_or_add(std::vector<T>& vec, const std::string& name) {
  T* ret = find_or_null(vec, name);
  if (ret)
    return *ret;
  vec.emplace_back(name);
  return vec.back();
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
  bool operator==(const SequenceItem& o) const {
    return num == o.num && mon == o.mon;
  }
};

using Sequence = std::vector<SequenceItem>;

struct Entity {
  EntityType type = EntityType::Unknown;
  PolymerType polymer_type = PolymerType::NA;
  Sequence sequence;

  std::string type_as_string() const {
    switch (type) {
      case EntityType::Polymer: return "polymer";
      case EntityType::NonPolymer: return "non-polymer";
      case EntityType::Water: return "water";
      default /*EntityType::Unknown*/: return "?";
    }
  }
  std::string polymer_type_as_string() const {
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

inline double calculate_dihedral_from_atoms(const Atom* a, const Atom* b,
                                            const Atom* c, const Atom* d) {
  if (a && b && c && d)
    return calculate_dihedral(a->pos, b->pos, c->pos, d->pos);
  return NAN;
}


struct ResidueId {
  struct OptionalNum {
    enum { None=-10000 };
    int value = None;
    bool has_value() const { return value != None; }
    std::string str() const {
      return has_value() ? std::to_string(value) : "?";
    }
    OptionalNum& operator=(int n) { value = n; return *this; }
    bool operator==(const OptionalNum& o) const { return value == o.value; }
    bool operator==(int n) const { return value == n; }
    explicit operator int() const { return value; }
    explicit operator bool() const { return has_value(); }
    // these are defined for partial compatibility with C++17 std::optional
    using value_type = int;
    int& operator*() { return value; }
    const int& operator*() const { return value; }
    int& emplace(int n) { value = n; return value; }
  };

  OptionalNum label_seq;  // mmCIF _atom_site.label_seq_id

  // traditional residue sequence numbers are coupled with insertion codes
  OptionalNum seq_num; // sequence number
  char icode = '\0';  // insertion code

  //bool in_main_conformer/is_point_mut
  //uint32_t segment_id; // number or 4 characters
  std::string segment; // normally up to 4 characters in the PDB file
  std::string name;

  char printable_icode() const { return icode ? icode : ' '; }
  bool same_seq_id(const ResidueId& o) const {
    return seq_num == o.seq_num && icode == o.icode;
  }
  int seq_num_for_pdb() const { return int(seq_num ? seq_num : label_seq); }
  std::string seq_id() const {
    std::string r = seq_num.str();
    if (icode)
      r += icode;
    return r;
  }
};

struct Residue : public ResidueId {
  bool is_cis = false;  // bond to the next residue marked as cis
  char het_flag = '\0';  // 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  std::vector<Atom> atoms;
  // Connection::name (from Model::connections)
  std::vector<std::string> conn;
  Chain* parent = nullptr;

  Residue() = default;
  explicit Residue(const ResidueId& rid) noexcept : ResidueId(rid) {}

  std::string ident() const;

  std::vector<Atom>& children() { return atoms; }
  const std::vector<Atom>& children() const { return atoms; }
  bool matches(const ResidueId& rid) const;

  const Atom* find_by_element(El el) const {
    for (const Atom& a : atoms)
      if (a.element == el)
        return &a;
    return nullptr;
  }

  // default values accept anything
  const Atom* find_atom(const std::string& atom_name, char altloc='*',
                        El el=El::X) const {
    for (const Atom& a : atoms)
      if (a.name == atom_name && (altloc == '*' || a.altloc == altloc)
          && (el == El::X || a.element == el))
        return &a;
    return nullptr;
  }
  Atom* find_atom(const std::string& atom_name, char altloc='*', El el=El::X) {
    const Residue* const_this = this;
    return const_cast<Atom*>(const_this->find_atom(atom_name, altloc, el));
  }

  const Atom* get_ca() const {
    static const std::string CA("CA");
    return find_atom(CA, '*', El::C);
  }

  const Residue* prev_bonded_aa() const;
  const Residue* next_bonded_aa() const;

  double calculate_omega(const Residue& next) const;
  bool calculate_phi_psi_omega(double* phi, double* psi, double* omega) const;
};

struct Chain {
  std::string name;
  std::string auth_name;
  std::vector<Residue> residues;
  std::string entity_id;
  Model* parent = nullptr;
  int force_pdb_serial = 0;

  explicit Chain(std::string cname) noexcept : name(cname) {}
  Residue* find_residue(const ResidueId& rid);
  Residue* find_or_add_residue(const ResidueId& rid);
  void append_residues(std::vector<Residue> new_resi);
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
  enum Type { Covale, Disulf, Hydrog, MetalC, None };
  std::string name;  // the id is refered by Residue::conn;
  Type type = None;
  SymmetryImage image = SymmetryImage::Unspecified;
  char altloc[2] = {'\0', '\0'};
  ResidueId res_id[2];
  // The pointers get invalidated by some changes to the model.
  Residue* res[2] = {nullptr, nullptr};
  // _struct_conn.ptnr[12]_label_atom_id, only for LINK not for SSBOND
  std::string atom[2];
};

inline const char* get_mmcif_connection_type_id(Connection::Type t) {
  static constexpr const char* type_ids[] = {
    "covale", "disulf", "hydrog", "metalc" };
  return type_ids[t];
}

struct Model {
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  std::vector<Connection> connections;
  Structure* parent = nullptr;
  explicit Model(std::string mname) noexcept : name(mname) {}

  Chain* find_chain(const std::string& chain_name) {
    return impl::find_or_null(chains, chain_name);
  }
  Chain& find_or_add_chain(const std::string& chain_name) {
    return impl::find_or_add(chains, chain_name);
  }
  Connection* find_connection_by_name(const std::string& conn_name) {
    return impl::find_or_null(connections, conn_name);
  }
  void invalidate_connection_pointers() {
    for (Connection& conn : connections)
      conn.res[0] = conn.res[1] = nullptr;
  }

  std::vector<Chain>& children() { return chains; }
  const std::vector<Chain>& children() const { return chains; }
};

struct NcsOp {
  std::string id;
  bool given;
  Transform tr;
  Position apply(const Position& p) const { return Position(tr.apply(p)); }
};

struct Structure {
  std::string name;
  gemmi::UnitCell cell;
  std::string sg_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::map<std::string, Entity> entities;

  // Store ORIGXn / _database_PDB_matrix.origx*
  Transform origx;

  // Minimal metadata with keys being mmcif tags: _entry.id, _exptl.method, ...
  std::map<std::string, std::string> info;

  const char* get_info(const std::string& tag, const char* def=nullptr) const {
    auto it = info.find(tag);
    return it != info.end() ? it->second.c_str() : def;
  }
  Model* find_model(const std::string& model_name) {
    return impl::find_or_null(models, model_name);
  }
  Model& find_or_add_model(const std::string& model_name) {
    return impl::find_or_add(models, model_name);
  }

  Entity* find_entity(const std::string& ent_id) {
    auto ent = entities.find(ent_id);
    return ent == entities.end() ? nullptr : &ent->second;
  }
  const Entity* find_entity(const std::string& ent_id) const {
    auto ent = entities.find(ent_id);
    return ent == entities.end() ? nullptr : &ent->second;
  }


  Entity& find_or_add_entity(const std::string& ent_id) {
    auto ent = entities.find(ent_id);
    if (ent == entities.end())
      ent = entities.emplace(ent_id, Entity()).first;
    return ent->second;
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

inline std::string Residue::ident() const {
  return (parent ? parent->name_for_pdb() + "/" : "") +
         (seq_num ? seq_id() : label_seq.str());
}

inline bool Residue::matches(const ResidueId& rid) const {
  if (!rid.label_seq && rid.seq_num && !same_seq_id(rid))
    return false;
  return label_seq == rid.label_seq &&
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
  const Atom* C = find_atom("C", '*', El::C);
  const Atom* nextN = next.find_atom("N", '*', El::N);
  return calculate_dihedral_from_atoms(get_ca(), C, nextN, next.get_ca());
}

inline bool Residue::calculate_phi_psi_omega(double* phi, double* psi,
                                             double* omega) const {
  const Atom* CA = get_ca();
  if (!CA)
    return false;
  const Residue* prev = prev_bonded_aa();
  const Residue* next = next_bonded_aa();
  if (!prev && !next)
    return false;
  const Atom* C = find_atom("C", '*', El::C);
  const Atom* N = find_atom("N", '*', El::N);
  const Atom* prevC = prev ? prev->find_atom("C", '*', El::C) : nullptr;
  const Atom* nextN = next ? next->find_atom("N", '*', El::N) : nullptr;
  if (phi)
    *phi = calculate_dihedral_from_atoms(prevC, N, CA, C);
  if (psi)
    *psi = calculate_dihedral_from_atoms(N, CA, C, nextN);
  if (omega)
    *omega = calculate_dihedral_from_atoms(CA, C, nextN,
                                           next ? next->get_ca() : nullptr);
  return true;
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

inline void Chain::append_residues(std::vector<Residue> new_resi) {
  size_t init_capacity = residues.capacity();
  int seqnum = 0;
  for (const Residue& res : residues)
    seqnum = std::max({seqnum, int(res.label_seq), int(res.seq_num)});
  for (Residue& res : new_resi) {
    res.label_seq = res.seq_num = ++seqnum;
    res.icode = '\0';
  }
  std::move(new_resi.begin(), new_resi.end(), std::back_inserter(residues));
  add_backlinks(*this);
  if (parent && residues.capacity() != init_capacity)
    parent->invalidate_connection_pointers();
}

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

  // TODO: if "entities" were not specifed, deduce them based on sequence

  // fill UnitCell::images
  if (const SpaceGroup* sg = find_spacegroup_by_name(sg_hm)) {
    for (Op op : sg->operations()) {
      // TODO strict NCS
      if (op == Op::identity())
        continue;
      Matrix33 rot = {
        double(op.rot[0][0]), double(op.rot[0][1]), double(op.rot[0][2]),
        double(op.rot[1][0]), double(op.rot[1][1]), double(op.rot[1][2]),
        double(op.rot[2][0]), double(op.rot[2][1]), double(op.rot[2][2]) };
      double mult = 1.0 / Op::TDEN;
      Vec3 tran(mult * op.tran[0], mult * op.tran[1], mult * op.tran[2]);
      cell.images.push_back(FTransform{rot, tran});
    }
  }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
