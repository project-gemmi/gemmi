// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if
#include <array>
#include <cstdlib>    // for strtol
#include <iterator>   // for back_inserter
#include <map>        // for map
#include <stdexcept>  // for out_of_range
#include <string>
#include <vector>

#include "elem.hpp"
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "iterator.hpp"

namespace gemmi {

namespace impl {

// Optional int value. N is a special value that means not-set.
template<int N> struct OptionalInt {
  enum { None=N };
  int value = None;
  bool has_value() const { return value != None; }
  std::string str() const {
    return has_value() ? std::to_string(value) : "?";
  }
  OptionalInt& operator=(int n) { value = n; return *this; }
  bool operator==(const OptionalInt& o) const { return value == o.value; }
  bool operator==(int n) const { return value == n; }
  explicit operator int() const { return value; }
  explicit operator bool() const { return has_value(); }
  // these are defined for partial compatibility with C++17 std::optional
  using value_type = int;
  int& operator*() { return value; }
  const int& operator*() const { return value; }
  int& emplace(int n) { value = n; return value; }
};

template<typename T>
T* find_or_null(std::vector<T>& vec, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(),
                         [&name](const T& m) { return m.name == name; });
  return it != vec.end() ? &*it : nullptr;
}

template<typename T>
T& find_or_add(std::vector<T>& vec, const std::string& name) {
  if (T* ret = find_or_null(vec, name))
    return *ret;
  vec.emplace_back(name);
  return vec.back();
}

template<typename T> typename std::vector<T>::iterator
find_iter(std::vector<T>& vec, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(),
                         [&name](const T& m) { return m.name == name; });
  if (it == vec.end())
    throw std::invalid_argument(
        T::what() + (" " + name) + " not found (only [" +
        join_str(vec, ' ', [](const T& x) { return x.name; }) +
        "])");
  return it;
}

} // namespace impl


// File format with macromolecular model.
enum class CoorFormat { Unknown, Pdb, Mmcif, Mmjson };

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
  Unknown, // unknown or not applicable
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
  EntityType entity_type = EntityType::Unknown;
  PolymerType polymer_type = PolymerType::Unknown;
  Sequence sequence;

  std::string type_as_string() const {
    switch (entity_type) {
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
      default /*PolymerType::Unknown*/: return "?";
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
  return PolymerType::Unknown;
}

// sometimes a name shorter than "polydeoxyribonucleotide" is more readable
inline const char* polymer_type_abbr(PolymerType ptype) {
  switch (ptype) {
    case PolymerType::PeptideL: return "AAL";
    case PolymerType::PeptideD: return "AAD";
    case PolymerType::Dna: return "DNA";
    case PolymerType::Rna: return "RNA";
    case PolymerType::DnaRnaHybrid: return "xNA";
    default: return "";
  }
}

inline bool is_same_conformer(char altloc1, char altloc2) {
  return altloc1 == '\0' || altloc2 == '\0' || altloc1 == altloc2;
}

struct Atom {
  static const char* what() { return "Atom"; }
  std::string name;
  char altloc; // 0 if not set
  signed char charge;  // [-8, +8]
  Element element = El::X;
  int custom = 0;
  Position pos;
  float occ;
  float b_iso;
  float u11=0, u22=0, u33=0, u12=0, u13=0, u23=0;

  bool same_conformer(const Atom& other) const {
    return is_same_conformer(altloc, other.altloc);
  }
  bool same_group(const Atom& other) const { return name == other.name; }
  // Name as a string left-padded like in the PDB format:
  // the first two characters make the element name.
  std::string padded_name() const {
    std::string s;
    if (element.uname()[1] == '\0' && name.size() < 4)
      s += ' ';
    s += name;
    return s;
  }
};


struct ResidueId {
  using OptionalNum = impl::OptionalInt<-999>;

  // traditional residue sequence numbers are coupled with insertion codes
  OptionalNum seq_num; // sequence number
  char icode = ' ';  // insertion code

  //bool in_main_conformer/is_point_mut
  //uint32_t segment_id; // number or 4 characters
  std::string segment; // normally up to 4 characters in the PDB file
  std::string name;

  char has_icode() const { return icode != ' '; }
  bool same_seq_id(const ResidueId& o) const {
    return seq_num == o.seq_num && icode == o.icode;
  }
  std::string seq_id() const {
    std::string r = seq_num.str();
    if (has_icode())
      r += icode;
    return r;
  }
  std::string str() const { return seq_id() + "(" + name + ")"; }
  bool matches(const ResidueId& rid) const;
  // for first_conformation iterators
  bool same_group(const ResidueId& o) const { return same_seq_id(o); }
};

struct Residue : public ResidueId {
  OptionalNum label_seq;  // mmCIF _atom_site.label_seq_id
  bool is_cis = false;  // bond to the next residue marked as cis
  char het_flag = '\0';  // 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  std::vector<Atom> atoms;

  Residue() = default;
  explicit Residue(const ResidueId& rid) noexcept : ResidueId(rid) {}

  std::vector<Atom>& children() { return atoms; }
  const std::vector<Atom>& children() const { return atoms; }

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

  void remove_atom(const std::string& atom_name) {
    atoms.erase(impl::find_iter(atoms, atom_name));
  }

  // short-cuts to access peptide backbone atoms
  const Atom* get_ca() const {
    static const std::string CA("CA");
    return find_atom(CA, '*', El::C);
  }
  const Atom* get_c() const {
    static const std::string C("C");
    return find_atom(C, '*', El::C);
  }
  const Atom* get_n() const {
    static const std::string N("N");
    return find_atom(N, '*', El::N);
  }

  // short-cuts to access nucleic acid atoms
  const Atom* get_p() const {
    static const std::string P("P");
    return find_atom(P, '*', El::P);
  }
  const Atom* get_o3prim() const {
    static const std::string P("O3'");
    return find_atom(P, '*', El::O);
  }

  bool same_conformer(const Residue& other) const {
    return atoms.empty() || other.atoms.empty() ||
           atoms[0].same_conformer(other.atoms[0]) ||
           other.find_atom(other.atoms[0].name, atoms[0].altloc) != nullptr;
  }

  bool has_peptide_bond_to(const Residue& next) const {
    const Atom* a1 = get_c();
    const Atom* a2 = next.get_n();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < 2.0 * 2.0;
  }

  // Iterators that in case of multiple conformations (alt. locations)
  // skip all but the first conformation.
  UniqProxy<Atom> first_conformer() { return {atoms}; }
  ConstUniqProxy<Atom> first_conformer() const { return {atoms}; }
};

// ResidueGroup represents residues with the same sequence number and insertion
// code, but different residue names. I.e. microheterogeneity.
// Usually, there is only one residue in the group.
// The residues must be consecutive.
struct ResidueGroup {
  using iterator = std::vector<Residue>::iterator;
  using const_iterator = std::vector<Residue>::const_iterator;

  iterator begin_, end_;

  const_iterator begin() const { return begin_; }
  const_iterator end() const { return end_; }
  iterator begin() { return begin_; }
  iterator end() { return end_; }
  int size() const { return end_ - begin_; }
  const Residue& operator[](int i) const { return *(begin_ + i); }
  Residue& operator[](int i) { return *(begin_ + i); }
  Residue& at(int i) {
    if (i >= size())
      throw std::out_of_range("ResidueGroup: no item " + std::to_string(i));
    return *(begin_ + i);
  }
  const Residue& at(int i) const {
    return const_cast<ResidueGroup*>(this)->at(i);
  }
  Residue& by_resname(const std::string& name) {
    for (auto it = begin_; it != end_; ++it)
      if (it->name == name)
        return *it;
    throw std::invalid_argument("ResidueGroup has no residue " + name);
  }
  bool empty() const { return begin_ == end_; }
  explicit operator bool() const { return begin_ != end_; }
};


struct Chain {
  static const char* what() { return "Chain"; }
  std::string name;
  std::string auth_name;
  std::vector<Residue> residues;
  std::string entity_id;

  explicit Chain(std::string cname) noexcept : name(cname) {}
  ResidueGroup find_residue_group(int seqnum, char icode=' ');
  Residue* find_residue(const ResidueId& rid);
  Residue* find_or_add_residue(const ResidueId& rid);
  void append_residues(std::vector<Residue> new_resi);
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }
  const std::string& name_for_pdb() const {
    return auth_name.empty() ? name : auth_name;
  }

  ResidueGroup find_by_seqid(const std::string& seqid) {
    char* endptr;
    int seqnum = std::strtol(seqid.c_str(), &endptr, 10);
    if (endptr == seqid.c_str() || (*endptr != '\0'  && endptr[1] != '\0'))
      throw std::invalid_argument("Not a seqid: " + seqid);
    return find_residue_group(seqnum, *endptr);
  }

  // Returns the previous residue or nullptr.
  // Got complicated by handling of multi-conformations / microheterogeneity.
  const Residue* previous_residue(const Residue& res) const {
    const Residue* start = residues.data();
    for (const Residue* p = &res; p != start; )
      if (!res.same_seq_id(*--p)) {
        while (p != start && p->same_seq_id(*(p-1)) && !res.same_conformer(*p))
          --p;
        return p;
      }
    return nullptr;
  }

  // Returns the next residue or nullptr.
  const Residue* next_residue(const Residue& res) const {
    const Residue* end = residues.data() + residues.size();
    for (const Residue* p = &res + 1; p != end; ++p)
      if (!res.same_seq_id(*p)) {
        while (p+1 != end && p->same_seq_id(*(p+1)) && !res.same_conformer(*p))
          ++p;
        return p;
      }
    return nullptr;
  }

  const Residue* prev_bonded_aa(const Residue& res) const {
    if (const Residue* prev = previous_residue(res))
      if (prev->has_peptide_bond_to(res))
        return prev;
    return nullptr;
  }

  const Residue* next_bonded_aa(const Residue& res) const {
    if (const Residue* next = next_residue(res))
      if (res.has_peptide_bond_to(*next))
        return next;
    return nullptr;
  }

  // Iterators that in case of multiple conformations (microheterogeneity)
  // skip all but the first conformation.
  UniqProxy<Residue> first_conformer() { return {residues}; }
  ConstUniqProxy<Residue> first_conformer() const { return {residues}; }
};

struct AtomAddress {
  std::string chain_name;
  ResidueId res_id;
  std::string atom_name;
  char altloc = '\0';
  bool use_auth_name = true;

  AtomAddress() = default;
  AtomAddress(const Chain& ch, const Residue& res, const Atom& at)
    : chain_name(ch.name), res_id(res), atom_name(at.name), altloc(at.altloc),
      use_auth_name(false) {}

  bool operator==(const AtomAddress& o) const {
    return chain_name == o.chain_name && res_id.matches(o.res_id) &&
           atom_name == o.atom_name && altloc == o.altloc &&
           use_auth_name == o.use_auth_name;
  }

  std::string str() const {
    std::string r = chain_name + "/" + res_id.name + " " +
                    res_id.seq_id() + "/" + atom_name;
    if (altloc) {
      r += '.';
      r += altloc;
    }
    return r;
  }
};

struct CRA {
  Chain* chain;
  Residue* residue;
  Atom* atom;
};

struct const_CRA {
  const Chain* chain;
  const Residue* residue;
  const Atom* atom;
};

// A connection. Corresponds to _struct_conn.
// Symmetry operators are not trusted and not stored.
// We assume that the nearest symmetry mate is connected.
struct Connection {
  enum Type { Covale, Disulf, Hydrog, MetalC, None };
  std::string name;
  Type type = None;
  SameAsu asu = SameAsu::Any;
  AtomAddress atom[2];
};

inline const char* get_mmcif_connection_type_id(Connection::Type t) {
  static constexpr const char* type_ids[] = {
    "covale", "disulf", "hydrog", "metalc", "." };
  return type_ids[t];
}

struct Model {
  static const char* what() { return "Model"; }
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  std::vector<Connection> connections;
  explicit Model(std::string mname) noexcept : name(mname) {}

  Chain* find_chain(const std::string& chain_name) {
    return impl::find_or_null(chains, chain_name);
  }

  Chain& add_chain(const std::string& chain_name) {
    if (find_chain(chain_name))
      throw std::runtime_error("The chain '" + chain_name + "' already exists");
    chains.emplace_back(chain_name);
    return chains.back();
  }
  Chain& find_or_add_chain(const std::string& chain_name) {
    return impl::find_or_add(chains, chain_name);
  }

  void remove_chain(const std::string& chain_name) {
    chains.erase(impl::find_iter(chains, chain_name));
  }

  ResidueGroup residues(const std::string& auth_chain, int resnum, char icode) {
    ResidueGroup rg;
    for (Chain& chain : chains)
      if (chain.auth_name == auth_chain) {
        rg = chain.find_residue_group(resnum, icode);
        if (rg)
          break;
      }
    return rg;
  }

  Residue& residue(const std::string& auth_chain, int resnum, char icode) {
    ResidueGroup rg = residues(auth_chain, resnum, icode);
    if (rg.size() != 1) {
      std::string err = rg.empty() ? "No residue " : "Multiple residues ";
      err += auth_chain + " " + std::to_string(resnum);
      if (icode != ' ')
        err += icode;
      throw std::runtime_error(err);
    }
    return rg[0];
  }

  Connection* find_connection_by_name(const std::string& conn_name) {
    return impl::find_or_null(connections, conn_name);
  }

  CRA find_cra(const AtomAddress& address) {
    for (Chain& chain : chains)
      if ((address.use_auth_name ? chain.auth_name : chain.name) ==
          address.chain_name)
        if (Residue* res = chain.find_residue(address.res_id)) {
          Atom *at = res->find_atom(address.atom_name, address.altloc);
          return {&chain, res, at};
        }
    return {nullptr, nullptr, nullptr};
  }

  const_CRA find_cra(const AtomAddress& address) const {
    CRA cra = const_cast<Model*>(this)->find_cra(address);
    return {cra.chain, cra.residue, cra.atom};
  }


  Atom* find_atom(const AtomAddress& address) { return find_cra(address).atom; }

  std::array<int, 3> get_indices(const Chain* c, const Residue* r,
                                 const Atom* a) const {
    return {{ c      ? static_cast<int>(c - chains.data()) : -1,
              c && r ? static_cast<int>(r - c->residues.data()) : -1,
              r && a ? static_cast<int>(a - r->atoms.data()) : -1 }};
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
  UnitCell cell;
  std::string spacegroup_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::map<std::string, Entity> entities;

  // Store ORIGXn / _database_PDB_matrix.origx*
  bool has_origx = false;
  Transform origx;

  // Minimal metadata with keys being mmcif tags: _entry.id, _exptl.method, ...
  std::map<std::string, std::string> info;
  // original REMARK records stored if the file was read from the PDB format
  std::vector<std::string> raw_remarks;
  // simplistic resolution value from/for REMARK 2
  double resolution = 0;

  CoorFormat input_format = CoorFormat::Unknown;

  const std::string& get_info(const std::string& tag) const {
    static const std::string empty;
    auto it = info.find(tag);
    return it != info.end() ? it->second : empty;
  }
  Model* find_model(const std::string& model_name) {
    return impl::find_or_null(models, model_name);
  }
  Model& find_or_add_model(const std::string& model_name) {
    return impl::find_or_add(models, model_name);
  }

  void remove_model(const std::string& model_name) {
    models.erase(impl::find_iter(models, model_name));
  }

  Entity* get_entity(const std::string& ent_id) {
    auto ent = entities.find(ent_id);
    return ent == entities.end() ? nullptr : &ent->second;
  }
  const Entity* get_entity(const std::string& ent_id) const {
    auto ent = entities.find(ent_id);
    return ent == entities.end() ? nullptr : &ent->second;
  }

  Entity* get_entity_of(const Chain& chain) {
    return get_entity(chain.entity_id);
  }
  const Entity* get_entity_of(const Chain& chain) const {
    return get_entity(chain.entity_id);
  }

  Entity& find_or_add_entity(const std::string& ent_id) {
    auto ent = entities.find(ent_id);
    if (ent == entities.end())
      ent = entities.emplace(ent_id, Entity()).first;
    return ent->second;
  }

  double get_ncs_multiplier() const {
    int given = std::count_if(ncs.begin(), ncs.end(),
                              [](const NcsOp& o) { return o.given; });
    return (ncs.size() + 1.0) / (given + 1.0);  // +1 b/c identity not included
  }

  std::vector<Model>& children() { return models; }
  const std::vector<Model>& children() const { return models; }
  void setup_cell_images();
};

inline bool ResidueId::matches(const ResidueId& rid) const {
  return same_seq_id(rid) && segment == rid.segment && name == rid.name;
}

inline ResidueGroup Chain::find_residue_group(int seqnum, char icode) {
  auto match = [&](const Residue& r) {
    return r.seq_num == seqnum && (r.icode | 0x20) == (icode | 0x20);
  };
  auto begin_ = std::find_if(residues.begin(), residues.end(), match);
  auto end_ = std::find_if_not(begin_, residues.end(), match);
  return ResidueGroup{begin_, end_};
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

inline void Chain::append_residues(std::vector<Residue> new_resi) {
  int seqnum = 0;
  for (const Residue& res : residues)
    seqnum = std::max({seqnum, int(res.label_seq), int(res.seq_num)});
  for (Residue& res : new_resi) {
    res.label_seq = res.seq_num = ++seqnum;
    res.icode = ' ';
  }
  std::move(new_resi.begin(), new_resi.end(), std::back_inserter(residues));
}

inline void Structure::setup_cell_images() {
  if (const SpaceGroup* sg = find_spacegroup_by_name(spacegroup_hm)) {
    for (Op op : sg->operations()) {
      if (op == Op::identity())
        continue;
      Mat33 rot(
        double(op.rot[0][0]), double(op.rot[0][1]), double(op.rot[0][2]),
        double(op.rot[1][0]), double(op.rot[1][1]), double(op.rot[1][2]),
        double(op.rot[2][0]), double(op.rot[2][1]), double(op.rot[2][2]));
      double mult = 1.0 / Op::TDEN;
      Vec3 tran(mult * op.tran[0], mult * op.tran[1], mult * op.tran[2]);
      cell.images.emplace_back(rot, tran);
    }
  }
  // Strict NCS from MTRIXn.
  size_t n = cell.images.size();
  for (const NcsOp& op : ncs)
    if (!op.given) {
      // We need it to operates on fractional, not orthogonal coordinates.
      FTransform f = cell.frac.combine(op.tr.combine(cell.orth));
      cell.images.emplace_back(f);
      for (size_t i = 0; i < n; ++i)
        cell.images.emplace_back(cell.images[i].combine(f));
    }
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
