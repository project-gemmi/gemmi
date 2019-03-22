// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if
#include <array>
#include <iterator>   // for back_inserter
#include <map>        // for map
#include <stdexcept>  // for out_of_range
#include <string>
#include <vector>

#include "elem.hpp"
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "metadata.hpp"
#include "iterator.hpp"
#include "seqid.hpp"

namespace gemmi {

namespace impl {

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

template<typename Iter, typename T = typename Iter::value_type>
Iter find_iter(Iter begin, Iter end, const std::string& name) {
  Iter i = std::find_if(begin, end, [&](const T& x) { return x.name == name; });
  if (i == end)
    throw std::invalid_argument(
        T::what() + (" " + name) + " not found (only [" +
        join_str(begin, end, ' ', [](const T& x) { return x.name; }) +
        "])");
  return i;
}

template<typename T>
typename std::vector<T>::iterator find_iter(std::vector<T>& vec,
                                            const std::string& name) {
  return find_iter(vec.begin(), vec.end(), name);
}

} // namespace impl


// File format with macromolecular model.
// Unknown = unknown coordinate format (not ChemComp)
// UnknownAny = any format (coordinate file for a monomer/ligand/chemcomp)
enum class CoorFormat { Unknown, UnknownAny, Pdb, Mmcif, Mmjson, ChemComp };

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

inline bool is_polypeptide(PolymerType pt) {
  return pt == PolymerType::PeptideL || pt == PolymerType::PeptideD;
}

inline bool is_polynucleotide(PolymerType pt) {
  return pt == PolymerType::Dna || pt == PolymerType::Rna ||
         pt == PolymerType::DnaRnaHybrid;
}

struct PolySeqItem {
  int num;  // _entity_poly_seq.num or -1 if from SEQRES
  std::string mon;  // _entity_poly_seq.mon_id or resName from SEQRES
  explicit PolySeqItem(std::string m) noexcept : num(-1), mon(m) {}
  PolySeqItem(int n, std::string m) noexcept : num(n), mon(m) {}
  bool operator==(const PolySeqItem& o) const {
    return num == o.num && mon == o.mon;
  }
};

struct Entity {
  std::string name;
  std::vector<std::string> subchains;
  EntityType entity_type = EntityType::Unknown;
  PolymerType polymer_type = PolymerType::Unknown;
  std::vector<PolySeqItem> poly_seq;  // SEQRES / entity_poly_seq

  explicit Entity(std::string name_) noexcept : name(name_) {}

  //ConstUniqProxy<PolySeqItem>
  //seq_first_conformer() const { return {poly_seq}; }

  // TODO: is it worth to use first_conformer UniqProxy
  bool is_seq_first_conformer(int idx) const {
    int num = poly_seq[idx].num;
    return num < 0 || idx == 0 || num != poly_seq[idx-1].num;
  }

  // handles point mutations, unlike poly_seq.size()
  int seq_length() const {
    int len = 0;
    for (size_t i = 0; i != poly_seq.size(); ++i)
      if (is_seq_first_conformer(i))
        ++len;
    return len;
  }
};

inline bool is_same_conformer(char altloc1, char altloc2) {
  return altloc1 == '\0' || altloc2 == '\0' || altloc1 == altloc2;
}

struct Atom {
  static const char* what() { return "Atom"; }
  std::string name;
  char altloc = '\0'; // 0 if not set
  signed char charge = 0;  // [-8, +8]
  Element element = El::X;
  char flag = '\0'; // custom flag
  int serial = 0;
  Position pos;
  float occ = 1.0f;
  float b_iso = 20.0f; // arbitrary default value
  float u11=0, u22=0, u33=0, u12=0, u13=0, u23=0;

  char altloc_or(char null_char) const { return altloc ? altloc : null_char; }
  bool same_conformer(const Atom& other) const {
    return is_same_conformer(altloc, other.altloc);
  }
  // same_group() is for use in UniqIter
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
  bool has_anisou() const {
    return u11 != 0.f || u22 != 0.f || u33 != 0.f ||
           u12 != 0.f || u13 != 0.f || u23 != 0.f;
  }
  bool is_hydrogen() const { return gemmi::is_hydrogen(element); }
};

// Sequence ID (sequence number + insertion code) + residue name + segment ID
struct ResidueId {
  SeqId seqid;
  std::string segment; // segid - up to 4 characters in the PDB file
  std::string name;

  // checks for equality of sequence ID; used for first_conformation iterators
  bool same_group(const ResidueId& o) const { return seqid == o.seqid; }

  bool matches(const ResidueId& o) const {
    return seqid == o.seqid && segment == o.segment && name == o.name;
  }
  std::string str() const { return seqid.str() + "(" + name + ")"; }
};

struct Residue : public ResidueId {
  using OptionalNum = SeqId::OptionalNum;
  static const char* what() { return "Residue"; }

  std::string subchain;   // mmCIF _atom_site.label_asym_id
  OptionalNum label_seq;  // mmCIF _atom_site.label_seq_id
  EntityType entity_type = EntityType::Unknown;
  char het_flag = '\0';   // 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  bool is_cis = false;    // bond to the next residue marked as cis
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
      if (a.name == atom_name
          && (altloc == '*' || a.altloc == '\0' || a.altloc == altloc)
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

  // convenience function that duplicates functionality from resinfo.hpp
  bool is_water() const {
    if (name.length() != 3)
      return false;
    int id = ialpha4_id(name.c_str());
    return id == ialpha4_id("HOH") || id == ialpha4_id("DOD") ||
           id == ialpha4_id("WAT") || id == ialpha4_id("H2O");
  }

  // Iterators that in case of multiple conformations (alt. locations)
  // skip all but the first conformation.
  UniqProxy<Atom> first_conformer() { return {atoms}; }
  ConstUniqProxy<Atom> first_conformer() const { return {atoms}; }
};

// ResidueSpan represents consecutive residues within the same chain.
// It's used as return value of get_polymer(), get_ligands(), get_waters()
// and get_subchain().
struct ResidueSpan {
  using iterator = std::vector<Residue>::iterator;
  using const_iterator = std::vector<Residue>::const_iterator;

  iterator begin_;
  std::size_t size_ = 0;
  std::vector<Residue>* residues_ = nullptr;  // for remove_residue()

  ResidueSpan() = default;
  ResidueSpan(std::vector<Residue>& v, iterator begin, std::size_t n)
    : begin_(begin), size_(n), residues_(&v) {}
  ResidueSpan(std::vector<Residue>& v)
    : begin_(v.begin()), size_(v.size()), residues_(&v) {}

  const_iterator begin() const { return begin_; }
  const_iterator end() const { return begin_ + size_; }
  iterator begin() { return begin_; }
  iterator end() { return begin_ + size_; }
  Residue& front() { return *begin_; }
  const Residue& front() const { return *begin_; }
  Residue& back() { return *(begin_ + size_ - 1); }
  const Residue& back() const { return *(begin_ + size_ - 1); }
  std::size_t size() const { return size_; }
  int length() const {
    int length = size_;
    for (int n = int(size_) - 1; n > 0; --n)
      if ((begin_ + n)->same_group(*(begin_ + n - 1)))
        --length;
    return length;
  }

  SeqId::OptionalNum extreme_num(bool label, bool min) const {
    SeqId::OptionalNum result;
    for (const Residue& r : *this) {
      if (auto num = label ? r.label_seq : r.seqid.num)
        if (!result || (min ? int(num) < int(result) : int(num) > int(result)))
          result = num;
    }
    return result;
  }
  SeqId::OptionalNum min_seqnum() const { return extreme_num(false, true); }
  SeqId::OptionalNum max_seqnum() const { return extreme_num(false, false); }
  SeqId::OptionalNum min_label_seq() const { return extreme_num(true, true); }
  SeqId::OptionalNum max_label_seq() const { return extreme_num(true, false); }

  const Residue& operator[](std::size_t i) const { return *(begin_ + i); }
  Residue& operator[](std::size_t i) { return *(begin_ + i); }
  Residue& at(std::size_t i) {
    if (i >= size())
      throw std::out_of_range("ResidueSpan: no item " + std::to_string(i));
    return *(begin_ + i);
  }
  const Residue& at(std::size_t i) const {
    return const_cast<ResidueSpan*>(this)->at(i);
  }
  bool empty() const { return size_ == 0; }
  explicit operator bool() const { return size_ != 0; }

  UniqProxy<Residue, ResidueSpan> first_conformer() { return {*this}; }
  ConstUniqProxy<Residue, ResidueSpan> first_conformer() const {return {*this};}

  const std::string& subchain_id() const {
    if (size_ == 0)
      throw std::out_of_range("No ResidueSpan::subchain_id() for empty span");
    return begin_->subchain;
  }
};

// ResidueGroup represents residues with the same sequence number and insertion
// code, but different residue names. I.e. microheterogeneity.
// Usually, there is only one residue in the group.
// The residues must be consecutive.
struct ResidueGroup : ResidueSpan {
  ResidueGroup() = default;
  ResidueGroup(ResidueSpan&& span) : ResidueSpan(std::move(span)) {}
  Residue& by_resname(const std::string& name) {
    return *impl::find_iter(begin_, begin_ + size_, name);
  }
  void remove_residue(const std::string& name) {
    residues_->erase(impl::find_iter(begin_, begin_ + size_, name));
  }
};

struct Chain {
  static const char* what() { return "Chain"; }
  std::string name;
  std::vector<Residue> residues;

  explicit Chain(std::string cname) noexcept : name(cname) {}

  template<typename T> ResidueSpan get_residue_span(T&& func) {
    auto begin = std::find_if(residues.begin(), residues.end(), func);
    auto size = std::find_if_not(begin, residues.end(), func) - begin;
    return ResidueSpan(residues, begin, size);
  }

  ResidueSpan whole() { return ResidueSpan(residues); }
  const ResidueSpan whole() const { return const_cast<Chain*>(this)->whole(); }

  ResidueSpan get_polymer() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::Polymer;
    });
  }
  const ResidueSpan get_polymer() const {
    return const_cast<Chain*>(this)->get_polymer();
  }

  ResidueSpan get_ligands() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::NonPolymer;
    });
  }
  const ResidueSpan get_ligands() const {
    return const_cast<Chain*>(this)->get_ligands();
  }

  ResidueSpan get_waters() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::Water;
    });
  }
  const ResidueSpan get_waters() const {
    return const_cast<Chain*>(this)->get_waters();
  }

  ResidueSpan get_subchain(const std::string& s) {
    return get_residue_span([&](const Residue& r) { return r.subchain == s; });
  }
  const ResidueSpan get_subchain(const std::string& s) const {
    return const_cast<Chain*>(this)->get_subchain(s);
  }

  std::vector<ResidueSpan> subchains() {
    std::vector<ResidueSpan> v;
    for (auto i = residues.begin(); i != residues.end(); i += v.back().size())
      v.emplace_back(get_residue_span([&](const Residue& r) {
            return r.subchain == i->subchain;
      }));
    return v;
  }

  ResidueGroup find_residue_group(SeqId id) {
    return get_residue_span([&](const Residue& r) { return r.seqid == id; });
  }

  Residue* find_residue(const ResidueId& rid);
  Residue* find_or_add_residue(const ResidueId& rid);
  void append_residues(std::vector<Residue> new_resi, int min_sep=0);
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }

  // Returns the previous residue or nullptr.
  // Got complicated by handling of multi-conformations / microheterogeneity.
  const Residue* previous_residue(const Residue& res) const {
    const Residue* start = residues.data();
    for (const Residue* p = &res; p != start; )
      if (!res.same_group(*--p)) {
        while (p != start && p->same_group(*(p-1)) &&
               (res.atoms.at(0).altloc == '\0' || !res.same_conformer(*p)))
          --p;
        return p;
      }
    return nullptr;
  }

  // Returns the next residue or nullptr.
  const Residue* next_residue(const Residue& res) const {
    const Residue* end = residues.data() + residues.size();
    for (const Residue* p = &res + 1; p != end; ++p)
      if (!res.same_group(*p)) {
        while (p+1 != end && p->same_group(*(p+1)) && !res.same_conformer(*p))
          ++p;
        return p;
      }
    return nullptr;
  }

  const Residue* previous_bonded_aa(const Residue& res) const {
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

inline std::string atom_str(const std::string& chain_name,
                            const ResidueId& res_id,
                            const std::string& atom_name,
                            char altloc) {
  std::string r = chain_name;
  r += '/';
  r += res_id.name;
  r += ' ';
  r += res_id.seqid.str();
  r += '/';
  r += atom_name;
  if (altloc) {
    r += '.';
    r += altloc;
  }
  return r;
}

inline std::string atom_str(const Chain& chain,
                            const ResidueId& res_id,
                            const Atom& atom) {
  return atom_str(chain.name, res_id, atom.name, atom.altloc);
}

struct AtomAddress {
  std::string chain_name;
  ResidueId res_id;
  std::string atom_name;
  char altloc = '\0';

  AtomAddress() = default;
  AtomAddress(const Chain& ch, const Residue& res, const Atom& at)
    : chain_name(ch.name), res_id(res), atom_name(at.name), altloc(at.altloc)
  {}

  bool operator==(const AtomAddress& o) const {
    return chain_name == o.chain_name && res_id.matches(o.res_id) &&
           atom_name == o.atom_name && altloc == o.altloc;
  }

  std::string str() const {
    return atom_str(chain_name, res_id, atom_name, altloc);
  }
};

struct const_CRA {
  const Chain* chain;
  const Residue* residue;
  const Atom* atom;
};

struct CRA {
  Chain* chain;
  Residue* residue;
  Atom* atom;
  operator const_CRA() const { return const_CRA{chain, residue, atom}; }
};

inline std::string atom_str(const const_CRA& cra) {
  static const ResidueId null_residue_id = {};
  return atom_str(cra.chain ? cra.chain->name : "null",
                  cra.residue ? *cra.residue : null_residue_id,
                  cra.atom ? cra.atom->name : "null",
                  cra.atom ? cra.atom->altloc : '\0');
}

// A connection. Corresponds to _struct_conn.
// Symmetry operators are not trusted and not stored.
// We assume that the nearest symmetry mate is connected.
struct Connection {
  enum Type { Covale, Disulf, Hydrog, MetalC, None };
  std::string name;
  Type type = None;
  SameAsu asu = SameAsu::Any;
  AtomAddress atom[2];
  double reported_distance = 0.0;
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

  // Returns the first chain with given name, or nullptr.
  Chain* find_chain(const std::string& chain_name) {
    return impl::find_or_null(chains, chain_name);
  }

  // Returns the last chain with given name, or nullptr.
  Chain* find_last_chain(const std::string& chain_name) {
    auto it = std::find_if(chains.rbegin(), chains.rend(),
                         [&](const Chain& c) { return c.name == chain_name; });
    return it != chains.rend() ? &*it : nullptr;
  }


  void remove_chain(const std::string& chain_name) {
    vector_remove_if(chains,
                     [&](const Chain& c) { return c.name == chain_name; });
  }

  void merge_chain_parts(int min_sep=0) {
    for (auto i = chains.begin(); i != chains.end(); ++i)
      for (auto j = i + 1; j != chains.end(); ++j)
        if (i->name == j->name) {
          i->append_residues(j->residues, min_sep);
          chains.erase(j--);
        }
  }

  ResidueSpan get_subchain(const std::string& sub_name) {
    for (Chain& chain : chains)
      if (ResidueSpan sub = chain.get_subchain(sub_name))
        return sub;
    return ResidueSpan();
  }
  const ResidueSpan get_subchain(const std::string& sub_name) const {
    return const_cast<Model*>(this)->get_subchain(sub_name);
  }

  std::vector<ResidueSpan> subchains() {
    std::vector<ResidueSpan> v;
    for (Chain& chain : chains)
      vector_move_extend(v, chain.subchains());
    return v;
  }

  Residue* find_residue(const std::string& chain_name, const ResidueId& rid) {
    for (Chain& chain : chains)
      if (chain.name == chain_name)
        if (Residue* residue = chain.find_residue(rid))
          return residue;
    return nullptr;
  }

  ResidueGroup find_residue_group(const std::string& chain_name, SeqId seqid) {
    for (Chain& chain : chains)
      if (chain.name == chain_name)
        if (ResidueGroup res_group = chain.find_residue_group(seqid))
          return res_group;
    fail("No such chain or residue: " + chain_name + " " + seqid.str());
  }

  Residue& sole_residue(const std::string& chain_name, SeqId seqid) {
    ResidueSpan rr = find_residue_group(chain_name, seqid);
    if (rr.size() != 1)
      fail("Multiple residues " + chain_name + " " + seqid.str());
    return rr[0];
  }

  Connection* find_connection_by_name(const std::string& conn_name) {
    return impl::find_or_null(connections, conn_name);
  }

  CRA find_cra(const AtomAddress& address) {
    for (Chain& chain : chains)
      if (chain.name == address.chain_name)
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

inline const Entity* get_entity_of(const ResidueSpan& sub,
                                   const std::vector<Entity>& entities) {
  if (sub && !sub.subchain_id().empty())
    for (const Entity& ent : entities)
      if (in_vector(sub.subchain_id(), ent.subchains))
        return &ent;
  return nullptr;
}

struct Structure {
  std::string name;
  UnitCell cell;
  std::string spacegroup_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::vector<Entity> entities;
  Metadata meta;

  // Store ORIGXn / _database_PDB_matrix.origx*
  bool has_origx = false;
  Transform origx;

  // Minimal metadata with keys being mmcif tags: _entry.id, _cell.Z_PDB, ...
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
    return impl::find_or_null(entities, ent_id);
  }
  const Entity* get_entity(const std::string& ent_id) const {
    return const_cast<Structure*>(this)->get_entity(ent_id);
  }

  const Entity* get_entity_of(const ResidueSpan& sub) const {
    return gemmi::get_entity_of(sub, entities);
  }
  Entity* get_entity_of(const ResidueSpan& sub) {
    return const_cast<Entity*>(gemmi::get_entity_of(sub, entities));
  }

  double get_ncs_multiplier() const {
    int given = std::count_if(ncs.begin(), ncs.end(),
                              [](const NcsOp& o) { return o.given; });
    return (ncs.size() + 1.0) / (given + 1.0);  // +1 b/c identity not included
  }

  void merge_chain_parts(int min_sep=0) {
    for (Model& model : models)
      model.merge_chain_parts(min_sep);
  }

  std::vector<Model>& children() { return models; }
  const std::vector<Model>& children() const { return models; }
  void setup_cell_images();
};

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

inline void Chain::append_residues(std::vector<Residue> new_resi, int min_sep) {
  if (min_sep > 0) {
    auto diff = ResidueSpan(new_resi).min_seqnum() - whole().max_seqnum();
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.seqid.num += min_sep - int(diff);
    diff = ResidueSpan(new_resi).min_label_seq() - whole().max_label_seq();
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.label_seq += min_sep - int(diff);
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
