// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if
#include <array>
#include <bitset>
#include <iterator>   // for back_inserter
#include <map>        // for map
#include <stdexcept>  // for out_of_range
#include <string>
#include <vector>

#include "elem.hpp"
#include "fail.hpp"      // for fail
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "metadata.hpp"
#include "iterator.hpp"
#include "span.hpp"      // for Span, MutableVectorSpan
#include "seqid.hpp"
#include "util.hpp"      // for join_str, vector_move_extend, in_vector

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

template<typename Span, typename T = typename Span::value_type>
typename Span::iterator find_iter(Span& span, const std::string& name) {
  auto i = std::find_if(span.begin(), span.end(),
                        [&](const T& x) { return x.name == name; });
  if (i == span.end())
    throw std::invalid_argument(
        T::what() + (" " + name) + " not found (only [" +
        join_str(span.begin(), span.end(), ' ',
                 [](const T& x) { return x.name; }) +
        "])");
  return i;
}

template<typename Group>
typename Group::element_type& get_by_altloc(Group& group, char alt) {
  for (auto& atom : group)
    if (atom.altloc == alt)
      return atom;
  fail("No such altloc");
}

template<typename T, typename M> std::vector<T> model_subchains(M* model) {
  std::vector<T> v;
  for (auto& chain : model->chains)
    vector_move_extend(v, chain.subchains());
  return v;
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

// values corresponding to mmCIF _entity_poly.type
enum class PolymerType : unsigned char {
  Unknown,       // unknown or not applicable
  PeptideL,      // polypeptide(L) in mmCIF (168923 values in the PDB in 2017)
  PeptideD,      // polypeptide(D) (57 values)
  Dna,           // polydeoxyribonucleotide (9905)
  Rna,           // polyribonucleotide (4559)
  DnaRnaHybrid,  // polydeoxyribonucleotide/polyribonucleotide hybrid (156)
  SaccharideD,   // polysaccharide(D) (18)
  SaccharideL,   // polysaccharide(L) (0)
  Pna,           // peptide nucleic acid (2)
  CyclicPseudoPeptide,  // cyclic-pseudo-peptide (1)
  Other,         // other (4)
};

inline bool is_polypeptide(PolymerType pt) {
  return pt == PolymerType::PeptideL || pt == PolymerType::PeptideD;
}

inline bool is_polynucleotide(PolymerType pt) {
  return pt == PolymerType::Dna || pt == PolymerType::Rna ||
         pt == PolymerType::DnaRnaHybrid;
}

struct Entity {
  std::string name;
  std::vector<std::string> subchains;
  EntityType entity_type = EntityType::Unknown;
  PolymerType polymer_type = PolymerType::Unknown;
  // SEQRES or entity_poly_seq with microheterogeneity as comma-separated names
  std::vector<std::string> full_sequence;

  explicit Entity(std::string name_) noexcept : name(name_) {}
  static std::string first_mon(const std::string& mon_list) {
    return mon_list.substr(0, mon_list.find(','));
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
  bool has_altloc() const { return altloc != '\0'; }
  bool has_anisou() const {
    return u11 != 0.f || u22 != 0.f || u33 != 0.f ||
           u12 != 0.f || u13 != 0.f || u23 != 0.f;
  }
  double b_iso_from_aniso() const {
    return 8 * pi() * pi() / 3. * (u11 + u22 + u33);
  }
  bool is_hydrogen() const { return gemmi::is_hydrogen(element); }
};

struct AtomGroup : MutableVectorSpan<Atom> {
  using MutableVectorSpan::MutableVectorSpan;
  std::string name() const { return size() != 0 ? front().name : ""; }
  Atom& by_altloc(char alt) { return impl::get_by_altloc(*this, alt); }
};

struct ConstAtomGroup : Span<const Atom> {
  ConstAtomGroup(const Atom* begin, size_t n) : Span(begin, n) {}
  ConstAtomGroup(const AtomGroup& o) : Span(o.begin(), o.size()) {}
  std::string name() const { return size() != 0 ? front().name : ""; }
  const Atom& by_altloc(char a) const { return impl::get_by_altloc(*this, a); }
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
  const Atom* find_atom(const std::string& atom_name, char altloc,
                        El el=El::X) const {
    for (const Atom& a : atoms)
      if (a.name == atom_name
          && (altloc == '*' || a.altloc == '\0' || a.altloc == altloc)
          && (el == El::X || a.element == el))
        return &a;
    return nullptr;
  }
  Atom* find_atom(const std::string& atom_name, char altloc, El el=El::X) {
    const Residue* const_this = this;
    return const_cast<Atom*>(const_this->find_atom(atom_name, altloc, el));
  }

  std::vector<Atom>::iterator find_atom_iter(const std::string& atom_name,
                                             char altloc, El el=El::X) {
    if (Atom* a = find_atom(atom_name, altloc, el))
      return atoms.begin() + (a - atoms.data());
    fail("Atom to be removed not found.");
  }

  AtomGroup get(const std::string& atom_name) {
    auto func = [&](const Atom& a) { return a.name == atom_name; };
    auto g_begin = std::find_if(atoms.begin(), atoms.end(), func);
    if (g_begin == atoms.end())
      fail("No such atom: " + atom_name);
    auto g_end = std::find_if_not(g_begin, atoms.end(), func);
    if (std::find_if(g_end, atoms.end(), func) != atoms.end())
      fail("Non-consecutive alternative location of the atom");
    return AtomGroup(atoms, &*g_begin, g_end - g_begin);
  }

  Atom& sole_atom(const std::string& atom_name) {
    AtomGroup aa = get(atom_name);
    if (aa.size() != 1)
      fail("Multiple alternative atoms " + atom_name);
    return aa[0];
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

struct ResidueGroup;
struct ConstResidueGroup;

struct ConstResidueSpan : Span<const Residue> {
  using Parent = Span<const Residue>;
  using Parent::Span;
  ConstResidueSpan(Parent&& span) : Parent(std::move(span)) {}

  int length() const {
    int length = (int) size();
    for (int n = length - 1; n > 0; --n)
      if ((begin() + n)->same_group(*(begin() + n - 1)))
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

  ConstUniqProxy<Residue, ConstResidueSpan> first_conformer() const {
    return {*this};
  }

  const std::string& subchain_id() const {
    if (this->empty())
      throw std::out_of_range("subchain_id(): empty span");
    if (this->size() > 1 && this->front().subchain != this->back().subchain)
      fail("subchain id varies");
    return this->begin()->subchain;
  }

  ConstResidueGroup find_residue_group(SeqId id) const;
};

// ResidueSpan represents consecutive residues within the same chain.
// It's used as return value of get_polymer(), get_ligands(), get_waters()
// and get_subchain().
struct ResidueSpan : MutableVectorSpan<Residue> {
  using Parent = MutableVectorSpan<Residue>;
  ResidueSpan() = default;
  ResidueSpan(Parent&& span) : Parent(std::move(span)) {}
  ResidueSpan(vector_type& v, iterator begin, std::size_t n)
    : Parent(v, begin, n) {}
  int length() const { return const_().length(); }
  SeqId::OptionalNum min_seqnum() const { return const_().min_seqnum(); }
  SeqId::OptionalNum max_seqnum() const { return const_().max_seqnum(); }
  SeqId::OptionalNum min_label_seq() const { return const_().min_label_seq(); }
  SeqId::OptionalNum max_label_seq() const { return const_().max_label_seq(); }
  UniqProxy<Residue, ResidueSpan> first_conformer() { return {*this}; }
  ConstUniqProxy<Residue, ResidueSpan> first_conformer() const { return {*this}; }
  const std::string& subchain_id() const { return const_().subchain_id(); }
  ResidueGroup find_residue_group(SeqId id);
  ConstResidueGroup find_residue_group(SeqId id) const;
private:
  ConstResidueSpan const_() const { return ConstResidueSpan(begin(), size()); }
};

// ResidueGroup represents residues with the same sequence number and insertion
// code, but different residue names. I.e. microheterogeneity.
// Usually, there is only one residue in the group.
// The residues must be consecutive.
struct ResidueGroup : ResidueSpan {
  ResidueGroup() = default;
  ResidueGroup(ResidueSpan&& span) : ResidueSpan(std::move(span)) {}
  Residue& by_resname(const std::string& name) {
    return *impl::find_iter(*this, name);
  }
  void remove_residue(const std::string& name) {
    erase(impl::find_iter(*this, name));
  }
};

struct ConstResidueGroup : ConstResidueSpan {
  ConstResidueGroup() = default;
  ConstResidueGroup(ConstResidueSpan&& sp) : ConstResidueSpan(std::move(sp)) {}
  const Residue& by_resname(const std::string& name) {
    return *impl::find_iter(*this, name);
  }
};

inline ResidueGroup ResidueSpan::find_residue_group(SeqId id) {
  return ResidueSpan(subspan([&](const Residue& r) { return r.seqid == id; }));
}


namespace impl {
template<typename T, typename Ch> std::vector<T> chain_subchains(Ch* ch) {
  std::vector<T> v;
  auto span_start = ch->residues.begin();
  for (auto i = span_start; i != ch->residues.end(); ++i)
    if (i->subchain != span_start->subchain) {
      v.push_back(ch->whole().sub(span_start, i));
      span_start = i;
    }
  v.push_back(ch->whole().sub(span_start, ch->residues.end()));
  return v;
}
} // namespace impl

struct Chain {
  static const char* what() { return "Chain"; }
  std::string name;
  std::vector<Residue> residues;

  explicit Chain(std::string cname) noexcept : name(cname) {}

  ResidueSpan whole() {
    auto begin = residues.empty() ? nullptr : &residues.at(0);
    return ResidueSpan(residues, begin, residues.size());
  }
  ConstResidueSpan whole() const {
    auto begin = residues.empty() ? nullptr : &residues.at(0);
    return ConstResidueSpan(begin, residues.size());
  }

  template<typename F> ResidueSpan get_residue_span(F&& func) {
    return whole().subspan(func);
  }
  template<typename F> ConstResidueSpan get_residue_span(F&& func) const {
    return whole().subspan(func);
  }

  ResidueSpan get_polymer() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::Polymer;
    });
  }
  ConstResidueSpan get_polymer() const {
    return const_cast<Chain*>(this)->get_polymer();
  }

  ResidueSpan get_ligands() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::NonPolymer;
    });
  }
  ConstResidueSpan get_ligands() const {
    return const_cast<Chain*>(this)->get_ligands();
  }

  ResidueSpan get_waters() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::Water;
    });
  }
  ConstResidueSpan get_waters() const {
    return const_cast<Chain*>(this)->get_waters();
  }

  ResidueSpan get_subchain(const std::string& s) {
    return get_residue_span([&](const Residue& r) { return r.subchain == s; });
  }
  ConstResidueSpan get_subchain(const std::string& s) const {
    return const_cast<Chain*>(this)->get_subchain(s);
  }

  std::vector<ResidueSpan> subchains() {
    return impl::chain_subchains<ResidueSpan>(this);
  }
  std::vector<ConstResidueSpan> subchains() const {
    return impl::chain_subchains<ConstResidueSpan>(this);
  }

  ResidueGroup find_residue_group(SeqId id) {
    return whole().find_residue_group(id);
  }
  ConstResidueGroup find_residue_group(SeqId id) const {
    return whole().find_residue_group(id);
  }

  Residue* find_residue(const ResidueId& rid);
  const Residue* find_residue(const ResidueId& rid) const {
    return const_cast<Chain*>(this)->find_residue(rid);
  }

  Residue* find_or_add_residue(const ResidueId& rid);
  void append_residues(std::vector<Residue> new_resi, int min_sep=0);
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }

  // Returns false only for alternative conformation (microheterogeneity).
  bool is_first_in_group(const Residue& res) const {
    return &res == residues.data() || !(&res - 1)->same_group(res);
  }

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
  AtomAddress(const std::string& ch, const SeqId& seqid, const std::string& res,
              const std::string& atom, char alt='\0')
    : chain_name(ch), res_id({seqid, "", res}), atom_name(atom), altloc(alt) {}
  AtomAddress(const Chain& ch, const Residue& res, const Atom& at)
    : chain_name(ch.name), res_id(res), atom_name(at.name), altloc(at.altloc) {}

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

inline bool atom_matches(const const_CRA& cra, const AtomAddress& addr) {
  return cra.chain && cra.chain->name == addr.chain_name &&
         cra.residue && cra.residue->matches(addr.res_id) &&
         cra.atom && cra.atom->name == addr.atom_name &&
         cra.atom->altloc == addr.altloc;
}

// A connection. Corresponds to _struct_conn.
// Symmetry operators are not trusted and not stored.
// We assume that the nearest symmetry mate is connected.
struct Connection {
  enum Type { Covale, Disulf, Hydrog, MetalC, None };
  std::string name;
  std::string link_id;  // _struct_conn.ccp4_link_id (== _chem_link.id)
  Type type = None;
  Asu asu = Asu::Any;
  AtomAddress partner1, partner2;
  double reported_distance = 0.0;
};

inline const char* get_mmcif_connection_type_id(Connection::Type t) {
  static constexpr const char* type_ids[] = {
    "covale", "disulf", "hydrog", "metalc", "." };
  return type_ids[t];
}

// Secondary structure. PDBx/mmCIF stores helices and sheets separately.

// mmCIF spec defines 32 possible values for _struct_conf.conf_type_id -
// "the type of the conformation of the backbone of the polymer (whether
// protein or nucleic acid)". But as of 2019 only HELX_P is used (not counting
// TURN_P that occurs in only 6 entries). The actual helix type is given
// by numeric value of _struct_conf.pdbx_PDB_helix_class, which corresponds
// to helixClass from the PDB HELIX record. These values are in the range 1-10.
// As of 2019 it's almost only type 1 and 5:
// 3116566 of  1 - right-handed alpha
//      16 of  2 - right-handed omega
//      84 of  3 - right-handed pi
//      79 of  4 - right-handed gamma
// 1063337 of  5 - right-handed 3-10
//      27 of  6 - left-handed alpha
//       5 of  7 - left-handed omega
//       2 of  8 - left-handed gamma
//       8 of  9 - 2-7 ribbon/helix
//      46 of 10 - polyproline
struct Helix {
  enum HelixClass {
    UnknownHelix, RAlpha, ROmega, RPi, RGamma, R310,
    LAlpha, LOmega, LGamma, Helix27, HelixPolyProlineNone
  };
  AtomAddress start, end;
  HelixClass pdb_helix_class = UnknownHelix;
  int length = -1;
  void set_helix_class_as_int(int n) {
    if (n >= 1 && n <= 10)
      pdb_helix_class = static_cast<HelixClass>(n);
  }
};

struct Sheet {
  struct Strand {
    AtomAddress start, end;
    AtomAddress hbond_atom2, hbond_atom1;
    int sense;  // 0 = first strand, 1 = parallel, -1 = anti-parallel.
    std::string name; // optional, _struct_sheet_range.id if from mmCIF
  };
  std::string name;
  std::vector<Strand> strands;
  explicit Sheet(std::string sheet_id) noexcept : name(sheet_id) {}
};


template<typename CraT>
class CraIterPolicy {
public:
  typedef CraT value_type;
  CraIterPolicy() : chains_end(nullptr), cra{nullptr, nullptr, nullptr} {}
  CraIterPolicy(const Chain* end, CraT cra_)
    : chains_end(end), cra(cra_) {}
  void increment() {
    if (cra.atom)
      if (++cra.atom == vector_end_ptr(cra.residue->atoms)) {
        if (++cra.residue == vector_end_ptr(cra.chain->residues)) {
          if (++cra.chain == chains_end) {
            cra.atom = nullptr;
            return;
          }
          cra.residue = &cra.chain->residues.at(0);
        }
        cra.atom = &cra.residue->atoms.at(0);
      }
  }
  void decrement() {
    if (cra.atom)
      if (cra.atom-- == cra.residue->atoms.data()) {
        if (cra.residue-- == cra.chain->residues.data())
          cra.residue = &(--cra.chain)->residues.back();
        cra.atom = &cra.residue->atoms.back();
      }
  }
  bool equal(const CraIterPolicy& o) const { return cra.atom == o.cra.atom; }
  CraT& dereference() { return cra; }
  using const_policy = CraIterPolicy<const_CRA>;
  operator const_policy() const { return const_policy(chains_end, cra); }
private:
  const Chain* chains_end;
  CraT cra;
};

struct ConstCraProxy {
  const std::vector<Chain>& chains;
  using iterator = BidirIterator<CraIterPolicy<const_CRA>>;
  iterator begin() {
    auto chain = &chains.at(0);
    auto residue = &chain->residues.at(0);
    auto atom = &residue->atoms.at(0);
    return CraIterPolicy<const_CRA>{vector_end_ptr(chains), {chain, residue, atom}};
  }
  iterator end() { return {}; }
};

struct Model {
  static const char* what() { return "Model"; }
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
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
  ConstResidueSpan get_subchain(const std::string& sub_name) const {
    return const_cast<Model*>(this)->get_subchain(sub_name);
  }

  std::vector<ResidueSpan> subchains() {
    return impl::model_subchains<ResidueSpan>(this);
  }
  std::vector<ConstResidueSpan> subchains() const {
    return impl::model_subchains<ConstResidueSpan>(this);
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

  std::vector<std::string> get_all_residue_names() const {
    std::vector<std::string> names;
    for (const Chain& chain : chains)
      for (const Residue& res : chain.residues)
        if (!in_vector(res.name, names))
          names.push_back(res.name);
    return names;
  }

  CRA find_cra(const AtomAddress& address) {
    for (Chain& chain : chains)
      if (chain.name == address.chain_name)
        if (Residue* res = chain.find_residue(address.res_id)) {
          Atom *at = nullptr;
          if (!address.atom_name.empty())
            at = res->find_atom(address.atom_name, address.altloc);
          return {&chain, res, at};
        }
    return {nullptr, nullptr, nullptr};
  }

  const_CRA find_cra(const AtomAddress& address) const {
    return const_cast<Model*>(this)->find_cra(address);
  }

  ConstCraProxy all() const { return {chains}; }

  Atom* find_atom(const AtomAddress& address) { return find_cra(address).atom; }

  std::array<int, 3> get_indices(const Chain* c, const Residue* r,
                                 const Atom* a) const {
    return {{ c      ? static_cast<int>(c - chains.data()) : -1,
              c && r ? static_cast<int>(r - c->residues.data()) : -1,
              r && a ? static_cast<int>(a - r->atoms.data()) : -1 }};
  }

  std::bitset<(size_t)El::END> present_elements() const {
    std::bitset<(size_t)El::END> table;
    for (const Chain& chain : chains)
      for (const Residue& res : chain.residues)
        for (const Atom& a : res.atoms)
          table.set((int)a.element.elem);
    return table;
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

// bioassembly / biomolecule
struct Assembly {
  struct Oper {
    std::string name; // optional
    std::string type; // optional (from mmCIF only)
    Transform transform;
  };
  struct Gen {
    std::vector<std::string> chains;
    std::vector<std::string> subchains;
    std::vector<Oper> opers;
  };
  enum class SpecialKind {
    NA, CompleteIcosahedral, RepresentativeHelical, CompletePoint
  };
  std::string name;
  bool author_determined = false;
  bool software_determined = false;
  SpecialKind special_kind = SpecialKind::NA;
  int oligomeric_count = 0;
  std::string oligomeric_details;
  std::string software_name;
  double absa = NAN; // TOTAL BURIED SURFACE AREA: ... ANGSTROM**2
  double ssa = NAN;  // SURFACE AREA OF THE COMPLEX: ... ANGSTROM**2
  double more = NAN; // CHANGE IN SOLVENT FREE ENERGY: ... KCAL/MOL
  std::vector<Gen> generators;
  Assembly(const std::string& name_) : name(name_) {}
};

inline const Entity* get_entity_of(const ConstResidueSpan& sub,
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
  std::vector<Connection> connections;
  std::vector<Helix> helices;
  std::vector<Sheet> sheets;
  std::vector<Assembly> assemblies;
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

  const SpaceGroup* find_spacegroup() const {
    return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma);
  }

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

  void renumber_models() {
    for (size_t i = 0; i != models.size(); ++i)
      models[i].name = std::to_string(i+1);
  }

  Entity* get_entity(const std::string& ent_id) {
    return impl::find_or_null(entities, ent_id);
  }
  const Entity* get_entity(const std::string& ent_id) const {
    return const_cast<Structure*>(this)->get_entity(ent_id);
  }

  const Entity* get_entity_of(const ConstResidueSpan& sub) const {
    return gemmi::get_entity_of(sub, entities);
  }
  Entity* get_entity_of(const ConstResidueSpan& sub) {
    return const_cast<Entity*>(gemmi::get_entity_of(sub, entities));
  }

  Connection* find_connection_by_name(const std::string& conn_name) {
    return impl::find_or_null(connections, conn_name);
  }

  Connection* find_connection_by_cra(const const_CRA& cra1,
                                     const const_CRA& cra2) {
    for (Connection& c : connections)
      if ((atom_matches(cra1, c.partner1) && atom_matches(cra2, c.partner2)) ||
          (atom_matches(cra1, c.partner2) && atom_matches(cra2, c.partner1)))
        return &c;
    return nullptr;
  }

  Connection* find_connection(const AtomAddress& a1, const AtomAddress& a2) {
    for (Connection& c : connections)
      if ((a1 == c.partner1 && a2 == c.partner2) ||
          (a1 == c.partner2 && a2 == c.partner1))
        return &c;
    return nullptr;
  }

  double get_ncs_multiplier() const {
    size_t given = std::count_if(ncs.begin(), ncs.end(),
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
  if (new_resi.empty())
    return;
  if (min_sep > 0) {
    ConstResidueSpan new_span(&new_resi[0], new_resi.size());
    // adjust sequence numbers if necessary
    auto diff = new_span.min_seqnum() - whole().max_seqnum();
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.seqid.num += min_sep - int(diff);
    // adjust label_seq_id if necessary
    diff = new_span.min_label_seq() - whole().max_label_seq();
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.label_seq += min_sep - int(diff);
  }
  std::move(new_resi.begin(), new_resi.end(), std::back_inserter(residues));
}

inline void Structure::setup_cell_images() {
  const SpaceGroup* sg = find_spacegroup();
  cell.set_cell_images_from_spacegroup(sg);

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
