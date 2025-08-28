// Copyright 2017 Global Phasing Ltd.
//
// Data structures to store macromolecular structure models.

#ifndef GEMMI_MODEL_HPP_
#define GEMMI_MODEL_HPP_

#include <algorithm>  // for find_if, count_if, lower_bound
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
auto get_id(const T& m) -> decltype(m.name) { return m.name; }
template<typename T>
auto get_id(const T& m) -> decltype(m.num) { return m.num; }


template<typename Vec, typename S>
auto find_iter_(Vec& vec, const S& name) {
  return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
}

template<typename T, typename S>
T* find_or_null(std::vector<T>& vec, const S& name) {
  auto it = find_iter_(vec, name);
  return it != vec.end() ? &*it : nullptr;
}

template<typename T, typename S>
T& find_or_add(std::vector<T>& vec, const S& name) {
  if (T* ret = find_or_null(vec, name))
    return *ret;
  vec.emplace_back(name);
  return vec.back();
}

template<typename Span, typename S>
typename Span::iterator find_iter(Span& span, const S& name) {
  using T = typename Span::value_type;
  auto it = find_iter_(span, name);
  if (it == span.end())
    throw std::invalid_argument(cat(
        T::what(), ' ', name, " not found (only [",
        join_str(span.begin(), span.end(), ' ', [](const T& item) { return cat(get_id(item)); }),
        "])"));
  return it;
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


/// File format of a macromolecular model. When passed to read_structure():
/// Unknown = guess format from the extension,
/// Detect = guess format from the content.
enum class CoorFormat : unsigned char {
  Unknown, Detect, Pdb, Mmcif, Mmjson, ChemComp
};

/// corresponds to _atom_site.calc_flag in mmCIF
enum class CalcFlag : signed char {
  // NoHydrogen is the same as NotSet; it's used internally to mark atoms
  // which should not have riding hydrogens.
  // NB: add_cif_atoms() relies on this order.
  NotSet=0, NoHydrogen, Determined, Calculated, Dummy
};

/// options affecting how pdb file is read
struct PdbReadOptions {
  int max_line_length = 0;
  bool ignore_ter = false; // ignores TER records completely
  bool split_chain_on_ter = false;
  bool skip_remarks = false;
};
// end of PdbReadOptions for mol.rst

// remove empty residues from chain, empty chains from model, etc
template<class T> void remove_empty_children(T& obj) {
  using Item = typename T::child_type;
  vector_remove_if(obj.children(), [](const Item& x) { return x.children().empty(); });
}

inline bool is_same_conformer(char altloc1, char altloc2) {
  return altloc1 == '\0' || altloc2 == '\0' || altloc1 == altloc2;
}

/// Represents atom site in macromolecular structure (~100 bytes).
struct Atom {
  static const char* what() { return "Atom"; }
  std::string name;
  char altloc = '\0'; // 0 if not set
  signed char charge = 0;  // [-8, +8]
  Element element = El::X;
  CalcFlag calc_flag = CalcFlag::NotSet;  // mmCIF _atom_site.calc_flag
  char flag = '\0';  // a custom flag
  short tls_group_id = -1;
  int serial = 0;
  float fraction = 0.f;  // custom value, one use is Refmac's ccp4_deuterium_fraction
  Position pos;
  float occ = 1.0f;
  // ADP - in MX it's usual to give isotropic ADP as B and anisotropic as U
  float b_iso = 20.0f; // arbitrary default value
  SMat33<float> aniso = {0, 0, 0, 0, 0, 0};

  char altloc_or(char null_char) const { return altloc ? altloc : null_char; }
  bool same_conformer(const Atom& other) const {
    return is_same_conformer(altloc, other.altloc);
  }
  bool altloc_matches(char request) const {
    return request == '*' || altloc == '\0' || altloc == request;
  }
  // group_key() is used in UniqIter and similar tools
  const std::string& group_key() const { return name; }
  bool has_altloc() const { return altloc != '\0'; }
  double b_eq() const { return u_to_b() / 3. * aniso.trace(); }
  bool is_hydrogen() const { return gemmi::is_hydrogen(element); }

  // Name as a string left-padded like in the PDB format:
  // the first two characters make the element name.
  std::string padded_name() const {
    std::string s;
    const char* el = element.uname();
    if (el[1] == '\0' &&
        (el[0] == alpha_up(name[0]) || (is_hydrogen() && alpha_up(name[0]) == 'H')) &&
        name.size() < 4)
      s += ' ';
    s += name;
    return s;
  }

  // a method present in Atom, Residue, ... Structure - used in templates
  Atom empty_copy() const { return Atom(*this); }
};

template<typename AtomType>
struct AtomGroup_ : ItemGroup<AtomType> {
  using ItemGroup<AtomType>::ItemGroup;
  std::string name() const { return !this->empty() ? this->front().name : ""; }
  AtomType& by_altloc(char alt) {
    for (int i = 0; i != this->extent(); ++i) {
      AtomType* a = &this->front() + i;
      if (a->altloc == alt && (a->name == this->front().name))
        return *a;
    }
    fail("No such altloc");
  }
};

using AtomGroup = AtomGroup_<Atom>;
using ConstAtomGroup = AtomGroup_<const Atom>;

struct Residue : public ResidueId {
  using OptionalNum = SeqId::OptionalNum;
  static const char* what() { return "Residue"; }

  std::string subchain;   // mmCIF _atom_site.label_asym_id
  std::string entity_id;  // mmCIF _atom_site.label_entity_id
  OptionalNum label_seq;  // mmCIF _atom_site.label_seq_id
  EntityType entity_type = EntityType::Unknown;
  char het_flag = '\0';   // 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  char flag = '\0';       // custom flag
  SiftsUnpResidue sifts_unp;  // UniProt reference from SIFTS
  short group_idx = 0;        // ignore - internal variable
  std::vector<Atom> atoms;

  Residue() = default;
  explicit Residue(const ResidueId& rid) noexcept : ResidueId(rid) {}

  // copy all but atoms (children) - for use in templates
  Residue empty_copy() const {
    Residue res((ResidueId&)*this);
    res.subchain = subchain;
    res.entity_id = entity_id;
    res.label_seq = label_seq;
    res.entity_type = entity_type;
    res.het_flag = het_flag;
    res.flag = flag;
    res.sifts_unp = sifts_unp;
    return res;
  }
  using child_type = Atom;
  std::vector<Atom>& children() { return atoms; }
  const std::vector<Atom>& children() const { return atoms; }

  const Atom* find_by_element(El el) const {
    for (const Atom& a : atoms)
      if (a.element == el)
        return &a;
    return nullptr;
  }

  /// El::X means anything;
  /// in strict_altloc mode, '*' = any altloc, otherwise it's \0
  Atom* find_atom(const std::string& atom_name, char altloc, El el=El::X,
                  bool strict_altloc=true) {
    if (!strict_altloc && altloc == '\0')
      altloc = '*';
    for (Atom& a : atoms)
      if (a.name == atom_name && a.altloc_matches(altloc) && (el == El::X || a.element == el))
        return &a;
    return nullptr;
  }
  const Atom* find_atom(const std::string& atom_name, char altloc, El el=El::X,
                        bool strict_altloc=true) const {
    return const_cast<Residue*>(this)->find_atom(atom_name, altloc, el, strict_altloc);
  }

  std::vector<Atom>::iterator find_atom_iter(const std::string& atom_name,
                                             char altloc, El el=El::X) {
    if (Atom* a = find_atom(atom_name, altloc, el))
      return atoms.begin() + (a - atoms.data());
    fail("Atom not found.");
  }

  AtomGroup get(const std::string& atom_name) {
    for (Atom& atom : atoms)
      if (atom.name == atom_name)
        return AtomGroup(&atom, atoms.data() + atoms.size());
    fail("No such atom: " + atom_name);
  }

  Atom& sole_atom(const std::string& atom_name) {
    AtomGroup aa = get(atom_name);
    if (aa.size() != 1)
      fail("Multiple alternative atoms " + atom_name);
    return aa.front();
  }

  // short-cuts to access peptide backbone atoms
  const Atom* get_ca() const { return find_atom("CA", '*', El::C); }
  const Atom* get_c() const { return find_atom("C", '*', El::C); }
  const Atom* get_n() const { return find_atom("N", '*', El::N); }
  const Atom* get_o() const { return find_atom("O", '*', El::O); }
  // short-cuts to access nucleic acid atoms
  const Atom* get_p() const { return find_atom("P", '*', El::P); }
  const Atom* get_o3prim() const { return find_atom("O3'", '*', El::O); }

  bool same_conformer(const Residue& other) const {
    return atoms.empty() || other.atoms.empty() ||
           atoms[0].same_conformer(other.atoms[0]) ||
           other.find_atom(other.atoms[0].name, atoms[0].altloc) != nullptr;
  }

  /// Convenience function that duplicates functionality from resinfo.hpp.
  /// Returns true for HOH and DOD (and old alternative names of HOH),
  /// but not for OH and H3O/D3O.
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

inline void add_distinct_altlocs(const Residue& res, std::string& altlocs) {
  for (const Atom& atom : res.atoms)
    if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
      altlocs += atom.altloc;
}

struct ResidueGroup;
struct ConstResidueGroup;

struct ConstResidueSpan : Span<const Residue> {
  using Parent = Span<const Residue>;
  using Parent::Span;
  ConstResidueSpan(Parent&& span) : Parent(std::move(span)) {}

  int length() const {
    int length = (int) size();
    for (int n = length - 1; n > 0; --n)
      if ((begin() + n)->group_key() == (begin() + n - 1)->group_key())
        --length;
    return length;
  }

  // sign=-1 for min, sign=1 for max
  SeqId::OptionalNum extreme_num(bool label, int sign) const {
    SeqId::OptionalNum result;
    for (const Residue& r : *this) {
      if (auto num = label ? r.label_seq : r.seqid.num)
        if (!result || sign * int(num) > sign * int(result))
          result = num;
    }
    return result;
  }

  ConstUniqProxy<Residue, ConstResidueSpan> first_conformer() const {
    return {*this};
  }

  const std::string& subchain_id() const {
    if (this->empty())
      throw std::out_of_range("subchain_id(): empty span");
    if (this->size() > 1 && this->front().subchain != this->back().subchain)
      fail("subchain id varies in a residue span: ", this->front().subchain,
           " vs ", this->back().subchain);
    return this->begin()->subchain;
  }

  ConstResidueGroup find_residue_group(SeqId id) const;

  std::vector<std::string> extract_sequence() const {
    std::vector<std::string> seq;
    for (const Residue& res : first_conformer())
      seq.push_back(res.name);
    return seq;
  }

  // We assume residues are ordered. It works (approximately) also with
  // missing numbers which can be present in DBREF.
  SeqId label_seq_id_to_auth(SeqId::OptionalNum label_seq_id) const {
    if (empty())
      throw std::out_of_range("label_seq_id_to_auth(): empty span");
    const auto* it = std::lower_bound(begin(), end(), label_seq_id,
        [](const Residue& r, SeqId::OptionalNum v){ return r.label_seq < v; });
    if (it == end())
      --it;
    else if (it->label_seq == label_seq_id)
      return it->seqid;
    else if (it != begin() &&
             label_seq_id - (it-1)->label_seq < it->label_seq - label_seq_id)
      --it;
    return {it->seqid.num + (label_seq_id - it->label_seq), ' '};
  }
  // The residue numbers (auth) are sometimes not ordered.
  // That is why we use this multi-step heuristic.
  SeqId::OptionalNum auth_seq_id_to_label(SeqId auth_seq_id) const {
    if (empty())
      throw std::out_of_range("auth_seq_id_to_label(): empty span");
    for (const Residue& r : *this)
      if (r.seqid == auth_seq_id)
        return r.label_seq;
    const_iterator it;
    if (auth_seq_id.num < front().seqid.num) {
      it = begin();
    } else if (back().seqid.num < auth_seq_id.num) {
      it = end() - 1;
    } else {
      it = std::lower_bound(begin(), end(), auth_seq_id.num,
        [](const Residue& r, SeqId::OptionalNum v){ return r.seqid.num < v; });
      while (it != end() && it->seqid.num == auth_seq_id.num &&
             it->seqid.icode != auth_seq_id.icode)
        ++it;
      if (it == end())
        --it;
    }
    return it->label_seq + (auth_seq_id.num - it->seqid.num);
  }
};

// ResidueSpan represents consecutive residues within the same chain.
// It's used as return value of get_polymer(), get_ligands(), get_waters()
// and get_subchain().
struct ResidueSpan : MutableVectorSpan<Residue> {
  using Parent = MutableVectorSpan<Residue>;
  struct GroupingProxy;
  ResidueSpan() = default;
  ResidueSpan(Parent&& span) : Parent(std::move(span)) {}
  ResidueSpan(vector_type& v, iterator begin, std::size_t n)
    : Parent(v, begin, n) {}
  int length() const { return const_().length(); }
  SeqId::OptionalNum extreme_num(bool label, int sign) const {
    return const_().extreme_num(label, sign);
  }
  UniqProxy<Residue, ResidueSpan> first_conformer() { return {*this}; }
  ConstUniqProxy<Residue, ResidueSpan> first_conformer() const { return {*this}; }
  GroupingProxy residue_groups();
  const std::string& subchain_id() const { return const_().subchain_id(); }
  ResidueGroup find_residue_group(SeqId id);
  std::vector<std::string> extract_sequence() const { return const_().extract_sequence(); }
  ConstResidueGroup find_residue_group(SeqId id) const;
  SeqId label_seq_id_to_auth(SeqId::OptionalNum label_seq_id) const {
    return const_().label_seq_id_to_auth(label_seq_id);
  }
  SeqId::OptionalNum auth_seq_id_to_label(SeqId auth_seq_id) const {
    return const_().auth_seq_id_to_label(auth_seq_id);
  }
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
inline ConstResidueGroup ResidueSpan::find_residue_group(SeqId id) const {
  return const_().find_residue_group(id);
}
inline ConstResidueGroup ConstResidueSpan::find_residue_group(SeqId id) const {
  return ConstResidueSpan(subspan([&](const Residue& r) { return r.seqid == id; }));
}

struct ResidueSpan::GroupingProxy {
  ResidueSpan& span;
  using iterator = GroupingIter<ResidueSpan, ResidueGroup>;
  iterator begin() {
    return ++iterator{ResidueGroup(span.sub(span.begin(), span.begin()))};
  }
  iterator end() {
    return iterator{ResidueGroup(span.sub(span.end(), span.end()))};
  }
};

inline ResidueSpan::GroupingProxy ResidueSpan::residue_groups() {
  return {*this};
}


namespace impl {
template<typename T, typename Ch> std::vector<T> chain_subchains(Ch* ch) {
  std::vector<T> v;
  for (auto start = ch->residues.begin(); start != ch->residues.end(); ) {
    auto end = start + 1;
    while (end != ch->residues.end() && end->subchain == start->subchain)
      ++end;
    v.push_back(ch->whole().sub(start, end));
    start = end;
  }
  return v;
}
} // namespace impl

struct Chain {
  static const char* what() { return "Chain"; }
  std::string name;
  std::vector<Residue> residues;

  Chain() = default;
  explicit Chain(const std::string& name_) noexcept : name(name_) {}

  ResidueSpan whole() {
    Residue* begin = residues.empty() ? nullptr : &residues[0];
    return ResidueSpan(residues, begin, residues.size());
  }
  ConstResidueSpan whole() const {
    const Residue* begin = residues.empty() ? nullptr : &residues[0];
    return ConstResidueSpan(begin, residues.size());
  }

  template<typename F> ResidueSpan get_residue_span(F&& func) {
    return whole().subspan(func);
  }
  template<typename F> ConstResidueSpan get_residue_span(F&& func) const {
    return whole().subspan(func);
  }

  ResidueSpan get_polymer() {
    auto begin = residues.begin();
    while (begin != residues.end() && begin->entity_type != EntityType::Polymer)
      ++begin;
    auto end = begin;
    while (end != residues.end() && end->entity_type == EntityType::Polymer
                                 && end->subchain == begin->subchain)
      ++end;
    return ResidueSpan(residues, &*begin, end - begin);
  }
  ConstResidueSpan get_polymer() const {
    return const_cast<Chain*>(this)->get_polymer();
  }

  ResidueSpan get_ligands() {
    return get_residue_span([](const Residue& r) {
        return r.entity_type == EntityType::NonPolymer ||
               r.entity_type == EntityType::Branched;
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

  // methods present in Structure, Model, ... - used in templates
  Chain empty_copy() const { return Chain(name); }
  using child_type = Residue;
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }

  // Returns false only for alternative conformation (microheterogeneity).
  bool is_first_in_group(const Residue& res) const {
    return &res == residues.data() || (&res - 1)->group_key() != res.group_key();
  }

  // Returns the previous residue or nullptr.
  // Got complicated by handling of multi-conformations / microheterogeneity.
  const Residue* previous_residue(const Residue& res) const {
    const Residue* start = residues.data();
    for (const Residue* p = &res; p-- != start; )
      if (res.group_key() != p->group_key()) {
        while (p != start && p->group_key() == (p-1)->group_key() &&
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
      if (res.group_key() != p->group_key()) {
        while (p+1 != end && p->group_key() == (p+1)->group_key() &&
               !res.same_conformer(*p))
          ++p;
        return p;
      }
    return nullptr;
  }

  // Iterators that in case of multiple conformations (microheterogeneity)
  // skip all but the first conformation.
  UniqProxy<Residue> first_conformer() { return {residues}; }
  ConstUniqProxy<Residue> first_conformer() const { return {residues}; }
};

inline std::string atom_str(const Chain& chain,
                            const ResidueId& res_id,
                            const Atom& atom) {
  return atom_str(chain.name, res_id, atom.name, atom.altloc);
}

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

inline bool atom_matches(const const_CRA& cra, const AtomAddress& addr, bool ignore_segment=false) {
  return cra.chain && cra.chain->name == addr.chain_name &&
         cra.residue && cra.residue->matches_noseg(addr.res_id) &&
         (ignore_segment || cra.residue->segment == addr.res_id.segment) &&
         cra.atom && cra.atom->name == addr.atom_name &&
         cra.atom->altloc == addr.altloc;
}

inline AtomAddress make_address(const Chain& ch, const Residue& res, const Atom& at) {
  return AtomAddress(ch.name, res, at.name, at.altloc);
}


template<typename CraT>
class CraIterPolicy {
public:
  using value_type = CraT;
  using reference = const CraT;
  CraIterPolicy() : chains_end(nullptr), cra{nullptr, nullptr, nullptr} {}
  CraIterPolicy(const Chain* end, CraT cra_) : chains_end(end), cra(cra_) {}
  void increment() {
    if (cra.atom == nullptr)
      return;
    if (++cra.atom == vector_end_ptr(cra.residue->atoms)) {
      do {
        if (++cra.residue == vector_end_ptr(cra.chain->residues)) {
          do {
            if (++cra.chain == chains_end) {
              cra.atom = nullptr;
              return;
            }
          } while (cra.chain->residues.empty());
          cra.residue = &cra.chain->residues[0];
        }
      } while (cra.residue->atoms.empty());
      cra.atom = &cra.residue->atoms[0];
    }
  }
  void decrement() {
    while (cra.atom == nullptr || cra.atom == cra.residue->atoms.data()) {
      while (cra.residue == nullptr || cra.residue == cra.chain->residues.data()) {
        // iterating backward beyond begin() will have undefined effects
        while ((--cra.chain)->residues.empty()) {}
        cra.residue = vector_end_ptr(cra.chain->residues);
      }
      --cra.residue;
      cra.atom = vector_end_ptr(cra.residue->atoms);
    }
    --cra.atom;
  }
  bool equal(const CraIterPolicy& o) const { return cra.atom == o.cra.atom; }
  CraT dereference() { return cra; }  // make copy b/c increment() modifies cra
  using const_policy = CraIterPolicy<const_CRA>;
  operator const_policy() const { return const_policy(chains_end, cra); }
private:
  const Chain* chains_end;
  CraT cra;
};

template<typename CraT, typename ChainsRefT>
struct CraProxy_ {
  ChainsRefT chains;
  using iterator = BidirIterator<CraIterPolicy<CraT>>;
  iterator begin() {
    for (auto& chain : chains)
      for (auto& residue : chain.residues)
        for (auto& atom : residue.atoms)
          return CraIterPolicy<CraT>{vector_end_ptr(chains), CraT{&chain, &residue, &atom}};
    return {};
  }
  iterator end() {
    auto* chains_end = vector_end_ptr(chains);
    return CraIterPolicy<CraT>{chains_end, CraT{chains_end, nullptr, nullptr}};
  }
};

using CraProxy = CraProxy_<CRA, std::vector<Chain>&>;
using ConstCraProxy = CraProxy_<const_CRA, const std::vector<Chain>&>;

struct Model {
  static const char* what() { return "Model"; }
  int num = 0;
  std::vector<Chain> chains;

  Model() = default;
  explicit Model(int num_) noexcept : num(num_) {}

  // Returns the first chain with given name, or nullptr.
  Chain* find_chain(const std::string& chain_name) {
    return impl::find_or_null(chains, chain_name);
  }
  const Chain* find_chain(const std::string& chain_name) const {
    return const_cast<Model*>(this)->find_chain(chain_name);
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

  std::map<std::string, std::string> subchain_to_chain() const {
    std::map<std::string, std::string> mapping;
    for (const Chain& chain : chains) {
      std::string prev;
      for (const Residue& res : chain.residues)
        if (!res.subchain.empty() && res.subchain != prev) {
          prev = res.subchain;
          mapping[res.subchain] = chain.name;
        }
    }
    return mapping;
  }

  Residue* find_residue(const std::string& chain_name, const ResidueId& rid) {
    for (Chain& chain : chains)
      if (chain.name == chain_name)
        if (Residue* residue = chain.find_residue(rid))
          return residue;
    return nullptr;
  }
  const Residue* find_residue(const std::string& chain_name, const ResidueId& rid) const {
    return const_cast<Model*>(this)->find_residue(chain_name, rid);
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

  CRA find_cra(const AtomAddress& address, bool ignore_segment=false) {
    for (Chain& chain : chains)
      if (chain.name == address.chain_name) {
        for (Residue& res : chain.residues)
          if (address.res_id.matches_noseg(res) &&
              (ignore_segment || address.res_id.segment == res.segment)) {
            Atom *at = nullptr;
            if (!address.atom_name.empty())
              at = res.find_atom(address.atom_name, address.altloc);
            return {&chain, &res, at};
          }
      }
    return {nullptr, nullptr, nullptr};
  }

  const_CRA find_cra(const AtomAddress& address, bool ignore_segment=false) const {
    return const_cast<Model*>(this)->find_cra(address, ignore_segment);
  }

  CraProxy all() { return {chains}; }
  ConstCraProxy all() const { return {chains}; }

  Atom* find_atom(const AtomAddress& address) { return find_cra(address).atom; }
  const Atom* find_atom(const AtomAddress& address) const { return find_cra(address).atom; }

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
          table.set(a.element.ordinal());
    return table;
  }

  // methods present in Structure, Model, ... - used in templates
  Model empty_copy() const { return Model(num); }
  using child_type = Chain;
  std::vector<Chain>& children() { return chains; }
  const std::vector<Chain>& children() const { return chains; }
};

inline Entity* find_entity_of_subchain(const std::string& subchain_id,
                                       std::vector<Entity>& entities) {
  if (!subchain_id.empty())
    for (Entity& ent : entities)
      if (in_vector(subchain_id, ent.subchains))
        return &ent;
  return nullptr;
}
inline const Entity* find_entity_of_subchain(const std::string& subchain_id,
                                             const std::vector<Entity>& entities) {
  return find_entity_of_subchain(subchain_id, const_cast<std::vector<Entity>&>(entities));
}

struct Structure {
  static const char* what() { return "Structure"; }
  std::string name;
  UnitCell cell;
  std::string spacegroup_hm;  // usually pdb symbol, cf. SpaceGroup::pdb_name()
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::vector<Entity> entities;
  std::vector<Connection> connections;
  std::vector<CisPep> cispeps;
  std::vector<ModRes> mod_residues;
  std::vector<Helix> helices;
  std::vector<Sheet> sheets;
  std::vector<Assembly> assemblies;
  std::map<int, std::vector<int>> conect_map;
  Metadata meta;

  CoorFormat input_format = CoorFormat::Unknown;
  bool has_d_fraction = false;  // uses Refmac's ccp4_deuterium_fraction
  /// in input PDB file: y = TER records were read, e = errors were detected
  char ter_status = '\0';

  /// Store ORIGXn / _database_PDB_matrix.origx*
  bool has_origx = false;
  Transform origx;

  /// Minimal metadata with keys being mmcif tags: _entry.id, _cell.Z_PDB, ...
  std::map<std::string, std::string> info;
  /// Mapping of long (4+) CCD codes (residue names) to PDB-compatible ones
  std::vector<std::pair<std::string,std::string>> shortened_ccd_codes;
  /// original REMARK records stored if the file was read from the PDB format
  std::vector<std::string> raw_remarks;
  /// simplistic resolution value from/for REMARK 2
  double resolution = 0;

  const SpaceGroup* find_spacegroup() const {
    return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma);
  }

  const std::string& get_info(const std::string& tag) const {
    static const std::string empty;
    auto it = info.find(tag);
    return it != info.end() ? it->second : empty;
  }

  Model& first_model() {
    if (models.empty())
      fail("no structural models");
    return models[0];
  }
  const Model& first_model() const {
    return const_cast<Structure*>(this)->first_model();
  }

  Model* find_model(int model_num) {
    return impl::find_or_null(models, model_num);
  }
  const Model* find_model(int model_num) const {
    return const_cast<Structure*>(this)->find_model(model_num);
  }
  Model& find_or_add_model(int model_num) {
    return impl::find_or_add(models, model_num);
  }

  void renumber_models() {
    for (size_t i = 0; i != models.size(); ++i)
      models[i].num = i + 1;
  }

  Entity* get_entity(const std::string& ent_id) {
    return impl::find_or_null(entities, ent_id);
  }
  const Entity* get_entity(const std::string& ent_id) const {
    return const_cast<Structure*>(this)->get_entity(ent_id);
  }

  Entity* get_entity_of(const ConstResidueSpan& sub) {
    return sub ? find_entity_of_subchain(sub.subchain_id(), entities) : nullptr;
  }
  const Entity* get_entity_of(const ConstResidueSpan& sub) const {
    return const_cast<Structure*>(this)->get_entity_of(sub);
  }

  Assembly* find_assembly(const std::string& assembly_id) {
    return impl::find_or_null(assemblies, assembly_id);
  }

  Connection* find_connection_by_name(const std::string& conn_name) {
    return impl::find_or_null(connections, conn_name);
  }
  const Connection* find_connection_by_name(const std::string& conn_name) const {
    return const_cast<Structure*>(this)->find_connection_by_name(conn_name);
  }

  Connection* find_connection_by_cra(const const_CRA& cra1,
                                     const const_CRA& cra2,
                                     bool ignore_segment=false) {
    for (Connection& c : connections)
      if ((atom_matches(cra1, c.partner1, ignore_segment) &&
           atom_matches(cra2, c.partner2, ignore_segment)) ||
          (atom_matches(cra1, c.partner2, ignore_segment) &&
           atom_matches(cra2, c.partner1, ignore_segment)))
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

  size_t ncs_given_count() const {
    return std::count_if(ncs.begin(), ncs.end(), [](const NcsOp& o) { return o.given; });
  }
  double get_ncs_multiplier() const {
    return (ncs.size() + 1.0) / (ncs_given_count() + 1.0);  // +1 b/c identity not included
  }
  bool ncs_not_expanded() const {
    return std::any_of(ncs.begin(), ncs.end(), [](const NcsOp& o) { return !o.given; });
  }

  void add_conect_one_way(int serial_a, int serial_b, int order) {
    auto& vec = conect_map[serial_a];
    for (int i = 0; i < order; ++i)
      vec.insert(std::upper_bound(vec.begin(), vec.end(), serial_b), serial_b);
  }
  void add_conect(int serial1, int serial2, int order) {
    add_conect_one_way(serial1, serial2, order);
    add_conect_one_way(serial2, serial1, order);
  }

  void merge_chain_parts(int min_sep=0) {
    for (Model& model : models)
      model.merge_chain_parts(min_sep);
  }

  void remove_empty_chains() {
    for (Model& model : models)
      remove_empty_children(model);
  }

  // copy all but models (in general, empty_copy copies all but children)
  Structure empty_copy() const {
    Structure st;
    st.name = name;
    st.cell = cell;
    st.spacegroup_hm = spacegroup_hm;
    st.ncs = ncs;
    st.entities = entities;
    st.connections = connections;
    st.cispeps = cispeps;
    st.mod_residues = mod_residues;
    st.helices = helices;
    st.sheets = sheets;
    st.assemblies = assemblies;
    st.meta = meta;
    st.input_format = input_format;
    st.has_origx = has_origx;
    st.origx = origx;
    st.info = info;
    st.raw_remarks = raw_remarks;
    st.resolution = resolution;
    return st;
  }

  using child_type = Model;
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
    auto diff = new_span.extreme_num(false, -1) - whole().extreme_num(false, 1);
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.seqid.num += min_sep - int(diff);
    // adjust label_seq_id if necessary
    diff = new_span.extreme_num(true, -1) - whole().extreme_num(true, 1);
    if (diff && int(diff) < min_sep)
      for (Residue& res : new_resi)
        res.label_seq += min_sep - int(diff);
  }
  std::move(new_resi.begin(), new_resi.end(), std::back_inserter(residues));
}

inline void Structure::setup_cell_images() {
  const SpaceGroup* sg = find_spacegroup();
  cell.set_cell_images_from_spacegroup(sg);
  cell.add_ncs_images_to_cs_images(ncs);
}

} // namespace gemmi

#endif
