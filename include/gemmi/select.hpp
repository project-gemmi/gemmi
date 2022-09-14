// Copyright 2018 Global Phasing Ltd.
//
// Selections.

#ifndef GEMMI_SELECT_HPP_
#define GEMMI_SELECT_HPP_

#include <string>
#include <cstdlib>   // for strtol
#include <cctype>    // for isalpha
#include <climits>   // for INT_MIN, INT_MAX
#include "fail.hpp"  // for fail
#include "util.hpp"  // for is_in_list
#include "model.hpp" // for Model, Chain, etc
#include "iterator.hpp" // for FilterProxy
#include "sprintf.hpp" // for to_str
#include "atof.hpp"  // for fast_from_chars

namespace gemmi {

// from http://www.ccp4.ac.uk/html/pdbcur.html
// Specification of the selection sets:
// either
//     /mdl/chn/s1.i1-s2.i2/at[el]:aloc
// or
//     /mdl/chn/*(res).ic/at[el]:aloc
//

struct Selection {
  struct List {
    bool all = true;
    bool inverted = false;
    std::string list;  // comma-separated

    std::string str() const {
      if (all)
        return "*";
      return inverted ? "!" + list : list;
    }

    bool has(const std::string& name) const {
      if (all)
        return true;
      bool found = is_in_list(name, list);
      return inverted ? !found : found;
    }
  };

  struct FlagList {
    std::string pattern;
    bool has(char flag) const {
      if (pattern.empty())
        return true;
      bool invert = (pattern[0] == '!');
      bool found = (pattern.find(flag, invert ? 1 : 0) != std::string::npos);
      return invert ? !found : found;
    }
  };

  struct SequenceId {
    int seqnum;
    char icode;

    bool empty() const {
      return seqnum == INT_MIN || seqnum == INT_MAX;
    }

    std::string str() const {
      std::string s;
      if (!empty()) {
        s = std::to_string(seqnum);
        if (icode != '*') {
          s += '.';
          if (icode != ' ')
            s += icode;
        }
      }
      return s;
    }

    int compare(const SeqId& seqid) const {
      if (seqnum != *seqid.num)
        return seqnum < *seqid.num ? -1 : 1;
      if (icode != '*' && icode != seqid.icode)
        return icode < seqid.icode ? -1 : 1;
      return 0;
    }
  };

  struct AtomInequality {
    char property;
    int relation;
    double value;

    bool matches(const Atom& a) const {
      double atom_value = 0.;
      if (property == 'q')
        atom_value = a.occ;
      else if (property == 'b')
        atom_value = a.b_iso;
      if (relation < 0)
        return atom_value < value;
      if (relation > 0)
        return atom_value > value;
      return atom_value == value;
    }

    std::string str() const {
      std::string r = ";";
      r += property;
      r += relation == 0 ? '=' : relation < 0 ? '<' : '>';
      r += to_str(value);
      return r;
    }
  };

  int mdl = 0;            // 0 = all
  List chain_ids;
  SequenceId from_seqid = {INT_MIN, '*'};
  SequenceId to_seqid = {INT_MAX, '*'};
  List residue_names;
  List entity_types;
  // array corresponding to enum EntityType
  std::array<char, 6> et_flags;
  List atom_names;
  std::vector<char> elements;
  List altlocs;
  FlagList residue_flags;
  FlagList atom_flags;
  std::vector<AtomInequality> atom_inequalities;

  Selection() = default;
  Selection(const std::string& cid);

  std::string str() const {
    std::string cid = "/";
    if (mdl != 0)
      cid += std::to_string(mdl);
    cid += '/';
    cid += chain_ids.str();
    cid += '/';
    cid += from_seqid.str();
    if (!residue_names.all) {
      cid += '(';
      cid += residue_names.str();
      cid += ')';
    }
    if (!from_seqid.empty() || !to_seqid.empty()) {
      cid += '-';
      cid += to_seqid.str();
    }
    cid += '/';
    if (!atom_names.all)
      cid += atom_names.str();
    if (!elements.empty()) {
      cid += '[';
      bool inv = (std::count(elements.begin(), elements.end(), 1) > 64);
      if (inv)
        cid += '!';
      for (size_t i = 0; i < elements.size(); ++i)
        if (elements[i] != char(inv)) {
          cid += element_name(static_cast<El>(i));
          cid += ',';
        }
      cid.back() = ']';
    }
    if (!altlocs.all) {
      cid += ':';
      cid += altlocs.str();
    }
    if (!entity_types.all) {
      cid += ';';
      cid += entity_types.str();
    }
    for (const AtomInequality& ai : atom_inequalities)
      cid += ai.str();
    return cid;
  }

  bool matches(const Structure&) const { return true; }
  bool matches(const Model& model) const {
    return mdl == 0 || std::to_string(mdl) == model.name;
  }
  bool matches(const Chain& chain) const {
    return chain_ids.has(chain.name);
  }
  bool matches(const Residue& res) const {
    return (entity_types.all || et_flags[(int)res.entity_type]) &&
           residue_names.has(res.name) &&
           from_seqid.compare(res.seqid) <= 0 &&
           to_seqid.compare(res.seqid) >= 0 &&
           residue_flags.has(res.flag);
  }
  bool matches(const Atom& a) const {
    return atom_names.has(a.name) &&
           (elements.empty() || elements[a.element.ordinal()]) &&
           (altlocs.all || altlocs.has(std::string(a.altloc ? 1 : 0, a.altloc))) &&
           atom_flags.has(a.flag) &&
           std::all_of(atom_inequalities.begin(), atom_inequalities.end(),
                       [&](const AtomInequality& i) { return i.matches(a); });
  }
  bool matches(const CRA& cra) const {
    return (cra.chain == nullptr || matches(*cra.chain)) &&
           (cra.residue == nullptr || matches(*cra.residue)) &&
           (cra.atom == nullptr || matches(*cra.atom));
  }

  FilterProxy<Selection, Model> models(Structure& st) const {
    return {*this, st.models};
  }
  FilterProxy<Selection, Chain> chains(Model& model) const {
    return {*this, model.chains};
  }
  FilterProxy<Selection, Residue> residues(Chain& chain) const {
    return {*this, chain.residues};
  }
  FilterProxy<Selection, Atom> atoms(Residue& residue) const {
    return {*this, residue.atoms};
  }

  CRA first_in_model(Model& model) const {
    if (matches(model))
      for (Chain& chain : model.chains) {
        if (matches(chain))
          for (Residue& res : chain.residues) {
            if (matches(res))
              for (Atom& atom : res.atoms) {
                if (matches(atom))
                  return {&chain, &res, &atom};
              }
          }
      }
    return {nullptr, nullptr, nullptr};
  }

  std::pair<Model*, CRA> first(Structure& st) const {
    for (Model& model : st.models) {
      CRA cra = first_in_model(model);
      if (cra.chain)
        return {&model, cra};
    }
    return {nullptr, {nullptr, nullptr, nullptr}};
  }

  template<typename T>
  void add_matching_children(const T& orig, T& target) const {
    for (const auto& orig_child : orig.children())
      if (matches(orig_child)) {
        target.children().push_back(orig_child.empty_copy());
        add_matching_children(orig_child, target.children().back());
      }
  }
  void add_matching_children(const Atom&, Atom&) const {}

  Selection& set_residue_flags(const std::string& pattern) {
    residue_flags.pattern = pattern;
    return *this;
  }
  Selection& set_atom_flags(const std::string& pattern) {
    atom_flags.pattern = pattern;
    return *this;
  }

  template<typename T>
  T copy_selection(const T& orig) const {
    T copied = orig.empty_copy();
    add_matching_children(orig, copied);
    return copied;
  }

  template<typename T>
  void remove_selected(T& t) const {
    for (auto& child : t.children())
      if (matches(child))
        remove_selected(child);
    vector_remove_if(t.children(),
                     [&](typename T::child_type& c) { return c.children().empty(); });
  }
  void remove_selected(Residue& res) const {
    if (atom_names.all && elements.empty() && altlocs.all &&
        atom_flags.pattern.empty() && atom_inequalities.empty())
      res.atoms.clear();
    else
      vector_remove_if(res.atoms, [&](Atom& c) { return matches(c); });
  }

  template<typename T>
  void remove_not_selected(T& t) const {
    vector_remove_if(t.children(), [&](typename T::child_type& c) { return !matches(c); });
    for (auto& child : t.children())
      remove_not_selected(child);
  }
  void remove_not_selected(Atom&) const {}
};

namespace impl {

[[noreturn]]
inline GEMMI_COLD void wrong_syntax(const std::string& cid, size_t pos,
                                    const char* info=nullptr) {
  std::string msg = "Invalid selection syntax";
  if (info)
    msg += info;
  if (pos != 0) {
    msg += " near \"";
    msg += cid.substr(pos, 8);
    msg += '"';
  }
  msg += ": ";
  msg += cid;
  fail(msg);
}

inline int determine_omitted_cid_fields(const std::string& cid) {
  if (cid[0] == '/')
    return 0; // model
  if (std::isdigit(cid[0]) || cid[0] == '.' || cid[0] == '(' || cid[0] == '-')
    return 2; // residue
  size_t sep = cid.find_first_of("/([:;");
  if (sep == std::string::npos || cid[sep] == '/' || cid[sep] == ';')
    return 1; // chain
  if (cid[sep] == '(')
    return 2; // residue
  return 3;  // atom
}

inline Selection::List make_cid_list(const std::string& cid, size_t pos, size_t end) {
  Selection::List list;
  list.all = (cid[pos] == '*');
  list.inverted = (cid[pos] == '!');
  if (list.all || list.inverted)
    ++pos;
  list.list = cid.substr(pos, end - pos);
  // if a list have punctation other than ',' something must be wrong
  size_t idx = list.list.find_first_of("[]()!/*-.:;");
  if (idx != std::string::npos)
    wrong_syntax(cid, pos + idx, " in a list");
  return list;
}

inline void parse_cid_elements(const std::string& cid, size_t pos,
                               std::vector<char>& elements) {
  elements.clear();  // just in case
  if (cid[pos] == '*')
    return;
  bool inverted = false;
  if (cid[pos] == '!') {
    inverted = true;
    ++pos;
  }
  elements.resize((size_t)El::END, char(inverted));
  for (;;) {
    size_t sep = cid.find_first_of(",]", pos);
    if (sep == pos || sep > pos + 2)
      wrong_syntax(cid, 0, "in [...]");
    char elem_str[2] = {cid[pos], sep > pos+1 ? cid[pos+1] : '\0'};
    Element el = find_element(elem_str);
    if (el == El::X && (alpha_up(elem_str[0]) != 'X' || elem_str[1] != '\0'))
      wrong_syntax(cid, 0, " (invalid element in [...])");
    elements[el.ordinal()] = char(!inverted);
    pos = sep + 1;
    if (cid[sep] == ']')
      break;
  }
}

inline Selection::SequenceId parse_cid_seqid(const std::string& cid, size_t& pos,
                                             int default_seqnum) {
  size_t initial_pos = pos;
  int seqnum = default_seqnum;
  char icode = ' ';
  if (cid[pos] == '*') {
    ++pos;
    icode = '*';
  } else if (std::isdigit(cid[pos])) {
    char* endptr;
    seqnum = std::strtol(&cid[pos], &endptr, 10);
    pos = endptr - &cid[0];
  }
  if (cid[pos] == '.')
    ++pos;
  if (initial_pos != pos && (std::isalpha(cid[pos]) || cid[pos] == '*'))
    icode = cid[pos++];
  return {seqnum, icode};
}

inline Selection::AtomInequality parse_atom_inequality(const std::string& cid,
                                                       size_t pos, size_t end) {
  Selection::AtomInequality r;
  if (cid[pos] != 'q' && cid[pos] != 'b')
    wrong_syntax(cid, pos);
  r.property = cid[pos];
  ++pos;
  while (cid[pos] == ' ')
    ++pos;
  if (cid[pos] == '<')
    r.relation = -1;
  else if (cid[pos] == '>')
    r.relation = 1;
  else if (cid[pos] == '=')
    r.relation = 0;
  else
    wrong_syntax(cid, pos);
  ++pos;
  auto result = fast_from_chars(cid.c_str() + pos, r.value);
  if (result.ec != std::errc())
    wrong_syntax(cid, pos, " (expected number)");
  pos = size_t(result.ptr - cid.c_str());
  while (cid[pos] == ' ')
    ++pos;
  if (pos != end)
    wrong_syntax(cid, pos);
  return r;
}

inline bool has_inequality(const std::string& cid, size_t start, size_t end) {
  for (size_t i = start; i < end; ++i)
    if (cid[i] == '<' || cid[i] == '=' || cid[i] == '>')
      return true;
  return false;
}

inline void parse_cid(const std::string& cid, Selection& sel) {
  if (cid.empty() || (cid.size() == 1 && cid[0] == '*'))
    return;
  int omit = determine_omitted_cid_fields(cid);
  size_t sep = 0;
  size_t semi = cid.find(';');
  // model
  if (omit == 0) {
    sep = std::min(cid.find('/', 1), semi);
    if (sep != 1 && cid[1] != '*') {
      char* endptr;
      sel.mdl = std::strtol(&cid[1], &endptr, 10);
      size_t end_pos = endptr - &cid[0];
      if (end_pos != sep && end_pos != cid.size())
        wrong_syntax(cid, 0, " (at model number)");
    }
  }

  // chain
  if (omit <= 1 && sep < semi) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    sep = std::min(cid.find('/', pos), semi);
    sel.chain_ids = make_cid_list(cid, pos, sep);
  }

  // residue; MMDB CID syntax: s1.i1-s2.i2 or *(res).ic
  // In gemmi both 14.a and 14a are accepted.
  // *(ALA). and *(ALA) and (ALA). can be used instead of (ALA) for
  // compatibility with MMDB.
  if (omit <= 2 && sep < semi) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    if (cid[pos] != '(')
      sel.from_seqid = parse_cid_seqid(cid, pos, INT_MIN);
    if (cid[pos] == '(') {
      ++pos;
      size_t right_br = cid.find(')', pos);
      sel.residue_names = make_cid_list(cid, pos, right_br);
      pos = right_br + 1;
    }
    // allow "(RES)." and "(RES).*" and "(RES)*"
    if (cid[pos] == '.')
      ++pos;
    if (cid[pos] == '*')
      ++pos;
    if (cid[pos] == '-') {
      ++pos;
      sel.to_seqid = parse_cid_seqid(cid, pos, INT_MAX);
    }
    sep = pos;
    if (cid[sep] != '/' && cid[sep] != ';' && cid[sep] != '\0')
      wrong_syntax(cid, 0);
  }

  // atom;  at[el]:aloc
  if (sep < std::min(cid.size(), semi)) {
    size_t pos = (sep == 0 ? 0 : sep + 1);
    size_t end = cid.find_first_of("[:;", pos);
    if (end != pos) {
      sel.atom_names = make_cid_list(cid, pos, end);
      // Chain name can be empty, but not atom name,
      // so we interpret empty atom name as *.
      if (!sel.atom_names.inverted && sel.atom_names.list.empty())
        sel.atom_names.all = true;
      if (end == std::string::npos)
        return;
    }
    if (cid[end] == '[') {
      pos = end + 1;
      end = cid.find(']', pos);
      if (end == std::string::npos)
        wrong_syntax(cid, 0, " (no matching ']')");
      parse_cid_elements(cid, pos, sel.elements);
      ++end;
    }
    if (cid[end] == ':') {
      pos = end + 1;
      sel.altlocs = make_cid_list(cid, pos, semi);
    }
  }

  // extensions after semicolon(s)
  while (semi < cid.size()) {
    size_t pos = semi + 1;
    while (cid[pos] == ' ')
      ++pos;
    semi = std::min(cid.find(';', pos), cid.size());
    size_t end = semi;
    while (end > pos && cid[end-1] == ' ')
      --end;
    if (has_inequality(cid, pos, end)) {
      sel.atom_inequalities.push_back(parse_atom_inequality(cid, pos, end));
    } else {
      sel.entity_types = make_cid_list(cid, pos, end);
      bool inv = sel.entity_types.inverted;
      std::fill(sel.et_flags.begin(), sel.et_flags.end(), char(inv));
      for (const std::string& item : split_str(sel.entity_types.list, ',')) {
        EntityType et = EntityType::Unknown;
        if (item == "polymer")
          et = EntityType::Polymer;
        else if (item == "solvent")
          et = EntityType::Water;
        else
          wrong_syntax(cid, 0, (" at " + item).c_str());
        sel.et_flags[(int)et] = char(!inv);
      }
    }
  }
}

} // namespace impl


inline Selection::Selection(const std::string& cid) {
  impl::parse_cid(cid, *this);
}


template<class T> size_t count_atom_sites(const T& obj, const Selection* sel=nullptr) {
  size_t sum = 0;
  if (!sel || sel->matches(obj))
    for (const auto& child : obj.children())
      sum += count_atom_sites(child, sel);
  return sum;
}
template<> inline size_t count_atom_sites(const Atom& atom, const Selection* sel) {
  return (!sel || sel->matches(atom)) ? 1 : 0;
}

template<class T> double count_occupancies(const T& obj, const Selection* sel=nullptr) {
  double sum = 0;
  if (!sel || sel->matches(obj))
    for (const auto& child : obj.children())
        sum += count_occupancies(child, sel);
  return sum;
}
template<> inline double count_occupancies(const Atom& atom, const Selection* sel) {
  return (!sel || sel->matches(atom)) ? atom.occ : 0;
}

} // namespace gemmi
#endif
