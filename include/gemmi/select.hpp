// Copyright 2018 Global Phasing Ltd.
//
// Selections.

#ifndef GEMMI_SELECT_HPP_
#define GEMMI_SELECT_HPP_

#include <climits>   // for INT_MIN, INT_MAX
#include "model.hpp" // for Model

namespace gemmi {

// from http://www.ccp4.ac.uk/html/pdbcur.html
// Specification of the selection sets:
// either
//     /mdl/chn/s1.i1-s2.i2/at[el]:aloc
// or
//     /mdl/chn/*(res).ic/at[el]:aloc
//

struct GEMMI_DLL Selection {
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

    std::string str() const;

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

    std::string str() const;
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

  std::string str() const;

  bool matches(const Structure&) const { return true; }
  bool matches(const Model& model) const {
    return mdl == 0 || mdl == model.num;
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

} // namespace gemmi
#endif
