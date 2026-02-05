// Copyright 2018 Global Phasing Ltd.
//
// Selections.

#ifndef GEMMI_SELECT_HPP_
#define GEMMI_SELECT_HPP_

#include <climits>   // for INT_MIN, INT_MAX
#include "model.hpp" // for Model

namespace gemmi {

//! @brief Atom/residue selection using CID (coordinate ID) syntax.
//!
//! Based on CCP4 pdbcur selection syntax:
//! - /mdl/chn/s1.i1-s2.i2/at[el]:aloc  (sequence range)
//! - /mdl/chn/*(res).ic/at[el]:aloc   (residue name)
//!
//! Where: mdl=model, chn=chain, s=seqnum, i=icode, res=residue name,
//!        at=atom name, el=element, aloc=altloc
//!
//! Example: "//A/10-20/CA" selects CA atoms in residues 10-20 of chain A.
struct GEMMI_DLL Selection {
  //! @brief Comma-separated list with wildcard and inversion support.
  //!
  //! Used for chain IDs, residue names, atom names, etc.
  //! Supports '*' for all, '!' prefix for inversion.
  struct List {
    bool all = true;  //!< True if wildcard '*' (match all)
    bool inverted = false;  //!< True if list is negated with '!'
    std::string list;  //!< Comma-separated values

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

  //! @brief Single-character flag pattern for filtering.
  //!
  //! Used for residue_flags and atom_flags.
  struct FlagList {
    std::string pattern;  //!< Pattern string (e.g., "PS" or "!W")

    //! @brief Check if a flag matches the pattern.
    //! @param flag Character flag to test
    //! @return True if flag matches (or pattern is empty)
    bool has(char flag) const {
      if (pattern.empty())
        return true;
      bool invert = (pattern[0] == '!');
      bool found = (pattern.find(flag, invert ? 1 : 0) != std::string::npos);
      return invert ? !found : found;
    }
  };

  //! @brief Sequence identifier for residue range selection.
  //!
  //! Combines sequence number and insertion code.
  struct SequenceId {
    int seqnum;  //!< Sequence number
    char icode;  //!< Insertion code

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

  //! @brief Inequality constraint on atom properties.
  //!
  //! Used for filtering atoms by occupancy (q) or B-factor (b).
  struct AtomInequality {
    char property;  //!< Property: 'q' (occupancy) or 'b' (B-factor)
    int relation;  //!< Comparison: <0 (less), >0 (greater), 0 (equal)
    double value;  //!< Threshold value

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

  int mdl = 0;  //!< Model number (0 = all models)
  List chain_ids;  //!< Chain ID filter
  SequenceId from_seqid = {INT_MIN, '*'};  //!< Start of sequence range
  SequenceId to_seqid = {INT_MAX, '*'};  //!< End of sequence range
  List residue_names;  //!< Residue name filter
  List entity_types;  //!< Entity type filter
  std::array<char, 6> et_flags;  //!< Array indexed by EntityType enum
  List atom_names;  //!< Atom name filter
  std::vector<char> elements;  //!< Element filter (indexed by atomic number)
  List altlocs;  //!< Alternate location filter
  FlagList residue_flags;  //!< Residue flag filter
  FlagList atom_flags;  //!< Atom flag filter
  std::vector<AtomInequality> atom_inequalities;  //!< Property constraints

  Selection() = default;

  //! @brief Create selection from CID string.
  //! @param cid Coordinate ID selection string
  Selection(const std::string& cid);

  std::string str() const;

  //! @brief Check if structure matches (always true).
  bool matches(const Structure&) const { return true; }

  //! @brief Check if model matches selection.
  //! @param model Model to test
  //! @return True if model number matches
  bool matches(const Model& model) const {
    return mdl == 0 || mdl == model.num;
  }

  //! @brief Check if chain matches selection.
  //! @param chain Chain to test
  //! @return True if chain ID matches
  bool matches(const Chain& chain) const {
    return chain_ids.has(chain.name);
  }

  //! @brief Check if residue matches selection.
  //! @param res Residue to test
  //! @return True if all residue criteria match
  bool matches(const Residue& res) const {
    return (entity_types.all || et_flags[(int)res.entity_type]) &&
           residue_names.has(res.name) &&
           from_seqid.compare(res.seqid) <= 0 &&
           to_seqid.compare(res.seqid) >= 0 &&
           residue_flags.has(res.flag);
  }
  //! @brief Check if atom matches selection.
  //! @param a Atom to test
  //! @return True if all atom criteria match
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
