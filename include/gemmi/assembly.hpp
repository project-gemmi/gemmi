// Copyright 2020 Global Phasing Ltd.
//
// Generating biological assemblies by applying operations
// from struct Assembly to a Model.
// Includes chain (re)naming utilities.

#ifndef GEMMI_ASSEMBLY_HPP_
#define GEMMI_ASSEMBLY_HPP_

#include <ostream>    // for ostream
#include "model.hpp"  // for Model
#include "util.hpp"   // for in_vector

namespace gemmi {

enum class HowToNameCopiedChain { Short, AddNumber, Dup };

struct ChainNameGenerator {
  using How = HowToNameCopiedChain;
  How how;
  std::vector<std::string> used_names;

  ChainNameGenerator(How how_) : how(how_) {}
  ChainNameGenerator(const Model& model, How how_) : how(how_) {
    if (how != How::Dup)
      for (const Chain& chain : model.chains)
        used_names.push_back(chain.name);
  }
  bool has(const std::string& name) const {
    return in_vector(name, used_names);
  }
  const std::string& added(const std::string& name) {
    used_names.push_back(name);
    return name;
  }

  std::string make_short_name(const std::string& preferred) {
    static const char symbols[] = {
      'A','B','C','D','E','F','G','H','I','J','K','L','M',
      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
      'a','b','c','d','e','f','g','h','i','j','k','l','m',
      'n','o','p','q','r','s','t','u','v','w','x','y','z',
      '0','1','2','3','4','5','6','7','8','9'
    };
    if (!has(preferred))
      return added(preferred);
    std::string name(1, 'A');
    for (char symbol : symbols) {
      name[0] = symbol;
      if (!has(name))
        return added(name);
    }
    name += 'A';
    for (char symbol1 : symbols) {
      name[0] = symbol1;
      for (char symbol2 : symbols) {
        name[1] = symbol2;
        if (!has(name))
          return added(name);
      }
    }
    fail("run out of 1- and 2-letter chain names");
  }

  std::string make_name_with_numeric_postfix(const std::string& base, int n) {
    std::string name = base;
    name += std::to_string(n);
    while (has(name)) {
      name.resize(base.size());
      name += std::to_string(++n);
    }
    return added(name);
  }

  std::string make_new_name(const std::string& old, int n) {
    switch (how) {
      case How::Short: return make_short_name(old);
      case How::AddNumber: return make_name_with_numeric_postfix(old, n);
      case How::Dup: return old;
    }
    unreachable();
  }
};

inline void ensure_unique_chain_name(const Model& model, Chain& chain) {
  ChainNameGenerator namegen(HowToNameCopiedChain::Short);
  for (const Chain& ch : model.chains)
    if (&ch != &chain && !namegen.has(ch.name))
      namegen.added(ch.name);
  chain.name = namegen.make_short_name(chain.name);
}

GEMMI_DLL Model make_assembly(const Assembly& assembly, const Model& model,
                              HowToNameCopiedChain how, std::ostream* out);

inline Assembly pseudo_assembly_for_unit_cell(const UnitCell& cell) {
  Assembly assembly("unit_cell");
  std::vector<Assembly::Operator> operators(cell.images.size() + 1);
  // operators[0] stays as identity
  for (size_t i = 1; i != operators.size(); ++i) {
    const FTransform& op = cell.images[i-1];
    operators[i].transform = cell.orth.combine(op.combine(cell.frac));
  }
  assembly.generators.push_back({{"(all)"}, {}, operators});
  return assembly;
}

/// If called with assembly_name="unit_cell" changes structure to unit cell (P1).
/// \par keep_spacegroup preserves space group and unit cell - is it needed?
GEMMI_DLL void transform_to_assembly(Structure& st, const std::string& assembly_name,
                                     HowToNameCopiedChain how, std::ostream* out,
                                     bool keep_spacegroup=false, double merge_dist=0.2);


GEMMI_DLL Model expand_ncs_model(const Model& model, const std::vector<NcsOp>& ncs,
                                 HowToNameCopiedChain how);

/// Searches and merges overlapping equivalent atoms from different chains.
/// To be used after expand_ncs() and make_assembly().
GEMMI_DLL void merge_atoms_in_expanded_model(Model& model, const UnitCell& cell,
                                             double max_dist=0.2, bool compare_serial=true);


GEMMI_DLL void shorten_chain_names(Structure& st);

GEMMI_DLL void expand_ncs(Structure& st, HowToNameCopiedChain how, double merge_dist=0.2);

/// HowToNameCopiedChain::Dup adds segment name to chain name
GEMMI_DLL void split_chains_by_segments(Model& model, HowToNameCopiedChain how);

} // namespace gemmi
#endif
