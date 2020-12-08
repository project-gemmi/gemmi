// Copyright 2020 Global Phasing Ltd.
//
// A function that generates biological assemblies by applying operations
// from struct Assembly to a Model.
// And chain (re)naming utilities.

#ifndef GEMMI_ASSEMBLY_HPP_
#define GEMMI_ASSEMBLY_HPP_

#include <ostream>      // std::ostream
#include "model.hpp"
#include "util.hpp"

namespace gemmi {

enum class HowToNameCopiedChains { Short, AddNumber, Dup };

struct ChainNameGenerator {
  using How = HowToNameCopiedChains;
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
    static const char symbols[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                  "abcdefghijklmnopqrstuvwxyz0123456789";
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

inline Model make_assembly(const Assembly& assembly, const Model& model,
                           HowToNameCopiedChains how, std::ostream* out) {
  Model new_model(model.name);
  ChainNameGenerator namegen(how);
  std::map<std::string, std::string> subs = model.subchain_to_chain();
  for (const Assembly::Gen& gen : assembly.generators)
    for (const Assembly::Operator& oper : gen.operators) {
      if (out) {
        *out << "Applying " << oper.name << " to";
        if (!gen.chains.empty())
          *out << " chains: " << join_str(gen.chains, ',');
        else if (!gen.subchains.empty())
          *out << " subchains: " << join_str(gen.subchains, ',');
        *out << std::endl;
        for (const std::string& chain_name : gen.chains)
          if (!model.find_chain(chain_name))
            *out << "Warning: no chain " << chain_name << std::endl;
        for (const std::string& subchain_name : gen.subchains)
          if (subs.find(subchain_name) == subs.end())
            *out << "Warning: no subchain " << subchain_name << std::endl;
      }
      // PDB files specify bioassemblies in terms of chains,
      // mmCIF files in terms of subchains. We handle the two cases separately.
      if (!gen.chains.empty()) {
        // chains are not merged here, multiple chains may have the same name
        std::map<std::string, std::string> new_names;
        for (size_t i = 0; i != model.chains.size(); ++i) {
          if (in_vector(model.chains[i].name, gen.chains)) {
            new_model.chains.push_back(model.chains[i]);
            Chain& new_chain = new_model.chains.back();
            auto name_iter = new_names.find(model.chains[i].name);
            if (name_iter == new_names.end()) {
              new_chain.name = namegen.make_new_name(new_chain.name, 1);
              new_names.emplace(model.chains[i].name, new_chain.name);
            } else {
              new_chain.name = name_iter->second;
            }
            for (Residue& res : new_chain.residues) {
              for (Atom& a : res.atoms)
                a.pos = Position(oper.transform.apply(a.pos));
              if (!res.subchain.empty())
                res.subchain = new_chain.name + ":" + res.subchain;
            }
          }
        }
      } else if (!gen.subchains.empty()) {
        std::map<std::string, std::string> new_names;
        for (const std::string& subchain_name : gen.subchains) {
          auto sub_iter = subs.find(subchain_name);
          if (sub_iter == subs.end())
            continue;
          auto name_iter = new_names.find(sub_iter->second);
          Chain* new_chain;
          if (name_iter == new_names.end()) {
            std::string new_name = namegen.make_new_name(sub_iter->second, 1);
            new_names.emplace(sub_iter->second, new_name);
            new_model.chains.emplace_back(new_name);
            new_chain = &new_model.chains.back();
          } else {
            new_chain = new_model.find_chain(name_iter->second);
          }
          for (const Residue& res : model.get_subchain(subchain_name)) {
            new_chain->residues.push_back(res);
            Residue& new_res = new_chain->residues.back();
            new_res.subchain = new_chain->name + ":" + res.subchain;
            for (Atom& a : new_res.atoms)
              a.pos = Position(oper.transform.apply(a.pos));
          }
        }
      }
    }
  return new_model;
}

inline void change_to_assembly(Structure& st, const std::string& assembly_name,
                               HowToNameCopiedChains how, std::ostream* out) {
  Assembly* assembly = st.find_assembly(assembly_name);
  if (!assembly) {
    if (st.assemblies.empty())
      fail("no bioassemblies are listed for this structure");
    fail("wrong assembly name, use one of: " +
        join_str(st.assemblies, ' ', [](const Assembly& a) { return a.name; }));
  }
  for (Model& model : st.models)
    model = make_assembly(*assembly, model, how, out);
  st.connections.clear();
}

// chain is assumed to be from st.models[0]
inline void rename_chain(Structure& st, Chain& chain,
                         const std::string& new_name) {
  auto rename_if_matches = [&](AtomAddress& aa) {
    if (aa.chain_name == chain.name)
      aa.chain_name = new_name;
  };
  for (Connection& con : st.connections) {
    rename_if_matches(con.partner1);
    rename_if_matches(con.partner2);
  }
  for (Helix& helix : st.helices) {
    rename_if_matches(helix.start);
    rename_if_matches(helix.end);
  }
  for (Sheet& sheet : st.sheets)
    for (Sheet::Strand& strand : sheet.strands) {
      rename_if_matches(strand.start);
      rename_if_matches(strand.end);
      rename_if_matches(strand.hbond_atom2);
      rename_if_matches(strand.hbond_atom1);
    }
  for (RefinementInfo& ri : st.meta.refinement)
    for (TlsGroup& tls : ri.tls_groups)
      for (TlsGroup::Selection& sel : tls.selections)
        if (sel.chain == chain.name)
          sel.chain = new_name;
  for (auto it = st.models.begin() + 1; it != st.models.end(); ++it)
    if (Chain* ch = it->find_chain(chain.name))
      ch->name = new_name;
  chain.name = new_name;
}

inline void shorten_chain_names(Structure& st) {
  ChainNameGenerator namegen(HowToNameCopiedChains::Short);
  Model& model0 = st.models[0];
  size_t max_len = model0.chains.size() < 63 ? 1 : 2;
  for (const Chain& chain : model0.chains)
    if (chain.name.length() <= max_len)
      namegen.used_names.push_back(chain.name);
  for (Chain& chain : model0.chains)
    if (chain.name.length() > max_len)
      rename_chain(st, chain,
                   namegen.make_short_name(chain.name.substr(0, max_len)));
}


inline void expand_ncs(gemmi::Structure& st, HowToNameCopiedChains how) {
  for (gemmi::Model& model : st.models) {
    size_t orig_size = model.chains.size();
    gemmi::ChainNameGenerator namegen(model, how);
    for (const gemmi::NcsOp& op : st.ncs)
      if (!op.given) {
        for (size_t i = 0; i != orig_size; ++i) {
          if (how == HowToNameCopiedChains::Dup)
            for (gemmi::Residue& res : model.chains[i].residues)
              res.segment = "0";

          model.chains.push_back(model.chains[i]);
          gemmi::Chain& new_chain = model.chains.back();
          new_chain.name = namegen.make_new_name(new_chain.name, (int)i+1);

          for (gemmi::Residue& res : new_chain.residues) {
            for (gemmi::Atom& a : res.atoms) {
              a.pos = op.apply(a.pos);
              if (a.aniso.nonzero()) {
                SMat33<double> u = a.aniso.transformed_by(op.tr.mat);
                a.aniso = {static_cast<float>(u.u11), static_cast<float>(u.u22), static_cast<float>(u.u33),
                           static_cast<float>(u.u12), static_cast<float>(u.u13), static_cast<float>(u.u23)};
              }
            }
            if (!res.subchain.empty())
              res.subchain = new_chain.name + ":" + res.subchain;
            if (how == HowToNameCopiedChains::Dup)
              res.segment = op.id;
          }
        }
      }
  }
  for (gemmi::NcsOp& op : st.ncs)
    op.given = true;
  st.setup_cell_images();
}


} // namespace gemmi
#endif
