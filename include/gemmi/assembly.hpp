// Copyright 2020 Global Phasing Ltd.
//
// Generating biological assemblies by applying operations
// from struct Assembly to a Model.
// Includes chain (re)naming utilities.

#ifndef GEMMI_ASSEMBLY_HPP_
#define GEMMI_ASSEMBLY_HPP_

#include <ostream>      // ostream
#include <memory>       // unique_ptr
#include "model.hpp"
#include "modify.hpp"   // transform_pos_and_adp
#include "neighbor.hpp" // merge_atoms_in_expanded_model
#include "util.hpp"

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

namespace impl {

inline bool any_subchain_matches(const Chain& chain, const Assembly::Gen& gen) {
  const std::string* prev_subchain = nullptr;
  for (const Residue& res : chain.residues)
    if (prev_subchain == nullptr || res.subchain != *prev_subchain) {
      if (in_vector(res.subchain, gen.subchains))
        return true;
      prev_subchain = &res.subchain;
    }
  return false;
}

struct AssemblyMapping {
  std::vector<std::string> sub;  // records subchain name correspondence
  std::vector<std::map<std::string, std::string>> chain_maps;
};

} // namespace impl

/// @par mapping is for internal use
inline Model make_assembly(const Assembly& assembly, const Model& model,
                           HowToNameCopiedChain how, std::ostream* out,
                           impl::AssemblyMapping* mapping=nullptr) {
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
      // chains are not merged here, multiple chains may have the same name
      std::map<std::string, std::string> new_names;
      bool all_chains = (!gen.chains.empty() && gen.chains[0] == "(all)");
      for (const Chain& chain : model.chains) {
        // PDB files specify bioassemblies in terms of chains,
        // mmCIF files in terms of subchains.
        bool whole_chain = (all_chains || in_vector(chain.name, gen.chains));
        if (whole_chain ||
            (!gen.subchains.empty() && impl::any_subchain_matches(chain, gen))) {
          // add a new empty chain, but first figure out the name for it
          auto result = new_names.emplace(chain.name, "");
          if (result.second)  // insertion happened - generate a new chain name
            result.first->second = namegen.make_new_name(chain.name, 1);
          new_model.chains.emplace_back(result.first->second);
          Chain& new_chain = new_model.chains.back();
          // add residues to the chain
          for (const Residue& res : chain.residues)
            if (whole_chain || in_vector(res.subchain, gen.subchains)) {
              new_chain.residues.push_back(res);
              Residue& new_res = new_chain.residues.back();
              transform_pos_and_adp(new_res, oper.transform);
              if (!new_res.subchain.empty()) {
                // change subchain name for the residue
                if (mapping && !mapping->sub.empty() &&
                    *(mapping->sub.end() - 2) == new_res.subchain) {
                  new_res.subchain = mapping->sub.back();
                  continue;
                }
                if (how == HowToNameCopiedChain::Short)
                  new_res.subchain = new_chain.name + ":" + new_res.subchain;
                else if (how == HowToNameCopiedChain::AddNumber)
                  new_res.subchain += new_chain.name.substr(chain.name.size());
                if (mapping) {
                  mapping->sub.push_back(res.subchain);
                  mapping->sub.push_back(new_res.subchain);
                }
              }
            }
        }
      }
      if (mapping)
        mapping->chain_maps.push_back(std::move(new_names));
    }
  return new_model;
}

inline Assembly expand_to_p1(const UnitCell& cell) {
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

inline void transform_to_assembly(Structure& st, const std::string& assembly_name,
                                  HowToNameCopiedChain how, std::ostream* out) {
  const Assembly* assembly = st.find_assembly(assembly_name);
  std::unique_ptr<Assembly> p1_assembly;
  if (!assembly) {
    if (assembly_name == "unit_cell") {
      p1_assembly.reset(new Assembly(expand_to_p1(st.cell)));
      assembly = p1_assembly.get();
    } else if (st.assemblies.empty()) {
      fail("no bioassemblies are listed for this structure");
    } else {
      fail("wrong assembly name, use one of: " +
           join_str(st.assemblies, ' ', [](const Assembly& a) { return a.name; }));
    }
  }
  impl::AssemblyMapping mapping;
  for (Model& model : st.models) {
    bool set_mapping = (&model == &st.models[0] && how != HowToNameCopiedChain::Dup);
    model = make_assembly(*assembly, model, how, out,
                          set_mapping ? &mapping : nullptr);
    merge_atoms_in_expanded_model(model, gemmi::UnitCell());
    assign_serial_numbers(model);
  }
  // update Entity::subchains
  if (!mapping.sub.empty())
    for (Entity& ent : st.entities) {
      std::vector<std::string> new_subchains;
      for (const std::string& s : ent.subchains)
        for (size_t i = 0; i < mapping.sub.size(); i += 2) {
          if (mapping.sub[i] == s)
            new_subchains.push_back(mapping.sub[i+1]);
        }
      ent.subchains = std::move(new_subchains);
    }

  // connections
  std::vector<Connection> new_connections;
  for (const Connection& conn : st.connections)
    if (conn.asu == Asu::Same) {
      int counter = 0;
      for (const std::map<std::string, std::string>& ch_map : mapping.chain_maps) {
        auto ch1 = ch_map.find(conn.partner1.chain_name);
        auto ch2 = ch_map.find(conn.partner2.chain_name);
        if (ch1 != ch_map.end() && ch2 != ch_map.end()) {
          Connection new_conn = conn;
          new_conn.partner1.chain_name = ch1->second;
          new_conn.partner2.chain_name = ch2->second;
          if (st.models[0].find_atom(new_conn.partner1) &&
              st.models[0].find_atom(new_conn.partner2)) {
            if (counter != 0)
              cat_to(new_conn.name, '.', counter);
            ++counter;
            new_connections.push_back(new_conn);
          }
        }
      }
    } else {
      // connections other than 1_555 are lost for now
      // if it's needed - get in touch
    }
  st.connections = std::move(new_connections);

  // Should Assembly instructions be kept or removed? Currently - removing.
  st.assemblies.clear();
  // Should st.spacegroup_hm and st.cell be kept? Here we remove only:
  st.cell.images.clear();
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
  ChainNameGenerator namegen(HowToNameCopiedChain::Short);
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


inline void expand_ncs(Structure& st, HowToNameCopiedChain how) {
  size_t orig_conn_size = st.connections.size();
  for (Model& model : st.models) {
    if (how == HowToNameCopiedChain::Dup) {
      // change segment of original chains to "0" - is this a good idea?
      for (Chain& chain : model.chains)
        for (Residue& res : chain.residues)
          res.segment = "0";
    }
    size_t orig_size = model.chains.size();
    ChainNameGenerator namegen(model, how);
    for (const NcsOp& op : st.ncs)
      if (!op.given) {
        std::map<std::string, std::string> chain_mapping;
        for (size_t i = 0; i != orig_size; ++i) {
          model.chains.push_back(model.chains[i]);
          Chain& new_chain = model.chains.back();
          const std::string& old_name = model.chains[i].name;
          auto it = chain_mapping.find(old_name);
          if (it == chain_mapping.end()) {
            new_chain.name = namegen.make_new_name(old_name, (int)i+1);
            chain_mapping.emplace(old_name, new_chain.name);
          } else {
            new_chain.name = it->second;
          }

          for (Residue& res : new_chain.residues) {
            transform_pos_and_adp(res, op.tr);
            if (!res.subchain.empty())
              res.subchain = new_chain.name + ":" + res.subchain;
            if (how == HowToNameCopiedChain::Dup)
              res.segment = op.id;
          }
        }
        // add connections when processing the first model
        if (&model == &st.models[0]) {
          for (size_t i = 0; i != orig_conn_size; ++i) {
            st.connections.push_back(st.connections[i]);
            Connection& c = st.connections.back();
            c.name += '-';
            c.name += op.id;
            for (int j = 0; j < 2; ++j) {
              AtomAddress& aa = j == 0 ? c.partner1 : c.partner2;
              if (how == HowToNameCopiedChain::Dup) {
                aa.res_id.segment = op.id;
              } else {
                auto it = chain_mapping.find(aa.chain_name);
                if (it != chain_mapping.end())
                  aa.chain_name = it->second;
                else
                  st.connections.pop_back();
              }
            }
          }
        }
      }
  }
  // adjust connections after changing segment of original chains to "0"
  if (how == HowToNameCopiedChain::Dup) {
    for (size_t i = 0; i != orig_conn_size; ++i) {
      st.connections[i].partner1.res_id.segment = "0";
      st.connections[i].partner2.res_id.segment = "0";
    }
  }
  for (NcsOp& op : st.ncs)
    op.given = true;
  st.setup_cell_images();
}

inline std::vector<Chain> split_chain_by_segments(Chain& orig, ChainNameGenerator& namegen) {
  std::vector<Chain> chains;
  std::vector<Residue> orig_res;
  orig_res.swap(orig.residues);
  int n = 0;
  for (auto start = orig_res.begin(); start != orig_res.end(); ) {
    const std::string& seg = start->segment;
    auto ch = std::find_if(chains.begin(), chains.end(), [&](Chain& c) {
                return !c.residues.empty() && c.residues[0].segment == seg; });
    if (ch == chains.end()) {
      chains.push_back(orig);
      ch = chains.end() - 1;
      // Naming. Here, Dup means chain name + segment.
      switch (namegen.how) {
        case HowToNameCopiedChain::Short:
          ch->name = namegen.make_short_name(ch->name + seg);
          break;
        case HowToNameCopiedChain::AddNumber:
          ch->name = namegen.make_name_with_numeric_postfix(ch->name, ++n);
          break;
        case HowToNameCopiedChain::Dup:
          ch->name += seg;
          break;
      }
    }
    auto end = std::find_if(start, orig_res.end(),
                            [&](Residue& r) { return r.segment != seg; });
    ch->residues.insert(ch->residues.end(), std::make_move_iterator(start),
                                            std::make_move_iterator(end));
    start = end;
  }
  return chains;
}

// HowToNameCopiedChain::Dup adds segment name to
inline void split_chains_by_segments(Model& model, HowToNameCopiedChain how) {
  ChainNameGenerator namegen(how);
  std::vector<Chain> new_chains;
  for (Chain& chain : model.chains)
    vector_move_extend(new_chains, split_chain_by_segments(chain, namegen));
  model.chains = std::move(new_chains);
}

} // namespace gemmi
#endif
