// Copyright 2020 Global Phasing Ltd.

#include "gemmi/assembly.hpp"

#include <memory>             // unique_ptr
#include "gemmi/modify.hpp"   // transform_pos_and_adp
#include "gemmi/neighbor.hpp" // NeighborSearch

namespace gemmi {

namespace {

bool any_subchain_matches(const Chain& chain, const Assembly::Gen& gen) {
  const std::string* prev_subchain = nullptr;
  for (const Residue& res : chain.residues)
    if (prev_subchain == nullptr || res.subchain != *prev_subchain) {
      if (in_vector(res.subchain, gen.subchains))
        return true;
      prev_subchain = &res.subchain;
    }
  return false;
}

// chain name mappings old->new for each Assembly::Operator used
struct ChainMap {
  bool uses_segments = false;
  std::string id;
  std::map<std::string, std::string> names;
};
struct AssemblyMapping {
  HowToNameCopiedChain how;
  // records subchain name correspondence (new->old), used for updating Entity::subchain
  std::map<std::string, std::string> sub;
  std::vector<ChainMap> chain_maps;
};

void update_address(AtomAddress& a, const ChainMap& chain_map) {
  auto it = chain_map.names.find(a.chain_name);
  if (it != chain_map.names.end())
    a.chain_name = it->second;
}

void remove_cras(Model& model, std::vector<CRA>& vec) {
  // sort in reverse order, so items can be erased without invalidating pointers
  std::sort(vec.begin(), vec.end(), [](const CRA& a, const CRA& b) {
      return std::tie(a.chain, a.residue, a.atom) > std::tie(b.chain, b.residue, b.atom);
  });
  const Atom* prev_a = nullptr;
  for (CRA& cra : vec) {
    if (cra.atom == prev_a)
      continue;
    prev_a = cra.atom;
    auto atom_idx = cra.atom - cra.residue->atoms.data();
    cra.residue->atoms.erase(cra.residue->atoms.begin() + atom_idx);
    if (cra.residue->atoms.empty()) {
      auto res_idx = cra.residue - cra.chain->residues.data();
      cra.chain->residues.erase(cra.chain->residues.begin() + res_idx);
      if (cra.chain->residues.empty()) {
        auto chain_idx = cra.chain - model.chains.data();
        model.chains.erase(model.chains.begin() + chain_idx);
      }
    }
  }
}

Model make_assembly_(const Assembly& assembly, const Model& model,
                     HowToNameCopiedChain how, std::ostream* out,
                     AssemblyMapping* mapping) {
  Model new_model(model.name);
  ChainNameGenerator namegen(how);
  std::map<std::string, std::string> subs = model.subchain_to_chain();
  int counter = 0;
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
      ChainMap chain_map;
      if (counter != 0) {
        chain_map.uses_segments = (how == HowToNameCopiedChain::Dup);
        chain_map.id = std::to_string(counter);
      }
      bool all_chains = (!gen.chains.empty() && gen.chains[0] == "(all)");
      for (const Chain& chain : model.chains) {
        // PDB files specify bioassemblies in terms of chains,
        // mmCIF files in terms of subchains.
        bool whole_chain = (all_chains || in_vector(chain.name, gen.chains));
        if (whole_chain ||
            (!gen.subchains.empty() && any_subchain_matches(chain, gen))) {
          // try add a new empty chain, but first figure out the name for it
          auto result = chain_map.names.emplace(chain.name, "");
          if (result.second)  // insertion happened - generate a new chain name
            result.first->second = namegen.make_new_name(chain.name, counter+1);
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
                if (how == HowToNameCopiedChain::Short)
                  new_res.subchain = new_chain.name + ":" + new_res.subchain;
                else if (how == HowToNameCopiedChain::AddNumber)
                  new_res.subchain += new_chain.name.substr(chain.name.size());
                if (mapping)
                  mapping->sub.emplace(new_res.subchain, res.subchain);
              }
              if (chain_map.uses_segments)
                new_res.segment = chain_map.id;
            }
        }
      }
      if (mapping)
        mapping->chain_maps.push_back(std::move(chain_map));
      ++counter;
    }
  return new_model;
}


void expand_ncs_model_(Model& model, const std::vector<NcsOp>& ncs,
                       HowToNameCopiedChain how, AssemblyMapping* mapping) {
  // For HowToNameCopiedChain::Dup we used to set segment="0" for original
  // residues. Now it's left blank. Not sure which is better.
  size_t orig_size = model.chains.size();
  ChainNameGenerator namegen(model, how);
  if (mapping) {  // add identities to AssemblyMapping
    mapping->chain_maps.emplace_back();
    ChainMap& chain_map = mapping->chain_maps.back();
    // chain_map.id is left empty
    for (const Chain& chain : model.chains) {
      chain_map.names.emplace(chain.name, chain.name);
      for (const ConstResidueSpan& span : chain.subchains()) {
        const std::string& sub_id = span.subchain_id();
        mapping->sub.emplace(sub_id, sub_id);
      }
    }
  }
  int counter = 0;
  for (const NcsOp& op : ncs)
    if (!op.given) {
      ChainMap chain_map;
      chain_map.uses_segments = (how == HowToNameCopiedChain::Dup);
      chain_map.id = op.id;
      ++counter;
      for (size_t i = 0; i != orig_size; ++i) {
        model.chains.push_back(model.chains[i]);
        Chain& new_chain = model.chains.back();
        const std::string& old_name = model.chains[i].name;
        auto result = chain_map.names.emplace(old_name, "");
        if (how != HowToNameCopiedChain::Dup) {
          if (result.second) // if insertion happened - generate a new chain name
            result.first->second = namegen.make_new_name(old_name, counter);
          new_chain.name = result.first->second;
        }
        for (Residue& new_res : new_chain.residues) {
          transform_pos_and_adp(new_res, op.tr);
          if (!new_res.subchain.empty()) {
            std::string old_subchain = new_res.subchain;
            new_res.subchain = new_chain.name + ":" + new_res.subchain;
            if (mapping)
              mapping->sub.emplace(new_res.subchain, old_subchain);
          }
          if (chain_map.uses_segments)
            new_res.segment = chain_map.id;
        }
      }
      if (mapping)
        mapping->chain_maps.push_back(std::move(chain_map));
    }
}


void finalize_expansion(Structure& st, const AssemblyMapping& mapping,
                        double merge_dist, bool expanding_ncs) {
  if (merge_dist > 0)
    for (Model& model : st.models) {
      merge_atoms_in_expanded_model(model, gemmi::UnitCell(), merge_dist);
      assign_serial_numbers(model);
    }

  // update Entity::subchains
  if (!mapping.sub.empty())
    for (Entity& ent : st.entities) {
      std::vector<std::string> new_subchains;
      for (const std::string& s : ent.subchains)
        for (const auto& new_old : mapping.sub)
          if (new_old.second == s)
            new_subchains.push_back(new_old.first);
      ent.subchains = std::move(new_subchains);
    }

  // connections
  std::vector<Connection> new_connections;
  for (const Connection& conn : st.connections)
    if (expanding_ncs || conn.asu == Asu::Same) {
      bool first = true;
      for (const ChainMap& chain_map : mapping.chain_maps) {
        auto ch1 = chain_map.names.find(conn.partner1.chain_name);
        auto ch2 = chain_map.names.find(conn.partner2.chain_name);
        if (ch1 != chain_map.names.end() &&
            ch2 != chain_map.names.end()) {
          Connection new_conn = conn;
          new_conn.partner1.chain_name = ch1->second;
          new_conn.partner2.chain_name = ch2->second;
          if (chain_map.uses_segments)
            new_conn.partner1.res_id.segment =
            new_conn.partner2.res_id.segment = chain_map.id;
          if (st.models[0].find_atom(new_conn.partner1) &&
              st.models[0].find_atom(new_conn.partner2)) {
            if (!first) {
              cat_to(new_conn.name, '.', chain_map.id);
              first = false;
            }
            new_connections.push_back(new_conn);
          }
        }
      }
    } else {
      // connections other than 1_555 are lost when making assembly
      // if it's needed - get in touch
    }
  st.connections = std::move(new_connections);

  if (mapping.how == HowToNameCopiedChain::Dup)
    return;

  // cispeps
  std::vector<CisPep> new_cispeps;
  new_cispeps.reserve(st.cispeps.size() * mapping.chain_maps.size());
  for (const CisPep& cispep : st.cispeps)
    for (const ChainMap& chain_map : mapping.chain_maps) {
      new_cispeps.push_back(cispep);
      update_address(new_cispeps.back().partner_c, chain_map);
      update_address(new_cispeps.back().partner_n, chain_map);
    }
  st.cispeps = std::move(new_cispeps);

  // secondary structure - helices
  std::vector<Helix> new_helices;
  new_helices.reserve(st.helices.size() * mapping.chain_maps.size());
  for (const Helix& helix : st.helices)
    for (const ChainMap& chain_map : mapping.chain_maps) {
      new_helices.push_back(helix);
      update_address(new_helices.back().start, chain_map);
      update_address(new_helices.back().end, chain_map);
    }
  st.helices = std::move(new_helices);

  // secondary structure - sheets
  std::vector<Sheet> new_sheets;
  new_sheets.reserve(st.sheets.size() * mapping.chain_maps.size());
  for (const Sheet& sheet : st.sheets)
    for (const ChainMap& chain_map : mapping.chain_maps) {
      new_sheets.push_back(sheet);
      for (Sheet::Strand& strand : new_sheets.back().strands) {
        update_address(strand.start, chain_map);
        update_address(strand.end, chain_map);
        update_address(strand.hbond_atom2, chain_map);
        update_address(strand.hbond_atom1, chain_map);
      }
    }
  st.sheets = new_sheets;
}

} // anonymous namespace

Model make_assembly(const Assembly& assembly, const Model& model,
                    HowToNameCopiedChain how, std::ostream* out) {
  return make_assembly_(assembly, model, how, out, nullptr);
}

void transform_to_assembly(Structure& st, const std::string& assembly_name,
                           HowToNameCopiedChain how, std::ostream* out,
                           bool keep_spacegroup, double merge_dist) {
  const Assembly* assembly = st.find_assembly(assembly_name);
  std::unique_ptr<Assembly> p1_assembly;
  if (!assembly) {
    if (assembly_name == "unit_cell") {
      p1_assembly.reset(new Assembly(pseudo_assembly_for_unit_cell(st.cell)));
      assembly = p1_assembly.get();
    } else if (st.assemblies.empty()) {
      fail("no bioassemblies are listed for this structure");
    } else {
      fail("wrong assembly name, use one of: " +
           join_str(st.assemblies, ' ', [](const Assembly& a) { return a.name; }));
    }
  }

  AssemblyMapping mapping;
  mapping.how = how;
  AssemblyMapping* mapping_ptr = &mapping;
  for (Model& model : st.models) {
    model = make_assembly_(*assembly, model, how, out, mapping_ptr);
    mapping_ptr = nullptr;  // AssemblyMapping is based only on the first model
  }
  finalize_expansion(st, mapping, merge_dist, false);

  // Should Assembly instructions be kept or removed? Currently - removing.
  st.assemblies.clear();

  if (!keep_spacegroup) {
    st.spacegroup_hm = "P 1";  // maybe it should be empty?
    if (assembly_name != "unit_cell")
      st.cell = UnitCell();
  } else {
    st.cell.images.clear();
  }
}

Model expand_ncs_model(const Model& model, const std::vector<NcsOp>& ncs,
                       HowToNameCopiedChain how) {
  Model model_copy = model;
  expand_ncs_model_(model_copy, ncs, how, nullptr);
  return model_copy;
}

void expand_ncs(Structure& st, HowToNameCopiedChain how, double merge_dist) {
  AssemblyMapping mapping;
  mapping.how = how;
  AssemblyMapping* mapping_ptr = &mapping;
  for (Model& model : st.models) {
    expand_ncs_model_(model, st.ncs, how, mapping_ptr);
    mapping_ptr = nullptr;  // AssemblyMapping is based only on the first model
  }
  finalize_expansion(st, mapping, merge_dist, true);

  for (NcsOp& op : st.ncs)
    op.given = true;

  st.setup_cell_images();
}


void merge_atoms_in_expanded_model(Model& model, const UnitCell& cell, double max_dist,
                                   bool compare_serial) {
  using Mark = NeighborSearch::Mark;
  NeighborSearch ns(model, cell, 4.0);
  ns.populate(true);
  std::vector<CRA> to_be_deleted;
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    Chain& chain = model.chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      Residue& res = chain.residues[n_res];
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        Atom& atom = res.atoms[n_atom];
        std::vector<std::pair<CRA, int>> equiv;
        ns.for_each_cell(atom.pos, [&](std::vector<Mark>& marks, const Fractional& fr) {
            for (Mark& m : marks) {
              // We look for the same atoms, but copied to a different chain.
              // First quick check that filters out most of non-matching pairs.
              if (m.altloc != atom.altloc || m.element != atom.element ||
                  m.chain_idx == n_ch || m.atom_idx != n_atom)
                continue;
              // Now check if everything else matches.
              CRA cra = m.to_cra(model);
              if (cra.atom &&
                  (!compare_serial || cra.atom->serial == atom.serial) &&
                  cra.atom->name == atom.name &&
                  cra.atom->b_iso == atom.b_iso &&
                  cra.residue->matches_noseg(res) &&
                  m.pos.dist_sq(ns.grid.unit_cell.orthogonalize(fr)) < sq(max_dist))
                equiv.emplace_back(cra, m.image_idx);
            }
        });
        if (!equiv.empty()) {
          Position pos_sum = atom.pos;
          for (auto& t : equiv) {
            CRA& cra = t.first;
            pos_sum += ns.grid.unit_cell.find_nearest_pbc_position(
                                          atom.pos, cra.atom->pos, t.second);
            // The atoms in equiv are to be discarded later.
            // Deleting now would invalidate indices in NeighborSearch.
            to_be_deleted.push_back(cra);
            // Modify the atoms to avoid processing them again.
            cra.atom->serial = -1;
            cra.atom->name.clear();
          }
          size_t n = 1 + equiv.size();
          atom.pos = pos_sum / double(n);
          atom.occ = std::min(1.f, n * atom.occ);
        }
      }
    }
  }
  remove_cras(model, to_be_deleted);
}


void shorten_chain_names(Structure& st) {
  ChainNameGenerator namegen(HowToNameCopiedChain::Short);
  Model& model0 = st.models[0];
  size_t max_len = model0.chains.size() < 63 ? 1 : 2;
  for (const Chain& chain : model0.chains)
    if (chain.name.length() <= max_len)
      namegen.used_names.push_back(chain.name);
  for (Chain& chain : model0.chains)
    if (chain.name.length() > max_len)
      rename_chain(st, chain.name,
                   namegen.make_short_name(chain.name.substr(0, max_len)));
}


static std::vector<Chain> split_chain_by_segments(Chain& orig, ChainNameGenerator& namegen) {
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

void split_chains_by_segments(Model& model, HowToNameCopiedChain how) {
  ChainNameGenerator namegen(how);
  std::vector<Chain> new_chains;
  for (Chain& chain : model.chains)
    vector_move_extend(new_chains, split_chain_by_segments(chain, namegen));
  model.chains = std::move(new_chains);
}

} // namespace gemmi
