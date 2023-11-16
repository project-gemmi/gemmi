// Copyright 2017-2022 Global Phasing Ltd.

#include <gemmi/polyheur.hpp>

namespace gemmi {

PolymerType check_polymer_type(const ConstResidueSpan& span, bool ignore_entity_type) {
  if (span.empty())
    return PolymerType::Unknown;
  size_t counts[(int)ResidueKind::ELS+1] = {0};
  size_t aa = 0;
  size_t na = 0;
  size_t total = 0;
  bool has_atom_record = false;
  for (const Residue& r : span)
    if (ignore_entity_type ||
        r.entity_type == EntityType::Unknown ||
        r.entity_type == EntityType::Polymer) {
      if (r.het_flag == 'A')
        has_atom_record = true;
      ResidueInfo info = find_tabulated_residue(r.name);
      if (info.found()) {
        // Exclude water and ions - it can make difference
        // if this function is called for the whole chain.
        // Components w/o hydrogens are often ions and always non-polymers
        // (and almost never are in a polymer - except PO4, PO2 and AZI (N3)
        // which in a few PDB entries are included in polymers - but it
        // doesn't matter here).
        if (info.kind == ResidueKind::HOH || info.hydrogen_count == 0)
          continue;
        if (info.is_peptide_linking())
          ++aa;
        if (info.is_na_linking())
          ++na;
        counts[(int)info.kind]++;
      } else if (r.get_ca()) {
        ++aa;
      } else if (r.get_p()) {
        ++na;
      }
      ++total;
    }
  if (total == 0)
    return PolymerType::Unknown;
  // One residue is not a polymer, but it may happen that only a single residue
  // of a chain is modelled. OTOH, a single non-standard residue is usually
  // a ligand.
  if (total == 1 && !has_atom_record)
    return PolymerType::Unknown;
  // ATOM records suggest a polymer, so weaken the condition for AA/NA polymers.
  size_t bonus = has_atom_record ? 1 : 0;
  if (2 * aa + bonus > total)
    return counts[(int)ResidueKind::AA] >= counts[(int)ResidueKind::AAD]
           ? PolymerType::PeptideL : PolymerType::PeptideD;
  if (2 * na + bonus > total) {
    if (counts[(int)ResidueKind::DNA] == 0)
      return PolymerType::Rna;
    if (counts[(int)ResidueKind::RNA] == 0)
      return PolymerType::Dna;
    return PolymerType::DnaRnaHybrid;
  }
  return PolymerType::Unknown;
}

static std::vector<Residue>::iterator infer_polymer_end(Chain& chain, PolymerType ptype) {
  auto b = chain.residues.begin();
  auto e = chain.residues.end();
  for (auto it = b; it != e; ++it) {
    ResidueInfo info = find_tabulated_residue(it->name);
    if (info.found()) {
      if (info.is_water()) {
        e = it;
        break;
      }
      bool maybe_linking = (is_polypeptide(ptype) && info.is_peptide_linking())
                        || (is_polynucleotide(ptype) && info.is_na_linking());
      // The first residue could be non-polymer.
      if (!maybe_linking && b != chain.residues.begin()) {
        e = it;
        break;
      }
      // If a standard residue is HETATM we assume that it is in the buffer.
      if (info.is_standard() && it->het_flag == 'H') {
        e = it;
        break;
      }
      if (maybe_linking && info.is_standard())
        b = it;
    }
  }
  if (b == e || b + 1 == e)
    return e;
  // Ligands are often separated by a significant gap in the sequence ID numeration.
  // But such gap can also mean that part of the chain is not modelled.
  auto last = std::min(e, chain.residues.end() - 1);
  for (auto it = b; it < last; ++it) {
    int gap = *(it+1)->seqid.num - *it->seqid.num;
    // The gap should be non-negative, but you can find exceptions in the PDB.
    if (gap < -1 || gap > 10)
      return it+1;
    // Usually polymers are longer than 1-2 residues, although there are
    // exceptions (example: 1-residue polymers in 5N22), so we can't be sure.
    // OTOH a protein can be capped with monomers different from amino-acid
    // and are_connected2() may return false negative. So if there is no gap
    // in numbering, it seems better to assume the polymer didn't end yet.
    if (gap == 1 && it - chain.residues.begin() < 2)
      continue;
    if (!are_connected2(*it, *(it+1), ptype))
      return it+1;
  }
  return e;
}

void add_entity_types(Chain& chain, bool overwrite) {
  if (!overwrite &&
      std::all_of(chain.residues.begin(), chain.residues.end(),
                  [](const Residue& r) { return r.entity_type != EntityType::Unknown; }))
    return;
  PolymerType ptype = check_polymer_type(chain.whole(), /*ignore_entity_type=*/overwrite);
  auto it = chain.residues.begin();
  if (ptype != PolymerType::Unknown) {
    auto polymer_end = infer_polymer_end(chain, ptype);
    for (; it != polymer_end; ++it)
      if (overwrite || it->entity_type == EntityType::Unknown)
        it->entity_type = EntityType::Polymer;
  }
  for (; it != chain.residues.end(); ++it)
    if (overwrite || it->entity_type == EntityType::Unknown)
      it->entity_type = it->is_water() ? EntityType::Water
                                       : EntityType::NonPolymer;
}

void add_entity_types(Structure& st, bool overwrite) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      add_entity_types(chain, overwrite);
}

void remove_entity_types(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        res.entity_type = EntityType::Unknown;
}

void add_entity_ids(Structure& st, bool overwrite) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (ResidueSpan& sub : chain.subchains()) {
        if (Entity* ent = st.get_entity_of(sub))
          for (Residue& res : sub) {
            if (overwrite || res.entity_id.empty())
              res.entity_id = ent->name;
          }
      }
}

void assign_subchain_names(Chain& chain, int& nonpolymer_counter) {
  for (Residue& res : chain.residues) {
    res.subchain = chain.name;
    // We'd use '-' as a separator (A-p or B-4 is more clear), but although
    // such names are valid in mmCIF, OneDep refuses to accept them.
    res.subchain += "x";
    switch (res.entity_type) {
      case EntityType::Polymer:
        res.subchain += 'p';
        break;
      case EntityType::NonPolymer:
        ++nonpolymer_counter;
        // to keep the name short use base36 for 2+ digit numbers:
        // 1, 2, ..., 9, 00, 01, ..., 09, 0A, 0B, ..., 0Z, 10, ...
        if (nonpolymer_counter < 10) {
          res.subchain += ('0' + nonpolymer_counter);
        } else {
          const char base36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
          int n = nonpolymer_counter - 10;
          if (n < 36)
            res.subchain += '0';
          size_t pos = res.subchain.size();
          while (n != 0) {
            res.subchain.insert(res.subchain.begin() + pos, base36[n % 36]);
            n /= 36;
          }
        }
        break;
      case EntityType::Water:
        res.subchain += 'w';
        break;
      // In the wwPDB branched are kept each in separate auth/label chain.
      // So we have one subchain in chain.
      case EntityType::Branched:
        res.subchain += 'b';
        break;
      case EntityType::Unknown:
        break;
    }
  }
}

static std::pair<bool,bool> has_entity_types_and_subchains(const Chain& chain) {
  bool has_entity_types = true;
  bool has_subchains = true;
  for (const Residue& res : chain.residues) {
    if (res.subchain.empty())
      has_subchains = false;
    if (res.entity_type == EntityType::Unknown)
      has_entity_types = false;
  }
  return {has_entity_types, has_subchains};
}

void assign_subchains(Structure& st, bool force, bool fail_if_unknown) {
  for (Model& model : st.models) {
    std::map<std::string, int> counters;
    for (Chain& chain : model.chains) {
      auto has = has_entity_types_and_subchains(chain);
      if (force || !has.second) {
        if (has.first)  // all chain's residues have known entity_type
          assign_subchain_names(chain, counters[chain.name]);
        else if (fail_if_unknown)
          fail("assign_subchains(): missing entity_type in chain " + chain.name);
      }
    }
  }
}

void ensure_entities(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (ResidueSpan& sub : chain.subchains()) {
        Entity* ent = st.get_entity_of(sub);
        if (!ent) {
          EntityType etype = sub[0].entity_type;
          std::string name = sub[0].entity_id;
          if (name.empty()) {
            if (etype == EntityType::Polymer || etype == EntityType::Branched)
              name = chain.name;
            else if (etype == EntityType::NonPolymer)
              name = sub[0].name + "!";
            else if (etype == EntityType::Water)
              name = "water";
          }
          if (!name.empty()) {
            ent = &impl::find_or_add(st.entities, name);
            ent->entity_type = etype;
            ent->subchains.push_back(sub.subchain_id());
          }
        }
        // ensure we have polymer_type set where needed
        if (ent && ent->entity_type == EntityType::Polymer &&
            ent->polymer_type == PolymerType::Unknown)
          ent->polymer_type = check_polymer_type(sub);
      }
}

static bool operator==(const Entity::DbRef& a, const Entity::DbRef& b) {
  return a.db_name == b.db_name &&
         a.id_code == b.id_code &&
         a.isoform == b.isoform &&
         a.seq_begin == b.seq_begin && a.seq_end == b.seq_end &&
         a.db_begin == b.db_begin && a.db_end == b.db_end;
}

void deduplicate_entities(Structure& st) {
  for (auto i = st.entities.begin(); i != st.entities.end(); ++i)
    if (!i->full_sequence.empty())
      for (auto j = i + 1; j != st.entities.end(); ++j)
        if (j->polymer_type == i->polymer_type &&
            j->full_sequence == i->full_sequence &&
            j->dbrefs == i->dbrefs) {
          vector_move_extend(i->subchains, std::move(j->subchains));
          st.entities.erase(j--);
        }
}

void change_ccd_code(Structure& st, const std::string& old, const std::string& new_) {
  auto process = [&](ResidueId& rid) {
    if (rid.name == old)
      rid.name = new_;
  };
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        process(res);
  for (Entity& ent : st.entities)
    for (std::string& mon_ids : ent.full_sequence)
      for (size_t start = 0;;) {
        size_t end = mon_ids.find(',', start);
        if (mon_ids.compare(start, end-start, old) == 0) {
          mon_ids.replace(start, end-start, new_);
          if (end != std::string::npos)
            end = start + new_.size();
        }
        if (end == std::string::npos)
          break;
        start = end + 1;
      }
  for (Connection& conn : st.connections) {
    process(conn.partner1.res_id);
    process(conn.partner2.res_id);
  }
  for (CisPep& cispep : st.cispeps) {
    process(cispep.partner_c.res_id);
    process(cispep.partner_n.res_id);
  }
  for (ModRes& modres : st.mod_residues)
    process(modres.res_id);
  for (Helix& helix : st.helices) {
    process(helix.start.res_id);
    process(helix.end.res_id);
  }
  for (Sheet& sheet : st.sheets)
    for (Sheet::Strand& strand : sheet.strands) {
      process(strand.start.res_id);
      process(strand.end.res_id);
      process(strand.hbond_atom2.res_id);
      process(strand.hbond_atom1.res_id);
    }
}

void shorten_ccd_codes(Structure& st) {
  // find all long residue names in both models and sequences
  auto already_in = [&](const std::string& name) -> bool {
    return in_vector_f([&](const OldToNew& x) { return x.old == name; },
                       st.shortened_ccd_codes);
  };
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        if (res.name.size() > 3 && !already_in(res.name))
          st.shortened_ccd_codes.push_back({res.name, {}});
  for (const Entity& ent : st.entities)
    for (const std::string& mon_ids : ent.full_sequence) {
      for (size_t start = 0;;) {
        size_t end = mon_ids.find(',', start);
        size_t len = std::min(end, mon_ids.size()) - start;
        if (len > 3) {
          std::string s(mon_ids, start, len);
          if (!already_in(s))
            st.shortened_ccd_codes.push_back({s, {}});
        }
        if (end == std::string::npos)
          break;
        start = end + 1;
      }
    }
  // pick a new residue name and call change_ccd_code()
  for (OldToNew& item : st.shortened_ccd_codes) {
    std::string short_code = item.old.substr(0, 3);
    short_code[2] = '~';
    // if short_code it's already used, try X[0-9]~ and then [0-9][0-9]~
    char c0 = '0', c1 = '0';
    while (in_vector_f([&](const OldToNew& x) { return x.new_ == short_code; },
                       st.shortened_ccd_codes)) {
      short_code[1] = c1++;
      if (c1 > '9') {
        short_code[0] = c0++;
        c1 = '0';
        if (c0 > 'Z')  // shouldn't happen
          break;
      }
    }
    item.new_ = short_code;
    change_ccd_code(st, item.old, item.new_);
  }
}

} // namespace gemmi
