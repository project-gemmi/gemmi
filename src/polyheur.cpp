// Copyright 2017-2022 Global Phasing Ltd.

#include <gemmi/polyheur.hpp>

namespace gemmi {

PolymerType check_polymer_type(const ConstResidueSpan& span) {
  if (span.empty())
    return PolymerType::Unknown;
  size_t counts[ResidueInfo::ELS+1] = {0};
  size_t aa = 0;
  size_t na = 0;
  size_t total = 0;
  bool has_atom_record = false;
  for (const Residue& r : span)
    if (r.entity_type == EntityType::Unknown ||
        r.entity_type == EntityType::Polymer) {
      if (r.het_flag == 'A')
        has_atom_record = true;
      ResidueInfo info = find_tabulated_residue(r.name);
      if (info.found()) {
        // Exclude water and ions - it can make difference
        // if function is called for the whole chain.
        if (info.kind == ResidueInfo::HOH)
          continue;
        // buffer molecules w/o hydrogens are mostly ions
        if (info.kind == ResidueInfo::BUF && info.hydrogen_count == 0)
          continue;
        if (info.is_peptide_linking())
          ++aa;
        if (info.is_na_linking())
          ++na;
        counts[info.kind]++;
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
    return counts[ResidueInfo::AA] >= counts[ResidueInfo::AAD]
           ? PolymerType::PeptideL : PolymerType::PeptideD;
  if (2 * na + bonus > total) {
    if (counts[ResidueInfo::DNA] == 0)
      return PolymerType::Rna;
    if (counts[ResidueInfo::RNA] == 0)
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
  PolymerType ptype = check_polymer_type(chain.whole());
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
          std::string name;
          if (etype == EntityType::Polymer)
            name = chain.name;
          else if (etype == EntityType::NonPolymer)
            name = sub[0].name + "!";
          else if (etype == EntityType::Water)
            name = "water";
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

} // namespace gemmi
