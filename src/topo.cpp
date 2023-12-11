// Copyright 2018-2022 Global Phasing Ltd.

#include <gemmi/topo.hpp>
#include <cmath>               // for sqrt, round
#include <map>                 // for multimap
#include <unordered_map>       // for unordered_multimap
#include <gemmi/polyheur.hpp>  // for get_or_check_polymer_type, ...
#include <gemmi/riding_h.hpp>  // for place_hydrogens_on_all_atoms, ...
#include <gemmi/modify.hpp>    // for remove_hydrogens

namespace gemmi {

std::unique_ptr<ChemComp> make_chemcomp_with_restraints(const Residue& res) {
  std::unique_ptr<ChemComp> cc(new ChemComp());
  cc->name = res.name;
  cc->group = ChemComp::Group::Null;
  cc->has_coordinates = true;
  // add atoms
  cc->atoms.reserve(res.atoms.size());
  for (const Atom& a : res.atoms) {
    if (!a.same_conformer(res.atoms[0]))
      continue;
    Element el = a.element;
    if (el == El::X)
      el = El::N;
    if (el == El::D)
      el = El::H;
    const std::string& chem_type = el.uname();
    cc->atoms.push_back(ChemComp::Atom{a.name, "", el, float(a.charge), chem_type, a.pos});
  }
  // prepare pairs of atoms
  struct Pair {
    size_t n1, n2;
    double dist;
  };
  std::vector<Pair> pairs;
  // first heavy atoms only
  for (size_t i = 0; i != res.atoms.size(); ++i) {
    const Atom& at1 = res.atoms[i];
    if (at1.is_hydrogen() || !at1.same_conformer(res.atoms[0]))
      continue;
    float r1 = at1.element.covalent_r();
    for (size_t j = i+1; j != res.atoms.size(); ++j) {
      const Atom& at2 = res.atoms[j];
      if (at2.is_hydrogen() || !at2.same_conformer(res.atoms[0]))
        continue;
      double d2 = at1.pos.dist_sq(at2.pos);
      float r2 = at2.element.covalent_r();
      double dmax = std::max(2.0, 1.3 * std::max(r1, r2));
      if (d2 < sq(dmax))
        pairs.push_back(Pair{i, j, std::sqrt(d2)});
    }
  }
  // now each hydrogen with the nearest heavy atom
  for (size_t i = 0; i != res.atoms.size(); ++i) {
    const Atom& at1 = res.atoms[i];
    if (at1.is_hydrogen() && at1.same_conformer(res.atoms[0])) {
      size_t nearest = (size_t)-1;
      double min_d2 = sq(2.5);
      for (size_t j = 0; j != res.atoms.size(); ++j) {
        const Atom& at2 = res.atoms[j];
        if (!at2.is_hydrogen() && at2.same_conformer(at1)) {
          double d2 = at1.pos.dist_sq(at2.pos);
          if (d2 < min_d2) {
            min_d2 = d2;
            nearest = j;
          }
        }
      }
      if (nearest != (size_t)-1) {
        pairs.push_back(Pair{nearest, i, std::sqrt(min_d2)});
      }
    }
  }

  // add bonds
  for (const Pair& p : pairs) {
    Restraints::Bond bond;
    bond.id1 = Restraints::AtomId{1, res.atoms[p.n1].name};
    bond.id2 = Restraints::AtomId{1, res.atoms[p.n2].name};
    bond.type = BondType::Unspec;
    bond.aromatic = false;
    double rounded_dist = 0.001 * std::round(1000 * p.dist);
    bond.value = bond.value_nucleus = rounded_dist;
    bond.esd = bond.esd_nucleus = 0.02;
    cc->rt.bonds.push_back(bond);
  }
  // add angles
  struct Triple {
    size_t n1, n2, n3;
  };
  std::vector<Triple> triples;
  for (size_t i = 0; i != pairs.size(); ++i)
    for (size_t j = i+1; j != pairs.size(); ++j) {
      if (pairs[i].n1 == pairs[j].n1)
        triples.push_back(Triple{pairs[i].n2, pairs[i].n1, pairs[j].n2});
      else if (pairs[i].n1 == pairs[j].n2)
        triples.push_back(Triple{pairs[i].n2, pairs[i].n1, pairs[j].n1});
      else if (pairs[i].n2 == pairs[j].n1)
        triples.push_back(Triple{pairs[i].n1, pairs[i].n2, pairs[j].n2});
      else if (pairs[i].n2 == pairs[j].n2)
        triples.push_back(Triple{pairs[i].n1, pairs[i].n2, pairs[j].n1});
    }
  for (const Triple& triple : triples) {
    Restraints::Angle angle;
    angle.id1 = Restraints::AtomId{1, res.atoms[triple.n1].name};
    angle.id2 = Restraints::AtomId{1, res.atoms[triple.n2].name};
    angle.id3 = Restraints::AtomId{1, res.atoms[triple.n3].name};
    double angle_rad = calculate_angle(res.atoms[triple.n1].pos,
                                       res.atoms[triple.n2].pos,
                                       res.atoms[triple.n3].pos);
    angle.value = 0.01 * std::round(100 * deg(angle_rad));
    angle.esd = 3.0;
    cc->rt.angles.push_back(angle);
  }
  return cc;
}

Topo::ChainInfo::ChainInfo(ResidueSpan& subchain, const Chain& chain, const Entity* ent)
  : chain_ref(chain) {
  subchain_name = subchain.at(0).subchain;
  res_infos.reserve(subchain.size());
  if (ent) {
    entity_id = ent->name;
    polymer = ent->entity_type == EntityType::Polymer;
    polymer_type = get_or_check_polymer_type(ent, subchain);
  } else {
    polymer = false;
    polymer_type = PolymerType::Unknown;
  }
  for (Residue& res : subchain)
    res_infos.emplace_back(&res);
}

// Add a ChemLink that restraints only bond length.
static
const std::string& add_auto_chemlink(MonLib& monlib,
                                     const std::string& resname1, const std::string& aname1,
                                     const std::string& resname2, const std::string& aname2,
                                     double ideal_dist, double esd) {
  ChemLink cl;
  cl.side1.comp = resname1;
  cl.side2.comp = resname2;
  cl.id = resname1 + resname2;
  cl.name = "auto-" + cl.id;
  cl.rt.bonds.push_back({Restraints::AtomId{1, aname1},
                         Restraints::AtomId{2, aname2},
                         BondType::Unspec, false,
                         ideal_dist, esd,
                         ideal_dist, esd});
  // ensure unique link id
  size_t orig_len = cl.id.size();
  for (int n = 0; monlib.get_link(cl.id) != nullptr; ++n)
    cl.id.replace(orig_len, cl.id.size(), std::to_string(n));

  auto it = monlib.links.emplace(cl.id, cl);
  return it.first->first;
}

static const ChemLink* setup_link_if_matches(Topo::Link& link, const MonLib& monlib,
                                             const std::string& atom1,
                                             const std::string& atom2) {
  bool invert;
  const ChemLink* match = nullptr;
  std::tie(match, invert, link.aliasing1, link.aliasing2) =
    monlib.match_link(*link.res1, atom1, link.alt1, *link.res2, atom2, link.alt2);
  if (match) {
    link.link_id = match->id;
    if (invert) {
      std::swap(link.res1, link.res2);
      std::swap(link.alt1, link.alt2);
      std::swap(link.aliasing1, link.aliasing2);
    }
  }
  return match;
}

static void add_polymer_links(PolymerType polymer_type,
                              const Topo::ResInfo& ri1,
                              Topo::ResInfo& ri2,
                              MonLib* monlib) {
  Topo::Link link;
  link.res1 = ri1.res;
  link.res2 = ri2.res;
  assert(&ri1 - &ri2 == link.res_distance());
  bool groups_ok = (ri1.orig_chemcomp != nullptr && ri2.orig_chemcomp != nullptr);

  if (is_polypeptide(polymer_type)) {
    std::string c = "C";
    std::string ca1 = "CA";
    std::string n = "N";
    std::string ca2 = "CA";
    if (ri1.orig_chemcomp && !ChemComp::is_peptide_group(ri1.orig_chemcomp->group)) {
      for (const ChemComp::Aliasing& aliasing : ri1.orig_chemcomp->aliases)
        if (ChemComp::is_peptide_group(aliasing.group)) {
          link.aliasing1 = &aliasing;
          if (const std::string* ptr = aliasing.name_from_alias(c))
            c = *ptr;
          if (const std::string* ptr = aliasing.name_from_alias(ca1))
            ca1 = *ptr;
          break;
        }
      if (!link.aliasing1)
        groups_ok = false;
    }
    ChemComp::Group n_terminus_group = ri2.orig_chemcomp ? ri2.orig_chemcomp->group
                                                         : ChemComp::Group::Null;
    if (ri2.orig_chemcomp && !ChemComp::is_peptide_group(ri2.orig_chemcomp->group)) {
      for (const ChemComp::Aliasing& aliasing : ri2.orig_chemcomp->aliases)
        if (ChemComp::is_peptide_group(aliasing.group)) {
          link.aliasing2 = &aliasing;
          n_terminus_group = aliasing.group;
          if (const std::string* ptr = aliasing.name_from_alias(n))
            n = *ptr;
          if (const std::string* ptr = aliasing.name_from_alias(ca2))
            ca2 = *ptr;
          break;
        }
      if (!link.aliasing2)
        groups_ok = false;
    }
    for (const Atom& a1 : ri1.res->atoms)
      if (a1.name == c && a1.element == El::C) {
        for (const Atom& a2 : ri2.res->atoms)
          if (a2.name == n && a2.element == El::N &&
              (a2.altloc == a1.altloc || a2.altloc == '\0' || a1.altloc == '\0') &&
              in_peptide_bond_distance(&a1, &a2)) {
            // One C-N pair of atoms can create here only one link.
            // Ignoring artificial configuration of no-altloc C and N atoms,
            // and CA atoms in 2+ conformations making both CIS and TRANS links.
            link.alt1 = a1.altloc;
            link.alt2 = a2.altloc;
            if (groups_ok) {
              // Deciding CIS/TRANS based on omega angle.
              char alt = a1.altloc_or(a2.altloc_or('*'));
              const Atom* ca1_atom = ri1.res->find_atom(ca1, alt, El::C);
              const Atom* ca2_atom = ri2.res->find_atom(ca2, alt, El::C);
              link.is_cis = is_peptide_bond_cis(ca1_atom, &a1, &a2, ca2_atom);
              if (n_terminus_group == ChemComp::Group::PPeptide)
                link.link_id = link.is_cis ? "PCIS" : "PTRANS";
              else if (n_terminus_group == ChemComp::Group::MPeptide)
                link.link_id = link.is_cis ? "NMCIS" : "NMTRANS";
              else
                link.link_id = link.is_cis ? "CIS" : "TRANS";
            } else if (monlib) {
              if (!setup_link_if_matches(link, *monlib, c, n))
                link.link_id = add_auto_chemlink(*monlib,
                                                 ri1.res->name, c,
                                                 ri2.res->name, n,
                                                 1.34, 0.04);
            }
            ri2.prev.push_back(link);
          }
      }

  } else if (is_polynucleotide(polymer_type)) {
    std::string o3p = "O3'";
    std::string p = "P";
    if (ri1.orig_chemcomp && !ChemComp::is_nucleotide_group(ri1.orig_chemcomp->group)) {
      for (const ChemComp::Aliasing& aliasing : ri1.orig_chemcomp->aliases)
        if (ChemComp::is_nucleotide_group(aliasing.group)) {
          link.aliasing1 = &aliasing;
          if (const std::string* o3p_ptr = aliasing.name_from_alias(o3p))
            o3p = *o3p_ptr;
          break;
        }
      if (!link.aliasing1)
        groups_ok = false;
    }
    if (ri2.orig_chemcomp && !ChemComp::is_nucleotide_group(ri2.orig_chemcomp->group)) {
      for (const ChemComp::Aliasing& aliasing : ri2.orig_chemcomp->aliases)
        if (ChemComp::is_nucleotide_group(aliasing.group)) {
          link.aliasing2 = &aliasing;
          if (const std::string* p_ptr = aliasing.name_from_alias(p))
            p = *p_ptr;
          break;
        }
      if (!link.aliasing2)
        groups_ok = false;
    }
    for (const Atom& a1 : ri1.res->atoms)
      if (a1.name == o3p && a1.element == El::O) {
        for (const Atom& a2 : ri2.res->atoms)
          if (a2.name == p && a2.element == El::P &&
              (a2.altloc == a1.altloc || a2.altloc == '\0' || a1.altloc == '\0') &&
              in_nucleotide_bond_distance(&a1, &a2)) {
            link.alt1 = a1.altloc;
            link.alt2 = a2.altloc;
            if (groups_ok) {
              link.link_id = "p";
            } else if (monlib) {
              if (!setup_link_if_matches(link, *monlib, o3p, p))
                link.link_id = add_auto_chemlink(*monlib,
                                                 ri1.res->name, o3p,
                                                 ri2.res->name, p,
                                                 1.606, 0.02);
            }
            ri2.prev.push_back(link);
          }
    }
  }

  if (ri2.prev.empty()) {
    link.link_id = "gap";
    ri2.prev.push_back(link);
  }
}

double Topo::ideal_chiral_abs_volume(const Chirality &ch) const {
  const Restraints::Bond* bond_c1 = take_bond(ch.atoms[0], ch.atoms[1]);
  const Restraints::Bond* bond_c2 = take_bond(ch.atoms[0], ch.atoms[2]);
  const Restraints::Bond* bond_c3 = take_bond(ch.atoms[0], ch.atoms[3]);
  const Restraints::Angle* angle_1c2 = take_angle(ch.atoms[1], ch.atoms[0], ch.atoms[2]);
  const Restraints::Angle* angle_2c3 = take_angle(ch.atoms[2], ch.atoms[0], ch.atoms[3]);
  const Restraints::Angle* angle_3c1 = take_angle(ch.atoms[3], ch.atoms[0], ch.atoms[1]);
  if (bond_c1 && bond_c2 && bond_c3 && angle_1c2 && angle_2c3 && angle_3c1)
    return chiral_abs_volume(bond_c1->value, bond_c2->value, bond_c3->value,
                             angle_1c2->value, angle_2c3->value, angle_3c1->value);
  return std::numeric_limits<double>::quiet_NaN();
}

std::vector<Topo::Rule> Topo::apply_restraints(const Restraints& rt,
                                               Residue& res, Residue* res2, Asu asu,
                                               char altloc1, char altloc2, bool require_alt) {
  std::string altlocs;
  if (altloc1 == '\0' && altloc2 == '\0') {
    add_distinct_altlocs(res, altlocs);
    if (res2)
      add_distinct_altlocs(*res2, altlocs);
  }
  if (altlocs.empty())
    altlocs += altloc1 != 0 ? altloc1 : altloc2;

  std::vector<Rule> rules;
  for (const Restraints::Bond& bond : rt.bonds)
    for (char alt : altlocs)
      if (Atom* at1 = bond.id1.get_from(res, res2, alt, altloc2))
        if (Atom* at2 = bond.id2.get_from(res, res2, alt, altloc2)) {
          bool with_alt = at1->altloc || at2->altloc;
          if (with_alt || !require_alt) {
            rules.push_back({RKind::Bond, bonds.size()});
            bonds.push_back({&bond, {{at1, at2}}, asu});
          }
          if (!with_alt)
            break;
        }

  if (only_bonds)
    return rules;

  for (const Restraints::Angle& angle : rt.angles)
    for (char alt : altlocs)
      if (Atom* at1 = angle.id1.get_from(res, res2, alt, altloc2))
        if (Atom* at2 = angle.id2.get_from(res, res2, alt, altloc2))
          if (Atom* at3 = angle.id3.get_from(res, res2, alt, altloc2)) {
            bool with_alt = at1->altloc || at2->altloc || at3->altloc;
            if (with_alt || !require_alt) {
              rules.push_back({RKind::Angle, angles.size()});
              angles.push_back({&angle, {{at1, at2, at3}}});
            }
            if (!with_alt)
              break;
          }
  for (const Restraints::Torsion& tor : rt.torsions)
    for (char alt : altlocs)
      if (Atom* at1 = tor.id1.get_from(res, res2, alt, altloc2))
        if (Atom* at2 = tor.id2.get_from(res, res2, alt, altloc2))
          if (Atom* at3 = tor.id3.get_from(res, res2, alt, altloc2))
            if (Atom* at4 = tor.id4.get_from(res, res2, alt, altloc2)) {
              bool with_alt = at1->altloc || at2->altloc || at3->altloc || at4->altloc;
              if (with_alt || !require_alt) {
                rules.push_back({RKind::Torsion, torsions.size()});
                torsions.push_back({&tor, {{at1, at2, at3, at4}}});
              }
              if (!with_alt)
                break;
        }
  for (const Restraints::Chirality& chir : rt.chirs)
    for (char alt : altlocs)
      if (Atom* at1 = chir.id_ctr.get_from(res, res2, alt, altloc2))
        if (Atom* at2 = chir.id1.get_from(res, res2, alt, altloc2))
          if (Atom* at3 = chir.id2.get_from(res, res2, alt, altloc2))
            if (Atom* at4 = chir.id3.get_from(res, res2, alt, altloc2)) {
              bool with_alt = at1->altloc || at2->altloc || at3->altloc || at4->altloc;
              if (with_alt || !require_alt) {
                rules.push_back({RKind::Chirality, chirs.size()});
                chirs.push_back({&chir, {{at1, at2, at3, at4}}});
              }
              if (!with_alt)
                break;
            }
  for (const Restraints::Plane& plane : rt.planes)
    for (char alt : altlocs) {
      std::vector<Atom*> atoms;
      bool with_alt = false;
      for (const Restraints::AtomId& id : plane.ids)
        if (Atom* atom = id.get_from(res, res2, alt, altloc2)) {
          with_alt = with_alt || atom->altloc;
          atoms.push_back(atom);
        }
      if (atoms.size() >= 4 && (with_alt || !require_alt)) {
        rules.push_back({RKind::Plane, planes.size()});
        planes.push_back({&plane, atoms});
      }
      if (!with_alt)
        break;
    }
  return rules;
}

void Topo::apply_restraints_from_link(Link& link, const MonLib& monlib) {
  if (link.link_id.empty())
    return;
  const ChemLink* chem_link = monlib.get_link(link.link_id);
  if (!chem_link) {
    err("ignoring link '" + link.link_id + "' as it is not in the monomer library");
    return;
  }
  const Restraints* rt = &chem_link->rt;
  if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
    err(cat("LINK between different conformers: ", link.alt1, " (in ",
            link.res1->name, ") and ", link.alt2, " (in " + link.res2->name, ")."));
  // aliases are a new feature - introduced in 2022
  if (link.aliasing1 || link.aliasing2) {
    std::unique_ptr<Restraints> rt_copy(new Restraints(*rt));
    if (link.aliasing1)
      for (const auto& p : link.aliasing1->related)
        rt_copy->rename_atom(Restraints::AtomId{1, p.second}, p.first);
    if (link.aliasing2)
      for (const auto& p : link.aliasing2->related)
        rt_copy->rename_atom(Restraints::AtomId{2, p.second}, p.first);
    rt = rt_copy.get();
    rt_storage.push_back(std::move(rt_copy));
  }
  link.link_rules = apply_restraints(*rt, *link.res1, link.res2, link.asu,
                                     link.alt1, link.alt2, false);
}

// see comments above the declaration
void Topo::initialize_refmac_topology(Structure& st, Model& model0,
                                      MonLib& monlib, bool ignore_unknown_links) {
  // initialize chain_infos
  for (Chain& chain : model0.chains)
    for (ResidueSpan& sub : chain.subchains()) {
      // set Residue::group_idx which is used in Restraints::AtomId::get_from()
      for (size_t i = 0; i != sub.size(); ++i) {
        sub[i].group_idx = 0;
        if (i != 0 && sub[i-1].seqid == sub[i].seqid)
          sub[i].group_idx = sub[i-1].group_idx + 1;
      }
      const Entity* ent = st.get_entity_of(sub);
      chain_infos.emplace_back(sub, chain, ent);
    }

  // add modifications from MODRES records (Refmac's extension)
  for (const ModRes& modres : st.mod_residues)
    if (!modres.mod_id.empty()) {
      for (Topo::ChainInfo& chain_info : chain_infos)
        if (chain_info.chain_ref.name == modres.chain_name) {
          for (Topo::ResInfo& res_info : chain_info.res_infos)
            if (*res_info.res == modres.res_id)
              res_info.add_mod(modres.mod_id, nullptr, '\0');
        }
    }

  // setup pointers to monomers and links in the polymer
  for (ChainInfo& ci : chain_infos) {
    for (ResInfo& ri : ci.res_infos) {
      auto it = monlib.monomers.find(ri.res->name);
      if (it != monlib.monomers.end())
        ri.orig_chemcomp = &it->second;
      // orig_chemcomp is left null for missing monomers
    }
    // setup polymer links
    if (ci.polymer && !ci.res_infos.empty()) {
      // handling of microheterogeneities makes it more complicated
      // and it should be even more complex to handle partial bonding
      auto prev_begin = ci.res_infos.begin();
      auto prev_end = ci.group_end(prev_begin);
      while (prev_end != ci.res_infos.end()) {
        auto group_begin = prev_end;
        auto group_end = ci.group_end(group_begin);
        for (auto ri = group_begin; ri != group_end; ++ri)
          for (auto prev_ri = prev_begin; prev_ri != prev_end; ++prev_ri) {
            MonLib* monlib_ptr = ignore_unknown_links ? nullptr : &monlib;
            add_polymer_links(ci.polymer_type, *prev_ri, *ri, monlib_ptr);
          }
        prev_begin = group_begin;
        prev_end = group_end;
      }
    }
  }

  // add extra links
  for (Connection& conn : st.connections)
    if (conn.type != Connection::Hydrog) // ignoring hydrogen bonds
      setup_connection(conn, model0, monlib, ignore_unknown_links);

  // Add modifications from standard links. We do it here b/c polymer links
  // could be disabled (link_id.clear()) in setup_connection().
  for (ChainInfo& ci : chain_infos)
    for (ResInfo& ri : ci.res_infos)
      for (Link& prev : ri.prev)
        if (const ChemLink* chem_link = monlib.get_link(prev.link_id)) {
          ResInfo* ri_prev = &ri + prev.res_distance();
          ri_prev->add_mod(chem_link->side1.mod, prev.aliasing1, prev.alt1);
          ri.add_mod(chem_link->side2.mod, prev.aliasing2, prev.alt2);
        }

  // Apply modifications to monomer restraints
  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& ri : chain_info.res_infos) {
      if (ri.orig_chemcomp) {
        // The final ChemComp restraints that we'll use are made from
        // original (_chem_comp) restraints with modifications (_chem_mod)
        // applied. There is a corner case in which different conformations
        // of the residue have different modifications applied.
        bool has_mod_altlocs = false;
        for (const Mod& mod : ri.mods)
          if (mod.altloc != '\0')
            has_mod_altlocs = true;
        std::string altlocs;
        if (has_mod_altlocs)
          add_distinct_altlocs(*ri.res, altlocs);  // cf. apply_restraints
        if (altlocs.empty())
          altlocs += '\0';
        for (char altloc : altlocs) {
          // key for caching in Topo::cc_cache: ChemComp::name + modifications
          std::string key = ri.orig_chemcomp->name;
          for (const Mod& mod : ri.mods)
            if (!mod.altloc || altloc == mod.altloc) {
              key += char(1 + static_cast<int>(mod.alias));
              key += mod.id;
            }
          auto it = cc_cache.find(key);
          if (it != cc_cache.end()) {
            ri.chemcomps.push_back({altloc, it->second.get()});
          } else {  // it's not in cache yet - we need to add it
            std::unique_ptr<ChemComp> cc_copy(new ChemComp(*ri.orig_chemcomp));
            // apply modifications
            for (const Mod& mod : ri.mods) {
              if (!mod.altloc || altloc == mod.altloc) {
                if (const ChemMod* chem_mod = monlib.get_mod(mod.id)) {
                  try {
                    chem_mod->apply_to(*cc_copy, mod.alias);
                  } catch(std::runtime_error& e) {
                    err(cat("failed to apply modification ", chem_mod->id,
                            " to ", ri.res->name, ": ", e.what()));
                  }
                } else {
                  err("modification not found: " + mod.id);
                }
              }
            }
            ri.chemcomps.push_back({altloc, cc_copy.get()});
            cc_cache.emplace(key, std::move(cc_copy));
          }
        }
        // Usually the same modifications are applied to all conformers.
        // In such case reduce chemcomps to a single value.
        if (ri.chemcomps.size() > 1 &&
            std::all_of(ri.chemcomps.begin() + 1, ri.chemcomps.end(),
              [&](const FinalChemComp& f) { return f.cc == ri.chemcomps[0].cc; })) {
          ri.chemcomps.resize(1);
          ri.chemcomps[0].altloc = '\0';
        }
      } else {  // orig_chemcomp not set - make ChemComp with ad-hoc restraints
        // no cache - ad-hoc restraints are separate for each residue
        cc_storage.emplace_back(make_chemcomp_with_restraints(*ri.res));
        ri.chemcomps.push_back({'\0', cc_storage.back().get()});
      }
    }
}

void Topo::apply_all_restraints(const MonLib& monlib) {
  bonds.clear();
  angles.clear();
  torsions.clear();
  chirs.clear();
  planes.clear();
  rt_storage.clear();
  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& ri : chain_info.res_infos) {
      // link restraints
      for (Link& link : ri.prev)
        apply_restraints_from_link(link, monlib);
      // monomer restraints
      auto it = ri.chemcomps.cbegin();
      ri.monomer_rules = apply_restraints(it->cc->rt, *ri.res, nullptr, Asu::Same,
                                          it->altloc, '\0', /*require_alt=*/false);
      while (++it != ri.chemcomps.end()) {
        auto rules = apply_restraints(it->cc->rt, *ri.res, nullptr, Asu::Same,
                                      it->altloc, '\0', /*require_alt=*/true);
        vector_move_extend(ri.monomer_rules, std::move(rules));
      }
    }
  for (Link& link : extras)
    apply_restraints_from_link(link, monlib);
}

void Topo::create_indices() {
  bond_index.clear();
  for (Bond& bond : bonds) {
    bond_index.emplace(bond.atoms[0], &bond);
    if (bond.atoms[1] != bond.atoms[0])
      bond_index.emplace(bond.atoms[1], &bond);
  }
  angle_index.clear();
  for (Angle& ang : angles)
    angle_index.emplace(ang.atoms[1], &ang);
  torsion_index.clear();
  for (Torsion& tor : torsions) {
    torsion_index.emplace(tor.atoms[1], &tor);
    if (tor.atoms[1] != tor.atoms[2])
      torsion_index.emplace(tor.atoms[2], &tor);
  }
  plane_index.clear();
  for (Plane& plane : planes)
    for (Atom* atom : plane.atoms)
      plane_index.emplace(atom, &plane);
}

// Tries to construct Topo::Link and append it to extras.
// Side-effects: it may modify conn.link_id and add ChemLink to monlib.links.
void Topo::setup_connection(Connection& conn, Model& model0, MonLib& monlib,
                            bool ignore_unknown_links) {
  if (conn.link_id == "gap") {
    Link* polymer_link = find_polymer_link(conn.partner1, conn.partner2);
    if (polymer_link)
      polymer_link->link_id.clear();  // disable polymer link
    return;
  }

  Link extra;
  CRA cra1 = model0.find_cra(conn.partner1, true);
  CRA cra2 = model0.find_cra(conn.partner2, true);
  if (!cra1.atom || !cra2.atom)
    return;
  extra.res1 = cra1.residue;
  extra.res2 = cra2.residue;
  extra.alt1 = conn.partner1.altloc;
  extra.alt2 = conn.partner2.altloc;
  extra.asu = conn.asu;

  const ChemLink* match = nullptr;

  // If we have link_id find ChemLink by name (and check if it matches).
  if (!conn.link_id.empty()) {
    match = monlib.get_link(conn.link_id);
    if (!match) {
      err("link not found in monomer library: " + conn.link_id);
      return;
    }
    if (match->rt.bonds.empty() ||
        !monlib.link_side_matches_residue(match->side1, extra.res1->name, &extra.aliasing1) ||
        !monlib.link_side_matches_residue(match->side2, extra.res2->name, &extra.aliasing2) ||
        !atom_match_with_alias(match->rt.bonds[0].id1.atom,
                               conn.partner1.atom_name, extra.aliasing1) ||
        !atom_match_with_alias(match->rt.bonds[0].id2.atom,
                               conn.partner2.atom_name, extra.aliasing2)) {
      err("link from the monomer library does not match: " + conn.link_id);
      return;
    }
  } else {
    // we don't have link_id - use the best matching link (if any)
    match = setup_link_if_matches(extra, monlib,
                                  conn.partner1.atom_name, conn.partner2.atom_name);
    if (match)
      conn.link_id = match->id;
  }

  // If a polymer link is also given in LINK/struct_conn,
  // use only one of them. If LINK has explicit name (ccp4_link_id),
  // or if it matches residue-specific link from monomer library, use it;
  // otherwise, LINK is repetition of TRANS/CIS, so ignore LINK.
  if (Link* polymer_link = find_polymer_link(conn.partner1, conn.partner2)) {
    if (conn.link_id.empty() && !cif::is_null(polymer_link->link_id) &&
        polymer_link->link_id != "gap" &&
        (!match || (match->side1.comp.empty() && match->side2.comp.empty())))
      return;
    polymer_link->link_id.clear();  // disable polymer link
  }

  if (match) {
    extra.link_id = match->id;
    // add modifications from the link
    find_resinfo(extra.res1)->add_mod(match->side1.mod, extra.aliasing1, extra.alt1);
    find_resinfo(extra.res2)->add_mod(match->side2.mod, extra.aliasing2, extra.alt2);
  } else {
    if (ignore_unknown_links)
      return;
    // create a new ChemLink and add it to the monomer library
    double ideal_dist = monlib.find_ideal_distance(cra1, cra2);
    extra.link_id = add_auto_chemlink(monlib,
                                      extra.res1->name, conn.partner1.atom_name,
                                      extra.res2->name, conn.partner2.atom_name,
                                      ideal_dist, 0.02);
  }
  extras.push_back(extra);
}

void Topo::set_cispeps_in_structure(Structure& st) {
  st.cispeps.clear();
  if (chain_infos.empty())
    return;
  // model is not stored in Topo, let's determine it from chain_infos[0]
  std::string model_str;
  for (const Model& model : st.models)
    if (!model.chains.empty() &&
        &model.chains[0] == &chain_infos[0].chain_ref)
      model_str = model.name;
  for (ChainInfo& chain_info : chain_infos)
    for (ResInfo& res_info : chain_info.res_infos)
      for (Link& link : res_info.prev)
        if (link.is_cis) {
          CisPep cp;
          cp.model_str = model_str;
          cp.partner_c = AtomAddress(chain_info.chain_ref.name, *link.res1, "", link.alt1);
          cp.partner_n = AtomAddress(chain_info.chain_ref.name, *link.res2, "", link.alt2);
          cp.only_altloc = link.alt1 ? link.alt1 : link.alt2;
          for (const Rule& rule : link.link_rules)
            if (rule.rkind == RKind::Torsion) {
              const Torsion& tor = torsions[rule.index];
              if (tor.restr->label == "omega") {
                cp.reported_angle = tor.calculate();
                break;
              }
            }
          st.cispeps.push_back(cp);
        }
}

namespace {

void remove_hydrogens_from_atom(Topo::ResInfo* ri,
                                       const std::string& atom_name, char alt) {
  if (!ri)
    return;
  std::vector<Atom>& atoms = ri->res->atoms;
  const Restraints& rt = ri->get_final_chemcomp(alt).rt;
  for (auto it = atoms.end(); it-- != atoms.begin(); ) {
    if (it->is_hydrogen() && is_same_conformer(it->altloc, alt)) {
      const Restraints::AtomId* heavy = rt.first_bonded_atom(it->name);
      if (heavy && heavy->atom == atom_name)
        atoms.erase(it);
    }
  }
}

void set_cis_in_link(Topo::Link& link, bool is_cis) {
  if (is_cis) {
    if (ends_with(link.link_id, "TRANS"))
      link.link_id.replace(link.link_id.size() - 5, 5, "CIS");
  } else {
    if (ends_with(link.link_id, "CIS"))
      link.link_id.replace(link.link_id.size() - 3, 3, "TRANS");
  }
  link.is_cis = is_cis;
}

void force_cispeps(Topo& topo, bool single_model, const Model& model,
                          const std::vector<CisPep>& cispeps,
                          std::ostream* warnings) {
  std::multimap<const Residue*, const CisPep*> cispep_index;
  for (const CisPep& cp : cispeps) {
    if (single_model || model.name == cp.model_str)
      if (const Residue* res_n = model.find_cra(cp.partner_n).residue)
        cispep_index.emplace(res_n, &cp);
  }
  for (Topo::ChainInfo& chain_info : topo.chain_infos)
    for (Topo::ResInfo& res_info : chain_info.res_infos) {
      auto range = cispep_index.equal_range(res_info.res);
      for (Topo::Link& link : res_info.prev) {
        bool is_cis = false;
        for (auto i = range.first; i != range.second; ++i) {
          const CisPep& cp = *i->second;
          if (cp.partner_c.res_id.matches_noseg(*link.res1) &&
              (cp.only_altloc == '\0' || cp.only_altloc == link.alt1
                                      || cp.only_altloc == link.alt2))
            is_cis = true;
        }
        if (is_cis != link.is_cis) {
          set_cis_in_link(link, is_cis);
          if (warnings)
            *warnings << "Link between "
                      << atom_str(chain_info.chain_ref.name, *link.res1, "", link.alt1)
                      << " and "
                      << atom_str(chain_info.chain_ref.name, *link.res2, "", link.alt2)
                      << " forced to " << link.link_id << std::endl;
        }
      }
    }
}

struct Neigh {
  char alt;
  float occ;
};
using NeighMap = std::unordered_multimap<int, Neigh>;

// Assumes no hydrogens in the residue.
// Position and serial number are not assigned for new atoms.
void add_hydrogens_without_positions(Topo::ResInfo& ri, const NeighMap& neighbor_altlocs) {
  std::vector<Atom>& atoms = ri.res->atoms;
  // Add H atom for each conformation (altloc) of the parent atom and its
  // first neighbors.
  for (size_t i = 0, size = atoms.size(); i != size; ++i) {
    if (atoms[i].calc_flag == CalcFlag::NoHydrogen)
      continue;
    char parent_alt = atoms[i].altloc;
    float parent_occ = atoms[i].occ;
    std::map<char, float> altlocs; // altloc + occupancy
    if (parent_alt == '\0') {
      // We could have neighbor 1 with altlocs A,B and neighbor 2 with C,D.
      // Or A,B and A,B,C. Or, neighbor 1 could have different occupancy
      // of altloc A than neighbor 2. We don't want to check for all possible
      // weird cases. Only ensure that the added H has total occupancy <= 1.
      float occ_sum = 0;
      auto range = neighbor_altlocs.equal_range(atoms[i].serial);
      for (auto it = range.first; it != range.second; ++it) {
        const Neigh& neigh = it->second;
        if (altlocs.count(neigh.alt) == 0 && occ_sum + neigh.occ <= 1.001f) {
          occ_sum += neigh.occ;
          altlocs.emplace(neigh.alt, neigh.occ);
        }
      }
      if (occ_sum - parent_occ > 0.001f)
        for (auto& it : altlocs)
          it.second *= parent_occ;
    }
    if (altlocs.empty())
      altlocs.emplace(parent_alt, parent_occ);
    const ChemComp& cc = ri.get_final_chemcomp(parent_alt);
    for (const Restraints::Bond& bond : cc.rt.bonds) {
      // atoms may get re-allocated, so we can't set parent earlier
      const Atom& parent = atoms[i];
      assert(!parent.is_hydrogen());
      const Restraints::AtomId* atom_id = bond.other(parent.name);
      if (!atom_id)
        continue;
      auto it = cc.find_atom(atom_id->atom);
      if (it == cc.atoms.end())
        fail("inconsistent _chem_comp " + cc.name);
      if (it->is_hydrogen()) {
        gemmi::Atom atom;
        atom.name = it->id;
        atom.element = it->el;
        // calc_flag will be changed to Calculated when the position is set
        atom.calc_flag = CalcFlag::Dummy;
        atom.b_iso = parent.b_iso;
        for (auto alt_occ : altlocs) {
          atom.altloc = alt_occ.first;
          atom.occ = alt_occ.second;
          atoms.push_back(atom);
        }
      }
    }
  }
}

NeighMap prepare_neighbor_altlocs(Topo& topo, const MonLib& monlib) {
  // disable warnings here, so they are not printed twice
  std::streambuf* warnings_orig = nullptr;
  if (topo.warnings)
    warnings_orig = topo.warnings->rdbuf(nullptr);
  // Prepare bonds. Fills topo.bonds, monomer_rules/link_rules and rt_storage,
  // but they are all reset when apply_all_restraints() is called again.
  topo.only_bonds = true;
  topo.apply_all_restraints(monlib);
  topo.only_bonds = false;
  // re-enable warnings
  if (warnings_orig)
    topo.warnings->rdbuf(warnings_orig);
  NeighMap neighbors;
  for (const Topo::Bond& bond : topo.bonds) {
    const Atom* a1 = bond.atoms[0];
    const Atom* a2 = bond.atoms[1];
    if (a2->altloc)
      neighbors.emplace(a1->serial, Neigh{a2->altloc, a2->occ});
    if (a1->altloc && a1 != a2)
      neighbors.emplace(a2->serial, Neigh{a1->altloc, a1->occ});
  }
  return neighbors;
}

}  // anonymous namespace

std::unique_ptr<Topo>
prepare_topology(Structure& st, MonLib& monlib, size_t model_index,
                 HydrogenChange h_change, bool reorder,
                 std::ostream* warnings, bool ignore_unknown_links, bool use_cispeps) {
  std::unique_ptr<Topo> topo(new Topo);
  topo->warnings = warnings;
  if (model_index >= st.models.size())
    fail("no such model index: " + std::to_string(model_index));
  Model& model = st.models[model_index];
  topo->initialize_refmac_topology(st, model, monlib, ignore_unknown_links);

  if (use_cispeps)
    force_cispeps(*topo, st.models.size() == 1, model, st.cispeps, warnings);

  // remove hydrogens, or change deuterium to fraction, or nothing
  // and then check atom names
  for (Topo::ChainInfo& chain_info : topo->chain_infos)
    for (Topo::ResInfo& ri : chain_info.res_infos) {
      Residue& res = *ri.res;
      if (h_change != HydrogenChange::NoChange && h_change != HydrogenChange::Shift
          // don't re-add H's if we don't have chemical component description
          && (ri.orig_chemcomp != nullptr || h_change == HydrogenChange::Remove)) {
        // remove/add hydrogens
        remove_hydrogens(res);
      } else {
        // Special handling of Deuterium - mostly for Refmac.
        // Note: if the model has deuterium, it gets modified.
        if (replace_deuterium_with_fraction(res)) {
          // deuterium names usually differ from the names in dictionary
          for (Atom& atom : res.atoms)
            if (atom.name[0] == 'D' && atom.fraction != 0) {
              const ChemComp& cc = ri.get_final_chemcomp(atom.altloc);
              if (cc.find_atom(atom.name) == cc.atoms.end())
                atom.name[0] = 'H';
            }
          st.has_d_fraction = true;
        }
      }
      // check atom names
      for (Atom& atom : res.atoms) {
        const ChemComp& cc = ri.get_final_chemcomp(atom.altloc);
        if (!cc.has_atom(atom.name)) {
          std::string msg = "definition not found for "
                          + atom_str(chain_info.chain_ref, *ri.res, atom);
          if (ri.orig_chemcomp && ri.orig_chemcomp->has_atom(atom.name)) {
            msg += " (linkage should remove this atom)";
          } else {
            auto it = cc.find_atom_by_old_name(atom.name);
            if (it != cc.atoms.end())
              cat_to(msg, " (replace ", atom.name, " with ", it->id, ')');
          }
          topo->err(msg);
        }
      }
    }

  // add hydrogens
  if (h_change == HydrogenChange::ReAdd ||
      h_change == HydrogenChange::ReAddButWater ||
      h_change == HydrogenChange::ReAddKnown) {
    NeighMap neighbor_altlocs = prepare_neighbor_altlocs(*topo, monlib);
    for (Topo::ChainInfo& chain_info : topo->chain_infos)
      for (Topo::ResInfo& ri : chain_info.res_infos) {
        Residue& res = *ri.res;
        if (ri.orig_chemcomp != nullptr &&
            (h_change == HydrogenChange::ReAdd || !res.is_water())) {
          add_hydrogens_without_positions(ri, neighbor_altlocs);

          // a special handling of HIS for compatibility with Refmac
          if (res.name == "HIS") {
            for (gemmi::Atom& atom : ri.res->atoms)
              if (atom.name == "HD1" || atom.name == "HE2")
                atom.occ = 0;
          }
        }
      }
  }

  // sort atoms in residues
  for (Topo::ChainInfo& chain_info : topo->chain_infos)
    for (Topo::ResInfo& ri : chain_info.res_infos)
      if (reorder && ri.orig_chemcomp) {
        Residue& res = *ri.res;
        const ChemComp& cc = *ri.orig_chemcomp;
        for (Atom& atom : res.atoms) {
          auto it = cc.find_atom(atom.name);
          // If atom.name is not found (b/c it was added in a modification),
          // the atom will be put after original atoms.
          atom.serial = int(it - cc.atoms.begin()); // temporary, for sorting only
        }
        std::sort(res.atoms.begin(), res.atoms.end(), [](const Atom& a, const Atom& b) {
                    return a.serial != b.serial ? a.serial < b.serial
                                                : a.altloc < b.altloc;
        });
        // check for missing altloc
        for (auto atom = res.atoms.begin(); atom + 1 < res.atoms.end(); ++atom)
          if (atom->name == (atom + 1)->name && atom->altloc == '\0')
            topo->err("missing altloc in " + atom_str(chain_info.chain_ref, *ri.res, *atom));
      }

  // for atoms with ad-hoc links, for now we don't want hydrogens
  if (!ignore_unknown_links && h_change != HydrogenChange::NoChange) {
    auto remove_h_from_auto_links = [&](const Topo::Link& link) {
      const ChemLink* cl = monlib.get_link(link.link_id);
      if (cl && starts_with(cl->name, "auto-")) {
        const Restraints::Bond& bond = cl->rt.bonds.at(0);
        remove_hydrogens_from_atom(topo->find_resinfo(link.res1), bond.id1.atom, link.alt1);
        remove_hydrogens_from_atom(topo->find_resinfo(link.res2), bond.id2.atom, link.alt2);
      }
    };
    for (const Topo::ChainInfo& chain_info : topo->chain_infos)
      for (const Topo::ResInfo& res_info : chain_info.res_infos)
        for (const Topo::Link& link : res_info.prev)
          remove_h_from_auto_links(link);
    for (const Topo::Link& link : topo->extras)
      remove_h_from_auto_links(link);
  }

  // fill Topo::bonds, angles, ... and ResInfo::monomer_rules, Links::link_rules
  topo->apply_all_restraints(monlib);
  // fill bond_index, angle_index, etc
  topo->create_indices();

  // the hydrogens added previously have positions not set
  if (h_change != HydrogenChange::NoChange)
    place_hydrogens_on_all_atoms(*topo);

  if (h_change == HydrogenChange::ReAddKnown) {
    // To leave only known hydrogens, we remove Hs with zero occupancy.
    // As a side-effect, it removes any H atoms on zero-occupancy parents.
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        vector_remove_if(res.atoms, [](Atom& a) { return a.is_hydrogen() && a.occ == 0; });

    // disable warnings here, so they are not printed twice
    std::streambuf *warnings_orig = nullptr;
    if (topo->warnings)
      warnings_orig = topo->warnings->rdbuf(nullptr);
    // re-set restraints and indices
    topo->apply_all_restraints(monlib);
    topo->create_indices();
    // re-enable warnings
    if (warnings_orig)
      topo->warnings->rdbuf(warnings_orig);
  }

  assign_serial_numbers(model);
  return topo;
}

std::vector<AtomAddress> find_missing_atoms(const Topo& topo, bool including_hydrogen) {
  std::vector<AtomAddress> ret;
  for (const Topo::ChainInfo& chain_info : topo.chain_infos)
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      // check only the first conformation
      const Topo::FinalChemComp& fcc = ri.chemcomps.at(0);
      char altloc = fcc.altloc != '\0' ? fcc.altloc : '*';
      for (const ChemComp::Atom& at : fcc.cc->atoms)
        if ((including_hydrogen || !at.is_hydrogen()) &&
            ri.res->find_atom(at.id, altloc) == nullptr)
          ret.emplace_back(chain_info.chain_ref.name, *ri.res, at.id, fcc.altloc);
    }
  return ret;
}

}
