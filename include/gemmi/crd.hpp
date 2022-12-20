// Copyright 2022 Global Phasing Ltd.
//
// Generate Refmac intermediate (prepared) files crd and rst

#ifndef GEMMI_CRD_HPP_
#define GEMMI_CRD_HPP_

#include <cmath>         // for sqrt
#include "topo.hpp"      // for Topo
#include "placeh.hpp"    // for prepare_topology
#include "enumstr.hpp"   // for entity_type_to_string
#include "to_mmcif.hpp"  // for write_struct_conn
#include "sprintf.hpp"   // for to_str, to_str_prec
#include "calculate.hpp" // for find_best_plane
#include "select.hpp"    // for count_atom_sites
#include "to_chemcomp.hpp" // for add_chemcomp_to_block
#include "contact.hpp"   // for ContactSearch
#include "version.hpp"   // for GEMMI_VERSION

namespace gemmi {

inline bool has_anisou(const Model& model) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& a : res.atoms)
        if (a.aniso.nonzero())
          return true;
  return false;
}

inline std::string refmac_calc_flag(const Atom& a) {
  switch (a.calc_flag) {
    case CalcFlag::NotSet: return ".";
    // Refmac seems to be using only flags R and M
    case CalcFlag::Calculated: return a.is_hydrogen() ? "R" : "M";
    // I think these are not used
    case CalcFlag::Determined: return "d";
    case CalcFlag::Dummy: return "dum";
  }
  unreachable();
}

inline void add_automatic_links(Model& model, Structure& st, const MonLib& monlib) {
  auto is_onsb = [](Element e) {
    return e == El::O || e == El::N || e == El::S || e == El::B;
  };
  NeighborSearch ns(model, st.cell, 5.0);
  ns.populate();
  ContactSearch contacts(3.1f);  // 3.1 > 130% of ZN-CYS bond (2.34)
  contacts.ignore = ContactSearch::Ignore::AdjacentResidues;
  int counter = 0;
  contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                    int image_idx, float dist_sq) {
    if (st.find_connection_by_cra(cra1, cra2))
      return;
    const ChemLink* link;
    bool invert;
    char altloc = cra1.atom->altloc_or(cra2.atom->altloc);
    double min_dist_sq = sq(1 / 1.4) * dist_sq;
    std::tie(link, invert) = monlib.match_link(*cra1.residue, cra1.atom->name,
                                               *cra2.residue, cra2.atom->name,
                                               altloc, min_dist_sq);
    if (!link) {
      // Similarly to "make link" in Refmac,
      // only search for links between metals and O,N,S,B.
      if (is_onsb(cra1.atom->element) && cra2.atom->element.is_metal())
        invert = false;
      else if (is_onsb(cra2.atom->element) && cra1.atom->element.is_metal())
        invert = true;
      else
        return;
      float rmax = std::max(cra1.atom->element.covalent_r(),
                            cra2.atom->element.covalent_r());
      if (dist_sq > sq(std::max(2.f, 1.3f * rmax)))
        return;
    }

    Connection conn;
    conn.name = "added" + std::to_string(++counter);
    if (link)
      conn.link_id = link->id;
    conn.type = Connection::Covale;
    conn.asu = (image_idx == 0 ? Asu::Same : Asu::Different);
    const CRA* c1 = &cra1;
    const CRA* c2 = &cra2;
    if (invert)
      std::swap(c1, c2);
    conn.partner1 = make_address(*c1->chain, *c1->residue, *c1->atom);
    conn.partner2 = make_address(*c2->chain, *c2->residue, *c2->atom);
    conn.reported_distance = std::sqrt(dist_sq);
    st.connections.push_back(conn);
  });
}


inline cif::Block prepare_crd(const Structure& st, const Topo& topo,
                              HydrogenChange h_change, const std::string& info_comment) {
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
  cif::Block block("structure_" + id);
  auto& items = block.items;

  if (!info_comment.empty())
    items.emplace_back(cif::CommentArg{info_comment});
  items.emplace_back("_entry.id", id);
  items.emplace_back("_database_2.code_PDB", id);
  auto keywords = st.info.find("_struct_keywords.pdbx_keywords");
  if (keywords != st.info.end())
    items.emplace_back("_struct_keywords.text", cif::quote(keywords->second));
  auto title = st.info.find("_struct.title");
  if (title != st.info.end())
    items.emplace_back(title->first, cif::quote(title->second));
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end())
    items.emplace_back("_audit.creation_date", initial_date->second);
  items.emplace_back("_software.name", "gemmi");

  // this corresponds to Refmac keyword "make hydr"
  const char* hydr = "A";  // appropriate for ReAdd*
  if (h_change == HydrogenChange::NoChange || h_change == HydrogenChange::Shift)
    hydr = "Y";
  else if (h_change == HydrogenChange::Remove)
    hydr = "N";
  items.emplace_back("_ccp4_refmac.hatom", hydr);

  items.emplace_back(cif::CommentArg{"############\n"
                                     "## ENTITY ##\n"
                                     "############"});
  cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
  for (const Entity& ent : st.entities)
    entity_loop.add_row({ent.name, entity_type_to_string(ent.entity_type)});

  std::string mod_info = "\n#### Applied modifications ####\n";
  for (const Topo::ChainInfo& ci : topo.chain_infos) {
    cat_to(mod_info, "# chain ", ci.chain_ref.name,
           " / ", ci.subchain_name, ", entity ", ci.entity_id);
    if (is_polypeptide(ci.polymer_type))
      mod_info += " (polypeptide)";
    if (is_polynucleotide(ci.polymer_type))
      mod_info += " (polynucleotide)";
    mod_info += '\n';
    for (const Topo::ResInfo& ri : ci.res_infos) {
      cat_to(mod_info, "#    ", ri.res->seqid.str(), ' ', ri.res->name, ':');
      if (ri.mods.empty())
        mod_info += " n/a";
      for (const Topo::Mod& mod : ri.mods) {
        cat_to(mod_info, ' ', mod.id);
        if (mod.altloc)
          cat_to(mod_info, ':', mod.altloc);
        if (mod.alias != ChemComp::Group::Null)
          cat_to(mod_info, "(alias ", ChemComp::group_str(mod.alias), ')');
      }
      mod_info += '\n';
    }
  }
  items.emplace_back(cif::CommentArg{mod_info});

  items.emplace_back(cif::CommentArg{"##########\n"
                                     "## CELL ##\n"
                                     "##########"});
  items.emplace_back("_cell.entry_id", id);
  items.emplace_back("_cell.length_a",    to_str(st.cell.a));
  items.emplace_back("_cell.length_b",    to_str(st.cell.b));
  items.emplace_back("_cell.length_c",    to_str(st.cell.c));
  items.emplace_back("_cell.angle_alpha", to_str(st.cell.alpha));
  items.emplace_back("_cell.angle_beta",  to_str(st.cell.beta));
  items.emplace_back("_cell.angle_gamma", to_str(st.cell.gamma));
  bool write_fract_matrix = true;
  if (write_fract_matrix) {
    items.emplace_back(cif::CommentArg{"##############################\n"
                                       "## FRACTIONALISATION MATRIX ##\n"
                                       "##############################"});
    std::string prefix = "_atom_sites.fract_transf_";
    for (int i = 0; i < 3; ++i) {
      std::string start = prefix + "matrix[" + std::to_string(i+1) + "][";
      const auto& row = st.cell.frac.mat[i];
      for (int j = 0; j < 3; ++j)
        items.emplace_back(start + std::to_string(j+1) + "]", to_str(row[j]));
    }
    for (int i = 0; i < 3; ++i)
      items.emplace_back(prefix + "vector[" + std::to_string(i + 1) + "]",
                         to_str(st.cell.frac.vec.at(i)));
  }

  items.emplace_back(cif::CommentArg{"##############\n"
                                     "## SYMMETRY ##\n"
                                     "##############"});
  items.emplace_back("_symmetry.entry_id", id);
  const std::string& hm = st.spacegroup_hm;
  items.emplace_back("_symmetry.space_group_name_H-M", cif::quote(hm));
  if (const SpaceGroup* sg = st.find_spacegroup())
    items.emplace_back("_symmetry.Int_Tables_number",
                       std::to_string(sg->number));
  const Model& model0 = st.first_model();
  items.emplace_back(cif::CommentArg{"#################\n"
                                     "## STRUCT_ASYM ##\n"
                                     "#################"});
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const Chain& chain : model0.chains)
    for (ConstResidueSpan sub : chain.subchains())
      if (!sub.subchain_id().empty()) {
        const Entity* ent = st.get_entity_of(sub);
        asym_loop.add_row({sub.subchain_id(), (ent ? ent->name : "?")});
      }
  if (!st.connections.empty()) {
    items.emplace_back(cif::CommentArg{"#################\n"
                                       "## STRUCT_CONN ##\n"
                                       "#################"});
    impl::write_struct_conn(st, block);
    // disable ptnrN_auth_asym_id - otherwise Refmac uses it
    // instead of ptnrN_label_asym_id
    auto it = block.items.end() - 2;
    assert(it->type == cif::ItemType::Loop && it->has_prefix("_struct_conn."));
    for (std::string& tag : it->loop.tags)
      if (ends_with(tag, "_auth_asym_id"))
        tag += "-disabled";
  }

  items.emplace_back(cif::CommentArg{"###############\n"
                                     "## ATOM_SITE ##\n"
                                     "###############"});
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
      "group_PDB",
      "id",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "auth_seq_id", // including insertion code (no pdbx_PDB_ins_code)
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "ccp4_deuterium_fraction",  // tags[11]
      "B_iso_or_equiv",
      "type_symbol",
      "calc_flag",
      "label_seg_id",
      "auth_atom_id",
      "label_chem_id"});
  bool write_anisou = has_anisou(model0);
  if (write_anisou)
    for (const char* idx : {"[1][1]", "[2][2]", "[3][3]",
                            "[1][2]", "[1][3]", "[2][3]"})
      atom_loop.tags.push_back(std::string("_atom_site.aniso_U") + idx);
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());


  for (const Topo::ChainInfo& chain_info : topo.chain_infos)
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      const Residue& res = *ri.res;
      std::string auth_seq_id = res.seqid.str();
      for (const Atom& a : res.atoms) {
        const ChemComp& cc = ri.get_final_chemcomp(a.altloc);
        const auto& cc_atom = cc.get_atom(a.name);
        vv.emplace_back("ATOM");
        vv.emplace_back(std::to_string(a.serial));
        vv.emplace_back(a.name);
        vv.emplace_back(1, a.altloc ? a.altloc : '.');
        vv.emplace_back(res.name);
        vv.emplace_back(cif::quote(res.subchain));
        vv.emplace_back(auth_seq_id);
        vv.emplace_back(to_str(a.pos.x));
        vv.emplace_back(to_str(a.pos.y));
        vv.emplace_back(to_str(a.pos.z));
        vv.emplace_back(to_str(a.occ));
        vv.emplace_back(st.has_d_fraction ? to_str(a.fraction) : "0");
        vv.emplace_back(to_str(a.b_iso));
        std::string type_symbol = cc_atom.el.uname();
        if (a.charge != 0) {
          if (a.charge > 0) type_symbol += '+';
          type_symbol += std::to_string(a.charge);
        }
        vv.emplace_back(type_symbol);
        vv.emplace_back(refmac_calc_flag(a));
        vv.emplace_back(1, '.'); // label_seg_id
        vv.emplace_back(a.name); // again
        vv.emplace_back(cc_atom.chem_type); // label_chem_id
        if (write_anisou) {
          if (a.aniso.nonzero()) {
            for (float u : {a.aniso.u11, a.aniso.u22, a.aniso.u33,
                            a.aniso.u12, a.aniso.u13, a.aniso.u23})
              vv.push_back(to_str(u));
          } else {
            vv.resize(vv.size() + 6, ".");
          }
        }
      }
    }
  return block;
}

template<int Prec>
std::string to_str_dot(double d) {
  static_assert(Prec >= 0 && Prec < 7, "unsupported precision");
  if (!std::isnan(d)) {
    char buf[16];
    int len = gstb_sprintf(buf, "%.*f", Prec, d);
    if (len > 0)
      return std::string(buf, len);
  }
  return ".";
}

inline void add_restraint_row(cif::Loop& restr_loop,
                              const char* record, int counter,
                              const std::string& label, const std::string& period,
                              std::initializer_list<const Atom*> atoms,
                              double value, double dev,
                              double value_nucleus, double dev_nucleus,
                              double obs) {
  // Ignore restraints with zero-occupancy atoms (NoZeroOccRestr)
  // Perhaps it could be done earlier, in Topo::apply_restraints().
  for (const Atom* a : atoms)
    if (a->occ == 0)
      return;

  auto& values = restr_loop.values;
  values.emplace_back(record);  // record
  values.emplace_back(std::to_string(counter));  // number
  values.emplace_back(label);  // label
  values.emplace_back(period);  // period
  for (const Atom* a : atoms)
    values.emplace_back(std::to_string(a->serial));  // atom_id_i
  for (size_t i = atoms.size(); i < 4; ++i)
    values.emplace_back(".");
  values.emplace_back(to_str_dot<4>(value));  // value
  values.emplace_back(to_str_dot<4>(dev));  // dev
  values.emplace_back(to_str_dot<4>(value_nucleus));  // value_nucleus
  values.emplace_back(to_str_dot<4>(dev_nucleus));  // dev_nucleus
  values.emplace_back(to_str_prec<3>(obs));  // val_obs
  std::string& last = values.back();
  last += " #";
  for (const Atom* a : atoms) {
    last += ' ';
    last += a->name;
    if (a->has_altloc()) {
      last += '.';
      last += a->altloc;
    }
  }
}

inline void add_restraints(const Topo::Rule rule, const Topo& topo,
                           cif::Loop& restr_loop, int (&counters)[6],
                           const UnitCell* cell=nullptr) {
  if (rule.rkind == Topo::RKind::Bond) {
    const Topo::Bond& t = topo.bonds[rule.index];
    if (cell == nullptr) {  // don't use symmetry
      add_restraint_row(restr_loop, "BOND", ++counters[0],
                        bond_type_to_string(t.restr->type), ".",
                        {t.atoms[0], t.atoms[1]},
                        t.restr->value, t.restr->esd,
                        t.restr->value_nucleus, t.restr->esd_nucleus,
                        t.calculate());
    } else {
      NearestImage im = cell->find_nearest_image(t.atoms[0]->pos, t.atoms[1]->pos,
                                                 Asu::Different);
      add_restraint_row(restr_loop, "BNDS", ++counters[5],
                        im.symmetry_code(true), ".",
                        {t.atoms[0], t.atoms[1]},
                        t.restr->value, t.restr->esd,
                        t.restr->value_nucleus, t.restr->esd_nucleus,
                        std::sqrt(im.dist_sq));
    }
  } else if (rule.rkind == Topo::RKind::Angle) {
    const Topo::Angle& t = topo.angles[rule.index];
    add_restraint_row(restr_loop, "ANGL", ++counters[1], ".", ".",
                      {t.atoms[0], t.atoms[1], t.atoms[2]},
                       t.restr->value, t.restr->esd, NAN, NAN,
                       deg(t.calculate()));
  } else if (rule.rkind == Topo::RKind::Torsion) {
    const Topo::Torsion& t = topo.torsions[rule.index];
    add_restraint_row(restr_loop, "TORS", ++counters[2],
                      t.restr->label, std::to_string(t.restr->period),
                      {t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]},
                      t.restr->value, t.restr->esd, NAN, NAN,
                      deg(t.calculate()));
  } else if (rule.rkind == Topo::RKind::Chirality) {
    const Topo::Chirality& t = topo.chirs[rule.index];
    add_restraint_row(restr_loop, "CHIR", ++counters[3],
                      chirality_to_string(t.restr->sign), ".",
                      {t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]},
                      topo.ideal_chiral_abs_volume(t), 0.02, NAN, NAN,
                      t.calculate());
  } else if (rule.rkind == Topo::RKind::Plane) {
    const Topo::Plane& t = topo.planes[rule.index];
    ++counters[4];
    auto coeff = find_best_plane(t.atoms);
    for (const Atom* atom : t.atoms)
      add_restraint_row(restr_loop, "PLAN", counters[4], t.restr->label, ".",
                        {atom},
                        t.restr->esd, NAN, NAN, NAN,
                        get_distance_from_plane(atom->pos, coeff));
  }
}

inline cif::Block prepare_rst(const Topo& topo, const MonLib& monlib, const UnitCell& cell) {
  cif::Block block("restraints");
  cif::Loop& restr_loop = block.init_mmcif_loop("_restr.", {
              "record", "number", "label", "period",
              "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
              "value", "dev", "value_nucleus", "dev_nucleus", "val_obs"});
  int counters[6] = {0, 0, 0, 0, 0, 0};
  for (const Topo::ChainInfo& chain_info : topo.chain_infos) {
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      // write link
      for (const Topo::Link& prev : ri.prev) {
        const ChemLink* link = monlib.get_link(prev.link_id);
        if (link && !prev.link_rules.empty()) {
          std::string comment = " link " + prev.link_id + " " +
                                prev.res1->seqid.str() + " " +
                                prev.res1->name + " - " +
                                ri.res->seqid.str() + " " + ri.res->name;
          restr_loop.add_comment_and_row({comment, "LINK", ".", cif::quote(prev.link_id), ".",
                                          ".", ".", ".", ".", ".", ".", ".", ".", "."});
          for (const Topo::Rule& rule : prev.link_rules)
            add_restraints(rule, topo, restr_loop, counters);
        }
      }
      // write monomer
      if (!ri.monomer_rules.empty()) {
        std::string res_info = " monomer " + chain_info.subchain_name + " " +
                               ri.res->seqid.str() + " " + ri.res->name;
        if (!ri.mods.empty()) {
          res_info += " modified by ";
          for (size_t i = 0; i != ri.mods.size(); ++i) {
            if (i != 0)
              res_info += ", ";
            res_info += ri.mods[i].id;
            if (ri.mods[i].altloc != '\0')
              cat_to(res_info, ':', ri.mods[i].altloc);
          }
        }

        std::string group_str = ".";
        if (ri.orig_chemcomp) {
          ChemComp::Group group = ri.orig_chemcomp->group;
          if (ChemComp::is_peptide_group(group))
            group_str = "L-peptid";  // we try to be compatible with Refmac
          else if (group != ChemComp::Group::NonPolymer)
            group_str = ChemComp::group_str(group);
        }

        restr_loop.add_comment_and_row({res_info, "MONO", ".", group_str, ".",
                                        ".", ".", ".", ".", ".", ".", ".", ".", "."});
        for (const Topo::Rule& rule : ri.monomer_rules)
          add_restraints(rule, topo, restr_loop, counters);
      }
    }
  }
  // explicit links
  for (const Topo::Link& extra : topo.extras) {
    if (extra.asu == Asu::Different) // symmetry links are left for later
      continue;
    const ChemLink* chem_link = monlib.get_link(extra.link_id);
    assert(chem_link);
    std::string comment = " link " + chem_link->id;
    restr_loop.add_comment_and_row({comment, "LINK", ".", cif::quote(chem_link->id), ".",
                                    ".", ".", ".", ".", ".", ".", ".", ".", "."});
    for (const Topo::Rule& rule : extra.link_rules)
      add_restraints(rule, topo, restr_loop, counters);
  }

  // special symmetry related links
  for (const Topo::Link& extra : topo.extras) {
    if (extra.asu != Asu::Different)
      continue;
    const ChemLink* chem_link = monlib.get_link(extra.link_id);
    assert(chem_link);
    std::string comment = " link (symmetry) " + chem_link->id;
    restr_loop.add_comment_and_row({comment, "LINK", ".", "symmetry", ".",
                                    ".", ".", ".", ".", ".", ".", ".", ".", "."});
    for (const Topo::Rule& rule : extra.link_rules)
      if (rule.rkind == Topo::RKind::Bond)
        add_restraints(rule, topo, restr_loop, counters, &cell);
  }

  return block;
}

inline cif::Document prepare_refmac_crd(const Structure& st, const Topo& topo,
                                        const MonLib& monlib, HydrogenChange h_change) {
  cif::Document doc;
  std::string info_comment = "# Refmac CRD file generated with gemmi " GEMMI_VERSION
                             "\n# Monomer library version: " + monlib.lib_version;
  doc.blocks.push_back(prepare_crd(st, topo, h_change, info_comment));
  doc.blocks.push_back(prepare_rst(topo, monlib, st.cell));

  doc.blocks.emplace_back("for_refmac_mmcif");
  std::vector<std::string> resnames = st.models.at(0).get_all_residue_names();
  for (const std::string& resname : resnames) {
    auto it = monlib.monomers.find(resname);
    if (it != monlib.monomers.end()) {
      const ChemComp& cc = it->second;
      doc.blocks.emplace_back(cc.name);
      cif::Block& block = doc.blocks.back();
      block.items.emplace_back("_chem_comp.id", cc.name);
      block.items.emplace_back("_chem_comp.group", ChemComp::group_str(cc.group));
      add_chemcomp_to_block(cc, block);
    }
  }

  // gather used links and mods
  std::vector<std::string> used_links;
  std::vector<std::string> used_mods;
  for (const Topo::ChainInfo& chain_info : topo.chain_infos)
    for (const Topo::ResInfo& res_info : chain_info.res_infos) {
      for (const Topo::Link& link : res_info.prev)
        if (!in_vector(link.link_id, used_links))
          used_links.push_back(link.link_id);
      for (const Topo::Mod& mod : res_info.mods)
        if (!in_vector(mod.id, used_mods))
          used_mods.push_back(mod.id);
    }
  for (const Topo::Link& extra : topo.extras)
    if (!in_vector(extra.link_id, used_links))
      used_links.push_back(extra.link_id);

  // add links and mods blocks to the document
  auto q = [](const std::string& s) { return s.empty() ? "?" : cif::quote(s); };
  for (const std::string& link_name : used_links) {
    const ChemLink* cl = monlib.get_link(link_name);
    // ignore ad-hoc links and dummy (empty) links such as "gap"
    if (cl && !starts_with(cl->name, "auto-") && !cl->block.items.empty()) {
      doc.blocks.push_back(cl->block);
      cif::Block& block = doc.blocks.back();
      block.init_mmcif_loop("_chem_link.", {"id", "name",
                                            "comp_id_1", "mod_id_1", "group_comp_1",
                                            "comp_id_2", "mod_id_2", "group_comp_2"})
        .add_row({q(cl->id), q(cl->name),
                  q(cl->side1.comp), q(cl->side1.mod), ChemComp::group_str(cl->side1.group),
                  q(cl->side2.comp), q(cl->side2.mod), ChemComp::group_str(cl->side2.group)});
    }
  }
  for (const std::string& mod_name : used_mods)
    if (const ChemMod* mod = monlib.get_mod(mod_name)) {
      doc.blocks.push_back(mod->block);
      cif::Block& block = doc.blocks.back();
      block.init_mmcif_loop("_chem_mod.", {"id", "name", "comp_id", "group_id"})
        .add_row({q(mod->id), q(mod->name), q(mod->comp_id), q(mod->group_id)});
    }

  return doc;
}

} // namespace gemmi
#endif
