// Copyright 2017 Global Phasing Ltd.
//
// mol::Structure -> cif::Document -> mmcif (PDBx/mmCIF) file

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include <string>
#include <utility>  // std::pair
#ifdef USE_STD_SNPRINTF
# include <cstdio>
# define stbsp_sprintf std::sprintf
#else
# include <stb_sprintf.h>
#endif
#include "cif.hpp"
#include "model.hpp"

namespace gemmi {
namespace mol {

inline std::string to_str(double d) {
  char buf[24];
  int len = stbsp_sprintf(buf, "%.9g", d);
  return std::string(buf, len > 0 ? len : 0);
}
inline std::string to_str(float d) {
  char buf[16];
  int len = stbsp_sprintf(buf, "%.6g", d);
  return std::string(buf, len > 0 ? len : 0);
}


void add_cif_atoms(const Structure& st, cif::Block& block) {
  // atom list
  cif::Loop& atom_loop = block.clear_or_add_loop("_atom_site.");
  atom_loop.tags = {cif::LoopTag("_atom_site.id"),
                    cif::LoopTag("_atom_site.type_symbol"),
                    cif::LoopTag("_atom_site.label_atom_id"),
                    cif::LoopTag("_atom_site.label_alt_id"),
                    cif::LoopTag("_atom_site.label_comp_id"),
                    cif::LoopTag("_atom_site.label_asym_id"),
                    cif::LoopTag("_atom_site.label_seq_id"),
                    cif::LoopTag("_atom_site.pdbx_PDB_ins_code"),
                    cif::LoopTag("_atom_site.Cartn_x"),
                    cif::LoopTag("_atom_site.Cartn_y"),
                    cif::LoopTag("_atom_site.Cartn_z"),
                    cif::LoopTag("_atom_site.occupancy"),
                    cif::LoopTag("_atom_site.B_iso_or_equiv"),
                    cif::LoopTag("_atom_site.pdbx_formal_charge"),
                    cif::LoopTag("_atom_site.auth_seq_id"),
                    cif::LoopTag("_atom_site.auth_asym_id"),
                    cif::LoopTag("_atom_site.pdbx_PDB_model_num")};
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());
  std::vector<std::pair<int, const Atom*>> aniso;
  int serial = 0;
  for (const Model& model : st.models) {
    for (const Chain& chain : model.chains) {
      for (const Residue& res : chain.residues) {
        std::string seq_id = std::to_string(res.seq_id);
        std::string auth_seq_id = std::to_string(res.auth_seq_id);
        std::string ins_code(1, res.ins_code ? res.ins_code : '?');
        for (const Atom& a : res.atoms) {
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.element.uname());
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc ? a.altloc : '.');
          vv.emplace_back(res.name);
          vv.emplace_back(chain.name);
          vv.emplace_back(seq_id);
          vv.emplace_back(ins_code);
          vv.emplace_back(to_str(a.pos.x));
          vv.emplace_back(to_str(a.pos.y));
          vv.emplace_back(to_str(a.pos.z));
          vv.emplace_back(to_str(a.occ));
          vv.emplace_back(to_str(a.b_iso));
          vv.emplace_back(a.charge == 0 ? "?" : std::to_string(a.charge));
          vv.emplace_back(auth_seq_id);
          vv.emplace_back(chain.auth_name);
          vv.emplace_back(model.name);
          if (a.u11 != 0.f)
            aniso.emplace_back(serial, &a);
        }
      }
    }
  }
  if (aniso.empty()) {
    block.delete_loop("_atom_site_anisotrop.id");
  } else {
    cif::Loop& aniso_loop = block.clear_or_add_loop("_atom_site_anisotrop.");
    aniso_loop.tags = {cif::LoopTag("_atom_site_anisotrop.id"),
                       cif::LoopTag("_atom_site_anisotrop.U[1][1]"),
                       cif::LoopTag("_atom_site_anisotrop.U[2][2]"),
                       cif::LoopTag("_atom_site_anisotrop.U[3][3]"),
                       cif::LoopTag("_atom_site_anisotrop.U[1][2]"),
                       cif::LoopTag("_atom_site_anisotrop.U[1][3]"),
                       cif::LoopTag("_atom_site_anisotrop.U[2][3]")};
    std::vector<std::string>& aniso_val = aniso_loop.values;
    aniso_val.reserve(aniso_loop.tags.size() * aniso.size());
    for (const auto& a : aniso) {
      aniso_val.emplace_back(std::to_string(a.first));
      aniso_val.emplace_back(to_str(a.second->u11));
      aniso_val.emplace_back(to_str(a.second->u22));
      aniso_val.emplace_back(to_str(a.second->u33));
      aniso_val.emplace_back(to_str(a.second->u12));
      aniso_val.emplace_back(to_str(a.second->u13));
      aniso_val.emplace_back(to_str(a.second->u23));
    }
  }
}

inline void update_cif_block(const Structure& st, cif::Block& block) {
  if (st.models.empty())
    return;
  block.name = st.name;
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
  block.update_value("_entry.id", id);
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end()) {
    block.update_value("_pdbx_database_status.entry_id", id);
    block.update_value(initial_date->first, initial_date->second);
  }

  // unit cell and symmetry
  block.update_value("_cell.entry_id", id);
  block.update_value("_cell.length_a",    to_str(st.cell.a));
  block.update_value("_cell.length_b",    to_str(st.cell.b));
  block.update_value("_cell.length_c",    to_str(st.cell.c));
  block.update_value("_cell.angle_alpha", to_str(st.cell.alpha));
  block.update_value("_cell.angle_beta",  to_str(st.cell.beta));
  block.update_value("_cell.angle_gamma", to_str(st.cell.gamma));
  auto z_pdb = st.info.find("_cell.Z_PDB");
  if (z_pdb != st.info.end())
    block.update_value(z_pdb->first, z_pdb->second);
  block.update_value("_symmetry.entry_id", id);
  block.update_value("_symmetry.space_group_name_H-M", cif::quote(st.sg_hm));

  // _entity
  cif::Loop& entity_loop = block.clear_or_add_loop("_entity.");
  entity_loop.tags = {cif::LoopTag("_entity.id"), cif::LoopTag("_entity.type")};
  for (const auto& ent : st.entities) {
    entity_loop.values.push_back(ent->id);
    entity_loop.values.push_back(ent->type_as_string());
  }

  // title, keywords, etc
  auto exptl_method = st.info.find("_exptl.method");
  if (exptl_method != st.info.end()) {
    block.update_value("_exptl.entry_id", id);
    block.update_value(exptl_method->first, cif::quote(exptl_method->second));
  }
  auto title = st.info.find("_struct.title");
  if (title != st.info.end()) {
    block.update_value("_struct.entry_id", id);
    block.update_value(title->first, cif::quote(title->second));
  }
  auto pdbx_keywords = st.info.find("_struct_keywords.pdbx_keywords");
  auto keywords = st.info.find("_struct_keywords.text");
  if (pdbx_keywords != st.info.end() || keywords != st.info.end())
    block.update_value("_struct_keywords.entry_id", id);
  if (pdbx_keywords != st.info.end())
    block.update_value(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
  if (keywords != st.info.end())
    block.update_value(keywords->first, cif::quote(keywords->second));

  // _struct_asym
  cif::Loop& asym_loop = block.clear_or_add_loop("_struct_asym.");
  asym_loop.tags.emplace_back("_struct_asym.id");
  asym_loop.tags.emplace_back("_struct_asym.entity_id");
  const std::vector<Chain>& asym_chains = st.get_chains();
  for (const auto& ch : asym_chains) {
    asym_loop.values.push_back(ch.name);
    asym_loop.values.push_back(ch.entity ? ch.entity->id : "?");
  }
  // SEQRES from PDB doesn't record microheterogeneity, so if the resulting
  // cif has unknown("?") _entity_poly_seq.num, it cannot be trusted.
  cif::Loop& poly_loop = block.clear_or_add_loop("_entity_poly_seq.");
  poly_loop.tags = {cif::LoopTag("_entity_poly_seq.entity_id"),
                    cif::LoopTag("_entity_poly_seq.num"),
                    cif::LoopTag("_entity_poly_seq.mon_id")};
  for (const auto& ent : st.entities)
    if (ent->type == EntityType::Polymer)
      for (const SequenceItem& si : ent->sequence) {
        poly_loop.values.emplace_back(ent->id);
        poly_loop.values.emplace_back(si.num >= 0 ? std::to_string(si.num)
                                                  : "?");
        poly_loop.values.emplace_back(si.mon);
      }

  // matrices (scaling, NCS, etc)
  // TODO

  add_cif_atoms(st, block);
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
