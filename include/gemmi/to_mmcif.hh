// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hh.

#ifndef GEMMI_TO_MMCIF_HH_
#define GEMMI_TO_MMCIF_HH_

#include <string>
#include "cif.hh"
#include "model.hh"

namespace gemmi {
namespace mol {

inline void update_block(const Structure& st, cif::Block& block) {
  if (st.models.empty())
    return;
  block.name = st.name;
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
  block.update_value("_entry.id", id);

  // unit cell and symmetry
  block.update_value("_cell.entry_id", id);
  block.update_value("_cell.length_a", std::to_string(st.cell.a));
  block.update_value("_cell.length_b", std::to_string(st.cell.b));
  block.update_value("_cell.length_c", std::to_string(st.cell.c));
  block.update_value("_cell.angle_alpha", std::to_string(st.cell.alpha));
  block.update_value("_cell.angle_beta",  std::to_string(st.cell.beta));
  block.update_value("_cell.angle_gamma", std::to_string(st.cell.gamma));
  auto z_pdb = st.info.find("_cell.Z_PDB");
  if (z_pdb != st.info.end())
    block.update_value(z_pdb->first, z_pdb->second);
  block.update_value("_cell.angle_gamma", std::to_string(st.cell.gamma));
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
  for (auto ch : st.models[0].chains) {
    asym_loop.values.push_back(ch.name);
    asym_loop.values.push_back(ch.entity ? ch.entity->id : "?");
  }

  // matrices (scaling, NCS, etc)
  // TODO

  // atom list

  // aniso U
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
