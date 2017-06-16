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
  cif::Loop& atom_loop = block.clear_or_add_loop("_atom_site.", {
      "id",
      "type_symbol",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "label_seq_id",
      "pdbx_PDB_ins_code",
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "B_iso_or_equiv",
      "pdbx_formal_charge",
      "auth_seq_id",
      "auth_asym_id",
      "pdbx_PDB_model_num"});
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
    cif::Loop& aniso_loop = block.clear_or_add_loop("_atom_site_anisotrop.", {
                                    "id", "U[1][1]", "U[2][2]", "U[3][3]",
                                    "U[1][2]", "U[1][3]", "U[2][3]"});
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
  cif::Loop& entity_loop = block.clear_or_add_loop("_entity.", {"id", "type"});
  for (const auto& ent : st.entities)
    entity_loop.append_row({ent->id, ent->type_as_string()});

  // _entity_poly
  cif::Loop& entity_poly_loop = block.clear_or_add_loop("_entity_poly.",
                                                        {"entity_id", "type"});
  for (const auto& ent : st.entities)
    if (ent->type == EntityType::Polymer)
      entity_poly_loop.append_row({ent->id, ent->polymer_type_as_string()});

  // _exptl
  cif::Loop& exptl_method_loop = block.clear_or_add_loop("_exptl.",
                                                        {"entry_id", "method"});
  auto exptl_method = st.info.find("_exptl.method");
  if (exptl_method != st.info.end())
    for (const std::string& m : gemmi::split_str(exptl_method->second, ','))
      exptl_method_loop.append_row({id, cif::quote(m)});

  // title, keywords
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

  // _struct_ncs_oper (MTRIX)
  if (!st.ncs.empty()) {
    cif::Loop& ncs_oper = block.clear_or_add_loop("_struct_ncs_oper.",
        {"id", "code",
         "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
         "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
         "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
    for (const NcsOp& op : st.ncs) {
      ncs_oper.values.emplace_back(op.id);
      ncs_oper.values.emplace_back(op.given ? "given" : "generate");
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
          ncs_oper.values.emplace_back(to_str(op.transform[j][i]));
    }
  }

  // _struct_asym
  cif::Loop& asym_loop = block.clear_or_add_loop("_struct_asym.",
                                                 {"id", "entity_id"});
  for (const auto& ch : st.get_chains())
    asym_loop.append_row({ch.name, (ch.entity ? ch.entity->id : "?")});

  // _database_PDB_matrix (ORIGX)
  if (st.origx != Mat4x4(linalg::identity)) {
    block.update_value("_database_PDB_matrix.entry_id", id);
    std::string prefix = "_database_PDB_matrix.origx";
    for (int i = 0; i < 3; ++i) {
      std::string s = "[" + std::to_string(i+1) + "]";
      block.update_value(prefix + s + "[1]", to_str(st.origx.x[i]));
      block.update_value(prefix + s + "[2]", to_str(st.origx.y[i]));
      block.update_value(prefix + s + "[3]", to_str(st.origx.z[i]));
      block.update_value(prefix + "_vector" + s, to_str(st.origx.w[i]));
    }
  }

  // _atom_sites (SCALE)
  if (st.cell.explicit_matrices) {
    block.update_value("_atom_sites.entry_id", id);
    std::string prefix = "_atom_sites.fract_transf_";
    block.update_value(prefix + "matrix[1][1]", to_str(st.cell.frac.a11));
    block.update_value(prefix + "matrix[1][2]", to_str(st.cell.frac.a12));
    block.update_value(prefix + "matrix[1][3]", to_str(st.cell.frac.a13));
    block.update_value(prefix + "matrix[2][1]", to_str(st.cell.frac.a21));
    block.update_value(prefix + "matrix[2][2]", to_str(st.cell.frac.a22));
    block.update_value(prefix + "matrix[2][3]", to_str(st.cell.frac.a23));
    block.update_value(prefix + "matrix[3][1]", to_str(st.cell.frac.a31));
    block.update_value(prefix + "matrix[3][2]", to_str(st.cell.frac.a32));
    block.update_value(prefix + "matrix[3][3]", to_str(st.cell.frac.a33));
    block.update_value(prefix + "vector[1]",    to_str(st.cell.shift.x));
    block.update_value(prefix + "vector[2]",    to_str(st.cell.shift.y));
    block.update_value(prefix + "vector[3]",    to_str(st.cell.shift.z));
  }

  // SEQRES from PDB doesn't record microheterogeneity, so if the resulting
  // cif has unknown("?") _entity_poly_seq.num, it cannot be trusted.
  cif::Loop& poly_loop = block.clear_or_add_loop("_entity_poly_seq.", {
                                      "entity_id", "num", "mon_id"});
  for (const auto& ent : st.entities)
    if (ent->type == EntityType::Polymer)
      for (const SequenceItem& si : ent->sequence) {
        poly_loop.append_row({ent->id,
                              (si.num >= 0 ? std::to_string(si.num) : "?"),
                              si.mon});
      }

  add_cif_atoms(st, block);
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
