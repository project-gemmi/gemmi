// Copyright 2017 Global Phasing Ltd.
//
// Structure -> cif::Document -> mmcif (PDBx/mmCIF) file

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include <string>
#include <utility>  // std::pair
#include "sprintf.hpp"
#include "cifdoc.hpp"
#include "model.hpp"
#include "calculate.hpp"  // for count_atom_sites

namespace gemmi {

namespace impl {

inline std::string pdbx_icode(const ResidueId& rid) {
  return std::string(1, rid.has_icode() ? rid.icode : '?');
}

inline void add_cif_atoms(const Structure& st, cif::Block& block) {
  // atom list
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
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
        std::string label_seq_id = res.label_seq.str();
        std::string auth_seq_id = res.seq_num.str();
        for (const Atom& a : res.atoms) {
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.element.uname());
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc ? a.altloc : '.');
          vv.emplace_back(res.name);
          vv.emplace_back(chain.name);
          vv.emplace_back(label_seq_id);
          vv.emplace_back(pdbx_icode(res));
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
    block.find_mmcif_category("_atom_site_anisotrop.").erase();
  } else {
    cif::Loop& aniso_loop = block.init_mmcif_loop("_atom_site_anisotrop.", {
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

inline void write_struct_conn(const Structure& st, cif::Block& block) {
  // example:
  // disulf1 disulf A CYS 3  SG ? 3 ? 1_555 A CYS 18 SG ? 18 ?  1_555 ? 2.045
  cif::Loop& conn_loop = block.init_mmcif_loop("_struct_conn.",
      {"id", "conn_type_id",
       "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id",
       "ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id",
       "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_symmetry",
       "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id",
       "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id",
       "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code", "ptnr2_symmetry",
       "details", "pdbx_dist_value"});
  for (const Connection& con : st.models.at(0).connections) {
    const_CRA cra1 = st.models[0].find_cra(con.atom[0]);
    const_CRA cra2 = st.models[0].find_cra(con.atom[1]);
    if (!cra1.atom || !cra2.atom)
      continue;
    SymImage im = st.cell.find_nearest_image(cra1.atom->pos,
                                             cra2.atom->pos, con.asu);
    conn_loop.add_row({
        con.name,                               // id
        get_mmcif_connection_type_id(con.type), // conn_type_id
        cra1.chain->name,                       // ptnr1_label_asym_id
        con.atom[0].res_id.name,                // ptnr1_label_comp_id
        ".",                                    // ptnr1_label_seq_id
        con.atom[0].atom_name,                  // ptnr1_label_atom_id
        std::string(1, con.atom[0].altloc ? con.atom[0].altloc : '?'),
        con.atom[0].res_id.seq_num.str(),       // ptnr1_auth_seq_id
        pdbx_icode(con.atom[0].res_id),         // ptnr1_PDB_ins_code
        "1_555",                                // ptnr1_symmetry
        cra2.chain->name,                       // ptnr2_label_asym_id
        con.atom[1].res_id.name,                // ptnr2_label_comp_id
        ".",                                    // ptnr2_label_seq_id
        con.atom[1].atom_name,                  // ptnr2_label_atom_id
        std::string(1, con.atom[1].altloc ? con.atom[1].altloc : '?'),
        con.atom[1].res_id.seq_num.str(),       // ptnr2_auth_seq_id
        pdbx_icode(con.atom[1].res_id),         // ptnr2_PDB_ins_code
        im.pdb_symbol(true),                    // ptnr2_symmetry
        "?",                                    // details
        to_str(im.dist())                       // pdbx_dist_value
    });
  }
}

} // namespace impl

inline void update_cif_block(const Structure& st, cif::Block& block) {
  using std::to_string;
  if (st.models.empty())
    return;
  block.name = st.name;
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
  block.set_pair("_entry.id", id);
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end()) {
    block.set_pair("_pdbx_database_status.entry_id", id);
    block.set_pair(initial_date->first, initial_date->second);
  }

  // unit cell and symmetry
  block.set_pair("_cell.entry_id", id);
  block.set_pair("_cell.length_a",    to_str(st.cell.a));
  block.set_pair("_cell.length_b",    to_str(st.cell.b));
  block.set_pair("_cell.length_c",    to_str(st.cell.c));
  block.set_pair("_cell.angle_alpha", to_str(st.cell.alpha));
  block.set_pair("_cell.angle_beta",  to_str(st.cell.beta));
  block.set_pair("_cell.angle_gamma", to_str(st.cell.gamma));
  auto z_pdb = st.info.find("_cell.Z_PDB");
  if (z_pdb != st.info.end())
    block.set_pair(z_pdb->first, z_pdb->second);
  block.set_pair("_symmetry.entry_id", id);
  block.set_pair("_symmetry.space_group_name_H-M",
                 cif::quote(st.spacegroup_hm));
  if (const SpaceGroup* sg = find_spacegroup_by_name(st.spacegroup_hm))
    block.set_pair("_symmetry.Int_Tables_number", to_string(sg->number));

  // _entity
  cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
  for (const auto& ent : st.entities)
    entity_loop.add_row({ent.first, ent.second.type_as_string()});

  // _entity_poly
  cif::Loop& ent_poly_loop = block.init_mmcif_loop("_entity_poly.",
                                                   {"entity_id", "type"});
  for (const auto& ent : st.entities)
    if (ent.second.entity_type == EntityType::Polymer)
      ent_poly_loop.add_row({ent.first, ent.second.polymer_type_as_string()});

  // _exptl
  cif::Loop& exptl_method_loop = block.init_mmcif_loop("_exptl.",
                                                       {"entry_id", "method"});
  auto exptl_method = st.info.find("_exptl.method");
  if (exptl_method != st.info.end())
    for (const std::string& m : gemmi::split_str(exptl_method->second, "; "))
      exptl_method_loop.add_row({id, cif::quote(m)});

  // title, keywords
  auto title = st.info.find("_struct.title");
  if (title != st.info.end()) {
    block.set_pair("_struct.entry_id", id);
    block.set_pair(title->first, cif::quote(title->second));
  }
  auto pdbx_keywords = st.info.find("_struct_keywords.pdbx_keywords");
  auto keywords = st.info.find("_struct_keywords.text");
  if (pdbx_keywords != st.info.end() || keywords != st.info.end())
    block.set_pair("_struct_keywords.entry_id", id);
  if (pdbx_keywords != st.info.end())
    block.set_pair(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
  if (keywords != st.info.end())
    block.set_pair(keywords->first, cif::quote(keywords->second));

  // _struct_ncs_oper (MTRIX)
  if (!st.ncs.empty()) {
    cif::Loop& ncs_oper = block.init_mmcif_loop("_struct_ncs_oper.",
        {"id", "code",
         "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
         "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
         "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
    for (const NcsOp& op : st.ncs) {
      ncs_oper.values.emplace_back(op.id);
      ncs_oper.values.emplace_back(op.given ? "given" : "generate");
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
        ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
      }
    }
  }

  // _struct_asym
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const auto& ch : st.models.at(0).chains)
    asym_loop.add_row({ch.name, (ch.entity_id.empty() ? "?" : ch.entity_id)});

  // _database_PDB_matrix (ORIGX)
  if (st.has_origx) {
    block.set_pair("_database_PDB_matrix.entry_id", id);
    std::string prefix = "_database_PDB_matrix.origx";
    for (int i = 0; i < 3; ++i) {
      std::string s = "[" + to_string(i+1) + "]";
      block.set_pair(prefix + s + "[1]", to_str(st.origx.mat[i][0]));
      block.set_pair(prefix + s + "[2]", to_str(st.origx.mat[i][1]));
      block.set_pair(prefix + s + "[3]", to_str(st.origx.mat[i][2]));
      block.set_pair(prefix + "_vector" + s, to_str(st.origx.vec.at(i)));
    }
  }

  impl::write_struct_conn(st, block);

  // _struct_mon_prot_cis
  cif::Loop& prot_cis_loop = block.init_mmcif_loop("_struct_mon_prot_cis.",
      {"pdbx_id", "pdbx_PDB_model_num", "label_asym_id", "label_seq_id",
       "auth_seq_id", "pdbx_PDB_ins_code",
       "label_comp_id", "auth_comp_id", "label_alt_id"});
  for (const Model& model : st.models)
    for (const Chain& chain : model.chains)
        for (const Residue& res : chain.residues)
          if (res.is_cis)
            prot_cis_loop.add_row({to_string(prot_cis_loop.length()+1),
                                   model.name, chain.name, res.label_seq.str(),
                                   res.seq_num.str(), impl::pdbx_icode(res),
                                   res.name, res.name, "."});

  // _atom_sites (SCALE)
  if (st.has_origx || st.cell.explicit_matrices) {
    block.set_pair("_atom_sites.entry_id", id);
    std::string prefix = "_atom_sites.fract_transf_";
    for (int i = 0; i < 3; ++i) {
      std::string idx = "[" + std::to_string(i + 1) + "]";
      const auto& frac = st.cell.frac;
      block.set_pair(prefix + "matrix" + idx + "[1]", to_str(frac.mat[i][0]));
      block.set_pair(prefix + "matrix" + idx + "[2]", to_str(frac.mat[i][1]));
      block.set_pair(prefix + "matrix" + idx + "[3]", to_str(frac.mat[i][2]));
      block.set_pair(prefix + "vector" + idx, to_str(frac.vec.at(i)));
    }
  }

  // SEQRES from PDB doesn't record microheterogeneity, so if the resulting
  // cif has unknown("?") _entity_poly_seq.num, it cannot be trusted.
  cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.",
                                         {"entity_id", "num", "mon_id"});
  for (const auto& ent : st.entities)
    if (ent.second.entity_type == EntityType::Polymer)
      for (const SequenceItem& si : ent.second.sequence) {
        std::string num = si.num >= 0 ? to_string(si.num) : "?";
        poly_loop.add_row({ent.first, num, si.mon});
      }

  impl::add_cif_atoms(st, block);
}

inline cif::Document make_mmcif_document(const Structure& st) {
  cif::Document doc;
  doc.blocks.resize(1);
  gemmi::update_cif_block(st, doc.blocks[0]);
  return doc;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
