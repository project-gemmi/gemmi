// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hpp.

#ifndef GEMMI_MMCIF_HPP_
#define GEMMI_MMCIF_HPP_

#include <array>
#include <string>
#include <unordered_map>
#include "cifdoc.hpp"
#include "numb.hpp"
#include "model.hpp"

namespace gemmi {

namespace impl {

inline std::unordered_map<std::string, std::array<float,6>>
get_anisotropic_u(const cif::Block& block) {
  cif::TableView aniso_tab = block.find("_atom_site_anisotrop.",
                                        {"id", "U[1][1]", "U[2][2]", "U[3][3]",
                                         "U[1][2]", "U[1][3]", "U[2][3]"});
  std::unordered_map<std::string, std::array<float,6>> aniso_map;
  for (auto ani : aniso_tab)
    aniso_map[ani[0]] = {{(float) cif::as_number(ani[1]),
                          (float) cif::as_number(ani[2]),
                          (float) cif::as_number(ani[3]),
                          (float) cif::as_number(ani[4]),
                          (float) cif::as_number(ani[5]),
                          (float) cif::as_number(ani[6])}};
  return aniso_map;
}

inline cif::TableView find_transform(const cif::Block& block,
                                     std::string category,
                                     std::string mstr="matrix",
                                     std::string vstr="vector") {
  return block.find(category, {
      mstr + "[1][1]", mstr + "[2][1]", mstr + "[3][1]",
      mstr + "[1][2]", mstr + "[2][2]", mstr + "[3][2]",
      mstr + "[1][3]", mstr + "[2][3]", mstr + "[3][3]",
      vstr + "[1]",    vstr + "[2]",    vstr + "[3]"});
}

inline Mat4x4 get_transform_matrix(const cif::TableView::Row& r) {
  auto num = [&r](int n) { return cif::as_number(r[n]); };
  return {{num(0), num(1),  num(2),  0},
          {num(3), num(4),  num(5),  0},
          {num(6), num(7),  num(8),  0},
          {num(9), num(10), num(11), 1}};
}

inline Mat4x4 Matrix33_to_Mat4x4(const Matrix33& m) {
  return {{m.a11, m.a21, m.a31, 0},
          {m.a12, m.a22, m.a32, 0},
          {m.a13, m.a23, m.a33, 0},
          {    0,     0,     0, 1}};
}

inline Structure structure_from_cif_block(const cif::Block& block) {
  using cif::as_number;
  using cif::as_string;
  Structure st;
  st.name = block.name;

  // unit cell and symmetry
  cif::TableView cell = block.find("_cell.",
                                   {"length_a", "length_b", "length_c",
                                   "angle_alpha", "angle_beta", "angle_gamma"});
  if (cell.ok()) {
    auto c = cell.one();
    st.cell.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
                as_number(c[3]), as_number(c[4]), as_number(c[5]));
  }
  st.sg_hm = as_string(block.find_value("_symmetry.space_group_name_H-M"));

  auto add_info = [&](std::string tag) {
    cif::TableView t = block.find(tag);
    if (t.length() >= 1) {
      bool first = true;
      for (const cif::TableView::Row &r : t)
        if (!cif::is_null(r[0])) {
          if (first)
            st.info[tag] = as_string(r[0]);
          else
            st.info[tag] += "," + as_string(r[0]);
          first = false;
        }
    }
  };
  add_info("_entry.id");
  add_info("_cell.Z_PDB");
  add_info("_exptl.method");
  add_info("_struct.title");
  // in pdbx/mmcif v5 date_original was replaced with a much longer tag
  std::string old_date_tag = "_database_PDB_rev.date_original";
  std::string new_date_tag
                      = "_pdbx_database_status.recvd_initial_deposition_date";
  add_info(old_date_tag);
  add_info(new_date_tag);
  if (st.info.count(old_date_tag) == 1 && st.info.count(new_date_tag) == 0)
    st.info[new_date_tag] = st.info[old_date_tag];
  add_info("_struct_keywords.pdbx_keywords");
  add_info("_struct_keywords.text");

  cif::TableView ncs_oper = find_transform(block, "_struct_ncs_oper.");
  int ncs_oper_id = block.add_field(ncs_oper, "_struct_ncs_oper.id");
  int ncs_code_idx = block.add_field(ncs_oper, "_struct_ncs_oper.code");
  for (auto op : ncs_oper) {
    bool given = (ncs_code_idx > 0 && op.str(ncs_code_idx) == "given");
    std::string id;
    if (ncs_oper_id)
      id = op.str(ncs_oper_id);
    Mat4x4 mat = get_transform_matrix(op);
    if (mat != Mat4x4(linalg::identity))
      st.ncs.push_back({id, given, mat});
  }

  // PDBx/mmcif spec defines both _database_PDB_matrix.scale* and
  // _atom_sites.fract_transf_* as equivalent of pdb SCALE, but the former
  // is not used, so we ignore it.
  cif::TableView fract_tv = find_transform(block, "_atom_sites.fract_transf_");
  if (fract_tv.length() > 0) {
    Mat4x4 fract = get_transform_matrix(fract_tv[0]);
    st.cell.set_matrices_from_fract(fract);
  }

  // We store origx just for completeness. It may never be useful
  // for anything but writing it back to a file.
  cif::TableView origx_tv = find_transform(block, "_database_PDB_matrix.origx",
                                           "", "_vector");
  if (origx_tv.length() > 0)
    st.origx = get_transform_matrix(origx_tv[0]);

  auto aniso_map = get_anisotropic_u(block);

  // atom list
  enum { kId=0, kSymbol, kAtomId, kAltId, kCompId, kAsymId, kSeqId, kInsCode,
         kX, kY, kZ, kOcc, kBiso, kCharge, kAuthSeqId, kAuthAsymId, kModelNum };
  cif::TableView atom_table = block.find("_atom_site.",
                                         {"id",
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
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  for (auto row : atom_table) {
    if (!model || row[kModelNum] != model->name) {
      model = st.find_or_add_model(row[kModelNum]);
      chain = nullptr;
    }
    if (!chain || row[kAsymId] != chain->name) {
      chain = model->find_or_add_chain(row.str(kAsymId));
      chain->auth_name = row.str(kAuthAsymId);
      resi = nullptr;
    }
    ResidueId rid(cif::as_int(row[kSeqId], Residue::UnknownId),
                  cif::as_int(row[kAuthSeqId], Residue::UnknownId),
                  as_string(row[kInsCode])[0],
                  as_string(row[kCompId]));
    if (!resi || !resi->matches(rid)) {
      // the insertion code happens to be always a single letter
      assert(row[kInsCode].size() == 1);
      resi = chain->find_or_add_residue(rid);
    } else {
      assert(resi->auth_seq_id == rid.auth_seq_id);
      assert(resi->ins_code == rid.ins_code);
    }
    Atom atom;
    atom.name = as_string(row[kAtomId]);
    atom.group = 0;
    atom.altloc = as_string(row[kAltId])[0];
    atom.charge = cif::is_null(row[kCharge]) ? 0 : cif::as_int(row[kCharge]);
    atom.element = gemmi::Element(as_string(row[kSymbol]));
    atom.pos.x = cif::as_number(row[kX]);
    atom.pos.y = cif::as_number(row[kY]);
    atom.pos.z = cif::as_number(row[kZ]);
    atom.occ = (float) cif::as_number(row[kOcc], 1.0);
    atom.b_iso = (float) cif::as_number(row[kBiso], 50.0);

    if (!aniso_map.empty()) {
      auto ani = aniso_map.find(row[kId]);
      if (ani != aniso_map.end()) {
        atom.u11 = ani->second[0];
        atom.u22 = ani->second[1];
        atom.u33 = ani->second[2];
        atom.u12 = ani->second[3];
        atom.u13 = ani->second[4];
        atom.u23 = ani->second[5];
      }
    }
    resi->atoms.emplace_back(atom);
  }

  cif::TableView polymer_types = block.find("_entity_poly.",
                                            {"entity_id", "type"});
  for (const auto& row : block.find("_entity.", {"id", "type"})) {
    std::string id = row.str(0);
    EntityType etype = entity_type_from_string(row.str(1));
    PolymerType ptype = PolymerType::NA;
    try {
      ptype = polymer_type_from_string(polymer_types.find_row(id).str(1));
    } catch (std::runtime_error&) {}
    st.entities.emplace_back(new Entity{id, etype, ptype, {}});
  }

  for (const auto& row : block.find("_entity_poly_seq.",
                                    {"entity_id", "num", "mon_id"})) {
    Entity *ent = st.find_or_add_entity(row.str(0));
    ent->sequence.push_back({cif::as_int(row[1], -1), row.str(2)});
  }

  auto chain_to_entity = block.find("_struct_asym.", {"id", "entity_id"});
  for (Model& mod : st.models)
    for (Chain& ch : mod.chains)
      try {
        std::string ent_id = chain_to_entity.find_row(ch.name).str(1);
        ch.entity = st.find_or_add_entity(ent_id);
      } catch (std::runtime_error&) {  // maybe _struct_asym is missing
        ch.entity = nullptr;
      }
  st.finish();
  return st;
}

} // namespace impl

inline Structure read_atoms(const cif::Document& doc) {
  return impl::structure_from_cif_block(doc.sole_block());
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
