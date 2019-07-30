// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hpp.

#ifndef GEMMI_MMCIF_HPP_
#define GEMMI_MMCIF_HPP_

#include <array>
#include <string>
#include <unordered_map>
#include "cifdoc.hpp"
#include "numb.hpp"   // for as_number
#include "atox.hpp"   // for string_to_int
#include "model.hpp"
#include "entstr.hpp" // for entity_type_from_string, polymer_type_from_string
#include "mmcif_impl.hpp" // for set_cell_from_mmcif

namespace gemmi {

namespace impl {

inline std::unordered_map<std::string, std::array<float,6>>
get_anisotropic_u(cif::Block& block) {
  cif::Table aniso_tab = block.find("_atom_site_anisotrop.",
                                    {"id", "U[1][1]", "U[2][2]", "U[3][3]",
                                     "U[1][2]", "U[1][3]", "U[2][3]"});
  std::unordered_map<std::string, std::array<float,6>> aniso_map;
  for (auto ani : aniso_tab)
    aniso_map.emplace(ani[0], std::array<float,6>{{
                                (float) cif::as_number(ani[1]),
                                (float) cif::as_number(ani[2]),
                                (float) cif::as_number(ani[3]),
                                (float) cif::as_number(ani[4]),
                                (float) cif::as_number(ani[5]),
                                (float) cif::as_number(ani[6])}});
  return aniso_map;
}

inline
std::vector<std::string> transform_tags(std::string mstr, std::string vstr) {
  return {mstr + "[1][1]", mstr + "[1][2]", mstr + "[1][3]", vstr + "[1]",
          mstr + "[2][1]", mstr + "[2][2]", mstr + "[2][3]", vstr + "[2]",
          mstr + "[3][1]", mstr + "[3][2]", mstr + "[3][3]", vstr + "[3]"};
}

inline Transform get_transform_matrix(const cif::Table::Row& r) {
  Transform t;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      t.mat[i][j] = cif::as_number(r[4*i+j]);
    t.vec.at(i) = cif::as_number(r[4*i+3]);
  }
  return t;
}

inline ResidueId make_resid(const std::string& name,
                            const std::string& seqid,
                            const std::string* icode) {
  ResidueId rid;
  rid.name = name;
  if (icode)
    // the insertion code happens to be always a single letter
    rid.seqid.icode = cif::as_char(*icode, ' ');
  if (!seqid.empty()) {
    // old mmCIF files have auth_seq_id as number + icode (e.g. 15A)
    if (seqid.back() >= 'A') {
      if (rid.seqid.icode == ' ')
        rid.seqid.icode = seqid.back();
      else if (rid.seqid.icode != seqid.back())
        fail("Inconsistent insertion code in " + seqid);
      rid.seqid.num = cif::as_int(seqid.substr(0, seqid.size() - 1));
    } else {
      rid.seqid.num = cif::as_int(seqid, Residue::OptionalNum::None);
    }
  }
  return rid;
}

inline void read_helices(cif::Block& block, Structure& st) {
  for (const auto row : block.find("_struct_conf.", {
        "conf_type_id",                                              // 0
        "beg_auth_asym_id", "beg_label_comp_id",                     // 1-2
        "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",                 // 3-4
        "end_auth_asym_id", "end_label_comp_id",                     // 5-6
        "end_auth_seq_id", "?pdbx_end_PDB_ins_code",                 // 7-8
        "?pdbx_PDB_helix_class", "?pdbx_PDB_helix_length"})) {       // 9-10
    if (alpha_up(row.str(0)[0]) != 'H')
      continue;
    Helix h;
    h.start.chain_name = row.str(1);
    h.start.res_id = make_resid(row.str(2), row.str(3), row.ptr_at(4));
    h.end.chain_name = row.str(5);
    h.end.res_id = make_resid(row.str(6), row.str(7), row.ptr_at(8));
    if (row.has(9))
      h.set_helix_class_as_int(cif::as_int(row[9], -1));
    if (row.has(10))
      h.length = cif::as_int(row[10], -1);
    st.helices.push_back(h);
  }
}

inline void read_connectivity(cif::Block& block, Structure& st) {
  // label_ identifiers are not sufficient for HOH:
  // waters have null label_seq_id so we need auth_seq_id+icode.
  // And since we need auth_seq_id, we also use auth_asym_id for consistency.
  for (const auto row : block.find("_struct_conn.", {
        "id", "conn_type_id", // 0-1
        "ptnr1_auth_asym_id", "ptnr2_auth_asym_id", // 2-3
        "ptnr1_label_comp_id", "ptnr2_label_comp_id", // 4-5
        "ptnr1_label_atom_id", "ptnr2_label_atom_id", // 6-7
        "?pdbx_ptnr1_label_alt_id", "?pdbx_ptnr2_label_alt_id", // 8-9
        "ptnr1_auth_seq_id", "ptnr2_auth_seq_id", // 10-11
        "?pdbx_ptnr1_PDB_ins_code", "?pdbx_ptnr2_PDB_ins_code", // 12-13
        "?ptnr1_symmetry", "?ptnr2_symmetry", "?pdbx_dist_value"/*14-16*/})) {
    Connection c;
    c.name = row.str(0);
    std::string type = row.str(1);
    for (int i = 0; i != Connection::None; ++i)
      if (get_mmcif_connection_type_id(Connection::Type(i)) == type) {
        c.type = Connection::Type(i);
        break;
      }
    if (row.has2(14) && row.has2(15)) {
      c.asu = (row.str(14) == row.str(15) ? Asu::Same : Asu::Different);
    }
    if (row.has2(16))
      c.reported_distance = cif::as_number(row[16]);
    for (int i = 0; i < 2; ++i) {
      AtomAddress& a = c.atom[i];
      a.chain_name = row.str(2+i);
      if (row.has2(12+i))
        a.res_id.seqid.icode = cif::as_char(row[12+i], ' ');
      a.res_id = make_resid(row.str(4+i), row.str(10+i), row.ptr_at(12+i));
      a.atom_name = row.str(6+i);
      a.altloc = row.has2(8+i) ? cif::as_char(row[8+i], '\0') : '\0';
    }
    for (Model& mdl : st.models)
      mdl.connections.emplace_back(c);
  }
}

inline void fill_residue_entity_type(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (ResidueSpan sub : chain.subchains()) {
        EntityType etype = EntityType::Unknown;
        if (const Entity* ent = st.get_entity_of(sub))
          etype = ent->entity_type;
        if (etype == EntityType::Unknown) {
          if (sub[0].is_water())
            etype = EntityType::Water;
          else if (sub.length() > 1)
            etype = EntityType::Polymer;
          else
            etype = EntityType::NonPolymer;
        }
        for (Residue& residue : sub)
          residue.entity_type = etype;
      }
}

inline Structure make_structure_from_block(const cif::Block& block_) {
  using cif::as_number;
  using cif::as_string;
  // find() and Table don't have const variants, but we don't change anything.
  cif::Block& block = const_cast<cif::Block&>(block_);
  Structure st;
  st.input_format = CoorFormat::Mmcif;
  st.name = block.name;
  set_cell_from_mmcif(block, st.cell);
  st.spacegroup_hm = as_string(impl::find_spacegroup_hm_value(block));

  auto add_info = [&](std::string tag) {
    bool first = true;
    for (const std::string& v : block.find_values(tag))
      if (!cif::is_null(v)) {
        if (first)
          st.info[tag] = as_string(v);
        else
          st.info[tag] += "; " + as_string(v);
        first = false;
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

  for (const std::string& d : block.find_values("_refine.ls_d_res_high")) {
    double resol = cif::as_number(d);
    if (resol > 0 && (st.resolution == 0 || resol < st.resolution))
      st.resolution = resol;
  }

  for (auto row : block.find("_exptl.", {"method", "?crystals_number"})) {
    st.meta.experiments.emplace_back();
    st.meta.experiments.back().method = row.str(0);
    if (row.has2(1))
      st.meta.experiments.back().number_of_crystals = cif::as_int(row[1]);
  }

  for (auto row : block.find("_software.", {"name",
                                            "?classification",
                                            "?version",
                                            "?pdbx_ordinal"})) {
    st.meta.software.emplace_back();
    SoftwareItem& item = st.meta.software.back();
    item.name = row.str(0);
    if (row.has2(1))
      item.classification = software_classification_from_string(row.str(1));
    if (row.has2(2))
      item.version = row.str(2);
    if (row.has2(3))
      item.pdbx_ordinal = cif::as_int(row[3]);
  }

  std::vector<std::string> ncs_oper_tags = transform_tags("matrix", "vector");
  ncs_oper_tags.emplace_back("id");  // 12
  ncs_oper_tags.emplace_back("?code");  // 13
  cif::Table ncs_oper = block.find("_struct_ncs_oper.", ncs_oper_tags);
  for (auto op : ncs_oper) {
    bool given = op.has(13) && op.str(13) == "given";
    Transform tr = get_transform_matrix(op);
    if (!tr.is_identity())
      st.ncs.push_back({op.str(12), given, tr});
  }

  // PDBx/mmcif spec defines both _database_PDB_matrix.scale* and
  // _atom_sites.fract_transf_* as equivalent of pdb SCALE, but the former
  // is not used, so we ignore it.
  cif::Table fract_tv = block.find("_atom_sites.fract_transf_",
                                   transform_tags("matrix", "vector"));
  if (fract_tv.length() > 0) {
    Transform fract = get_transform_matrix(fract_tv[0]);
    st.cell.set_matrices_from_fract(fract);
  }

  // We read/write origx just for completeness, it's not used anywhere.
  cif::Table origx_tv = block.find("_database_PDB_matrix.",
                                   transform_tags("origx", "origx_vector"));
  if (origx_tv.length() > 0) {
    st.has_origx = true;
    st.origx = get_transform_matrix(origx_tv[0]);
  }

  auto aniso_map = get_anisotropic_u(block);

  // atom list
  enum { kId=0, kGroupPdb, kSymbol, kLabelAtomId, kAltId, kLabelCompId,
         kLabelAsymId, kLabelSeqId, kInsCode, kX, kY, kZ, kOcc, kBiso, kCharge,
         kAuthSeqId, kAuthCompId, kAuthAsymId, kAuthAtomId, kModelNum };
  cif::Table atom_table = block.find("_atom_site.",
                                     {"id",
                                      "?group_PDB",
                                      "type_symbol",
                                      "label_atom_id",
                                      "label_alt_id",
                                      "label_comp_id",
                                      "label_asym_id",
                                      "?label_seq_id",
                                      "?pdbx_PDB_ins_code",
                                      "Cartn_x",
                                      "Cartn_y",
                                      "Cartn_z",
                                      "occupancy",
                                      "B_iso_or_equiv",
                                      "?pdbx_formal_charge",
                                      "auth_seq_id",
                                      "?auth_comp_id",
                                      "?auth_asym_id",
                                      "?auth_atom_id",
                                      "?pdbx_PDB_model_num"});
  const int kCompId = atom_table.has_column(kAuthCompId) ? kAuthCompId
                                                         : kLabelCompId;
  const int kAsymId = atom_table.has_column(kAuthAsymId) ? kAuthAsymId
                                                         : kLabelAsymId;
  const int kAtomId = atom_table.has_column(kAuthAtomId) ? kAuthAtomId
                                                         : kLabelAtomId;
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  if (atom_table.length() != 0) {
    if (atom_table.has_column(kModelNum))
      model = &st.find_or_add_model(atom_table[0].str(kModelNum));
    else
      model = &st.find_or_add_model("1");
  }
  for (auto row : atom_table) {
    if (row.has(kModelNum) && row[kModelNum] != model->name) {
      model = &st.find_or_add_model(row.str(kModelNum));
      chain = nullptr;
    }
    if (!chain || as_string(row[kAsymId]) != chain->name) {
      model->chains.emplace_back(as_string(row[kAsymId]));
      chain = &model->chains.back();
      resi = nullptr;
    }
    ResidueId rid = make_resid(as_string(row[kCompId]),
                               as_string(row[kAuthSeqId]),
                               row.has(kInsCode) ? &row[kInsCode] : nullptr);
    if (!resi || !resi->matches(rid)) {
      resi = chain->find_or_add_residue(rid);
      if (resi->atoms.empty()) {
        if (row.has2(kLabelSeqId))
          resi->label_seq = cif::as_int(row[kLabelSeqId]);
        resi->subchain = row.str(kLabelAsymId);
        // don't check if group_PDB is consistent, it's not that important
        if (row.has2(kGroupPdb))
          for (int i = 0; i < 2; ++i) { // first character could be " or '
            const char c = alpha_up(row[kGroupPdb][i]);
            if (c == 'A' || c == 'H' || c == '\0')
              resi->het_flag = c;
          }
      }
    } else if (resi->seqid != rid.seqid) {
      fail("Inconsistent sequence ID: " + resi->str() + " / " + rid.str());
    }
    Atom atom;
    atom.name = as_string(row[kAtomId]);
    // altloc is always a single letter (not guaranteed by the mmCIF spec)
    atom.altloc = cif::as_char(row[kAltId], '\0');
    atom.charge = row.has2(kCharge) ? cif::as_int(row[kCharge]) : 0;
    atom.element = gemmi::Element(as_string(row[kSymbol]));
    // According to the PDBx/mmCIF spec _atom_site.id can be a string,
    // but in all the files it is a serial number; its value is not essential,
    // so we just ignore non-integer ids.
    atom.serial = string_to_int(row[kId], false);
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

  cif::Table polymer_types = block.find("_entity_poly.", {"entity_id", "type"});
  for (auto row : block.find("_entity.", {"id", "type"})) {
    Entity ent(row.str(0));
    ent.entity_type = entity_type_from_string(row.str(1));
    ent.polymer_type = PolymerType::Unknown;
    if (polymer_types.ok()) {
      try {
        std::string poly_type = polymer_types.find_row(ent.name).str(1);
        ent.polymer_type = polymer_type_from_string(poly_type);
      } catch (std::runtime_error&) {}
    }
    st.entities.push_back(ent);
  }

  for (auto row : block.find("_entity_poly_seq.",
                             {"entity_id", "num", "mon_id"}))
    if (Entity* ent = st.get_entity(row.str(0)))
      ent->poly_seq.push_back({cif::as_int(row[1], -1), row.str(2)});

  for (auto row : block.find("_struct_asym.", {"id", "entity_id"}))
    if (Entity* ent = st.get_entity(row.str(1)))
      ent->subchains.push_back(row.str(0));

  fill_residue_entity_type(st);

  st.setup_cell_images();

  // CISPEP
  for (auto row : block.find("_struct_mon_prot_cis.",
                             {"pdbx_PDB_model_num", "auth_asym_id",  // 0-1
                              "auth_seq_id", "?pdbx_PDB_ins_code",   // 2-3
                              "?label_comp_id", "?auth_comp_id"})) { // 4-5
    if (row.has2(0) && row.has2(1) && row.has2(2) &&
        (row.has2(4) || row.has2(5)))
      if (Model* mdl = st.find_model(row[0])) {
        std::string comp = row.str(row.has2(4) ? 4 : 5);
        ResidueId rid = make_resid(comp, row.str(2), row.ptr_at(3));
        if (Residue* res = mdl->find_residue(row[1], rid))
          res->is_cis = true;
      }
  }

  read_helices(block, st);

  read_connectivity(block, st);

  return st;
}

} // namespace impl

inline Structure make_structure_from_block(const cif::Block& block) {
  return impl::make_structure_from_block(block);
}


} // namespace gemmi
#endif
