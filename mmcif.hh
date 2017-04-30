// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hh.

#ifndef GEMMI_MMCIF_HH_
#define GEMMI_MMCIF_HH_

#include <string>
#include <iostream> // temporary
#include "cif.hh"
#include "numb.hh"
#include "model.hh"

namespace gemmi {
namespace mol {

template<typename T>
T* find_or_add(std::vector<T>& vec, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(), [&name](const T& m) {
      return m.name == name;
  });
  if (it != vec.end())
    return &*it;
  vec.emplace_back(name);
  return &vec.back();
}

inline Residue* find_or_add(std::vector<Residue>& vec, int seq_id,
                            char ins_code, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(), [&](const Residue& r) {
      return r.seq_id == seq_id && r.ins_code == ins_code && r.name == name;
  });
  if (it != vec.end())
    return &*it;
  vec.emplace_back(seq_id, ins_code, name);
  return &vec.back();
}

inline Structure structure_from_cif_block(const cif::Block& block) {
  Structure st;

  // unit cell and symmetry
  cif::TableView cell = block.find("_cell.", {"length_a", "length_b",
                      "length_c", "angle_alpha", "angle_beta", "angle_gamma"});
  if (cell.ok()) {
    auto c = cell.one();
    st.cell.a = c.as_num(0);
    st.cell.b = c.as_num(1);
    st.cell.c = c.as_num(2);
    st.cell.alpha = c.as_num(3);
    st.cell.beta = c.as_num(4);
    st.cell.gamma = c.as_num(5);
  }
  st.sg_hm = block.find("_symmetry.space_group_name_H-M").one().as_str(0);

  auto add_info = [&](std::string tag) {
    cif::TableView t = block.find(tag);
    if (t.length() >= 1)
      st.info[tag] = t[0].as_str(0);
  };
  add_info("_entry.id");
  add_info("_cell.Z_PDB");
  add_info("_struct.title");
  add_info("_exptl.method");
  add_info("_database_PDB_rev.date_original");
  add_info("_struct_keywords.pdbx_keywords");
  add_info("_struct_keywords.text");

  // sequence
  // TODO

  // atom list
  enum { kSymbol=0, kAtomId, kAltId, kCompId, kAsymId, kSeqId, kInsCode,
         kX, kY, kZ, kOcc, kBiso, kCharge, kModelNum };
  cif::TableView atom_table = block.find("_atom_site.",
                                         {"type_symbol",
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
                                          "pdbx_PDB_model_num"});
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  for (auto row : atom_table) {
    if (!model || row[kModelNum] != model->name) {
      model = find_or_add(st.models, row[kModelNum]);
      chain = nullptr;
    }
    if (!chain || row[kAsymId] != chain->name) {
      chain = find_or_add(model->chains, row[kAsymId]);
      resi = nullptr;
    }
    int seq_id = row.is_null(kSeqId) ? 0 : row.as_int(kSeqId);
    // let us assume that the insertion code is a single letter
    assert(row[kInsCode].size() == 1);
    char ins_code = row[kInsCode][0];
    if (ins_code == '?' || ins_code == '.')
      ins_code = '\0';
    if (!resi || seq_id != resi->seq_id || row[kCompId] != resi->name) {
      resi = find_or_add(chain->residues,
                         seq_id, ins_code, row.as_str(kCompId));
    }
    Atom atom;
    atom.name = cif::as_string(row[kAtomId]);
    atom.altloc = cif::as_string(row[kAltId])[0];
    atom.charge = cif::is_null(row[kCharge]) ? 0 : cif::as_int(row[kCharge]);
    atom.element = Element(cif::as_string(row[kSymbol]));
    atom.x = cif::as_number(row[kX]);
    atom.y = cif::as_number(row[kY]);
    atom.z = cif::as_number(row[kZ]);
    atom.occ = cif::as_number(row[kOcc], 1.0);
    atom.b_iso = cif::as_number(row[kBiso], 50.0);
    resi->atoms.emplace_back(atom);
  }
  return st;
}

inline Structure read_atoms(const cif::Document& doc) {
  return structure_from_cif_block(doc.sole_block());
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
