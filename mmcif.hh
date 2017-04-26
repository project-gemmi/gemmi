// Copyright 2017 Global Phasing Ltd.
//
// Work in progress...

#ifndef GEMMI_MMCIF_HH_
#define GEMMI_MMCIF_HH_

#include <string>
#include <iostream> // temporary
#include "cif.hh"
#include "numb.hh"
#include "elem.hh"

namespace gemmi {
namespace mol {

struct Structure;
struct Model;
struct Chain;
struct Residue;

enum class EntityType { Unknown, Polymer, NonPolymer, Water };


struct Atom {
  std::string name;
  char altloc;
  signed char charge;  // [-8, +8]
  Element element = El::X;
  double x, y, z;
  float occ;
  float b_iso;
  Residue* parent = nullptr;
};

struct Residue {
  int seq_id = -1000;
  char ins_code = '\0';
  std::string name;
  std::vector<Atom> atoms;
  Chain* parent = nullptr;
  Residue(int id, char ins, std::string rname) noexcept
    : seq_id(id), ins_code(ins), name(rname) {}
};

struct Chain {
  std::string name;
  EntityType entity_type = EntityType::Unknown;
  std::vector<Residue> residues;
  Model* parent = nullptr;
  explicit Chain(std::string cname) noexcept : name(cname) {}
};

struct Model {
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  Structure* parent = nullptr;
  explicit Model(std::string mname) noexcept : name(mname) {}
};

struct UnitCell {
  double lengths[3];
  double angles[3];
};

struct Structure {
  std::string entry_id;
  UnitCell cell;
  int z;
  std::string sg_hm;
  std::vector<Model> models;
  // std::vector<Ops> ncs;
};

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
    atom.y = cif::as_number(row[kZ]);
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
