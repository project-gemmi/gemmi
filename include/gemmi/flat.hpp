// Copyright Global Phasing Ltd.
//
// FlatStructure, FlatAtom

#ifndef GEMMI_FLAT_HPP_
#define GEMMI_FLAT_HPP_

#include <vector>
#include "model.hpp"
#include "calculate.hpp"

namespace gemmi {

struct FlatAtom {
  char atom_name[8] {};
  char residue_name[8] = {};
  char chain_id[8] = {};
  char subchain[8] = {};
  char entity_id[8] = {};
  SeqId seq_id;
  Position pos;
  float occ = 1.0f;
  float b_iso = 20.0f; // arbitrary default value
  char altloc = '\0'; // 0 if not set
  char het_flag = '\0';   // 'A' = ATOM, 'H' = HETATM, 0 = unspecified
  EntityType entity_type = EntityType::Unknown;
  Element element = El::X;
  signed char charge = 0;  // [-8, +8]
  SMat33<float> aniso = {0, 0, 0, 0, 0, 0};
  int model_num;
  bool selected = false;
};
struct FlatStructure {
  Structure empty_st;
  std::vector<FlatAtom> table;
  FlatStructure(const Structure& st) {
    empty_st = st.empty_copy();
    size_t n = count_atom_sites(st);
    table.reserve(n);
    FlatAtom fa;
    for (const Model& model : st.models) {
      fa.model_num = model.num;
      for (const Chain& chain : model.chains) {
        std::strcpy(fa.chain_id, chain.name.c_str());
        for (const Residue& res : chain.residues) {
          std::strcpy(fa.residue_name, res.name.c_str());
          std::strcpy(fa.subchain, res.subchain.c_str());
          std::strcpy(fa.entity_id, res.entity_id.c_str());
          fa.seq_id = res.seqid;
          fa.het_flag = res.het_flag;
          fa.entity_type = res.entity_type;
          for (const Atom& atom : res.atoms) {
            std::strcpy(fa.atom_name, atom.name.c_str());
            fa.pos = atom.pos;
            fa.occ = atom.occ;
            fa.b_iso = atom.b_iso;
            fa.altloc = atom.altloc;
            fa.element = atom.element;
            fa.charge = atom.charge;
            fa.aniso = atom.aniso;
            table.push_back(fa);
          }
        }
      }
    }
  }
  Structure generate_structure() {
    Structure st(empty_st);
    for (const FlatAtom& fa : table) {
      Model& model = st.find_or_add_model(fa.model_num);
      Chain* chain = model.find_chain(fa.chain_id);
      if (!chain) {
        model.chains.emplace_back(fa.chain_id);
        chain = &model.chains.back();
      }
      ResidueId rid{fa.seq_id, {}, fa.residue_name};
      Residue* residue = chain->find_or_add_residue(rid);
      residue->het_flag = fa.het_flag;
      residue->subchain = fa.subchain;
      residue->entity_id = fa.entity_id;
      residue->entity_type = fa.entity_type;
      
      Atom atom;
      atom.name = fa.atom_name;
      atom.pos = fa.pos;
      atom.occ = fa.occ;
      atom.b_iso = fa.b_iso;
      atom.altloc = fa.altloc;
      atom.element = fa.element;
      atom.charge = fa.charge;
      atom.aniso = fa.aniso;
      
      residue->atoms.emplace_back(atom);
    }
    return st;
  }
};


namespace impl {
}
} // namespace gemmi

#endif
