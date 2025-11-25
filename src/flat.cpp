// Copyright Global Phasing Ltd.

#include "gemmi/flat.hpp"
#include <cstring>  // for strcpy

namespace gemmi {

FlatStructure::FlatStructure(const Structure& st) {
  empty_st = st.empty_copy();
  size_t n = count_atom_sites(st);
  table.reserve(n);
  FlatAtom fa;
  for (const Model& model : st.models) {
    fa.model_num = model.num;
    for (const Chain& chain : model.chains) {
      if (chain.name.size() > 7)
        fail("FlatStructure doesn't support 8+ char subchain names: ", chain.name);
      std::strcpy(fa.chain_id, chain.name.c_str());
      for (const Residue& res : chain.residues) {
        if (res.name.size() > 7)
          fail("FlatStructure doesn't support 8+ char residue names: ", res.name);
        std::strcpy(fa.residue_name, res.name.c_str());
        if (res.subchain.size() > 7)
          fail("FlatStructure doesn't support 8+ char subchain names: ", res.subchain);
        std::strcpy(fa.subchain, res.subchain.c_str());
        if (res.entity_id.size() > 7)
          fail("FlatStructure doesn't support 8+ char entity IDs: ", res.entity_id);
        std::strcpy(fa.entity_id, res.entity_id.c_str());
        fa.seq_id = res.seqid;
        fa.het_flag = res.het_flag;
        fa.entity_type = res.entity_type;
        for (const Atom& atom : res.atoms) {
          if (atom.name.size() > 7)
            fail("FlatStructure doesn't support 8+ char atom names: ", atom.name);
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

Structure FlatStructure::generate_structure() {
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

} // namespace gemmi
