// Copyright Global Phasing Ltd.

#include "gemmi/flat.hpp"
#include "gemmi/calculate.hpp"  // for count_atom_sites
#include <cstring>  // for memcpy, memset

namespace gemmi {

namespace {
template<size_t N>
void copy_padded(char (&dest)[N], const std::string& src) {
  std::memcpy(dest, src.c_str(), src.size());
  std::memset(dest + src.size(), 0, N - src.size());
}
}  // anonymous namespace

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
      copy_padded(fa.chain_id, chain.name);
      for (const Residue& res : chain.residues) {
        if (res.name.size() > 7)
          fail("FlatStructure doesn't support 8+ char residue names: ", res.name);
        copy_padded(fa.residue_name, res.name);
        if (res.subchain.size() > 7)
          fail("FlatStructure doesn't support 8+ char subchain names: ", res.subchain);
        copy_padded(fa.subchain, res.subchain);
        if (res.entity_id.size() > 7)
          fail("FlatStructure doesn't support 8+ char entity IDs: ", res.entity_id);
        copy_padded(fa.entity_id, res.entity_id);
        fa.seq_id = res.seqid;
        fa.het_flag = res.het_flag;
        fa.entity_type = res.entity_type;
        for (const Atom& atom : res.atoms) {
          if (atom.name.size() > 7)
            fail("FlatStructure doesn't support 8+ char atom names: ", atom.name);
          copy_padded(fa.atom_name, atom.name);
          fa.pos = atom.pos;
          fa.occ = atom.occ;
          fa.b_iso = atom.b_iso;
          fa.altloc = atom.altloc;
          fa.element = atom.element;
          fa.charge = atom.charge;
          fa.aniso = atom.aniso;
          fa.serial = atom.serial;
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
    atom.serial = fa.serial;
    residue->atoms.emplace_back(atom);
  }
  return st;
}

} // namespace gemmi
