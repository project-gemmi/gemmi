// Copyright Global Phasing Ltd.
//
// FlatStructure, FlatAtom

#ifndef GEMMI_FLAT_HPP_
#define GEMMI_FLAT_HPP_

#include <vector>
#include "model.hpp"

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
  int serial = 0;
  bool selected = false;

  std::string atom_str() const {
    ResidueId resid{seq_id, "", residue_name};
    return gemmi::atom_str(chain_id, resid, atom_name, altloc);
  }
};
struct GEMMI_DLL FlatStructure {
  Structure empty_st;
  std::vector<FlatAtom> table;

  // Structure <-> FlatStructure
  FlatStructure(const Structure& st);
  Structure generate_structure();
};

} // namespace gemmi

#endif
