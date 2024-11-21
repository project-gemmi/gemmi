// Copyright 2017 Global Phasing Ltd.
//
// Read mmCIF (PDBx/mmCIF) file into a Structure from model.hpp.

#ifndef GEMMI_MMCIF_HPP_
#define GEMMI_MMCIF_HPP_

#include <string>
#include "cifdoc.hpp"      // for Block, etc
#include "fail.hpp"        // for fail
#include "model.hpp"       // for Structure

namespace gemmi {

/// structure from a coordinate mmCIF block
GEMMI_DLL Structure make_structure_from_block(const cif::Block& block);

/// structure from a coordinate mmCIF document
inline Structure make_structure(cif::Document&& doc, cif::Document* save_doc=nullptr) {
  // mmCIF files for deposition may have more than one block:
  // coordinates in the first block and restraints in the others.
  for (size_t i = 1; i < doc.blocks.size(); ++i)
    if (doc.blocks[i].has_tag("_atom_site.id"))
      fail("2+ blocks are ok if only the first one has coordinates;\n"
           "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
  Structure st = make_structure_from_block(doc.blocks.at(0));
  if (save_doc)
    *save_doc = std::move(doc);
  return st;
}

// Reading chemical component as a coordinate file.
enum class ChemCompModel {
  Xyz      = 1, // _chem_comp_atom.x, etc
  Example  = 2, // _chem_comp_atom.model_Cartn_x
  Ideal    = 4  // _chem_comp_atom.pdbx_model_Cartn_x_ideal
};

constexpr int operator|(ChemCompModel a, ChemCompModel b) { return (int)a | (int)b; }

/// make_residue_from_chemcomp_block
GEMMI_DLL Residue make_residue_from_chemcomp_block(const cif::Block& block, ChemCompModel kind);

inline Model make_model_from_chemcomp_block(const cif::Block& block, ChemCompModel kind) {
  Model model;
  model.chains.emplace_back("");
  model.chains[0].residues.push_back(make_residue_from_chemcomp_block(block, kind));
  return model;
}

// For CCD input - returns a structure with two single-residue models:
// example (model_Cartn_x) and ideal (pdbx_model_Cartn_x_ideal).
// For Refmac dictionary (monomer library) files returns structure with
// a single model.
inline Structure make_structure_from_chemcomp_block(const cif::Block& block, int which=7) {
  Structure st;
  st.input_format = CoorFormat::ChemComp;
  if (const std::string* name = block.find_value("_chem_comp.id"))
    st.name = *name;
  auto ok = [which](ChemCompModel x) { return which & static_cast<int>(x); };
  if (ok(ChemCompModel::Xyz) && block.has_any_value("_chem_comp_atom.x"))
    st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Xyz));
  if (ok(ChemCompModel::Example) && block.has_any_value("_chem_comp_atom.model_Cartn_x"))
    st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Example));
  if (ok(ChemCompModel::Ideal) && block.has_any_value("_chem_comp_atom.pdbx_model_Cartn_x_ideal"))
    st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Ideal));
  st.renumber_models();
  return st;
}

// a helper function for use with make_structure_from_chemcomp_block():
//   int n = check_chemcomp_block_number(doc);
//   if (n != -1)
//     Structure st = make_structure_from_chemcomp_block(doc.blocks[n]);
inline int check_chemcomp_block_number(const cif::Document& doc) {
  // monomer library file without global_
  if (doc.blocks.size() == 2 && doc.blocks[0].name == "comp_list")
    return 1;
  // monomer library file with global_
  if (doc.blocks.size() == 3 && doc.blocks[0].name.empty() &&
      doc.blocks[1].name == "comp_list")
    return 2;
  // CCD file
  if (doc.blocks.size() == 1 &&
      !doc.blocks[0].has_tag("_atom_site.id") &&
      !doc.blocks[0].has_tag("_cell.length_a") &&
      doc.blocks[0].has_tag("_chem_comp_atom.atom_id"))
    return 0;
  return -1;
}

inline Structure make_structure_from_chemcomp_doc(const cif::Document& doc,
                                                  cif::Document* save_doc=nullptr,
                                                  int which=7) {
  int n = check_chemcomp_block_number(doc);
  if (n == -1)
    fail("Not a chem_comp format.");
  Structure st = make_structure_from_chemcomp_block(doc.blocks[n], which);
  if (save_doc)
    *save_doc = std::move(doc);
  return st;
}

} // namespace gemmi
#endif
