//! @file
//! @brief Read mmCIF (PDBx/mmCIF) files into Structure objects.
//!
//! Read mmCIF (PDBx/mmCIF) file into a Structure from model.hpp.

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

//! @brief Fill an existing Structure from an mmCIF coordinate block.
//! @param block_ CIF block containing coordinate data
//! @param st Structure object to populate
//!
//! Parses _atom_site, _cell, _entity, and other mmCIF categories into the structure.
GEMMI_DLL void populate_structure_from_block(const cif::Block& block_, Structure& st);

//! @brief Create Structure from an mmCIF coordinate block.
//! @param block_ CIF block containing coordinate data
//! @return New Structure object
//!
//! Convenience wrapper around populate_structure_from_block().
inline Structure make_structure_from_block(const cif::Block& block_) {
  gemmi::Structure st;
  populate_structure_from_block(block_, st);
  return st;
}

//! @brief Create Structure from an mmCIF document.
//! @param doc CIF document (moved, not copied)
//! @param save_doc Optional pointer to store the original document
//! @return Structure from the first block with coordinates
//!
//! Expects coordinates in the first block only. Additional blocks may
//! contain restraints but must not have _atom_site data.
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

//! @brief Chemical component coordinate model types.
//!
//! Used when reading CCD (Chemical Component Dictionary) or monomer library files.
//! Multiple models can be combined using bitwise OR.
enum class ChemCompModel {
  Xyz      = 1, //!< _chem_comp_atom.x, y, z coordinates
  Example  = 2, //!< _chem_comp_atom.model_Cartn_x (example coordinates)
  Ideal    = 4  //!< _chem_comp_atom.pdbx_model_Cartn_x_ideal (ideal coordinates)
};

constexpr int operator|(ChemCompModel a, ChemCompModel b) { return (int)a | (int)b; }

//! @brief Create Residue from a chemical component block.
//! @param block CIF block containing _chem_comp_atom data
//! @param kind Which coordinate model to use (Xyz, Example, or Ideal)
//! @return Residue with atoms from the chemical component
GEMMI_DLL Residue make_residue_from_chemcomp_block(const cif::Block& block, ChemCompModel kind);

//! @brief Create Model containing a single residue from a chemical component block.
//! @param block CIF block with _chem_comp_atom data
//! @param kind Which coordinate model to use
//! @return Model with one chain containing one residue
inline Model make_model_from_chemcomp_block(const cif::Block& block, ChemCompModel kind) {
  Model model;
  model.chains.emplace_back("");
  model.chains[0].residues.push_back(make_residue_from_chemcomp_block(block, kind));
  return model;
}

//! @brief Create Structure from a chemical component block.
//! @param block CIF block with _chem_comp_atom data
//! @param which Bitmask of ChemCompModel values (default 7 = all models)
//! @return Structure with models for each available coordinate set
//!
//! For CCD files: creates models for example and ideal coordinates.
//! For Refmac dictionary files: creates a single model.
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

//! @brief Identify which block in a document contains chemical component data.
//! @param doc CIF document to check
//! @return Block index if chemical component detected, -1 otherwise
//!
//! Detects CCD files and Refmac monomer library files by structure.
//! Usage: int n = check_chemcomp_block_number(doc);
//!        if (n != -1) st = make_structure_from_chemcomp_block(doc.blocks[n]);
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

//! @brief Create Structure from a chemical component document.
//! @param doc CIF document (CCD or Refmac dictionary)
//! @param save_doc Optional pointer to store the original document
//! @param which Bitmask of ChemCompModel values (default 7 = all models)
//! @return Structure with chemical component coordinate models
//!
//! Automatically detects and reads CCD and monomer library formats.
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
