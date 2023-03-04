// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hpp.

#ifndef GEMMI_MMCIF_HPP_
#define GEMMI_MMCIF_HPP_

#include <string>
#include "cifdoc.hpp"
#include "fail.hpp"        // for fail
#include "model.hpp"       // for Structure

namespace gemmi {

GEMMI_DLL Structure make_structure_from_block(const cif::Block& block);

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

} // namespace gemmi
#endif
