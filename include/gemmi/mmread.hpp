// Copyright 2017 Global Phasing Ltd.
//
// Read any supported coordinate file.

#ifndef GEMMI_MMREAD_HPP_
#define GEMMI_MMREAD_HPP_

#include "chemcomp_xyz.hpp" // for make_structure_from_chemcomp_block
#include "cif.hpp"       // for cif::read
#include "fail.hpp"      // for fail
#include "input.hpp"     // for BasicInput
#include "json.hpp"      // for read_mmjson
#include "mmcif.hpp"     // for make_structure_from_block
#include "model.hpp"     // for Structure
#include "pdb.hpp"       // for read_pdb
#include "util.hpp"      // for iends_with

namespace gemmi {

inline CoorFormat coor_format_from_ext(const std::string& path) {
  if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
    return CoorFormat::Pdb;
  if (iends_with(path, ".cif") || iends_with(path, ".mmcif"))
    return CoorFormat::Mmcif;
  if (iends_with(path, ".json"))
    return CoorFormat::Mmjson;
  return CoorFormat::Unknown;
}

inline Structure make_structure(const cif::Document& doc) {
  // mmCIF files for deposition may have more than one block:
  // coordinates in the first block and restraints in the others.
  for (size_t i = 1; i < doc.blocks.size(); ++i)
    if (doc.blocks[i].has_tag("_atom_site.id"))
      fail("2+ blocks are ok if only the first one has coordinates;\n"
           "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
  return make_structure_from_block(doc.blocks.at(0));
}

template<typename T>
Structure read_structure(T&& input, CoorFormat format=CoorFormat::Unknown) {
  bool any = (format == CoorFormat::UnknownAny);
  if (format == CoorFormat::Unknown || any)
    format = coor_format_from_ext(input.basepath());
  switch (format) {
    case CoorFormat::Pdb:
      return read_pdb(input);
    case CoorFormat::Mmcif: {
      cif::Document doc = cif::read(input);
      if (any) {
        int n = check_chemcomp_block_number(doc);
        // first handle special case - refmac dictionary or CCD file
        if (n != -1)
          return make_structure_from_chemcomp_block(doc.blocks[n]);
      }
      return make_structure(doc);
    }
    case CoorFormat::Mmjson:
      return make_structure_from_block(cif::read_mmjson(input).sole_block());
    case CoorFormat::ChemComp:
      return make_structure_from_chemcomp_doc(cif::read(input));
    case CoorFormat::Unknown:
    case CoorFormat::UnknownAny:
      fail("Unknown format of " +
           (input.path().empty() ? "coordinate file" : input.path()) + ".");
  }
  unreachable();
}

inline Structure read_structure_file(const std::string& path,
                                     CoorFormat format=CoorFormat::Unknown) {
  return read_structure(BasicInput(path), format);
}

} // namespace gemmi
#endif
