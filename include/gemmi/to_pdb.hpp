// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format (Structure -> pdb file).

#ifndef GEMMI_TO_PDB_HPP_
#define GEMMI_TO_PDB_HPP_

#include "model.hpp"
#include <ostream>

namespace gemmi {

struct PdbWriteOptions {
  bool seqres_records = true;
  bool ssbond_records = true;
  bool cryst1_record = true;
  bool link_records = true;
  bool cispep_records = true;
  bool ter_records = true;
  bool numbered_ter = true;
  bool ter_ignores_type = false;
  bool use_linkr = false;
};

GEMMI_DLL void write_pdb(const Structure& st, std::ostream& os,
                         PdbWriteOptions opt=PdbWriteOptions());
GEMMI_DLL void write_minimal_pdb(const Structure& st, std::ostream& os,
                                 PdbWriteOptions opt=PdbWriteOptions());
GEMMI_DLL std::string make_pdb_headers(const Structure& st);

} // namespace gemmi

#endif
