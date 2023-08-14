// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format (Structure -> pdb file).

#ifndef GEMMI_TO_PDB_HPP_
#define GEMMI_TO_PDB_HPP_

#include "model.hpp"
#include <ostream>

namespace gemmi {

struct PdbWriteOptions {
  bool minimal_file = false;    // disable many records not listed below (HEADER, TITLE, ...)
  bool atom_records = true;     // write atomic models (set to false for headers only)
  bool seqres_records = true;   // write SEQRES
  bool ssbond_records = true;   // write SSBOND
  bool link_records = true;     // write LINK
  bool cispep_records = true;   // write CISPEP
  bool cryst1_record = true;    // write CRYST1
  bool ter_records = true;      // write TER records
  bool conect_records = false;  // write CONECT - matters only if add_conect() was used
  bool end_record = true;       // write END
  bool numbered_ter = true;     // TER record gets own serial number
  bool ter_ignores_type = false; // put TER after last atom in Chain (even if it's water)
  bool use_linkr = false;       // use non-standard Refmac LINKR record instead of LINK
  bool preserve_serial = false; // use serial numbers from Atom.serial
  // end of snippet for mol.rst

  static PdbWriteOptions minimal() {
    PdbWriteOptions opt;
    opt.minimal_file = true;
    opt.seqres_records = false;
    opt.ssbond_records = false;
    opt.link_records = false;
    opt.cispep_records = false;
    opt.end_record = false;
    return opt;
  }
  static PdbWriteOptions headers_only() {
    PdbWriteOptions opt;
    opt.atom_records = false;
    opt.end_record = false;
    return opt;
  }
};

GEMMI_DLL void write_pdb(const Structure& st, std::ostream& os,
                         PdbWriteOptions opt=PdbWriteOptions());
GEMMI_DLL std::string make_pdb_string(const Structure& st,
                                      PdbWriteOptions opt=PdbWriteOptions());

// deprecated
inline void write_minimal_pdb(const Structure& st, std::ostream& os) {
  write_pdb(st, os, PdbWriteOptions::minimal());
}
// deprecated
inline std::string make_pdb_headers(const Structure& st) {
  return make_pdb_string(st, PdbWriteOptions::headers_only());
}

} // namespace gemmi

#endif
