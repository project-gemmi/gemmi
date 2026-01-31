//! @file
//! @brief Write Structure to PDB file format.

// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format (Structure -> pdb file).

#ifndef GEMMI_TO_PDB_HPP_
#define GEMMI_TO_PDB_HPP_

#include "model.hpp"
#include <ostream>

namespace gemmi {

//! @brief Options for PDB file output.
//!
//! Controls which record types are written and formatting details.
struct PdbWriteOptions {
  bool minimal_file = false;       //!< Disable many records (HEADER, TITLE, etc.)
  bool atom_records = true;        //!< Write ATOM/HETATM records
  bool seqres_records = true;      //!< Write SEQRES records
  bool ssbond_records = true;      //!< Write SSBOND records
  bool link_records = true;        //!< Write LINK records
  bool cispep_records = true;      //!< Write CISPEP records
  bool cryst1_record = true;       //!< Write CRYST1 record
  bool ter_records = true;         //!< Write TER records
  bool conect_records = false;     //!< Write CONECT records (if add_conect() used)
  bool end_record = true;          //!< Write END record
  bool numbered_ter = true;        //!< TER gets its own serial number
  bool ter_ignores_type = false;   //!< Put TER after last atom in chain (even water)
  bool use_linkr = false;          //!< Use Refmac LINKR instead of LINK
  bool use_link_id = false;        //!< Write Link::link_id in LINK (implied by use_linkr)
  bool preserve_serial = false;    //!< Use Atom.serial for serial numbers
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
  //! @brief Create options for headers-only PDB output.
  //! @return Options with atom_records disabled
  static PdbWriteOptions headers_only() {
    PdbWriteOptions opt;
    opt.atom_records = false;
    opt.end_record = false;
    return opt;
  }
};

//! @brief Determine if residue should use HETATM record.
//! @param res Residue to check
//! @return true if HETATM, false if ATOM
GEMMI_DLL bool use_hetatm(const Residue& res);

//! @brief Write Structure to stream in PDB format.
//! @param st Structure to write
//! @param os Output stream
//! @param opt Write options (default: all records)
GEMMI_DLL void write_pdb(const Structure& st, std::ostream& os, PdbWriteOptions opt={});

//! @brief Convert Structure to PDB format string.
//! @param st Structure to convert
//! @param opt Write options (default: all records)
//! @return PDB-formatted string
GEMMI_DLL std::string make_pdb_string(const Structure& st, PdbWriteOptions opt={});

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
