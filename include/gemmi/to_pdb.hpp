// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format (Structure -> pdb file).

/// @file
/// @brief Serialization of Structure to PDB fixed-column format.

#ifndef GEMMI_TO_PDB_HPP_
#define GEMMI_TO_PDB_HPP_

#include "model.hpp"
#include <ostream>

namespace gemmi {

/// @brief Control record types and fields written in PDB output.
///
/// Options to selectively enable/disable PDB record types and special
/// formatting modes when writing Structure to PDB format.
struct PdbWriteOptions {
  /// @brief Minimize output by omitting ancillary records (HEADER, TITLE, REMARK, etc.).
  bool minimal_file = false;
  /// @brief Write ATOM/HETATM records with coordinates. If false, omit atomic models.
  bool atom_records = true;
  /// @brief Write SEQRES records with polymer sequences.
  bool seqres_records = true;
  /// @brief Write SSBOND records (disulfide bonds between cysteines).
  bool ssbond_records = true;
  /// @brief Write LINK records for non-disulfide chemical bonds.
  bool link_records = true;
  /// @brief Write CISPEP records for cis-peptide bonds.
  bool cispep_records = true;
  /// @brief Write CRYST1 record with unit cell and space group information.
  bool cryst1_record = true;
  /// @brief Write TER records to mark chain terminations.
  bool ter_records = true;
  /// @brief Write CONECT records for inter-residue connectivity (requires add_conect() preprocessing).
  bool conect_records = false;
  /// @brief Write END record at end of file.
  bool end_record = true;
  /// @brief Assign own serial numbers to TER records if true; reuse ATOM serial if false.
  bool numbered_ter = true;
  /// @brief If true, place TER after the last atom in each chain regardless of residue type.
  /// If false, omit TER after water/heteroatom chains.
  bool ter_ignores_type = false;
  /// @brief Use non-standard Refmac LINKR record instead of standard LINK record.
  bool use_linkr = false;
  /// @brief Write Link::link_id field in LINK record instead of distance (if link_id is non-empty).
  /// Automatically enabled by use_linkr.
  bool use_link_id = false;
  /// @brief Use atom serial numbers from Atom::serial. If false, generate new serial numbers.
  bool preserve_serial = false;
  // end of snippet for mol.rst

  /// @brief Factory method: create options for minimal output.
  /// @return Options with only ATOM records and minimal header information.
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

  /// @brief Factory method: create options to write header records only (no atomic coordinates).
  /// @return Options with atom_records disabled and no END record.
  static PdbWriteOptions headers_only() {
    PdbWriteOptions opt;
    opt.atom_records = false;
    opt.end_record = false;
    return opt;
  }
};

/// @brief Determine the PDB record type (ATOM or HETATM) for a given residue.
///
/// Uses residue flags and entity type to determine appropriate record type:
/// - Residues marked 'H' (het flag) use HETATM.
/// - Standard polymeric residues use ATOM.
/// - Branched, non-polymer, and water residues use HETATM.
///
/// @param res The residue to classify.
/// @return true if HETATM should be used, false for ATOM.
GEMMI_DLL bool use_hetatm(const Residue& res);

/// @brief Write a Structure to PDB format on an output stream.
/// @param st The Structure to serialize.
/// @param os The output stream to write to.
/// @param opt Options controlling record types and formatting.
GEMMI_DLL void write_pdb(const Structure& st, std::ostream& os, PdbWriteOptions opt={});

/// @brief Serialize a Structure to a PDB format string.
/// @param st The Structure to serialize.
/// @param opt Options controlling record types and formatting.
/// @return A std::string containing the complete PDB format output.
GEMMI_DLL std::string make_pdb_string(const Structure& st, PdbWriteOptions opt={});

/// @brief Write minimal PDB output (atomic records only, no headers).
/// @deprecated Use write_pdb(st, os, PdbWriteOptions::minimal()) instead.
/// @param st The Structure to serialize.
/// @param os The output stream to write to.
inline void write_minimal_pdb(const Structure& st, std::ostream& os) {
  write_pdb(st, os, PdbWriteOptions::minimal());
}

/// @brief Generate PDB header records only (no atomic coordinates).
/// @deprecated Use make_pdb_string(st, PdbWriteOptions::headers_only()) instead.
/// @param st The Structure whose headers to serialize.
/// @return A std::string containing PDB header records.
inline std::string make_pdb_headers(const Structure& st) {
  return make_pdb_string(st, PdbWriteOptions::headers_only());
}

} // namespace gemmi

#endif
