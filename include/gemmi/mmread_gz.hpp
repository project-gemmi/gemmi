/// @file
/// @brief Read coordinate files with optional gzip decompression.

// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped coordinate files.

#ifndef GEMMI_MMREAD_GZ_HPP_
#define GEMMI_MMREAD_GZ_HPP_

#include "model.hpp"  // for Structure

namespace gemmi {

namespace cif { struct Document; }

/// Read a Structure from a potentially gzip-compressed coordinate file.
/// Detects gzip format from filename (.gz extension) or file content,
/// decompresses if needed, then parses using the appropriate handler.
/// @param path     Path to coordinate file (may have .gz suffix)
/// @param format   File format (Unknown = auto-detect from filename or content;
///                 Detect = force content-based detection after decompression)
/// @param save_doc Optional pointer to receive the parsed CIF document
/// @return A Structure parsed from the file
/// @throws Throws on I/O or parse errors
GEMMI_DLL Structure read_structure_gz(const std::string& path,
                                      CoorFormat format=CoorFormat::Unknown,
                                      cif::Document* save_doc=nullptr);

/// Read a PDB-format Structure from a potentially gzip-compressed file.
/// @param path    Path to PDB file (may have .gz suffix)
/// @param options Parsing options (max line length, handling of TER records, etc.)
/// @return A Structure parsed from the PDB file
/// @throws Throws on I/O or parse errors
GEMMI_DLL Structure read_pdb_gz(const std::string& path,
                                PdbReadOptions options=PdbReadOptions());

/// Read a chemical component Structure from a potentially gzip-compressed file.
/// @param path     Path to chemical component file (may have .gz suffix)
/// @param save_doc Optional pointer to receive the parsed CIF document
/// @param which    Bitmask of ChemCompModel values to include (default 7 = all)
/// @return A Structure with the requested coordinate models
/// @throws Throws on I/O or parse errors, or if not a valid chemical component file
GEMMI_DLL Structure read_structure_from_chemcomp_gz(const std::string& path,
                                                    cif::Document* save_doc=nullptr,
                                                    int which=7);

/// Detect file format from a filename, accounting for gzip compression.
/// Strips .gz suffix if present before checking file extension.
/// @param path File path or name
/// @return Detected format (Pdb, Mmcif, Mmjson) or Unknown if no match
GEMMI_DLL CoorFormat coor_format_from_ext_gz(const std::string& path);

} // namespace gemmi

#endif
