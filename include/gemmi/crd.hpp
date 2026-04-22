// Copyright 2022 Global Phasing Ltd.
//
// Generate Refmac intermediate (prepared) files crd and rst

/// @file
/// @brief Preparation of molecular structures for Refmac refinement (CRD format generation).

#ifndef GEMMI_CRD_HPP_
#define GEMMI_CRD_HPP_

#include "topo.hpp"      // for Topo

namespace gemmi {

/// @brief Prepare a Structure for CRD file output.
///
/// Sets up entities, assigns subchains, and ensures proper formatting for Refmac input:
/// - Adds entity type and ID information.
/// - Forces subchain assignment.
/// - Adjusts subchain names to Refmac conventions (using '_' as separator).
/// - Normalizes all water residues to "HOH".
/// - Deduplicates and validates entity definitions.
///
/// @param st The Structure to prepare (modified in-place).
void setup_for_crd(Structure& st);

/// @brief Identify and add missing chemical bonds using a monomer library.
///
/// Searches for atoms within bonding distance and, when not already defined
/// as connections in the structure, attempts to match them against the monomer
/// library. Automatically adds Connection records for metals coordinated to O, N, S, B.
///
/// @param model The Model containing atoms to search.
/// @param st The Structure whose connections will be updated.
/// @param monlib The MonLib (monomer library) for matching link definitions.
void add_automatic_links(Model& model, Structure& st, const MonLib& monlib);

/// @brief Add chemical component dictionary blocks to a CIF document.
///
/// Appends _chem_comp blocks for specified residue types extracted from
/// the topology and monomer library. Used to embed component definitions
/// in Refmac CRD files.
///
/// @param doc The cif::Document to append blocks to.
/// @param resnames List of residue names (component IDs) to include.
/// @param topo The Topo containing topology information.
/// @param monlib The MonLib containing chemical component definitions.
void add_dictionary_blocks(cif::Document& doc, const std::vector<std::string>& resnames,
                                     const Topo& topo, const MonLib& monlib);

/// @brief Generate a complete Refmac CRD (coordinate) file as a CIF document.
///
/// Creates a comprehensive CIF document suitable for Refmac input, including:
/// - Structure metadata (entry ID, cell, space group, symmetry).
/// - Entity and sequence definitions.
/// - Atomic coordinates with topology-derived labels.
/// - NCS operations, structural connectivity, cis-peptides.
/// - Optional anisotropic displacement parameters.
/// - Hydrogen handling flags based on h_change parameter.
///
/// The Structure is modified: shorten_ccd_codes() is called to abbreviate
/// component IDs for Refmac compatibility.
///
/// @param st The Structure to serialize (modified: CCD codes shortened).
/// @param topo The Topo with topology and restraint information.
/// @param monlib The MonLib for component definitions and link types.
/// @param h_change How to handle hydrogen atoms (remove, shift, or keep).
/// @return A cif::Document containing the prepared CRD file structure.
GEMMI_DLL cif::Document prepare_refmac_crd(Structure& st, const Topo& topo,
                                           const MonLib& monlib, HydrogenChange h_change);

} // namespace gemmi
#endif
