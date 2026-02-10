//! @file
//! @brief Generate Refmac intermediate (prepared) files crd and rst.

// Copyright 2022 Global Phasing Ltd.
//
// Generate Refmac intermediate (prepared) files crd and rst

#ifndef GEMMI_CRD_HPP_
#define GEMMI_CRD_HPP_

#include "topo.hpp"      // for Topo

namespace gemmi {

//! @brief Prepare structure for CRD file generation.
//! @param st Structure to prepare
GEMMI_DLL void setup_for_crd(Structure& st);

//! @brief Add automatic links between residues.
//! @param model Model to process
//! @param st Structure containing the model
//! @param monlib Monomer library
GEMMI_DLL void add_automatic_links(Model& model, Structure& st, const MonLib& monlib);

//! @brief Add dictionary blocks for residues to CIF document.
//! @param doc CIF document to add blocks to
//! @param resnames List of residue names
//! @param topo Topology information
//! @param monlib Monomer library
GEMMI_DLL void add_dictionary_blocks(cif::Document& doc, const std::vector<std::string>& resnames,
                                     const Topo& topo, const MonLib& monlib);

//! @brief Prepare Refmac CRD file from structure.
//! @param st Structure (not const because this function calls shorten_ccd_codes(st))
//! @param topo Topology information
//! @param monlib Monomer library
//! @param h_change Hydrogen change policy
//! @return CIF document containing CRD data
GEMMI_DLL cif::Document prepare_refmac_crd(Structure& st, const Topo& topo,
                                           const MonLib& monlib, HydrogenChange h_change);

} // namespace gemmi
#endif
