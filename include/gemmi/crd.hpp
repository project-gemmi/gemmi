// Copyright 2022 Global Phasing Ltd.
//
// Generate Refmac intermediate (prepared) files crd and rst

#ifndef GEMMI_CRD_HPP_
#define GEMMI_CRD_HPP_

#include "topo.hpp"      // for Topo

namespace gemmi {

GEMMI_DLL void setup_for_crd(Structure& st);

GEMMI_DLL void add_automatic_links(Model& model, Structure& st, const MonLib& monlib);

GEMMI_DLL void add_dictionary_blocks(cif::Document& doc, const std::vector<std::string>& resnames,
                                     const Topo& topo, const MonLib& monlib);

// Structure is not const b/c this function calls shorten_ccd_codes(st)
GEMMI_DLL cif::Document prepare_refmac_crd(Structure& st, const Topo& topo,
                                           const MonLib& monlib, HydrogenChange h_change);

} // namespace gemmi
#endif
