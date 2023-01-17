// Copyright 2022 Global Phasing Ltd.
//
// Generate Refmac intermediate (prepared) files crd and rst

#ifndef GEMMI_CRD_HPP_
#define GEMMI_CRD_HPP_

#include "topo.hpp"      // for Topo
#include "riding_h.hpp"  // for HydrogenChange

namespace gemmi {

void setup_for_crd(Structure& st);

void add_automatic_links(Model& model, Structure& st, const MonLib& monlib);

cif::Block prepare_crd(const Structure& st, const Topo& topo,
                       HydrogenChange h_change, const std::string& info_comment);

cif::Block prepare_rst(const Topo& topo, const MonLib& monlib, const UnitCell& cell);

cif::Document prepare_refmac_crd(const Structure& st, const Topo& topo,
                                 const MonLib& monlib, HydrogenChange h_change);

} // namespace gemmi
#endif
