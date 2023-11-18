// Copyright 2017 Global Phasing Ltd.
//
// Create cif::Document (for PDBx/mmCIF file) from Structure.

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include "model.hpp"
#include "cifdoc.hpp"

namespace gemmi {

struct MmcifOutputGroups {
  bool atoms:1;
  bool block_name:1;
  bool entry:1;
  bool database_status:1;
  bool author:1;
  bool cell:1;
  bool symmetry:1;
  bool entity:1;
  bool entity_poly:1;
  bool struct_ref:1;
  bool chem_comp:1;
  bool exptl:1;
  bool diffrn:1;
  bool reflns:1;
  bool refine:1;
  bool title_keywords:1;
  bool ncs:1;
  bool struct_asym:1;
  bool origx:1;
  bool struct_conf:1;
  bool struct_sheet:1;
  bool struct_biol:1;
  bool assembly:1;
  bool conn:1;
  bool cis:1;
  bool modres:1;
  bool scale:1;
  bool atom_type:1;
  bool entity_poly_seq:1;
  bool tls:1;
  bool software:1;
  bool group_pdb:1;  // include _atom_site.group_PDB
  bool auth_all:1;   // include _atom_site.auth_atom_id and auth_comp_id

  explicit MmcifOutputGroups(bool all)
    : atoms(all), block_name(all), entry(all), database_status(all),
      author(all), cell(all), symmetry(all), entity(all), entity_poly(all),
      struct_ref(all), chem_comp(all), exptl(all), diffrn(all),
      reflns(all), refine(all), title_keywords(all), ncs(all),
      struct_asym(all), origx(all), struct_conf(all), struct_sheet(all),
      struct_biol(all), assembly(all), conn(all), cis(all), modres(all),
      scale(all), atom_type(all), entity_poly_seq(all), tls(all),
      software(all), group_pdb(all), auth_all(false) {}
};

GEMMI_DLL void update_mmcif_block(const Structure& st, cif::Block& block,
                                  MmcifOutputGroups groups=MmcifOutputGroups(true));
GEMMI_DLL cif::Document make_mmcif_document(const Structure& st,
                                            MmcifOutputGroups groups=MmcifOutputGroups(true));
GEMMI_DLL cif::Block make_mmcif_block(const Structure& st,
                                      MmcifOutputGroups groups=MmcifOutputGroups(true));
GEMMI_DLL cif::Block make_mmcif_headers(const Structure& st);
GEMMI_DLL void add_minimal_mmcif_data(const Structure& st, cif::Block& block);

// temporarily we use it in crd.cpp
GEMMI_DLL void write_struct_conn(const Structure& st, cif::Block& block);
GEMMI_DLL void write_cispeps(const Structure& st, cif::Block& block);

} // namespace gemmi

#endif
