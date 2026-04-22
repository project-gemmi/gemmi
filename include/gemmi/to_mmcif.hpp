// Copyright 2017 Global Phasing Ltd.
//
// Create cif::Document (for PDBx/mmCIF file) from Structure.

/// @file
/// @brief Serialization of Structure to mmCIF/PDBx format.

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include "model.hpp"
#include "cifdoc.hpp"

namespace gemmi {

/// @brief Control which mmCIF data groups to write to the output.
///
/// This struct uses bit-fields to selectively enable/disable writing of
/// individual data blocks and categories in mmCIF output. Each member
/// corresponds to a major category in the PDBx/mmCIF dictionary.
struct MmcifOutputGroups {
  /// @brief Write atom site coordinates (_atom_site category).
  bool atoms:1;
  /// @brief Write data block name (_entry.id).
  bool block_name:1;
  /// @brief Write entry metadata (database code, PDB ID).
  bool entry:1;
  /// @brief Write database status information.
  bool database_status:1;
  /// @brief Write author information (_audit_author category).
  bool author:1;
  /// @brief Write unit cell parameters (_cell category).
  bool cell:1;
  /// @brief Write crystal symmetry (_symmetry category).
  bool symmetry:1;
  /// @brief Write entity definitions (_entity category).
  bool entity:1;
  /// @brief Write polymer entity descriptions (_entity_poly category).
  bool entity_poly:1;
  /// @brief Write structural reference information (_struct_ref category).
  bool struct_ref:1;
  /// @brief Write chemical component information (_chem_comp category).
  bool chem_comp:1;
  /// @brief Write experimental method details (_exptl category).
  bool exptl:1;
  /// @brief Write diffraction data (_diffrn category).
  bool diffrn:1;
  /// @brief Write reflection statistics (_reflns category).
  bool reflns:1;
  /// @brief Write refinement details (_refine category).
  bool refine:1;
  /// @brief Write title and keywords.
  bool title_keywords:1;
  /// @brief Write NCS (non-crystallographic symmetry) operations.
  bool ncs:1;
  /// @brief Write structural asymmetric unit information (_struct_asym category).
  bool struct_asym:1;
  /// @brief Write crystal to fractional coordinate transformation (ORIGX).
  bool origx:1;
  /// @brief Write secondary structure definitions (_struct_conf category).
  bool struct_conf:1;
  /// @brief Write beta sheet information (_struct_sheet category).
  bool struct_sheet:1;
  /// @brief Write biological assembly metadata (_struct_biol category).
  bool struct_biol:1;
  /// @brief Write pdbx_struct_assembly transformations.
  bool assembly:1;
  /// @brief Write chemical bond connectivity (_struct_conn category).
  bool conn:1;
  /// @brief Write cis-peptide information (_struct_mon_prot_cis category).
  bool cis:1;
  /// @brief Write modified residues (_pdbx_chem_comp_atom_feature for non-standard residues).
  bool modres:1;
  /// @brief Write binding site information (_struct_site category).
  bool struct_site:1;
  /// @brief Write matrix scale transformations (_atom_sites_fract_tran category).
  bool scale:1;
  /// @brief Write atom type scattering info (_atom_type category).
  bool atom_type:1;
  /// @brief Write polymer sequence information (_entity_poly_seq category).
  bool entity_poly_seq:1;
  /// @brief Write TLS tensor information (thermal tensor correction).
  bool tls:1;
  /// @brief Write software used in structure processing (_software category).
  bool software:1;
  /// @brief @brief Write _atom_site.group_PDB field (ATOM/HETATM classification).
  bool group_pdb:1;
  /// @brief Write authentic atom and component IDs (_atom_site.auth_atom_id and auth_comp_id).
  bool auth_all:1;

  /// @brief Constructor initializing all groups to the same value.
  /// @param all Set all groups to true (write everything) or false (write nothing).
  explicit MmcifOutputGroups(bool all)
    : atoms(all), block_name(all), entry(all), database_status(all),
      author(all), cell(all), symmetry(all), entity(all), entity_poly(all),
      struct_ref(all), chem_comp(all), exptl(all), diffrn(all),
      reflns(all), refine(all), title_keywords(all), ncs(all),
      struct_asym(all), origx(all), struct_conf(all), struct_sheet(all),
      struct_biol(all), assembly(all), conn(all), cis(all), modres(all),
      struct_site(all),
      scale(all), atom_type(all), entity_poly_seq(all), tls(all),
      software(all), group_pdb(all), auth_all(false) {}
};

/// @brief Populate an existing mmCIF block with data from a Structure.
/// @param st The Structure to serialize.
/// @param block The mmCIF block to populate (typically obtained from make_mmcif_block()).
/// @param groups Selectively enable/disable categories to write.
GEMMI_DLL void update_mmcif_block(const Structure& st, cif::Block& block,
                                  MmcifOutputGroups groups=MmcifOutputGroups(true));

/// @brief Create a complete mmCIF/PDBx document from a Structure.
/// @param st The Structure to serialize.
/// @param groups Selectively enable/disable categories to write.
/// @return A cif::Document containing a single block with all structural information.
GEMMI_DLL cif::Document make_mmcif_document(const Structure& st,
                                            MmcifOutputGroups groups=MmcifOutputGroups(true));

/// @brief Create an mmCIF/PDBx block from a Structure.
/// @param st The Structure to serialize.
/// @param groups Selectively enable/disable categories to write.
/// @return A cif::Block containing all selected categories.
GEMMI_DLL cif::Block make_mmcif_block(const Structure& st,
                                      MmcifOutputGroups groups=MmcifOutputGroups(true));

/// @brief Create an mmCIF block with structural metadata and headers only.
/// @param st The Structure whose metadata to serialize.
/// @return A cif::Block containing metadata without atomic coordinates (_atom_site).
GEMMI_DLL cif::Block make_mmcif_headers(const Structure& st);

/// @brief Add minimal mmCIF category items required for valid output.
/// @param st The Structure to serialize.
/// @param block The mmCIF block to populate with essential items.
GEMMI_DLL void add_minimal_mmcif_data(const Structure& st, cif::Block& block);

/// @brief Write NCS (non-crystallographic symmetry) operations to an mmCIF block.
/// @param st The Structure containing NCS operations.
/// @param block The mmCIF block to populate.
GEMMI_DLL void write_ncs_oper(const Structure& st, cif::Block& block);

/// @brief Write chemical connectivity information to an mmCIF block.
/// @param st The Structure containing connection definitions.
/// @param block The mmCIF block to populate with _struct_conn items.
GEMMI_DLL void write_struct_conn(const Structure& st, cif::Block& block);

/// @brief Write cis-peptide information to an mmCIF block.
/// @param st The Structure containing cis-peptide annotations.
/// @param block The mmCIF block to populate with _struct_mon_prot_cis items.
GEMMI_DLL void write_cispeps(const Structure& st, cif::Block& block);

/// @brief Write binding site information to an mmCIF block.
/// @param st The Structure containing structural site definitions.
/// @param block The mmCIF block to populate with _struct_site items.
GEMMI_DLL void write_struct_sites(const Structure& st, cif::Block& block);

} // namespace gemmi

#endif
