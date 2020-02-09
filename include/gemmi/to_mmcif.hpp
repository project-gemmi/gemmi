// Copyright 2017 Global Phasing Ltd.
//
// Structure -> cif::Document -> mmcif (PDBx/mmCIF) file

#ifndef GEMMI_TO_MMCIF_HPP_
#define GEMMI_TO_MMCIF_HPP_

#include "model.hpp"
#include "cifdoc.hpp"

namespace gemmi {

void update_cif_block(const Structure& st, cif::Block& block, bool with_atoms);
cif::Document make_mmcif_document(const Structure& st);
cif::Block make_mmcif_headers(const Structure& st);

// temporarily we use it in crdrst.cpp
namespace impl {
void write_struct_conn(const Structure& st, cif::Block& block);
}

} // namespace gemmi

#ifdef GEMMI_WRITE_IMPLEMENTATION

#include <cassert>
#include <cmath>  // for isnan
#include <string>
#include <utility>  // std::pair
#include "sprintf.hpp"
#include "entstr.hpp" // for entity_type_to_string, polymer_type_to_string
#include "calculate.hpp"  // for count_atom_sites

namespace gemmi {

namespace impl {

inline std::string pdbx_icode(const ResidueId& rid) {
  return std::string(1, rid.seqid.has_icode() ? rid.seqid.icode : '?');
}

inline std::string subchain_or_dot(const Residue& res) {
  return res.subchain.empty() ? "." : cif::quote(res.subchain);
}

inline std::string number_or_dot(double d) {
  return std::isnan(d) ? "." : to_str(d);
}
inline std::string number_or_qmark(double d) {
  return std::isnan(d) ? "?" : to_str(d);
}

// for use with non-negative Metadata fields that use -1 for N/A
inline std::string int_or_dot(int n) {
  return n == -1 ? "." : std::to_string(n);
}
inline std::string int_or_qmark(int n) {
  return n == -1 ? "?" : std::to_string(n);
}

inline std::string string_or_dot(const std::string& s) {
  return s.empty() ? "." : cif::quote(s);
}
inline std::string string_or_qmark(const std::string& s) {
  return s.empty() ? "?" : cif::quote(s);
}


inline void add_cif_atoms(const Structure& st, cif::Block& block) {
  // atom list
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
      "id",
      "type_symbol",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "label_seq_id",
      "pdbx_PDB_ins_code",
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "B_iso_or_equiv",
      "pdbx_formal_charge",
      "auth_seq_id",
      "auth_asym_id",
      "pdbx_PDB_model_num"});
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());
  std::vector<std::pair<int, const Atom*>> aniso;
  int serial = 0;
  for (const Model& model : st.models) {
    for (const Chain& chain : model.chains) {
      for (const Residue& res : chain.residues) {
        std::string label_seq_id = res.label_seq.str('.');
        std::string auth_seq_id = res.seqid.num.str();
        for (const Atom& a : res.atoms) {
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.element.uname());
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc_or('.'));
          vv.emplace_back(res.name);
          vv.emplace_back(subchain_or_dot(res));
          vv.emplace_back(label_seq_id);
          vv.emplace_back(pdbx_icode(res));
          vv.emplace_back(to_str(a.pos.x));
          vv.emplace_back(to_str(a.pos.y));
          vv.emplace_back(to_str(a.pos.z));
          vv.emplace_back(to_str(a.occ));
          vv.emplace_back(to_str(a.b_iso));
          vv.emplace_back(a.charge == 0 ? "?" : std::to_string(a.charge));
          vv.emplace_back(auth_seq_id);
          vv.emplace_back(cif::quote(chain.name));
          vv.emplace_back(model.name);
          if (a.u11 != 0.f)
            aniso.emplace_back(serial, &a);
        }
      }
    }
  }
  if (aniso.empty()) {
    block.find_mmcif_category("_atom_site_anisotrop.").erase();
  } else {
    cif::Loop& aniso_loop = block.init_mmcif_loop("_atom_site_anisotrop.", {
                                    "id", "U[1][1]", "U[2][2]", "U[3][3]",
                                    "U[1][2]", "U[1][3]", "U[2][3]"});
    std::vector<std::string>& aniso_val = aniso_loop.values;
    aniso_val.reserve(aniso_loop.tags.size() * aniso.size());
    for (const auto& a : aniso) {
      aniso_val.emplace_back(std::to_string(a.first));
      aniso_val.emplace_back(to_str(a.second->u11));
      aniso_val.emplace_back(to_str(a.second->u22));
      aniso_val.emplace_back(to_str(a.second->u33));
      aniso_val.emplace_back(to_str(a.second->u12));
      aniso_val.emplace_back(to_str(a.second->u13));
      aniso_val.emplace_back(to_str(a.second->u23));
    }
  }
}

// the names are: monomeric, dimeric, ...meric, 21-meric, 22-meric, ...
int xmeric_to_number(const std::string& oligomeric) {
  static const char names[20][10] = {
    "mono", "di", "tri", "tetra", "penta",
    "hexa", "hepta", "octa", "nona", "deca",
    "undeca", "dodeca", "trideca", "tetradeca", "pentadeca",
    "hexadeca", "heptadeca", "octadeca", "nonadeca", "eicosa"
  };
  size_t len = oligomeric.length();
  const char* p = oligomeric.c_str();
  for (int i = 0; i != 20; ++i)
    if (len == std::strlen(names[i]) + 5 && strncmp(p, names[i], len-5) == 0)
      return i + 1;
  return no_sign_atoi(p);
}

void write_assemblies(const Structure& st, cif::Block& block) {
  block.items.reserve(block.items.size() + 4); // avoid re-allocation
  cif::Loop& a_loop = block.init_mmcif_loop("_pdbx_struct_assembly.",
      {"id", "details", "method_details",
       "oligomeric_details", "oligomeric_count"});
  cif::Loop& prop_loop = block.init_mmcif_loop("_pdbx_struct_assembly_prop.",
      {"biol_id", "type", "value"});
  cif::Loop& gen_loop = block.init_mmcif_loop("_pdbx_struct_assembly_gen.",
      {"assembly_id", "oper_expression", "asym_id_list"});
  cif::Loop& oper_loop = block.init_mmcif_loop("_pdbx_struct_oper_list.",
      {"id", "type",
       "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
       "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
       "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
  std::vector<const Assembly::Oper*> distinct_oper;
  for (const Assembly& as : st.assemblies) {
    std::string how_defined = "?";
    if (as.author_determined && as.software_determined)
      how_defined = "author_and_software_defined_assembly";
    else if (as.author_determined)
      how_defined = "author_defined_assembly";
    else if (as.software_determined)
      how_defined = "software_defined_assembly";
    else if (as.special_kind == Assembly::SpecialKind::CompleteIcosahedral)
      how_defined = "complete icosahedral assembly";
    else if (as.special_kind == Assembly::SpecialKind::RepresentativeHelical)
      how_defined = "representative helical assembly";
    else if (as.special_kind == Assembly::SpecialKind::CompletePoint)
      how_defined = "complete point assembly";
    std::string oligomer = to_lower(as.oligomeric_details);
    int nmer = as.oligomeric_count != 0 ? as.oligomeric_count
                                        : xmeric_to_number(oligomer);
    // _pdbx_struct_assembly
    a_loop.add_row({as.name,
                    how_defined,
                    impl::string_or_qmark(as.software_name),
                    cif::quote(oligomer),
                    nmer == 0 ? "?" : std::to_string(nmer)});

    // _pdbx_struct_assembly_prop
    if (!std::isnan(as.absa))
      prop_loop.add_row({as.name, "'ABSA (A^2)'", to_str(as.absa)});
    if (!std::isnan(as.ssa))
      prop_loop.add_row({as.name, "'SSA (A^2)'", to_str(as.ssa)});
    if (!std::isnan(as.more))
      prop_loop.add_row({as.name, "MORE", to_str(as.more)});

    // _pdbx_struct_assembly_gen and _pdbx_struct_oper_list
    for (const Assembly::Gen& gen : as.generators) {
      std::string subchain_str;
      for (const std::string& name : gen.subchains)
        string_append_sep(subchain_str, ',', name);
      if (subchain_str.empty()) // chain names to subchain names
        for (const Chain& chain : st.models[0].chains)
          if (in_vector(chain.name, gen.chains))
            for (const auto& sub : chain.subchains())
              string_append_sep(subchain_str, ',', sub.front().subchain);
      std::string oper_str;
      for (const Assembly::Oper& oper : gen.opers) {
        size_t k = 0;
        for (; k != distinct_oper.size(); ++k)
          if (distinct_oper[k]->transform.approx(oper.transform, 1e-9))
            break;
        string_append_sep(oper_str, ',', std::to_string(k+1));
        if (k != distinct_oper.size())
          continue;
        distinct_oper.emplace_back(&oper);
        oper_loop.values.emplace_back(std::to_string(k+1));
        if (!oper.type.empty()) {
          oper_loop.values.emplace_back(cif::quote(oper.type));
        } else if (oper.transform.is_identity()) {
          oper_loop.values.emplace_back("'identity operation'");
        } else if (as.author_determined || as.software_determined) {
          oper_loop.values.emplace_back("'crystal symmetry operation'");
        } else {
          oper_loop.values.emplace_back(".");
        }
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j)
            oper_loop.values.emplace_back(to_str(oper.transform.mat[i][j]));
          oper_loop.values.emplace_back(to_str(oper.transform.vec.at(i)));
        }
      }
      gen_loop.add_row({as.name,
                        oper_str.empty() ? "." : oper_str,
                        subchain_str.empty() ? "?" : subchain_str});
    }
  }
}

void write_struct_conn(const Structure& st, cif::Block& block) {
  // example:
  // disulf1 disulf A CYS 3  SG ? 3 ? 1_555 A CYS 18 SG ? 18 ?  1_555 ? 2.045
  cif::Loop& conn_loop = block.init_mmcif_loop("_struct_conn.",
      {"id", "conn_type_id",
       "ptnr1_auth_asym_id", "ptnr1_label_asym_id", "ptnr1_label_comp_id",
       "ptnr1_label_seq_id", "ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id",
       "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_symmetry",
       "ptnr2_auth_asym_id", "ptnr2_label_asym_id", "ptnr2_label_comp_id",
       "ptnr2_label_seq_id", "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id",
       "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code", "ptnr2_symmetry",
       "details", "pdbx_dist_value", "ccp4_link_id"});
  for (const Connection& con : st.connections) {
    const_CRA cra1 = st.models[0].find_cra(con.partner1);
    const_CRA cra2 = st.models[0].find_cra(con.partner2);
    if (!cra1.atom || !cra2.atom)
      continue;
    SymImage im = st.cell.find_nearest_image(cra1.atom->pos,
                                             cra2.atom->pos, con.asu);
    conn_loop.add_row({
        con.name,                                  // id
        get_mmcif_connection_type_id(con.type),    // conn_type_id
        cra1.chain->name,                          // ptnr1_auth_asym_id
        subchain_or_dot(*cra1.residue),            // ptnr1_label_asym_id
        cra1.residue->name,                        // ptnr1_label_comp_id
        cra1.residue->label_seq.str('.'),          // ptnr1_label_seq_id
        cra1.atom->name,                           // ptnr1_label_atom_id
        std::string(1, cra1.atom->altloc_or('?')), // pdbx_ptnr1_label_alt_id
        cra1.residue->seqid.num.str(),             // ptnr1_auth_seq_id
        pdbx_icode(con.partner1.res_id),           // ptnr1_PDB_ins_code
        "1_555",                                   // ptnr1_symmetry
        cra2.chain->name,                          // ptnr2_auth_asym_id
        subchain_or_dot(*cra2.residue),            // ptnr2_label_asym_id
        cra2.residue->name,                        // ptnr2_label_comp_id
        cra2.residue->label_seq.str('.'),          // ptnr2_label_seq_id
        cra2.atom->name,                           // ptnr2_label_atom_id
        std::string(1, cra2.atom->altloc_or('?')), // pdbx_ptnr2_label_alt_id
        cra2.residue->seqid.num.str(),             // ptnr2_auth_seq_id
        pdbx_icode(con.partner2.res_id),           // ptnr2_PDB_ins_code
        im.pdb_symbol(true),                       // ptnr2_symmetry
        "?",                                       // details
        to_str(im.dist()),                         // pdbx_dist_value
        impl::string_or_qmark(con.link_id)         // ccp4_link_id
    });
  }
}

} // namespace impl

void update_cif_block(const Structure& st, cif::Block& block, bool with_atoms) {
  using std::to_string;
  if (st.models.empty())
    return;
  block.name = st.name;
  auto e_id = st.info.find("_entry.id");
  std::string id = impl::string_or_dot(e_id != st.info.end() ? e_id->second : st.name);
  block.set_pair("_entry.id", id);
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end() && !initial_date->second.empty()) {
    block.set_pair("_pdbx_database_status.entry_id", id);
    block.set_pair(initial_date->first, initial_date->second);
  }

  // unit cell and symmetry
  block.set_pair("_cell.entry_id", id);
  block.set_pair("_cell.length_a",    to_str(st.cell.a));
  block.set_pair("_cell.length_b",    to_str(st.cell.b));
  block.set_pair("_cell.length_c",    to_str(st.cell.c));
  block.set_pair("_cell.angle_alpha", to_str(st.cell.alpha));
  block.set_pair("_cell.angle_beta",  to_str(st.cell.beta));
  block.set_pair("_cell.angle_gamma", to_str(st.cell.gamma));
  auto z_pdb = st.info.find("_cell.Z_PDB");
  if (z_pdb != st.info.end())
    block.set_pair(z_pdb->first, z_pdb->second);
  block.set_pair("_symmetry.entry_id", id);
  block.set_pair("_symmetry.space_group_name_H-M",
                 cif::quote(st.spacegroup_hm));
  if (const SpaceGroup* sg = st.find_spacegroup())
    block.set_pair("_symmetry.Int_Tables_number", to_string(sg->number));

  // _entity
  cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
  for (const Entity& ent : st.entities)
    entity_loop.add_row({ent.name, entity_type_to_string(ent.entity_type)});

  // _entity_poly
  cif::Loop& ent_poly_loop = block.init_mmcif_loop("_entity_poly.", {"entity_id", "type"});
  for (const Entity& ent : st.entities)
    if (ent.entity_type == EntityType::Polymer)
      ent_poly_loop.add_row({ent.name, polymer_type_to_qstring(ent.polymer_type)});

  // _exptl
  if (!st.meta.experiments.empty()) {
    cif::Loop& loop = block.init_mmcif_loop("_exptl.",
                                            {"entry_id", "method", "crystals_number"});
    for (const ExperimentInfo& exper : st.meta.experiments)
      loop.add_row({id, cif::quote(exper.method),
                    impl::int_or_qmark(exper.number_of_crystals)});
  } else {
    auto exptl_method = st.info.find("_exptl.method");
    if (exptl_method != st.info.end()) {
      cif::Loop& loop = block.init_mmcif_loop("_exptl.", {"entry_id", "method"});
      for (const std::string& m : gemmi::split_str(exptl_method->second, "; "))
        loop.add_row({id, cif::quote(m)});
    }
  }

  // _exptl_crystal
  if (!st.meta.crystals.empty()) {
    cif::Loop& loop = block.init_mmcif_loop("_exptl_crystal.",
                                            {"id", "description"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      loop.add_row({cryst.id, impl::string_or_qmark(cryst.description)});
  }

  // _exptl_crystal_grow
  if (std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
            [](const CrystalInfo& c) { return !c.ph_range.empty() || !std::isnan(c.ph); })) {
    cif::Loop& grow_loop = block.init_mmcif_loop("_exptl_crystal_grow.",
                                                 {"crystal_id", "pH", "pdbx_pH_range"});
    for (const CrystalInfo& crystal : st.meta.crystals)
      grow_loop.add_row({cif::quote(crystal.id),
                         impl::number_or_qmark(crystal.ph),
                         impl::string_or_qmark(crystal.ph_range)});
  }

  // _diffrn
  if (std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
                  [](const CrystalInfo& c) { return !c.diffractions.empty(); })) {

    cif::Loop& loop = block.init_mmcif_loop("_diffrn.", {"id", "crystal_id", "ambient_temp"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        loop.add_row({diffr.id, cryst.id, impl::number_or_qmark(diffr.temperature)});
    // _diffrn_detector
    cif::Loop& det_loop = block.init_mmcif_loop("_diffrn_detector.",
                                                {"diffrn_id",
                                                 "pdbx_collection_date",
                                                 "detector",
                                                 "type",
                                                 "details"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        det_loop.add_row({diffr.id,
                          impl::string_or_qmark(diffr.collection_date),
                          impl::string_or_qmark(diffr.detector),
                          impl::string_or_qmark(diffr.detector_make),
                          impl::string_or_qmark(diffr.optics)});

    // _diffrn_radiation
    cif::Loop& rad_loop = block.init_mmcif_loop("_diffrn_radiation.",
                                                {"diffrn_id",
                                                 "pdbx_scattering_type",
                                                 "pdbx_monochromatic_or_laue_m_l",
                                                 "monochromator"});
    for (const CrystalInfo& cryst : st.meta.crystals)
      for (const DiffractionInfo& diffr : cryst.diffractions)
        rad_loop.add_row({diffr.id,
                          impl::string_or_qmark(diffr.scattering_type),
                          std::string(1, diffr.mono_or_laue ? diffr.mono_or_laue : '?'),
                          impl::string_or_qmark(diffr.monochromator)});
    // _diffrn_source
    cif::Loop& source_loop = block.init_mmcif_loop("_diffrn_source.",
                                                   {"diffrn_id",
                                                    "source",
                                                    "type",
                                                    "pdbx_synchrotron_site",
                                                    "pdbx_synchrotron_beamline",
                                                    "pdbx_wavelength_list"});
    for (const CrystalInfo& crystal : st.meta.crystals)
      for (const DiffractionInfo& diffr : crystal.diffractions)
        source_loop.add_row({diffr.id,
                             impl::string_or_qmark(diffr.source),
                             impl::string_or_qmark(diffr.source_type),
                             impl::string_or_qmark(diffr.synchrotron),
                             impl::string_or_qmark(diffr.beamline),
                             impl::string_or_qmark(diffr.wavelengths)});
  }

  // _reflns
  if (!st.meta.experiments.empty()) {
    cif::Loop& loop = block.init_mmcif_loop("_reflns.", {
        "entry_id",
        "pdbx_ordinal",
        "pdbx_diffrn_id",
        "number_obs",
        "d_resolution_high",
        "d_resolution_low",
        "percent_possible_obs",
        "pdbx_redundancy",
        "pdbx_Rmerge_I_obs",
        "pdbx_Rsym_value",
        "pdbx_netI_over_sigmaI",
        /*"B_iso_Wilson_estimate"*/});
    int n = 0;
    for (const ExperimentInfo& exper : st.meta.experiments)
      loop.add_row({id,
                    std::to_string(++n),
                    impl::string_or_dot(join_str(exper.diffraction_ids, ",")),
                    impl::int_or_qmark(exper.unique_reflections),
                    impl::number_or_qmark(exper.reflections.resolution_high),
                    impl::number_or_qmark(exper.reflections.resolution_low),
                    impl::number_or_qmark(exper.reflections.completeness),
                    impl::number_or_qmark(exper.reflections.redundancy),
                    impl::number_or_qmark(exper.reflections.r_merge),
                    impl::number_or_qmark(exper.reflections.r_sym),
                    impl::number_or_qmark(exper.reflections.mean_I_over_sigma),
                    /*impl::number_or_qmark(exper.b_wilson)*/});
    // _reflns_shell
    cif::Loop& shell_loop = block.init_mmcif_loop("_reflns_shell.", {
        "pdbx_ordinal",
        "pdbx_diffrn_id",
        "d_res_high",
        "d_res_low",
        "percent_possible_all",
        "pdbx_redundancy",
        "Rmerge_I_obs",
        "pdbx_Rsym_value",
        "meanI_over_sigI_obs"});
    n = 0;
    for (const ExperimentInfo& exper : st.meta.experiments) {
      std::string diffrn_id =
        impl::string_or_dot(join_str(exper.diffraction_ids, ","));
      for (const ReflectionsInfo& shell : exper.shells)
        shell_loop.add_row({std::to_string(++n),
                            diffrn_id,
                            impl::number_or_qmark(shell.resolution_high),
                            impl::number_or_qmark(shell.resolution_low),
                            impl::number_or_qmark(shell.completeness),
                            impl::number_or_qmark(shell.redundancy),
                            impl::number_or_qmark(shell.r_merge),
                            impl::number_or_qmark(shell.r_sym),
                            impl::number_or_qmark(shell.mean_I_over_sigma)});
    }
  }

  // _refine
  if (!st.meta.refinement.empty()) {
    block.items.reserve(block.items.size() + 4);
    cif::Loop& loop = block.init_mmcif_loop("_refine.", {
        "entry_id",
        "pdbx_refine_id",
        "ls_d_res_high",
        "ls_d_res_low",
        "ls_percent_reflns_obs",
        "ls_number_reflns_obs"});
    cif::Loop& analyze_loop = block.init_mmcif_loop("_refine_analyze.", {
        "entry_id",
        "pdbx_refine_id",
        "Luzzati_coordinate_error_obs"});
    cif::Loop& restr_loop = block.init_mmcif_loop("_refine_ls_restr.", {
        "pdbx_refine_id", "type",
        "number", "weight", "pdbx_restraint_function", "dev_ideal"});
    cif::Loop& shell_loop = block.init_mmcif_loop("_refine_ls_shell.", {
        "pdbx_refine_id",
        "d_res_high",
        "d_res_low",
        "percent_reflns_obs",
        "number_reflns_obs",
        "number_reflns_R_free",
        "R_factor_obs",
        "R_factor_R_work",
        "R_factor_R_free"});
    for (size_t i = 0; i != st.meta.refinement.size(); ++i) {
      const RefinementInfo& ref = st.meta.refinement[i];
      loop.values.push_back(id);
      loop.values.push_back(cif::quote(ref.id));
      loop.values.push_back(impl::number_or_dot(ref.resolution_high));
      loop.values.push_back(impl::number_or_dot(ref.resolution_low));
      loop.values.push_back(impl::number_or_dot(ref.completeness));
      loop.values.push_back(impl::int_or_dot(ref.reflection_count));
      auto add = [&](const std::string& tag, const std::string& val) {
        if (i == 0)
          loop.tags.push_back("_refine." + tag);
        loop.values.push_back(val);
      };
      if (st.meta.has(&RefinementInfo::rfree_set_count))
        add("ls_number_reflns_R_free", impl::int_or_dot(ref.rfree_set_count));
      if (st.meta.has(&RefinementInfo::r_all))
        add("ls_R_factor_obs", impl::number_or_qmark(ref.r_all));
      if (st.meta.has(&RefinementInfo::r_work))
        add("ls_R_factor_R_work", impl::number_or_qmark(ref.r_work));
      if (st.meta.has(&RefinementInfo::r_free))
        add("ls_R_factor_R_free", impl::number_or_qmark(ref.r_free));
      if (st.meta.has(&RefinementInfo::cross_validation_method))
        add("pdbx_ls_cross_valid_method",
            impl::string_or_qmark(ref.cross_validation_method));
      if (st.meta.has(&RefinementInfo::rfree_selection_method))
        add("pdbx_R_Free_selection_details",
            impl::string_or_qmark(ref.rfree_selection_method));
      if (st.meta.has(&RefinementInfo::mean_b))
        add("B_iso_mean", impl::number_or_qmark(ref.mean_b));
      if (st.meta.has(&RefinementInfo::aniso_b)) {
        if (i == 0)
          for (const char* index : {"[1][1]", "[2][2]", "[3][3]",
                                    "[1][2]", "[1][3]", "[2][3]"})
            loop.tags.push_back(std::string("_refine.aniso_B") + index);
        const Mat33& t = ref.aniso_b;
        for (double d : {t[0][0], t[1][1], t[2][2], t[0][1], t[0][2], t[1][2]})
          loop.values.push_back(impl::number_or_qmark(d));
      }
      if (st.meta.has(&RefinementInfo::dpi_blow_r))
        add("pdbx_overall_SU_R_Blow_DPI",
            impl::number_or_qmark(ref.dpi_blow_r));
      if (st.meta.has(&RefinementInfo::dpi_blow_rfree))
        add("pdbx_overall_SU_R_free_Blow_DPI",
            impl::number_or_qmark(ref.dpi_blow_rfree));
      if (st.meta.has(&RefinementInfo::dpi_cruickshank_r))
        add("overall_SU_R_Cruickshank_DPI",
            impl::number_or_qmark(ref.dpi_cruickshank_r));
      if (st.meta.has(&RefinementInfo::dpi_cruickshank_rfree))
        add("pdbx_overall_SU_R_free_Cruickshank_DPI",
            impl::number_or_qmark(ref.dpi_cruickshank_rfree));
      if (st.meta.has(&RefinementInfo::cc_fo_fc))
        add("correlation_coeff_Fo_to_Fc", impl::number_or_qmark(ref.cc_fo_fc));
      if (st.meta.has(&RefinementInfo::cc_fo_fc_free))
        add("correlation_coeff_Fo_to_Fc_free",
            impl::number_or_qmark(ref.cc_fo_fc_free));
      if (!st.meta.solved_by.empty())
        add("pdbx_method_to_determine_struct", impl::string_or_qmark(st.meta.solved_by));
      if (!st.meta.starting_model.empty())
        add("pdbx_starting_model", impl::string_or_qmark(st.meta.starting_model));
      if (!std::isnan(ref.luzzati_error))
        analyze_loop.add_row({id,
                              cif::quote(ref.id),
                              impl::number_or_qmark(ref.luzzati_error)});
      for (const RefinementInfo::Restr& restr : ref.restr_stats)
        restr_loop.add_row({cif::quote(ref.id),
                            cif::quote(restr.name),
                            impl::int_or_qmark(restr.count),
                            impl::number_or_qmark(restr.weight),
                            impl::string_or_qmark(restr.function),
                            impl::number_or_qmark(restr.dev_ideal)});
      for (const BasicRefinementInfo& bin : ref.bins)
        shell_loop.add_row({cif::quote(ref.id),
                            impl::number_or_dot(bin.resolution_high),
                            impl::number_or_qmark(bin.resolution_low),
                            impl::number_or_qmark(bin.completeness),
                            impl::int_or_qmark(bin.reflection_count),
                            impl::int_or_qmark(bin.rfree_set_count),
                            impl::number_or_qmark(bin.r_all),
                            impl::number_or_qmark(bin.r_work),
                            impl::number_or_qmark(bin.r_free)});
    }
    assert(loop.values.size() % loop.tags.size() == 0);
  }

  // title, keywords
  auto title = st.info.find("_struct.title");
  if (title != st.info.end()) {
    block.set_pair("_struct.entry_id", id);
    block.set_pair(title->first, cif::quote(title->second));
  }
  auto pdbx_keywords = st.info.find("_struct_keywords.pdbx_keywords");
  auto keywords = st.info.find("_struct_keywords.text");
  if (pdbx_keywords != st.info.end() || keywords != st.info.end())
    block.set_pair("_struct_keywords.entry_id", id);
  if (pdbx_keywords != st.info.end())
    block.set_pair(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
  if (keywords != st.info.end())
    block.set_pair(keywords->first, cif::quote(keywords->second));

  // _struct_ncs_oper (MTRIX)
  if (!st.ncs.empty()) {
    cif::Loop& ncs_oper = block.init_mmcif_loop("_struct_ncs_oper.",
        {"id", "code",
         "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
         "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
         "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
    for (const NcsOp& op : st.ncs) {
      ncs_oper.values.emplace_back(op.id);
      ncs_oper.values.emplace_back(op.given ? "given" : "generate");
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
          ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
        ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
      }
    }
  }

  // _struct_asym
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const Chain& chain : st.models[0].chains)
    for (ConstResidueSpan& sub : chain.subchains())
      if (!sub.subchain_id().empty()) {
        const Entity* ent = st.get_entity_of(sub);
        asym_loop.add_row({sub.subchain_id(), (ent ? ent->name : "?")});
      }

  // _database_PDB_matrix (ORIGX)
  if (st.has_origx && !st.origx.is_identity()) {
    block.set_pair("_database_PDB_matrix.entry_id", id);
    std::string prefix = "_database_PDB_matrix.origx";
    for (int i = 0; i < 3; ++i) {
      std::string s = "[" + to_string(i+1) + "]";
      block.set_pair(prefix + s + "[1]", to_str(st.origx.mat[i][0]));
      block.set_pair(prefix + s + "[2]", to_str(st.origx.mat[i][1]));
      block.set_pair(prefix + s + "[3]", to_str(st.origx.mat[i][2]));
      block.set_pair(prefix + "_vector" + s, to_str(st.origx.vec.at(i)));
    }
  }

  if (!st.helices.empty()) {
    cif::Loop& struct_conf_loop = block.init_mmcif_loop("_struct_conf.",
        {"conf_type_id", "id",
         "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
         "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
         "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
         "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code",
         "pdbx_PDB_helix_class", "pdbx_PDB_helix_length"});
    int count = 0;
    for (const Helix& helix : st.helices) {
      const_CRA cra1 = st.models[0].find_cra(helix.start);
      const_CRA cra2 = st.models[0].find_cra(helix.end);
      if (!cra1.residue || !cra2.residue)
        continue;
      struct_conf_loop.add_row({
        "HELX_P",                                     // conf_type_id
        "H" + std::to_string(++count),                // id
        cra1.chain->name,                             // beg_auth_asym_id
        impl::subchain_or_dot(*cra1.residue),         // beg_label_asym_id
        cra1.residue->name,                           // beg_label_comp_id
        cra1.residue->label_seq.str(),                // beg_label_seq_id
        cra1.residue->seqid.num.str(),                // beg_auth_seq_id
        impl::pdbx_icode(*cra1.residue),              // beg_PDB_ins_code
        cra2.chain->name,                             // end_auth_asym_id
        impl::subchain_or_dot(*cra2.residue),         // end_label_asym_id
        cra2.residue->name,                           // end_label_comp_id
        cra2.residue->label_seq.str(),                // end_label_seq_id
        cra2.residue->seqid.num.str(),                // end_auth_seq_id
        impl::pdbx_icode(*cra2.residue),              // end_PDB_ins_code
        std::to_string((int)helix.pdb_helix_class),   // pdbx_PDB_helix_class
        impl::int_or_qmark(helix.length)              // pdbx_PDB_helix_length
      });
    }
    if (count != 0)
      block.set_pair("_struct_conf_type.id", "HELX_P");
  }

  // _struct_sheet*
  if (!st.sheets.empty()) {
    cif::Loop& sheet_loop = block.init_mmcif_loop("_struct_sheet.",
                                                  {"id", "number_strands"});
    for (const Sheet& sheet : st.sheets)
      sheet_loop.add_row({sheet.name, std::to_string(sheet.strands.size())});

    cif::Loop& order_loop = block.init_mmcif_loop("_struct_sheet_order.",
                    {"sheet_id", "range_id_1", "range_id_2", "sense"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 1; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        if (strand.sense != 0)
          order_loop.add_row({sheet.name,
                              std::to_string(i), std::to_string(i+1),
                              strand.sense > 0 ? "parallel" : "anti-parallel"});
      }

    cif::Loop& range_loop = block.init_mmcif_loop("_struct_sheet_range.",
        {"sheet_id", "id",
         "beg_auth_asym_id", "beg_label_asym_id", "beg_label_comp_id",
         "beg_label_seq_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code",
         "end_auth_asym_id", "end_label_asym_id", "end_label_comp_id",
         "end_label_seq_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 0; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        const_CRA cra1 = st.models[0].find_cra(strand.start);
        const_CRA cra2 = st.models[0].find_cra(strand.end);
        if (!cra1.residue || !cra2.residue)
          continue;
        range_loop.add_row({
          sheet.name,                                 // sheet_id
          std::to_string(i+1),                        // id
          cra1.chain->name,                           // beg_auth_asym_id
          impl::subchain_or_dot(*cra1.residue),       // beg_label_asym_id
          cra1.residue->name,                         // beg_label_comp_id
          cra1.residue->label_seq.str(),              // beg_label_seq_id
          cra1.residue->seqid.num.str(),              // beg_auth_seq_id
          impl::pdbx_icode(*cra1.residue),            // beg_PDB_ins_code
          cra2.chain->name,                           // end_auth_asym_id
          impl::subchain_or_dot(*cra2.residue),       // end_label_asym_id
          cra2.residue->name,                         // end_label_comp_id
          cra2.residue->label_seq.str(),              // end_label_seq_id
          cra2.residue->seqid.num.str(),              // end_auth_seq_id
          impl::pdbx_icode(*cra2.residue)             // end_PDB_ins_code
        });
    }

    cif::Loop& hbond_loop = block.init_mmcif_loop("_pdbx_struct_sheet_hbond.",
        {"sheet_id", "range_id_1", "range_id_2",
         "range_1_auth_asym_id", "range_1_label_asym_id",
         "range_1_label_comp_id", "range_1_label_seq_id", "range_1_auth_seq_id",
         "range_1_PDB_ins_code", "range_1_label_atom_id",
         "range_2_auth_asym_id", "range_2_label_asym_id",
         "range_2_label_comp_id", "range_2_label_seq_id", "range_2_auth_seq_id",
         "range_2_PDB_ins_code", "range_2_label_atom_id"});
    for (const Sheet& sheet : st.sheets)
      for (size_t i = 1; i < sheet.strands.size(); ++i) {
        const Sheet::Strand& strand = sheet.strands[i];
        if (strand.hbond_atom2.atom_name.empty())
          continue;
        // hbond_atomN is not a full atom "address": altloc is missing
        const_CRA cra1 = st.models[0].find_cra(strand.hbond_atom1);
        const_CRA cra2 = st.models[0].find_cra(strand.hbond_atom2);
        if (!cra1.residue || !cra2.residue)
          continue;
        hbond_loop.add_row({
          sheet.name,                                 // sheet_id
          std::to_string(i),                          // range_id_1
          std::to_string(i+1),                        // range_id_2
          cra1.chain->name,                           // range_1_auth_asym_id
          impl::subchain_or_dot(*cra1.residue),       // range_1_label_asym_id
          cra1.residue->name,                         // range_1_label_comp_id
          cra1.residue->label_seq.str(),              // range_1_label_seq_id
          cra1.residue->seqid.num.str(),              // range_1_auth_seq_id
          impl::pdbx_icode(*cra1.residue),            // range_1_PDB_ins_code
          strand.hbond_atom1.atom_name.c_str(),       // range_1_label_atom_id
          cra2.chain->name,                           // range_2_auth_asym_id
          impl::subchain_or_dot(*cra2.residue),       // range_2_label_asym_id
          cra2.residue->name,                         // range_2_label_comp_id
          cra2.residue->label_seq.str(),              // range_2_label_seq_id
          cra2.residue->seqid.num.str(),              // range_2_auth_seq_id
          impl::pdbx_icode(*cra2.residue),            // range_2_PDB_ins_code
          strand.hbond_atom2.atom_name.c_str()        // range_2_label_atom_id
        });
    }
  }

  // _pdbx_struct_assembly* and _struct_biol (REMARK 300/350)
  if (!st.meta.remark_300_detail.empty()) {
    block.set_pair("_struct_biol.id", "1");
    block.set_pair("_struct_biol.details", cif::quote(st.meta.remark_300_detail));
  }
  if (!st.assemblies.empty())
    impl::write_assemblies(st, block);

  impl::write_struct_conn(st, block);

  // _struct_mon_prot_cis
  cif::Loop& prot_cis_loop = block.init_mmcif_loop("_struct_mon_prot_cis.",
      {"pdbx_id", "pdbx_PDB_model_num", "label_asym_id", "label_seq_id",
       "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
       "label_comp_id", "label_alt_id"});
  for (const Model& model : st.models)
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        if (res.is_cis)
          prot_cis_loop.add_row({to_string(prot_cis_loop.length()+1),
                                 model.name, impl::subchain_or_dot(res),
                                 res.label_seq.str(), chain.name,
                                 res.seqid.num.str(), impl::pdbx_icode(res),
                                 res.name, "."});

  // _atom_sites (SCALE)
  if (st.has_origx || st.cell.explicit_matrices) {
    block.set_pair("_atom_sites.entry_id", id);
    std::string prefix = "_atom_sites.fract_transf_";
    for (int i = 0; i < 3; ++i) {
      std::string idx = "[" + std::to_string(i + 1) + "]";
      const auto& frac = st.cell.frac;
      block.set_pair(prefix + "matrix" + idx + "[1]", to_str(frac.mat[i][0]));
      block.set_pair(prefix + "matrix" + idx + "[2]", to_str(frac.mat[i][1]));
      block.set_pair(prefix + "matrix" + idx + "[3]", to_str(frac.mat[i][2]));
      block.set_pair(prefix + "vector" + idx, to_str(frac.vec.at(i)));
    }
  }

  // SEQRES from PDB doesn't record microheterogeneity, so if the resulting
  // cif has unknown("?") _entity_poly_seq.num, it cannot be trusted.
  cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.",
                                         {"entity_id", "num", "mon_id"});
  for (const Entity& ent : st.entities)
    if (ent.entity_type == EntityType::Polymer)
      for (size_t i = 0; i != ent.full_sequence.size(); ++i) {
        const std::string& mon_ids = ent.full_sequence[i];
        std::string num = std::to_string(i+1);
        size_t start = 0, end;
        while ((end = mon_ids.find(',', start)) != std::string::npos) {
          poly_loop.add_row({ent.name, num, mon_ids.substr(start, end-start)});
          start = end + 1;
        }
        poly_loop.add_row({ent.name, num, mon_ids.substr(start)});
      }

  if (with_atoms)
    impl::add_cif_atoms(st, block);

  if (st.meta.has_tls()) {
    cif::Loop& loop = block.init_mmcif_loop("_pdbx_refine_tls.", {
        "pdbx_refine_id", "id",
        "T[1][1]", "T[2][2]", "T[3][3]", "T[1][2]", "T[1][3]", "T[2][3]",
        "L[1][1]", "L[2][2]", "L[3][3]", "L[1][2]", "L[1][3]", "L[2][3]",
        "S[1][1]", "S[1][2]", "S[1][3]",
        "S[2][1]", "S[2][2]", "S[2][3]",
        "S[3][1]", "S[3][2]", "S[3][3]",
        "origin_x", "origin_y", "origin_z"});
    for (const RefinementInfo& ref : st.meta.refinement)
      for (const TlsGroup& tls : ref.tls_groups) {
        const Mat33& T = tls.T;
        const Mat33& L = tls.L;
        const Mat33& S = tls.S;
        auto q = impl::number_or_qmark;
        loop.add_row({cif::quote(ref.id), impl::string_or_dot(tls.id),
                      q(T[0][0]), q(T[1][1]), q(T[2][2]),
                      q(T[0][1]), q(T[0][2]), q(T[1][2]),
                      q(L[0][0]), q(L[1][1]), q(L[2][2]),
                      q(L[0][1]), q(L[0][2]), q(L[1][2]),
                      q(S[0][0]), q(S[0][1]), q(S[0][2]),
                      q(S[1][0]), q(S[1][1]), q(S[1][2]),
                      q(S[2][0]), q(S[2][1]), q(S[2][2]),
                      q(tls.origin.x), q(tls.origin.y), q(tls.origin.z)});
      }
    cif::Loop& group_loop = block.init_mmcif_loop("_pdbx_refine_tls_group.", {
        "id", "refine_tls_id", "beg_auth_asym_id", "beg_auth_seq_id",
        "end_auth_asym_id", "end_auth_seq_id", "selection_details"});
    int counter = 1;
    for (const RefinementInfo& ref : st.meta.refinement)
      for (const TlsGroup& tls : ref.tls_groups)
        for (const TlsGroup::Selection& sel : tls.selections)
          group_loop.add_row({std::to_string(counter++),
                              impl::string_or_dot(tls.id),
                              impl::string_or_qmark(sel.chain),
                              sel.res_begin.num ? sel.res_begin.str() : "?",
                              impl::string_or_qmark(sel.chain),
                              sel.res_end.num ? sel.res_end.str() : "?",
                              impl::string_or_qmark(sel.details)});
  }

  if (!st.meta.software.empty()) {
    cif::Loop& loop = block.init_mmcif_loop("_software.",
                 {"pdbx_ordinal", "classification", "name", "version", "date"});
    for (const SoftwareItem& item : st.meta.software)
      loop.add_row({
          std::to_string(item.pdbx_ordinal),
          cif::quote(software_classification_to_string(item.classification)),
          cif::quote(item.name),
          impl::string_or_dot(item.version),
          impl::string_or_qmark(item.date)});
  }
}

cif::Document make_mmcif_document(const Structure& st) {
  cif::Document doc;
  doc.blocks.resize(1);
  gemmi::update_cif_block(st, doc.blocks[0], true);
  return doc;
}

cif::Block make_mmcif_headers(const Structure& st) {
  cif::Block block;
  gemmi::update_cif_block(st, block, false);
  return block;
}

} // namespace gemmi
#endif // GEMMI_WRITE_IMPLEMENTATION

#endif
