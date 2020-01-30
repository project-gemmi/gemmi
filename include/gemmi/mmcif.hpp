// Copyright 2017 Global Phasing Ltd.
//
// Read mmcif (PDBx/mmCIF) file into a Structure from model.hpp.

#ifndef GEMMI_MMCIF_HPP_
#define GEMMI_MMCIF_HPP_

#include <array>
#include <string>
#include <unordered_map>
#include "cifdoc.hpp"
#include "fail.hpp"   // for fail
#include "numb.hpp"   // for as_number
#include "atox.hpp"   // for string_to_int
#include "model.hpp"
#include "entstr.hpp" // for entity_type_from_string, polymer_type_from_string
#include "mmcif_impl.hpp" // for set_cell_from_mmcif

namespace gemmi {

namespace impl {

void copy_int(const cif::Table::Row& row, int n, int& dest) {
  if (row.has2(n))
    dest = cif::as_int(row[n]);
}
void copy_double(const cif::Table::Row& row, int n, double& dest) {
  if (row.has2(n))
    dest = cif::as_number(row[n]);
}
void copy_string(const cif::Table::Row& row, int n, std::string& dest) {
  if (row.has2(n))
    dest = cif::as_string(row[n]);
}

inline std::unordered_map<std::string, std::array<float,6>>
get_anisotropic_u(cif::Block& block) {
  cif::Table aniso_tab = block.find("_atom_site_anisotrop.",
                                    {"id", "U[1][1]", "U[2][2]", "U[3][3]",
                                     "U[1][2]", "U[1][3]", "U[2][3]"});
  std::unordered_map<std::string, std::array<float,6>> aniso_map;
  for (auto ani : aniso_tab)
    aniso_map.emplace(ani[0], std::array<float,6>{{
                                (float) cif::as_number(ani[1]),
                                (float) cif::as_number(ani[2]),
                                (float) cif::as_number(ani[3]),
                                (float) cif::as_number(ani[4]),
                                (float) cif::as_number(ani[5]),
                                (float) cif::as_number(ani[6])}});
  return aniso_map;
}

inline
std::vector<std::string> transform_tags(std::string mstr, std::string vstr) {
  return {mstr + "[1][1]", mstr + "[1][2]", mstr + "[1][3]", vstr + "[1]",
          mstr + "[2][1]", mstr + "[2][2]", mstr + "[2][3]", vstr + "[2]",
          mstr + "[3][1]", mstr + "[3][2]", mstr + "[3][3]", vstr + "[3]"};
}

inline Transform get_transform_matrix(const cif::Table::Row& r) {
  Transform t;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      t.mat[i][j] = cif::as_number(r[4*i+j]);
    t.vec.at(i) = cif::as_number(r[4*i+3]);
  }
  return t;
}

inline ResidueId make_resid(const std::string& name,
                            const std::string& seqid,
                            const std::string* icode) {
  ResidueId rid;
  rid.name = name;
  if (icode)
    // the insertion code happens to be always a single letter
    rid.seqid.icode = cif::as_char(*icode, ' ');
  if (!seqid.empty()) {
    // old mmCIF files have auth_seq_id as number + icode (e.g. 15A)
    if (seqid.back() >= 'A') {
      if (rid.seqid.icode == ' ')
        rid.seqid.icode = seqid.back();
      else if (rid.seqid.icode != seqid.back())
        fail("Inconsistent insertion code in " + seqid);
      rid.seqid.num = cif::as_int(seqid.substr(0, seqid.size() - 1));
    } else {
      rid.seqid.num = cif::as_int(seqid, Residue::OptionalNum::None);
    }
  }
  return rid;
}

inline std::vector<Helix> read_helices(cif::Block& block) {
  std::vector<Helix> helices;
  for (const auto row : block.find("_struct_conf.", {
        "conf_type_id",                                          // 0
        "beg_auth_asym_id", "beg_label_comp_id",                 // 1-2
        "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",             // 3-4
        "end_auth_asym_id", "end_label_comp_id",                 // 5-6
        "end_auth_seq_id", "?pdbx_end_PDB_ins_code",             // 7-8
        "?pdbx_PDB_helix_class", "?pdbx_PDB_helix_length"})) {   // 9-10
    if (alpha_up(row.str(0)[0]) != 'H')
      continue;
    Helix h;
    h.start.chain_name = row.str(1);
    h.start.res_id = make_resid(row.str(2), row.str(3), row.ptr_at(4));
    h.end.chain_name = row.str(5);
    h.end.res_id = make_resid(row.str(6), row.str(7), row.ptr_at(8));
    if (row.has(9))
      h.set_helix_class_as_int(cif::as_int(row[9], -1));
    if (row.has(10))
      h.length = cif::as_int(row[10], -1);
    helices.push_back(h);
  }
  return helices;
}

inline std::vector<Sheet> read_sheets(cif::Block& block) {
  std::vector<Sheet> sheets;
  for (const std::string& sheet_id : block.find_values("_struct_sheet.id"))
    sheets.emplace_back(sheet_id);
  for (const auto row : block.find("_struct_sheet_range.", {
        "sheet_id", "id",                                        // 0-1
        "beg_auth_asym_id", "beg_label_comp_id",                 // 2-3
        "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",             // 4-5
        "end_auth_asym_id", "end_label_comp_id",                 // 6-7
        "end_auth_seq_id", "?pdbx_end_PDB_ins_code"})) {         // 8-9
    Sheet& sheet = impl::find_or_add(sheets, row.str(0));
    sheet.strands.emplace_back();
    Sheet::Strand& strand = sheet.strands.back();
    strand.name = row.str(1);
    strand.start.chain_name = row.str(2);
    strand.start.res_id = make_resid(row.str(3), row.str(4), row.ptr_at(5));
    strand.end.chain_name = row.str(6);
    strand.end.res_id = make_resid(row.str(7), row.str(8), row.ptr_at(9));
  }

  // below we assume that range_id_1 is the strand preceding range_id_2
  for (const auto row : block.find("_struct_sheet_order.", {
        "sheet_id", "range_id_2", "sense"}))
    if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
      if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1)))
        switch (alpha_up(row.str(2)[0])) {
          case 'P': ss->sense = 1; break; // parallel
          case 'A': ss->sense = -1; break; // anti-parallel
        }

  for (const auto row : block.find("_pdbx_struct_sheet_hbond.", {
        "sheet_id", "range_id_2",
        "range_1_auth_asym_id", "range_1_label_comp_id",
        "range_1_auth_seq_id", "?range_1_PDB_ins_code",
        "range_1_label_atom_id",
        "range_2_auth_asym_id", "range_2_label_comp_id",
        "range_2_auth_seq_id", "?range_2_PDB_ins_code",
        "range_2_label_atom_id"}))
    if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
      if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1))) {
        ss->hbond_atom1.chain_name = row.str(2);
        ss->hbond_atom1.res_id = make_resid(row.str(3), row.str(4),
                                            row.ptr_at(5));
        ss->hbond_atom1.atom_name = row.str(6);
        ss->hbond_atom2.chain_name = row.str(7);
        ss->hbond_atom2.res_id = make_resid(row.str(8), row.str(9),
                                            row.ptr_at(10));
        ss->hbond_atom2.atom_name = row.str(11);
      }

  return sheets;
}

inline void read_connectivity(cif::Block& block, Structure& st) {
  // label_ identifiers are not sufficient for HOH:
  // waters have null label_seq_id so we need auth_seq_id+icode.
  // And since we need auth_seq_id, we also use auth_asym_id for consistency.
  for (const auto row : block.find("_struct_conn.", {
        "id", "conn_type_id", // 0-1
        "ptnr1_auth_asym_id", "ptnr2_auth_asym_id", // 2-3
        "ptnr1_label_comp_id", "ptnr2_label_comp_id", // 4-5
        "ptnr1_label_atom_id", "ptnr2_label_atom_id", // 6-7
        "?pdbx_ptnr1_label_alt_id", "?pdbx_ptnr2_label_alt_id", // 8-9
        "ptnr1_auth_seq_id", "ptnr2_auth_seq_id", // 10-11
        "?pdbx_ptnr1_PDB_ins_code", "?pdbx_ptnr2_PDB_ins_code", // 12-13
        "?ptnr1_symmetry", "?ptnr2_symmetry", "?pdbx_dist_value", // 14-16
        "?ccp4_link_id"})) {
    Connection c;
    c.name = row.str(0);
    copy_string(row, 17, c.link_id);
    std::string type = row.str(1);
    for (int i = 0; i != Connection::None; ++i)
      if (get_mmcif_connection_type_id(Connection::Type(i)) == type) {
        c.type = Connection::Type(i);
        break;
      }
    if (row.has2(14) && row.has2(15)) {
      c.asu = (row.str(14) == row.str(15) ? Asu::Same : Asu::Different);
    }
    copy_double(row, 16, c.reported_distance);
    for (int i = 0; i < 2; ++i) {
      AtomAddress& a = (i == 0 ? c.partner1 : c.partner2);
      a.chain_name = row.str(2+i);
      if (row.has2(12+i))
        a.res_id.seqid.icode = cif::as_char(row[12+i], ' ');
      a.res_id = make_resid(row.str(4+i), row.str(10+i), row.ptr_at(12+i));
      a.atom_name = row.str(6+i);
      a.altloc = row.has2(8+i) ? cif::as_char(row[8+i], '\0') : '\0';
    }
    st.connections.emplace_back(c);
  }
}

// Operation expression is an item type used for *.oper_expression.
// Here, to keep it simple, we ignore products such as "(2)(3)".
// We parse "3", "1,3,5", "one,two", "(3)", "(a)", "(1-60)", "(2,3-8,XY)", etc
inline std::vector<std::string> parse_operation_expr(const std::string& expr) {
  std::vector<std::string> result;
  std::size_t start = 0;
  std::size_t close_br = std::string::npos;
  if (expr[0] == '(') {
    start = 1;
    close_br = expr.find(')');
  }
  for (;;) {
    std::size_t sep = std::min(expr.find(',', start), close_br);
    std::size_t minus = expr.find('-', start);
    if (minus < sep) {
      int n_min = no_sign_atoi(expr.c_str() + start);
      int n_max = no_sign_atoi(expr.c_str() + minus + 1);
      for (int n = n_min; n <= n_max; ++n)
        result.push_back(std::to_string(n));
    } else {
      result.emplace_back(expr, start, sep - start);
    }
    if (sep == close_br)
      break;
    start = sep + 1;
  }
  return result;
}

inline std::vector<Assembly> read_assemblies(cif::Block& block) {
  std::vector<Assembly> assemblies;
  cif::Table prop_tab = block.find("_pdbx_struct_assembly_prop.",
                                   {"biol_id", "type", "value"});
  cif::Table gen_tab = block.find("_pdbx_struct_assembly_gen.",
                          {"assembly_id", "oper_expression", "asym_id_list"});
  std::vector<Assembly::Oper> oper_list;
  std::vector<std::string> oper_list_tags = transform_tags("matrix", "vector");
  oper_list_tags.emplace_back("id");  // 12
  oper_list_tags.emplace_back("type");  // 13
  for (const auto row : block.find("_pdbx_struct_oper_list.", oper_list_tags)) {
    oper_list.emplace_back();
    oper_list.back().name = row.str(12);
    oper_list.back().type = row.str(13);
    oper_list.back().transform = get_transform_matrix(row);
  }
  for (const auto row : block.find("_pdbx_struct_assembly.", {
        "id", "details", "method_details",
        "oligomeric_details", "oligomeric_count"})) {
    assemblies.emplace_back(row.str(0));
    Assembly& a = assemblies.back();
    std::string detail = row.str(1);
    if (detail == "author_and_software_defined_assembly")
      a.author_determined = a.software_determined = true;
    else if (detail == "author_defined_assembly")
      a.author_determined = true;
    else if (detail == "software_defined_assembly")
      a.software_determined = true;
    else if (detail == "complete icosahedral assembly")
      a.special_kind = Assembly::SpecialKind::CompleteIcosahedral;
    else if (detail == "representative helical assembly")
      a.special_kind = Assembly::SpecialKind::RepresentativeHelical;
    else if (detail == "complete point assembly")
      a.special_kind = Assembly::SpecialKind::CompletePoint;

    if (!a.author_determined && !a.software_determined &&
        a.special_kind == Assembly::SpecialKind::NA && !detail.empty()) {
      assemblies.pop_back();
      continue;
    }
    if (a.software_determined && !cif::is_null(row[2]))
      a.software_name = row.str(2);  // method_details
    a.oligomeric_details = row.str(3);
    a.oligomeric_count = cif::as_int(row[4], 0);
    for (const auto row_p : prop_tab)
      if (row_p.str(0) == a.name) {
        std::string type = row_p.str(1);
        double value = cif::as_number(row_p[2]);
        if (type == "ABSA (A^2)")
          a.absa = value;
        else if (type == "SSA (A^2)")
          a.ssa = value;
        else if (type == "MORE")
          a.more = value;
      }
    for (const auto row_g : gen_tab)
      if (row_g.str(0) == a.name) {
        a.generators.emplace_back();
        Assembly::Gen& gen = a.generators.back();
        split_str_into(row_g.str(2), ',', gen.subchains);
        for (const std::string& name : parse_operation_expr(row_g.str(1)))
          if (const Assembly::Oper* oper = impl::find_or_null(oper_list, name))
            gen.opers.push_back(*oper);
      }
  }
  return assemblies;
}

inline void fill_residue_entity_type(Structure& st) {
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (ResidueSpan& sub : chain.subchains()) {
        EntityType etype = EntityType::Unknown;
        if (const Entity* ent = st.get_entity_of(sub))
          etype = ent->entity_type;
        if (etype == EntityType::Unknown) {
          if (sub[0].is_water())
            etype = EntityType::Water;
          else if (sub.length() > 1)
            etype = EntityType::Polymer;
          else
            etype = EntityType::NonPolymer;
        }
        for (Residue& residue : sub)
          residue.entity_type = etype;
      }
}

DiffractionInfo* find_diffrn(Metadata& meta, const std::string& diffrn_id) {
  for (CrystalInfo& crystal_info : meta.crystals)
    for (DiffractionInfo& diffr_info : crystal_info.diffractions)
      if (diffr_info.id == diffrn_id)
        return &diffr_info;
  return nullptr;
}

inline Structure make_structure_from_block(const cif::Block& block_) {
  using cif::as_number;
  using cif::as_string;
  // find() and Table don't have const variants, but we don't change anything.
  cif::Block& block = const_cast<cif::Block&>(block_);
  Structure st;
  st.input_format = CoorFormat::Mmcif;
  st.name = block.name;
  set_cell_from_mmcif(block, st.cell);
  st.spacegroup_hm = as_string(impl::find_spacegroup_hm_value(block));

  auto add_info = [&](std::string tag) {
    bool first = true;
    for (const std::string& v : block.find_values(tag))
      if (!cif::is_null(v)) {
        if (first)
          st.info[tag] = as_string(v);
        else
          st.info[tag] += "; " + as_string(v);
        first = false;
      }
  };
  add_info("_entry.id");
  add_info("_cell.Z_PDB");
  add_info("_exptl.method");
  add_info("_struct.title");
  // in pdbx/mmcif v5 date_original was replaced with a much longer tag
  std::string old_date_tag = "_database_PDB_rev.date_original";
  std::string new_date_tag
                      = "_pdbx_database_status.recvd_initial_deposition_date";
  add_info(old_date_tag);
  add_info(new_date_tag);
  if (st.info.count(old_date_tag) == 1 && st.info.count(new_date_tag) == 0)
    st.info[new_date_tag] = st.info[old_date_tag];
  add_info("_struct_keywords.pdbx_keywords");
  add_info("_struct_keywords.text");

  for (auto row : block.find("_refine.", {"pdbx_refine_id",           // 0
                                          "?ls_d_res_high",           // 1
                                          "?ls_d_res_low",            // 2
                                          "?ls_percent_reflns_obs",   // 3
                                          "?ls_number_reflns_obs"})) {
    st.meta.refinement.emplace_back();
    RefinementInfo& ref = st.meta.refinement.back();
    ref.id = row.str(0);
    if (row.has(1)) {
      ref.resolution_high = cif::as_number(row[1]);
      if (ref.resolution_high > 0 &&
          (st.resolution == 0 || ref.resolution_high < st.resolution))
        st.resolution = ref.resolution_high;
    }
    copy_double(row, 2, ref.resolution_low);
    copy_double(row, 3, ref.completeness);
    copy_int(row, 4, ref.reflection_count);
  }

  for (auto row : block.find("_exptl.", {"method", "?crystals_number"})) {
    st.meta.experiments.emplace_back();
    st.meta.experiments.back().method = row.str(0);
    copy_int(row, 1, st.meta.experiments.back().number_of_crystals);
  }

  for (auto row : block.find("_exptl_crystal.", {"id", "?description"})) {
    st.meta.crystals.emplace_back();
    st.meta.crystals.back().id = row.str(0);
    copy_string(row, 1, st.meta.crystals.back().description);
  }

  for (auto row : block.find("_diffrn.",
                             {"id", "crystal_id", "?ambient_temp"})) {
    std::string id = row.str(1);
    auto cryst = std::find_if(st.meta.crystals.begin(), st.meta.crystals.end(),
                              [&](const CrystalInfo& c) { return c.id == id; });
    if (cryst != st.meta.crystals.end()) {
      cryst->diffractions.emplace_back();
      cryst->diffractions.back().id = row.str(0);
      copy_double(row, 2, cryst->diffractions.back().temperature);
    }
  }
  for (auto row : block.find("_diffrn_detector.", {"diffrn_id",
                                                   "?pdbx_collection_date",
                                                   "?detector",
                                                   "?type",
                                                   "?details"}))
    if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
      copy_string(row, 1, di->collection_date);
      copy_string(row, 2, di->detector);
      copy_string(row, 3, di->detector_make);
      copy_string(row, 4, di->optics);
    }
  for (auto row : block.find("_diffrn_radiation.",
                             {"diffrn_id",
                              "?pdbx_scattering_type",
                              "?pdbx_monochromatic_or_laue_m_l",
                              "?monochromator"}))
    if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
      copy_string(row, 1, di->scattering_type);
      if (row.has2(2))
        di->mono_or_laue = row.str(2)[0];
      copy_string(row, 3, di->monochromator);
    }
  for (auto row : block.find("_diffrn_source.", {"diffrn_id",
                                                 "?source",
                                                 "?type",
                                                 "?pdbx_synchrotron_site",
                                                 "?pdbx_synchrotron_beamline",
                                                 "?pdbx_wavelength_list"}))
    if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
      copy_string(row, 1, di->source);
      copy_string(row, 2, di->source_type);
      copy_string(row, 3, di->synchrotron);
      copy_string(row, 4, di->beamline);
      copy_string(row, 5, di->wavelengths);
    }

  size_t n = 0;
  for (auto row : block.find("_reflns.", {"pdbx_diffrn_id",        // 0
                                          "?number_obs",           // 1
                                          "?d_resolution_high",    // 2
                                          "?d_resolution_low",     // 3
                                          "?percent_possible_obs", // 4
                                          "?pdbx_redundancy",      // 5
                                          "?pdbx_Rmerge_I_obs",    // 6
                                          "?pdbx_Rsym_value",      // 7
                                          "?pdbx_netI_over_sigmaI"})) {
    // In the case of multiple experiments (_exptl), which is rare,
    // it is not explicit to which experiment which data statistics
    // (_reflns) corresponds to. We assume they are in the same order.
    if (n >= st.meta.experiments.size())
      break;
    ExperimentInfo& exper = st.meta.experiments[n++];
    split_str_into(row.str(0), ',', exper.diffraction_ids);
    copy_int(row, 1, exper.unique_reflections);
    copy_double(row, 2, exper.reflections.resolution_high);
    copy_double(row, 3, exper.reflections.resolution_low);
    copy_double(row, 4, exper.reflections.completeness);
    copy_double(row, 5, exper.reflections.redundancy);
    copy_double(row, 6, exper.reflections.r_merge);
    copy_double(row, 7, exper.reflections.r_sym);
    copy_double(row, 8, exper.reflections.mean_I_over_sigma);
  }

  for (auto row : block.find("_software.", {"name",
                                            "?classification",
                                            "?version",
                                            "?date",
                                            "?pdbx_ordinal"})) {
    st.meta.software.emplace_back();
    SoftwareItem& item = st.meta.software.back();
    item.name = row.str(0);
    if (row.has2(1))
      item.classification = software_classification_from_string(row.str(1));
    copy_string(row, 2, item.version);
    copy_string(row, 3, item.date);
    copy_int(row, 4, item.pdbx_ordinal);
  }

  std::vector<std::string> ncs_oper_tags = transform_tags("matrix", "vector");
  ncs_oper_tags.emplace_back("id");  // 12
  ncs_oper_tags.emplace_back("?code");  // 13
  cif::Table ncs_oper = block.find("_struct_ncs_oper.", ncs_oper_tags);
  for (auto op : ncs_oper) {
    bool given = op.has(13) && op.str(13) == "given";
    Transform tr = get_transform_matrix(op);
    if (!tr.is_identity())
      st.ncs.push_back({op.str(12), given, tr});
  }

  // PDBx/mmcif spec defines both _database_PDB_matrix.scale* and
  // _atom_sites.fract_transf_* as equivalent of pdb SCALE, but the former
  // is not used, so we ignore it.
  cif::Table fract_tv = block.find("_atom_sites.fract_transf_",
                                   transform_tags("matrix", "vector"));
  if (fract_tv.length() > 0) {
    Transform fract = get_transform_matrix(fract_tv[0]);
    st.cell.set_matrices_from_fract(fract);
  }

  // We read/write origx just for completeness, it's not used anywhere.
  cif::Table origx_tv = block.find("_database_PDB_matrix.",
                                   transform_tags("origx", "origx_vector"));
  if (origx_tv.length() > 0) {
    st.has_origx = true;
    st.origx = get_transform_matrix(origx_tv[0]);
  }

  auto aniso_map = get_anisotropic_u(block);

  // atom list
  enum { kId=0, kGroupPdb, kSymbol, kLabelAtomId, kAltId, kLabelCompId,
         kLabelAsymId, kLabelSeqId, kInsCode, kX, kY, kZ, kOcc, kBiso, kCharge,
         kAuthSeqId, kAuthCompId, kAuthAsymId, kAuthAtomId, kModelNum };
  cif::Table atom_table = block.find("_atom_site.",
                                     {"id",
                                      "?group_PDB",
                                      "type_symbol",
                                      "label_atom_id",
                                      "label_alt_id",
                                      "label_comp_id",
                                      "label_asym_id",
                                      "?label_seq_id",
                                      "?pdbx_PDB_ins_code",
                                      "Cartn_x",
                                      "Cartn_y",
                                      "Cartn_z",
                                      "occupancy",
                                      "B_iso_or_equiv",
                                      "?pdbx_formal_charge",
                                      "auth_seq_id",
                                      "?auth_comp_id",
                                      "?auth_asym_id",
                                      "?auth_atom_id",
                                      "?pdbx_PDB_model_num"});
  const int kCompId = atom_table.has_column(kAuthCompId) ? kAuthCompId
                                                         : kLabelCompId;
  const int kAsymId = atom_table.has_column(kAuthAsymId) ? kAuthAsymId
                                                         : kLabelAsymId;
  const int kAtomId = atom_table.has_column(kAuthAtomId) ? kAuthAtomId
                                                         : kLabelAtomId;
  Model *model = nullptr;
  Chain *chain = nullptr;
  Residue *resi = nullptr;
  if (atom_table.length() != 0) {
    if (atom_table.has_column(kModelNum))
      model = &st.find_or_add_model(atom_table[0].str(kModelNum));
    else
      model = &st.find_or_add_model("1");
  }
  for (auto row : atom_table) {
    if (row.has(kModelNum) && row[kModelNum] != model->name) {
      model = &st.find_or_add_model(row.str(kModelNum));
      chain = nullptr;
    }
    if (!chain || as_string(row[kAsymId]) != chain->name) {
      model->chains.emplace_back(as_string(row[kAsymId]));
      chain = &model->chains.back();
      resi = nullptr;
    }
    ResidueId rid = make_resid(as_string(row[kCompId]),
                               as_string(row[kAuthSeqId]),
                               row.has(kInsCode) ? &row[kInsCode] : nullptr);
    if (!resi || !resi->matches(rid)) {
      resi = chain->find_or_add_residue(rid);
      if (resi->atoms.empty()) {
        if (row.has2(kLabelSeqId))
          resi->label_seq = cif::as_int(row[kLabelSeqId]);
        resi->subchain = row.str(kLabelAsymId);
        // don't check if group_PDB is consistent, it's not that important
        if (row.has2(kGroupPdb))
          for (int i = 0; i < 2; ++i) { // first character could be " or '
            const char c = alpha_up(row[kGroupPdb][i]);
            if (c == 'A' || c == 'H' || c == '\0')
              resi->het_flag = c;
          }
      }
    } else if (resi->seqid != rid.seqid) {
      fail("Inconsistent sequence ID: " + resi->str() + " / " + rid.str());
    }
    Atom atom;
    atom.name = as_string(row[kAtomId]);
    // altloc is always a single letter (not guaranteed by the mmCIF spec)
    atom.altloc = cif::as_char(row[kAltId], '\0');
    atom.charge = row.has2(kCharge) ? cif::as_int(row[kCharge]) : 0;
    atom.element = gemmi::Element(as_string(row[kSymbol]));
    // According to the PDBx/mmCIF spec _atom_site.id can be a string,
    // but in all the files it is a serial number; its value is not essential,
    // so we just ignore non-integer ids.
    atom.serial = string_to_int(row[kId], false);
    atom.pos.x = cif::as_number(row[kX]);
    atom.pos.y = cif::as_number(row[kY]);
    atom.pos.z = cif::as_number(row[kZ]);
    atom.occ = (float) cif::as_number(row[kOcc], 1.0);
    atom.b_iso = (float) cif::as_number(row[kBiso], 50.0);

    if (!aniso_map.empty()) {
      auto ani = aniso_map.find(row[kId]);
      if (ani != aniso_map.end()) {
        atom.u11 = ani->second[0];
        atom.u22 = ani->second[1];
        atom.u33 = ani->second[2];
        atom.u12 = ani->second[3];
        atom.u13 = ani->second[4];
        atom.u23 = ani->second[5];
      }
    }
    resi->atoms.emplace_back(atom);
  }

  cif::Table polymer_types = block.find("_entity_poly.", {"entity_id", "type"});
  for (auto row : block.find("_entity.", {"id", "type"})) {
    Entity ent(row.str(0));
    ent.entity_type = entity_type_from_string(row.str(1));
    ent.polymer_type = PolymerType::Unknown;
    if (polymer_types.ok()) {
      try {
        std::string poly_type = polymer_types.find_row(ent.name).str(1);
        ent.polymer_type = polymer_type_from_string(poly_type);
      } catch (std::runtime_error&) {}
    }
    st.entities.push_back(ent);
  }

  for (auto row : block.find("_entity_poly_seq.",
                             {"entity_id", "num", "mon_id"}))
    if (Entity* ent = st.get_entity(row.str(0))) {
      // According to the spec, num must be >= 1.
      int pos = cif::as_int(row[1], 0) - 1;
      if (pos == (int) ent->full_sequence.size())
        ent->full_sequence.push_back(row.str(2));
      else if (pos >= 0 && pos < (int) ent->full_sequence.size())
        ent->full_sequence[pos] += "," + row.str(2);
    }
  for (auto row : block.find("_struct_asym.", {"id", "entity_id"}))
    if (Entity* ent = st.get_entity(row.str(1)))
      ent->subchains.push_back(row.str(0));

  fill_residue_entity_type(st);

  st.setup_cell_images();

  // CISPEP
  for (auto row : block.find("_struct_mon_prot_cis.",
                             {"pdbx_PDB_model_num", "auth_asym_id",  // 0-1
                              "auth_seq_id", "?pdbx_PDB_ins_code",   // 2-3
                              "?label_comp_id", "?auth_comp_id"})) { // 4-5
    if (row.has2(0) && row.has2(1) && row.has2(2) &&
        (row.has2(4) || row.has2(5)))
      if (Model* mdl = st.find_model(row[0])) {
        std::string comp = row.str(row.has2(4) ? 4 : 5);
        ResidueId rid = make_resid(comp, row.str(2), row.ptr_at(3));
        if (Residue* res = mdl->find_residue(row[1], rid))
          res->is_cis = true;
      }
  }

  st.helices = read_helices(block);
  st.sheets = read_sheets(block);
  read_connectivity(block, st);
  st.assemblies = read_assemblies(block);

  return st;
}

} // namespace impl

inline Structure make_structure_from_block(const cif::Block& block) {
  return impl::make_structure_from_block(block);
}


} // namespace gemmi
#endif
