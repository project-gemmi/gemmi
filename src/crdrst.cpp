// Copyright 2017 Global Phasing Ltd.

#include <stdio.h>
#include "input.h"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"

#define GEMMI_PROG crdrst
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { Verbose=3 };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nMake intermediate files from one of PDB, mmCIF or mmJSON formats."
    "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { 0, 0, 0, 0, 0, 0 }
};

static std::string get_link_type(const gemmi::Residue& res,
                                 const gemmi::Residue* prev,
                                 gemmi::PolymerType ptype) {
  using gemmi::PolymerType;
  if (!prev)
    return ".";
  if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD) {
    std::string link = res.is_cis ? "CIS" : "TRANS";
    if (res.name == "PRO")
      link = "P" + link;
    return link;
  }
  if (ptype == PolymerType::Dna || ptype == PolymerType::Rna)
    return "p";
  return "?";
}

static std::string get_modification(const gemmi::Chain& chain,
                                    const gemmi::Residue& res,
                                    gemmi::PolymerType ptype) {
  using gemmi::PolymerType;
  if (&res == &chain.residues.back())
    return "TERMINUS";
  if (&res == &chain.residues.front()) {
    if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD)
      return "NH3";
    if (ptype == PolymerType::Dna || ptype == PolymerType::Rna)
      return "5*END";
  }
  return ".";
}

static cif::Document make_crd(const gemmi::Structure& st) {
  // consider update_cif_block()
  using gemmi::to_str;
  cif::Document crd;
  if (st.models.empty())
    return crd;
  auto e_id = st.info.find("_entry.id");
  std::string id = (e_id != st.info.end() ? e_id->second : st.name);
  crd.blocks.emplace_back("structure_" + id);
  cif::Block& block = crd.blocks[0];
  auto& items = block.items;

  items.emplace_back("_entry.id", id);
  items.emplace_back("_database_2.code_PDB", id);
  auto keywords = st.info.find("_struct_keywords.pdbx_keywords");
  if (keywords != st.info.end())
    items.emplace_back("_struct_keywords.text", cif::quote(keywords->second));
  auto title = st.info.find("_struct.title");
  if (title != st.info.end())
    items.emplace_back(title->first, cif::quote(title->second));
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end())
    items.emplace_back("_audit.creation_date", initial_date->second);
  items.emplace_back("_software.name", "gemmi");

  items.emplace_back(cif::CommentArg{"############\n"
                                     "## ENTITY ##\n"
                                     "############"});
  cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
  for (const auto& ent : st.entities) {
    entity_loop.values.push_back(ent.first);
    entity_loop.values.push_back(ent.second.type_as_string());
  }
  items.emplace_back(cif::CommentArg{"#####################\n"
                                     "## ENTITY_POLY_SEQ ##\n"
                                     "#####################"});
  cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.", {
              "mon_id", "ccp4_auth_seq_id", "entity_id",
              "ccp4_back_connect_type", "ccp4_num_mon_back", "ccp4_mod_id"});
  for (const gemmi::Chain& chain : st.models[0].chains) {
    const gemmi::Entity* ent = st.get_entity_of(chain);
    if (!ent || ent->entity_type != gemmi::EntityType::Polymer)
      continue;
    const gemmi::Residue* prev = nullptr;
    for (const gemmi::Residue& res : chain.residues) {
      poly_loop.add_row({res.name,
                         res.seq_id(),
                         chain.entity_id,
                         get_link_type(res, prev, ent->polymer_type),
                         prev ? prev->seq_id() : "n/a",
                         get_modification(chain, res, ent->polymer_type)});
      prev = &res;
    }
  }
  items.emplace_back(cif::CommentArg{"##########\n"
                                     "## CELL ##\n"
                                     "##########"});
  items.emplace_back("_cell.entry_id", id);
  items.emplace_back("_cell.length_a",    to_str(st.cell.a));
  items.emplace_back("_cell.length_b",    to_str(st.cell.b));
  items.emplace_back("_cell.length_c",    to_str(st.cell.c));
  items.emplace_back("_cell.angle_alpha", to_str(st.cell.alpha));
  items.emplace_back("_cell.angle_beta",  to_str(st.cell.beta));
  items.emplace_back("_cell.angle_gamma", to_str(st.cell.gamma));
  bool write_fract_matrix = true;
  if (write_fract_matrix) {
    items.emplace_back(cif::CommentArg{"##############################\n"
                                       "## FRACTIONALISATION MATRIX ##\n"
                                       "##############################"});
    std::string prefix = "_atom_sites.fract_transf_";
    for (int i = 0; i < 3; ++i) {
      std::string start = prefix + "matrix[" + std::to_string(i+1) + "][";
      const auto& row = st.cell.frac.mat[i];
      for (int j = 0; j < 3; ++j)
        items.emplace_back(start + std::to_string(j+1) + "]", to_str(row[j]));
    }
    for (int i = 0; i < 3; ++i)
      items.emplace_back(prefix + "vector[" + std::to_string(i + 1) + "]",
                         to_str(st.cell.frac.vec.at(i)));
  }

  items.emplace_back(cif::CommentArg{"##############\n"
                                     "## SYMMETRY ##\n"
                                     "##############"});
  items.emplace_back("_symmetry.entry_id", id);
  items.emplace_back("_symmetry.space_group_name_H-M", cif::quote(st.sg_hm));
  if (const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name(st.sg_hm))
    items.emplace_back("_symmetry.Int_Tables_number",
                       std::to_string(sg->number));
  items.emplace_back(cif::CommentArg{"#################\n"
                                     "## STRUCT_ASYM ##\n"
                                     "#################"});
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const auto& ch : st.models.at(0).chains) {
    asym_loop.values.push_back(ch.name);
    asym_loop.values.push_back(ch.entity_id.empty() ? "?" : ch.entity_id);
  }
  const auto& connections = st.models.at(0).connections;
  if (!connections.empty()) {
    items.emplace_back(cif::CommentArg{"#################\n"
                                       "## STRUCT_CONN ##\n"
                                       "#################"});
    gemmi::impl::write_struct_conn(st, block);
  }

  items.emplace_back(cif::CommentArg{"###############\n"
                                     "## ATOM_SITE ##\n"
                                     "###############"});
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
      "group_PDB",
      "id",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "auth_seq_id",
      //"pdbx_PDB_ins_code",
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "B_iso_or_equiv",
      "type_symbol",
      "calc_flag",
      "label_seg_id",
      "auth_atom_id",
      "label_chem_id"});
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());
  int serial = 0;
  for (const gemmi::Model& model : st.models) {
    for (const gemmi::Chain& chain : model.chains) {
      for (const gemmi::Residue& res : chain.residues) {
        //std::string label_seq = res.label_seq.str();
        std::string auth_seq_id = res.seq_num.str();
        //std::string ins_code(1, res.icode ? res.icode : '?');
        for (const gemmi::Atom& a : res.atoms) {
          vv.emplace_back("ATOM");
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc ? a.altloc : '.');
          vv.emplace_back(res.name);
          vv.emplace_back(chain.name);
          vv.emplace_back(auth_seq_id);
          //vv.emplace_back(ins_code);
          vv.emplace_back(to_str(a.pos.x));
          vv.emplace_back(to_str(a.pos.y));
          vv.emplace_back(to_str(a.pos.z));
          vv.emplace_back(to_str(a.occ));
          vv.emplace_back(to_str(a.b_iso));
          vv.emplace_back(a.element.uname());
          vv.emplace_back("."); // calc_flag
          vv.emplace_back("."); // label_seg_id
          vv.emplace_back(a.name); // again
          vv.emplace_back("?"); // label_chem_id
        }
      }
    }
  }
  return crd;
}


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  std::string input = p.nonOption(0);
  std::string output = p.nonOption(1);
  try {
    gemmi::Structure st = read_structure(input);
    cif::Document doc = make_crd(st);
    write_to_file(doc, output, cif::Style::NoBlankLines);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
