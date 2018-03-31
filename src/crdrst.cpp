// Copyright 2017 Global Phasing Ltd.

#include <stdio.h>
#include <cstdlib> // for getenv
#include <cctype>  // for tolower
#include <set>
#include <stdexcept>
#include "input.h"
#include "gemmi/chemcomp.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"
#include "gemmi/pdb.hpp"  // for split_nonpolymers
#include "gemmi/calculate.hpp"  // for calculate_angle

#define GEMMI_PROG crdrst
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { Verbose=3, Monomers };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_BASENAME"
    "\n\nMake intermediate files from one of PDB, mmCIF or mmJSON formats."
    "\n\nOptions:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { 0, 0, 0, 0, 0, 0 }
};


struct MonLib {
  cif::Document mon_lib_list;
  std::map<std::string, gemmi::ChemComp> monomers;
};

inline MonLib read_monomers(std::string monomer_dir,
                            const std::set<std::string>& resnames) {
  assert(!monomer_dir.empty());
  if (monomer_dir.back() != '/' && monomer_dir.back() != '\\')
    monomer_dir += '/';
  cif::Document doc = cif_read_any(monomer_dir + "list/mon_lib_list.cif");
  std::map<std::string, gemmi::ChemComp> monomers;
  for (const std::string& name : resnames) {
    std::string path = monomer_dir;
    path += std::tolower(name[0]);
    path += '/';
    path += name + ".cif";
    monomers.emplace(name,
                     gemmi::make_chemcomp_from_cif(name, cif_read_any(path)));
  }
  return {doc, monomers};
}

static double sq(double x) { return x * x; }

static bool are_connected(const gemmi::Residue& r1, const gemmi::Residue& r2,
                          gemmi::PolymerType ptype) {
  using gemmi::PolymerType;
  if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD) {
    // similar to has_peptide_bond_to()
    const gemmi::Atom* a1 = r1.get_c();
    const gemmi::Atom* a2 = r2.get_n();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.341 * 1.5);
  }
  if (ptype == PolymerType::Dna || ptype == PolymerType::Rna) {
    const gemmi::Atom* a1 = r1.get_o3prim();
    const gemmi::Atom* a2 = r2.get_p();
    return a1 && a2 && a1->pos.dist_sq(a2->pos) < sq(1.6 * 1.5);
  }
  return false;
}

static std::string get_link_type(const gemmi::Residue& res,
                                 const gemmi::Residue* prev,
                                 gemmi::PolymerType ptype) {
  using gemmi::PolymerType;
  if (!prev)
    return ".";
  if (!are_connected(*prev, res, ptype))
    return "gap";
  if (ptype == PolymerType::PeptideL || ptype == PolymerType::PeptideD) {
    std::string link = res.is_cis ? "CIS" : "TRANS";
    if (res.name == "PRO") {
      link = "P" + link;
    } else if (false /*check if is mpeptide*/) {
      link = "NM" + link;
    }
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

static cif::Document make_crd(const gemmi::Structure& st, MonLib& monlib) {
  using gemmi::to_str;
  cif::Document crd;
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
    // For now it won't work with microheterogeneity.
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
  for (const gemmi::Model& model : st.models) {
    for (const gemmi::Chain& chain : model.chains) {
      for (const gemmi::Residue& res : chain.residues) {
        //std::string label_seq = res.label_seq.str();
        std::string auth_seq_id = res.seq_num.str();
        //std::string ins_code(1, res.icode ? res.icode : '?');
        gemmi::ChemComp& cc = monlib.monomers.at(res.name);
        for (const gemmi::Atom& a : res.atoms) {
          vv.emplace_back("ATOM");
          vv.emplace_back(std::to_string(a.custom));
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
          vv.emplace_back(cc.get_atom(a.name).chem_type); // label_chem_id
        }
      }
    }
  }
  return crd;
}

static std::string bond_type_to_string(gemmi::ChemComp::BondType btype) {
  switch (btype) {
    case gemmi::ChemComp::Single: return "single";
    case gemmi::ChemComp::Double: return "double";
    case gemmi::ChemComp::Triple: return "triple";
    case gemmi::ChemComp::Aromatic: return "aromatic";
    case gemmi::ChemComp::Deloc: return "deloc";
    case gemmi::ChemComp::Metal: return "metal";
    default: return "???";
  }
}

static cif::Document make_rst(const gemmi::Structure& st, MonLib& monlib) {
  using gemmi::to_str;
  cif::Document doc;
  doc.blocks.emplace_back("restraints");
  cif::Block& block = doc.blocks[0];
  cif::Loop& restr_loop = block.init_mmcif_loop("_restr.", {
              "record", "number", "label", "period",
              "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
              "value", "dev", "val_obs"});
  int bond_cnt = 0;
  int angle_cnt = 0;
  for (const gemmi::Chain& chain : st.models[0].chains)
    for (const gemmi::Residue& res : chain.residues) {
      gemmi::ChemComp &cc = monlib.monomers.at(res.name);
      std::string comment = "# monomer " + chain.name + " " +
                            res.seq_id() + " " + res.name;
      restr_loop.add_row({comment + "\nMONO", ".", cif::quote(cc.group), ".",
                          ".", ".", ".", ".", ".", ".", "."});
      for (const gemmi::ChemComp::Bond& bond : cc.bonds)
        if (const gemmi::Atom* at1 = res.find_atom(bond.id1))
          if (const gemmi::Atom* at2 = res.find_atom(bond.id2)) {
            std::string obs = gemmi::to_str_prec<3>(at1->pos.dist(at2->pos));
            obs += " # " + at1->name + " " + at2->name;
            restr_loop.add_row({"BOND", std::to_string(++bond_cnt),
                                bond_type_to_string(bond.type), ".",
                                std::to_string(at1->custom),
                                std::to_string(at2->custom),
                                ".", ".",
                                to_str(bond.value), to_str(bond.esd), obs});
          }
      for (const gemmi::ChemComp::Angle& angle : cc.angles)
        if (const gemmi::Atom* at1 = res.find_atom(angle.id1))
          if (const gemmi::Atom* at2 = res.find_atom(angle.id2))
            if (const gemmi::Atom* at3 = res.find_atom(angle.id3)) {
              double a = gemmi::calculate_angle(at1->pos, at2->pos, at3->pos);
              std::string obs = gemmi::to_str_prec<3>(gemmi::deg(a));
              obs += " # " + at1->name + " " + at2->name + " " + at3->name;
              restr_loop.add_row({"ANGL", std::to_string(++angle_cnt),
                                  ".", ".",
                                  std::to_string(at1->custom),
                                  std::to_string(at2->custom),
                                  std::to_string(at3->custom),
                                  ".",
                                  to_str(angle.value), to_str(angle.esd), obs});
          }
    }
  return doc;
}


int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  if (monomer_dir == nullptr || *monomer_dir == '\0') {
    fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    return 1;
  }
  std::string input = p.nonOption(0);
  std::string output = p.nonOption(1);
  try {
    gemmi::Structure st = read_structure(input);
    if (st.input_format == gemmi::CoorFormat::Pdb)
      gemmi::split_nonpolymers(st);
    if (st.models.empty())
      return 1;
    std::set<std::string> resnames;
    for (const gemmi::Chain& chain : st.models[0].chains)
      for (const gemmi::Residue& res : chain.residues)
        resnames.insert(res.name);

    int serial = 0;
    for (gemmi::Model& model : st.models)
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
          for (gemmi::Atom& atom : res.atoms)
            atom.custom = ++serial;

    MonLib monlib = read_monomers(monomer_dir, resnames);
    cif::Document crd = make_crd(st, monlib);
    if (p.options[Verbose])
      printf("Writing coordinates to: %s.crd\n", output.c_str());
    write_to_file(crd, output + ".crd", cif::Style::NoBlankLines);
    cif::Document rst = make_rst(st, monlib);
    if (p.options[Verbose])
      printf("Writing restraints to: %s.rst\n", output.c_str());
    write_to_file(rst, output + ".rst", cif::Style::NoBlankLines);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
