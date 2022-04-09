// Copyright 2017-2022 Global Phasing Ltd.

#include <stdio.h>
#include <cstdlib>   // for getenv
#include <algorithm> // for count_if
#include <iostream>  // for cerr
#include <stdexcept>
#include "gemmi/chemcomp.hpp"  // for ChemComp
#include "gemmi/to_cif.hpp"    // for write_cif_to_stream
#include "gemmi/fstream.hpp"   // for Ofstream
#include "gemmi/enumstr.hpp"   // for entity_type_to_string
#include "gemmi/to_mmcif.hpp"  // for write_struct_conn
#include "gemmi/sprintf.hpp"   // for to_str, to_str_prec
#include "gemmi/calculate.hpp" // for find_best_plane
#include "gemmi/modify.hpp"    // for remove_hydrogens
#include "gemmi/polyheur.hpp"  // for setup_entities
#include "gemmi/monlib.hpp"    // for MonLib, read_monomer_lib
#include "gemmi/topo.hpp"      // for Topo
#include "gemmi/placeh.hpp"    // for place_hydrogens
#include "gemmi/read_cif.hpp"  // for read_cif_gz
#include "gemmi/read_coor.hpp" // for read_structure_gz

#define GEMMI_PROG prep
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Topo;
using gemmi::Restraints;

namespace {

enum OptionIndex {
  Split=4, Monomers, Libin, NoHydrogens, KeepHydrogens, NoZeroOccRestr
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n " EXE_NAME " --split [options] INPUT_FILE OUTPUT_BASENAME"
    "\n\nPrepare intermediate Refmac files."
    "\nINPUT_FILE can be in PDB, mmCIF or mmJSON format."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Split, 0, "", "split", Arg::None,
    "  --split  \tSplit output into two files: crd and rst." },
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { Libin, 0, "", "libin", Arg::Required,
    "  --libin=CIF  \tCustom additions to the monomer library." },
  { NoHydrogens, 0, "H", "no-hydrogens", Arg::None,
    "  -H, --no-hydrogens  \tRemove or do not add hydrogens." },
  { KeepHydrogens, 0, "", "keep-hydrogens", Arg::None,
    "  --keep-hydrogens  \tPreserve hydrogens from the input file." },
  //{ NoZeroOccRestr, 0, "", "no-zero-occ", Arg::None,
  //  "  --no-zero-occ  \tNo restraints for zero-occupancy atoms." },
  { 0, 0, 0, 0, 0, 0 }
};


bool has_anisou(const gemmi::Model& model) {
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& res : chain.residues)
      for (const gemmi::Atom& a : res.atoms)
        if (a.aniso.nonzero())
          return true;
  return false;
}

std::string refmac_calc_flag(const gemmi::Atom& a) {
  switch (a.calc_flag) {
    case gemmi::CalcFlag::NotSet: return ".";
    // Refmac seems to be using only flags R and M
    case gemmi::CalcFlag::Calculated: return a.is_hydrogen() ? "R" : "M";
    // I think these are not used
    case gemmi::CalcFlag::Determined: return "d";
    case gemmi::CalcFlag::Dummy: return "dum";
  }
  gemmi::unreachable();
}

// Get value of _entity_poly_seq.ccp4_mod_id compatible with makecif.
std::string get_ccp4_mod_id(const std::vector<std::string>& mods) {
  for (const std::string& m : mods)
    if (!gemmi::starts_with(m, "DEL-OXT") &&
        !gemmi::starts_with(m, "DEL-HN") && m != "DEL-NMH")
      return m;
  return ".";
}

cif::Block make_crd(const gemmi::Structure& st, const Topo& topo) {
  using gemmi::to_str;
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
  cif::Block block("structure_" + id);
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
  for (const gemmi::Entity& ent : st.entities)
    entity_loop.add_row({ent.name, entity_type_to_string(ent.entity_type)});
  items.emplace_back(cif::CommentArg{"#####################\n"
                                     "## ENTITY_POLY_SEQ ##\n"
                                     "#####################"});
  cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.", {
              "mon_id", "ccp4_auth_seq_id", "entity_id",
              "ccp4_back_connect_type", "ccp4_num_mon_back", "ccp4_mod_id"});
  for (const Topo::ChainInfo& chain_info : topo.chain_infos) {
    if (!chain_info.polymer)
      continue;
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      const Topo::Link* prev = ri.prev.empty() ? nullptr : &ri.prev[0];
      poly_loop.add_row({ri.res->name,
                         ri.res->seqid.str(),
                         chain_info.entity_id,
                         prev ? prev->link_id : ".",
                         prev ? prev->res1->seqid.str() : "n/a",
                         get_ccp4_mod_id(ri.mods)});
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
  const std::string& hm = st.spacegroup_hm;
  items.emplace_back("_symmetry.space_group_name_H-M", cif::quote(hm));
  if (const gemmi::SpaceGroup* sg = st.find_spacegroup())
    items.emplace_back("_symmetry.Int_Tables_number",
                       std::to_string(sg->number));
  const gemmi::Model& model0 = st.first_model();
  items.emplace_back(cif::CommentArg{"#################\n"
                                     "## STRUCT_ASYM ##\n"
                                     "#################"});
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const gemmi::Chain& chain : model0.chains)
    for (gemmi::ConstResidueSpan sub : chain.subchains())
      if (!sub.subchain_id().empty()) {
        const gemmi::Entity* ent = st.get_entity_of(sub);
        asym_loop.add_row({sub.subchain_id(), (ent ? ent->name : "?")});
      }
  if (!st.connections.empty()) {
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
  bool write_anisou = has_anisou(model0);
  if (write_anisou)
    for (const char* idx : {"[1][1]", "[2][2]", "[3][3]",
                            "[1][2]", "[1][3]", "[2][3]"})
      atom_loop.tags.push_back(std::string("_atom_site.aniso_U") + idx);
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());


  for (const Topo::ChainInfo& chain_info : topo.chain_infos)
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      const gemmi::ChemComp& cc = ri.chemcomp;
      const gemmi::Residue& res = *ri.res;
      std::string auth_seq_id = res.seqid.num.str();
      //std::string ins_code(1, res.icode != ' ' ? res.icode : '?');
      for (const gemmi::Atom& a : res.atoms) {
        vv.emplace_back("ATOM");
        vv.emplace_back(std::to_string(a.serial));
        vv.emplace_back(a.name);
        vv.emplace_back(1, a.altloc ? a.altloc : '.');
        vv.emplace_back(res.name);
        vv.emplace_back(cif::quote(chain_info.chain_ref.name));
        vv.emplace_back(auth_seq_id);
        //vv.emplace_back(ins_code);
        vv.emplace_back(to_str(a.pos.x));
        vv.emplace_back(to_str(a.pos.y));
        vv.emplace_back(to_str(a.pos.z));
        vv.emplace_back(to_str(a.occ));
        vv.emplace_back(to_str(a.b_iso));
        vv.emplace_back(a.element.uname());
        vv.emplace_back(refmac_calc_flag(a));
        vv.emplace_back(1, '.'); // label_seg_id
        vv.emplace_back(a.name); // again
        vv.emplace_back(cc.get_atom(a.name).chem_type); // label_chem_id
        if (write_anisou) {
          if (a.aniso.nonzero()) {
            for (float u : {a.aniso.u11, a.aniso.u22, a.aniso.u33,
                            a.aniso.u12, a.aniso.u13, a.aniso.u23})
              vv.push_back(to_str(u));
          } else {
            vv.resize(vv.size() + 6, ".");
          }
        }
      }
    }
  return block;
}

void add_restraint_row(cif::Loop& restr_loop, const char* record, int counter,
                       const std::string& label, const std::string& period,
                       std::initializer_list<const gemmi::Atom*> atoms,
                       double value, double dev,
                       double value_nucleus, double dev_nucleus,
                       double obs) {
  using namespace gemmi;
  // Ignore restraints with zero-occupancy atoms (NoZeroOccRestr)
  // Perhaps it could be done earlier, in Topo::apply_restraints().
  for (const Atom* a : atoms)
    if (a->occ == 0)
      return;

  auto& values = restr_loop.values;
  auto to_str_dot = [&](double x) { return std::isnan(x) ? "." : to_str_prec<3>(x); };
  values.emplace_back(record);  // record
  values.emplace_back(std::to_string(counter));  // number
  values.emplace_back(label);  // label
  values.emplace_back(period);  // period
  for (const Atom* a : atoms)
    values.emplace_back(std::to_string(a->serial));  // atom_id_i
  for (size_t i = atoms.size(); i < 4; ++i)
    values.emplace_back(".");
  values.emplace_back(to_str_dot(value));  // value
  values.emplace_back(to_str_dot(dev));  // dev
  values.emplace_back(to_str_dot(value_nucleus));  // value_nucleus
  values.emplace_back(to_str_dot(dev_nucleus));  // dev_nucleus
  values.emplace_back(to_str_prec<3>(obs));  // val_obs
  std::string& last = values.back();
  last += " #";
  for (const Atom* a : atoms) {
    last += ' ';
    last += a->name;
  }
}

void add_restraints(const Topo::Rule rule, const Topo& topo,
                    cif::Loop& restr_loop, int (&counters)[6],
                    const gemmi::UnitCell* cell=nullptr) {
  using namespace gemmi;
  if (rule.rkind == Topo::RKind::Bond) {
    const Topo::Bond& t = topo.bonds[rule.index];
    if (cell == nullptr) {  // don't use symmetry
      add_restraint_row(restr_loop, "BOND", ++counters[0],
                        bond_type_to_string(t.restr->type), ".",
                        {t.atoms[0], t.atoms[1]},
                        t.restr->value, t.restr->esd,
                        t.restr->value_nucleus, t.restr->esd_nucleus,
                        t.calculate());
    } else {
      NearestImage im = cell->find_nearest_image(t.atoms[0]->pos, t.atoms[1]->pos,
                                                 Asu::Different);
      add_restraint_row(restr_loop, "BNDS", ++counters[5],
                        im.symmetry_code(true), ".",
                        {t.atoms[0], t.atoms[1]},
                        t.restr->value, t.restr->esd,
                        t.restr->value_nucleus, t.restr->esd_nucleus,
                        std::sqrt(im.dist_sq));
    }
  } else if (rule.rkind == Topo::RKind::Angle) {
    const Topo::Angle& t = topo.angles[rule.index];
    add_restraint_row(restr_loop, "ANGL", ++counters[1], ".", ".",
                      {t.atoms[0], t.atoms[1], t.atoms[2]},
                       t.restr->value, t.restr->esd, NAN, NAN,
                       deg(t.calculate()));
  } else if (rule.rkind == Topo::RKind::Torsion) {
    const Topo::Torsion& t = topo.torsions[rule.index];
    add_restraint_row(restr_loop, "TORS", ++counters[2],
                      t.restr->label, std::to_string(t.restr->period),
                      {t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]},
                      t.restr->value, t.restr->esd, NAN, NAN,
                      deg(t.calculate()));
  } else if (rule.rkind == Topo::RKind::Chirality) {
    const Topo::Chirality& t = topo.chirs[rule.index];
    add_restraint_row(restr_loop, "CHIR", ++counters[3],
                      chirality_to_string(t.restr->sign), ".",
                      {t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]},
                      topo.ideal_chiral_abs_volume(t), 0.02, NAN, NAN,
                      t.calculate());
  } else if (rule.rkind == Topo::RKind::Plane) {
    const Topo::Plane& t = topo.planes[rule.index];
    ++counters[4];
    auto coeff = find_best_plane(t.atoms);
    for (const Atom* atom : t.atoms)
      add_restraint_row(restr_loop, "PLAN", counters[4], t.restr->label, ".",
                        {atom},
                        t.restr->esd, NAN, NAN, NAN,
                        get_distance_from_plane(atom->pos, coeff));
  }
}

cif::Block make_rst(const Topo& topo, const gemmi::MonLib& monlib,
                    const gemmi::UnitCell& cell) {
  cif::Block block("restraints");
  cif::Loop& restr_loop = block.init_mmcif_loop("_restr.", {
              "record", "number", "label", "period",
              "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
              "value", "dev", "value_nucleus", "dev_nucleus", "val_obs"});
  int counters[6] = {0, 0, 0, 0, 0, 0};
  for (const Topo::ChainInfo& chain_info : topo.chain_infos) {
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      // write link
      for (const Topo::Link& prev : ri.prev) {
        const gemmi::ChemLink* link = monlib.get_link(prev.link_id);
        if (link && !prev.link_rules.empty()) {
          std::string comment = " link " + prev.link_id + " " +
                                prev.res1->seqid.str() + " " +
                                prev.res1->name + " - " +
                                ri.res->seqid.str() + " " + ri.res->name;
          restr_loop.add_comment_and_row({comment, "LINK", ".", cif::quote(prev.link_id), ".",
                                          ".", ".", ".", ".", ".", ".", ".", ".", "."});
          for (const Topo::Rule& rule : prev.link_rules)
            add_restraints(rule, topo, restr_loop, counters);
        }
      }
      // write monomer
      if (!ri.monomer_rules.empty()) {
        std::string res_info = " monomer " + chain_info.subchain_name + " " +
                               ri.res->seqid.str() + " " + ri.res->name;
        if (!ri.mods.empty())
          res_info += " modified by " + gemmi::join_str(ri.mods, ", ");

        // need to revisit it later on
        std::string group = cif::quote(ri.chemcomp.group.substr(0, 8));
        if (group == "peptide" || group == "P-peptid" || group == "M-peptid")
          group = "L-peptid";
        else if (group == "NON-POLY")
          group = ".";

        restr_loop.add_comment_and_row({res_info, "MONO", ".", group, ".",
                                        ".", ".", ".", ".", ".", ".", ".", ".", "."});
        for (const Topo::Rule& rule : ri.monomer_rules)
          add_restraints(rule, topo, restr_loop, counters);
      }
    }
  }
  // explicit links
  for (const Topo::Link& extra : topo.extras) {
    if (extra.asu == gemmi::Asu::Different) // symmetry links are left for later
      continue;
    const gemmi::ChemLink* chem_link = monlib.get_link(extra.link_id);
    assert(chem_link);
    std::string comment = " link " + chem_link->id;
    restr_loop.add_comment_and_row({comment, "LINK", ".", cif::quote(chem_link->id), ".",
                                    ".", ".", ".", ".", ".", ".", ".", ".", "."});
    for (const Topo::Rule& rule : extra.link_rules)
      add_restraints(rule, topo, restr_loop, counters);
  }

  // special symmetry related links
  for (const Topo::Link& extra : topo.extras) {
    if (extra.asu != gemmi::Asu::Different)
      continue;
    const gemmi::ChemLink* chem_link = monlib.get_link(extra.link_id);
    assert(chem_link);
    std::string comment = " link (symmetry) " + chem_link->id;
    restr_loop.add_comment_and_row({comment, "LINK", ".", "symmetry", ".",
                                    ".", ".", ".", ".", ".", ".", ".", ".", "."});
    for (const Topo::Rule& rule : extra.link_rules)
      if (rule.rkind == Topo::RKind::Bond)
        add_restraints(rule, topo, restr_loop, counters, &cell);
  }

  return block;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  p.check_exclusive_pair(KeepHydrogens, NoHydrogens);
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  if (monomer_dir == nullptr || *monomer_dir == '\0') {
    fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    return 1;
  }
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  bool verbose = p.options[Verbose];
  try {
    if (verbose)
      printf("Reading %s ...\n", input.c_str());
    gemmi::Structure st = gemmi::read_structure_gz(input, gemmi::CoorFormat::Detect);
    if (st.input_format == gemmi::CoorFormat::Pdb ||
        st.input_format == gemmi::CoorFormat::ChemComp)
      gemmi::setup_entities(st);

    if (st.models.empty()) {
      fprintf(stderr, "No models found in the input file.\n");
      return 1;
    }
    gemmi::Model& model0 = st.models[0];

    std::string libin;
    if (p.options[Libin])
      libin = p.options[Libin].arg;
    if (verbose)
      printf("Reading monomer library...\n");
    gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir,
                                                   model0.get_all_residue_names(),
                                                   gemmi::read_cif_gz,
                                                   libin);

    if (verbose)
      printf("Preparing topology, hydrogens, restraints...\n");
    bool reorder = true;
    bool ignore_unknown_links = false;
    gemmi::HydrogenChange h_change;
    if (p.options[NoHydrogens])
      h_change = gemmi::HydrogenChange::Remove;
    else if (p.options[KeepHydrogens])
      h_change = gemmi::HydrogenChange::NoChange;
    else
      h_change = gemmi::HydrogenChange::ReAddButWater;
    auto topo = gemmi::prepare_topology(st, monlib, 0, h_change, reorder,
                                        &std::cerr, ignore_unknown_links);

    if (verbose)
      printf("Preparing structure data...\n");
    cif::Block crd = make_crd(st, *topo);
    if (p.options[Split])
      output += ".crd";
    if (verbose)
      printf("Writing coordinates to: %s\n", output.c_str());
    gemmi::Ofstream os(output);
    write_cif_block_to_stream(os.ref(), crd, cif::Style::NoBlankLines);

    if (verbose)
      printf("Preparing restraint data...\n");
    cif::Block rst = make_rst(*topo, monlib, st.cell);
    if (p.options[Split])
      output.replace(output.size()-3, 3, "rst");
    if (verbose)
      printf("Writing restraints to: %s\n", output.c_str());
    if (p.options[Split])
      os = gemmi::Ofstream(output);
    else
      os->write("\n\n", 2);
    write_cif_block_to_stream(os.ref(), rst, cif::Style::NoBlankLines);
  } catch (std::exception& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
