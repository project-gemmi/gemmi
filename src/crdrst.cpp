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
#include "gemmi/calculate.hpp"  // for calculate_angle, find_best_plane
#include "gemmi/polyheur.hpp"  // for are_connected

#define GEMMI_PROG crdrst
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Restraints;

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
  std::map<std::string, gemmi::ChemLink> links;
  std::map<std::string, gemmi::ChemMod> modifications;

  const gemmi::ChemLink* find_link(const std::string& name) const {
    auto link = links.find(name);
    return link != links.end() ? &link->second : nullptr;
  }
  const gemmi::ChemMod* find_mod(const std::string& name) const {
    auto modif = modifications.find(name);
    return modif != modifications.end() ? &modif->second : nullptr;
  }
};

inline MonLib read_monomers(std::string monomer_dir,
                            const std::set<std::string>& resnames) {
  MonLib monlib;
  assert(!monomer_dir.empty());
  if (monomer_dir.back() != '/' && monomer_dir.back() != '\\')
    monomer_dir += '/';
  monlib.mon_lib_list = cif_read_any(monomer_dir + "list/mon_lib_list.cif");
  for (const std::string& name : resnames) {
    std::string path = monomer_dir;
    path += std::tolower(name[0]);
    path += '/';
    path += name + ".cif";
    auto cc = gemmi::make_chemcomp_from_cif(name, cif_read_any(path));
    monlib.monomers.emplace(name, cc);
  }
  monlib.links = gemmi::read_chemlinks(monlib.mon_lib_list);
  monlib.modifications = gemmi::read_chemmods(monlib.mon_lib_list);
  return monlib;
}

struct Linkage {
  struct ResInfo {
    const gemmi::Residue* res;
    std::string prev_link;
    ResInfo* prev;
    std::vector<std::string> mods;
  };
  struct ChainInfo {
    std::string name;
    std::string entity_id;
    bool polymer;
    std::vector<ResInfo> residues;
  };

  std::vector<ChainInfo> chains;
  std::vector<gemmi::Connection> extra;
};

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
  return "";
}

Linkage::ChainInfo determine_linkage(const gemmi::Chain& chain,
                                    const gemmi::Entity* ent) {
  Linkage::ChainInfo lc;
  lc.residues.reserve(chain.residues.size());
  lc.name = chain.name;
  lc.entity_id = chain.entity_id;
  lc.polymer = ent && ent->entity_type == gemmi::EntityType::Polymer;
  if (lc.polymer) {
    // For now we ignore microheterogeneity.
    Linkage::ResInfo* prev = nullptr;
    for (const gemmi::Residue& res : chain.residues) {
      Linkage::ResInfo lr;
      lr.res = &res;
      lr.prev_link = get_link_type(res, (prev ? prev->res : nullptr),
                                   ent->polymer_type);
      lr.prev = prev;
      lr.mods.push_back(get_modification(chain, res, ent->polymer_type));
      lc.residues.push_back(lr);
      prev = &lc.residues.back();
    }
  } else {
    for (const gemmi::Residue& res : chain.residues)
      lc.residues.push_back({&res, "", nullptr, {}});
  }
  return lc;
}

static cif::Document make_crd(const gemmi::Structure& st, MonLib& monlib,
                              const Linkage& linkage) {
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
  for (const Linkage::ChainInfo& chain_info : linkage.chains) {
    if (!chain_info.polymer)
      continue;
    for (const Linkage::ResInfo& res_info : chain_info.residues) {
      std::string prev = res_info.prev ? res_info.prev->res->seq_id() : "n/a";
      std::string mod = res_info.mods.at(0);
      if (mod.empty())
        mod += '.';
      poly_loop.add_row({res_info.res->name, res_info.res->seq_id(),
                         chain_info.entity_id, res_info.prev_link, prev, mod});
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
  if (const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_name(hm))
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
        std::string auth_seq_id = res.seq_num.str();
        //std::string ins_code(1, res.icode != ' ' ? res.icode : '?');
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

static std::string chirality_to_string(Restraints::Chirality::Type ctype) {
  switch (ctype) {
    case Restraints::Chirality::Positive: return "positive";
    case Restraints::Chirality::Negative: return "negative";
    case Restraints::Chirality::Both: return "both";
    default: return "???";
  }
}

static
void add_restraints(const Restraints& rt,
                    const gemmi::Residue& res, const gemmi::Residue* res2,
                    cif::Loop& restr_loop, int (&counters)[5]) {
  //using gemmi::to_str;
  const auto& to_str = gemmi::to_str_prec<3>; // to make comparisons easier
  const auto& to_str3 = gemmi::to_str_prec<3>;

  for (const Restraints::Bond& bond : rt.bonds)
    if (const gemmi::Atom* at1 = bond.id1.get_from(res, res2))
      if (const gemmi::Atom* at2 = bond.id2.get_from(res, res2)) {
        std::string obs = to_str3(at1->pos.dist(at2->pos));
        obs += " # " + at1->name + " " + at2->name;
        restr_loop.add_row({"BOND", std::to_string(++counters[0]),
                            bond_type_to_string(bond.type), ".",
                            std::to_string(at1->custom),
                            std::to_string(at2->custom),
                            ".", ".",
                            to_str(bond.value), to_str(bond.esd), obs});
      }
  for (const Restraints::Angle& angle : rt.angles)
    if (const gemmi::Atom* at1 = angle.id1.get_from(res, res2))
      if (const gemmi::Atom* at2 = angle.id2.get_from(res, res2))
        if (const gemmi::Atom* at3 = angle.id3.get_from(res, res2)) {
          double a = gemmi::calculate_angle(at1->pos, at2->pos, at3->pos);
          std::string obs = to_str3(gemmi::deg(a));
          obs += " # " + at1->name + " " + at2->name + " " + at3->name;
          restr_loop.add_row({"ANGL", std::to_string(++counters[1]),
                              ".", ".",
                              std::to_string(at1->custom),
                              std::to_string(at2->custom),
                              std::to_string(at3->custom),
                              ".",
                              to_str(angle.value), to_str(angle.esd), obs});
      }
  for (const Restraints::Torsion& tor : rt.torsions)
    if (const gemmi::Atom* at1 = tor.id1.get_from(res, res2))
      if (const gemmi::Atom* at2 = tor.id2.get_from(res, res2))
        if (const gemmi::Atom* at3 = tor.id3.get_from(res, res2))
          if (const gemmi::Atom* at4 = tor.id4.get_from(res, res2)) {
            double d = gemmi::calculate_dihedral(at1->pos, at2->pos,
                                                 at3->pos, at4->pos);
            std::string obs = to_str3(gemmi::deg(d));
            obs += " # " + at1->name + " " + at2->name +
                   " " + at3->name + " " + at4->name;
            restr_loop.add_row({"TORS", std::to_string(++counters[2]),
                                tor.label, std::to_string(tor.period),
                                std::to_string(at1->custom),
                                std::to_string(at2->custom),
                                std::to_string(at3->custom),
                                std::to_string(at4->custom),
                                to_str(tor.value), to_str(tor.esd), obs});
      }
  for (const Restraints::Chirality& chir : rt.chirs)
    if (const gemmi::Atom* at1 = chir.id_ctr.get_from(res, res2))
      if (const gemmi::Atom* at2 = chir.id1.get_from(res, res2))
        if (const gemmi::Atom* at3 = chir.id2.get_from(res, res2))
          if (const gemmi::Atom* at4 = chir.id3.get_from(res, res2)) {
            double vol = rt.chiral_abs_volume(chir);
            double obs_vol = gemmi::calculate_chiral_volume(
                                  at1->pos, at2->pos, at3->pos, at4->pos);
            std::string obs = to_str3(obs_vol)
                              + " # " + at1->name + " " + at2->name
                              + " " + at3->name + " " + at4->name;
            restr_loop.add_row({"CHIR", std::to_string(++counters[3]),
                                chirality_to_string(chir.chir), ".",
                                std::to_string(at1->custom),
                                std::to_string(at2->custom),
                                std::to_string(at3->custom),
                                std::to_string(at4->custom),
                                to_str3(vol), "0.020", obs});
          }
  for (const auto& item : rt.planes) {
    const Restraints::Plane& plane = item.second;
    std::vector<const gemmi::Atom*> atoms;
    for (const Restraints::AtomId& id : plane.ids)
      if (const gemmi::Atom* atom = id.get_from(res, res2))
        atoms.push_back(atom);
    if (atoms.size() < 4)
      continue;
    ++counters[4];
    auto coeff = find_best_plane(atoms);
    for (const gemmi::Atom* atom : atoms) {
      double dist = coeff[0] * atom->pos.x + coeff[1] * atom->pos.y +
                    coeff[2] * atom->pos.z + coeff[3];
      std::string obs = to_str3(dist) + " # " + atom->name;
      restr_loop.add_row({"PLAN", std::to_string(counters[4]), item.first,
                          ".", std::to_string(atom->custom), ".", ".", ".",
                          to_str(plane.esd), ".", obs});
    }
  }
}

static cif::Document make_rst(const Linkage& linkage, MonLib& monlib) {
  cif::Document doc;
  doc.blocks.emplace_back("restraints");
  cif::Block& block = doc.blocks[0];
  cif::Loop& restr_loop = block.init_mmcif_loop("_restr.", {
              "record", "number", "label", "period",
              "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
              "value", "dev", "val_obs"});
  int counters[5] = {0, 0, 0, 0, 0};
  for (const Linkage::ChainInfo& chain_info : linkage.chains) {
    for (const Linkage::ResInfo& ri : chain_info.residues) {
      // write link
      if (ri.prev) {
        const gemmi::ChemLink* link = monlib.find_link(ri.prev_link);
        if (link && !link->rt.empty()) {
          const gemmi::Residue* prev = ri.prev->res;
          std::string comment = "# link " + ri.prev_link + " " +
                                 prev->seq_id() + " " + prev->name + " - " +
                                 ri.res->seq_id() + " " + ri.res->name;
          restr_loop.add_row({comment + "\nLINK", ".", cif::quote(ri.prev_link),
                              ".", ".", ".", ".", ".", ".", ".", "."});
          add_restraints(link->rt, *prev, ri.res, restr_loop, counters);
        }
      }

      gemmi::ChemComp chem_comp = monlib.monomers.at(ri.res->name);
      // apply modifications
      for (const std::string& modif : ri.mods) {
        if (!modif.empty()) {
          if (const gemmi::ChemMod* chem_mod = monlib.find_mod(modif))
            try {
              chem_mod->apply_to(chem_comp);
            } catch(std::runtime_error& e) {
              printf("Failed to apply modification %s to %s: %s\n",
                     chem_mod->id.c_str(), ri.res->name.c_str(), e.what());
            }
          else
            printf("Modification not found: %s\n", modif.c_str());
        }
      }
      // write monomer
      if (!chem_comp.rt.empty()) {
        // comments are added relying on how cif writing works
        std::string res_info = "# monomer " + chain_info.name + " " +
                               ri.res->seq_id() + " " + ri.res->name;
        restr_loop.add_row({res_info + "\nMONO", ".",
                            cif::quote(chem_comp.group.substr(0, 8)),
                            ".", ".", ".", ".", ".", ".", ".", "."});
        add_restraints(chem_comp.rt, *ri.res, nullptr, restr_loop, counters);
      }
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

    MonLib monlib = read_monomers(monomer_dir, resnames);
    int serial = 0;
    for (gemmi::Model& model : st.models)
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues) {
          const gemmi::ChemComp &cc = monlib.monomers.at(res.name);
          cc.reorder_atoms(res.atoms);
          for (gemmi::Atom& atom : res.atoms)
            atom.custom = ++serial;
        }

    Linkage linkage;
    linkage.chains.reserve(st.models[0].chains.size());
    for (gemmi::Chain& chain : st.models[0].chains) {
      const gemmi::Entity* ent = st.get_entity_of(chain);
      linkage.chains.push_back(determine_linkage(chain, ent));
    }
    // add modifications from links
    for (Linkage::ChainInfo& chain_info : linkage.chains)
      for (Linkage::ResInfo& ri : chain_info.residues)
        if (const gemmi::ChemLink* link = monlib.find_link(ri.prev_link)) {
          if (!link->mod[0].empty())
            ri.prev->mods.push_back(link->mod[0]);
          if (!link->mod[1].empty())
            ri.mods.push_back(link->mod[1]);
        }

    cif::Document crd = make_crd(st, monlib, linkage);
    if (p.options[Verbose])
      printf("Writing coordinates to: %s.crd\n", output.c_str());
    write_cif_to_file(crd, output + ".crd", cif::Style::NoBlankLines);
    cif::Document rst = make_rst(linkage, monlib);
    if (p.options[Verbose])
      printf("Writing restraints to: %s.rst\n", output.c_str());
    write_cif_to_file(rst, output + ".rst", cif::Style::NoBlankLines);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (std::out_of_range& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
