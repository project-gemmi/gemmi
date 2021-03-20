// Copyright 2017 Global Phasing Ltd.

#include <stdio.h>
#include <cstdlib> // for getenv
#include <algorithm> // for count_if
#include <stdexcept>
#include "gemmi/chemcomp.hpp"  // for ChemComp
#include "gemmi/to_cif.hpp"    // for write_cif_to_stream
#include "gemmi/fstream.hpp"   // for Ofstream
#include "gemmi/enumstr.hpp"   // for entity_type_to_string
#include "gemmi/to_mmcif.hpp"  // for write_struct_conn
#include "gemmi/sprintf.hpp"   // for to_str, to_str_prec
#include "gemmi/calculate.hpp" // for find_best_plane
#include "gemmi/polyheur.hpp"  // for remove_hydrogens
#include "gemmi/monlib.hpp"    // for MonLib, read_monomer_lib
#include "gemmi/topo.hpp"      // for Topo
#include "gemmi/placeh.hpp"    // for place_hydrogens
#include "gemmi/read_cif.hpp"  // for read_cif_gz
#include "gemmi/read_coor.hpp" // for read_structure_gz

#define GEMMI_PROG crdrst
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Topo;
using gemmi::Restraints;

namespace {

enum OptionIndex { Monomers=4, NoHydrogens, KeepHydrogens, NoZeroOccRestr };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_BASENAME"
    "\n\nMake intermediate files from one of PDB, mmCIF or mmJSON formats."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { NoHydrogens, 0, "H", "no-hydrogens", Arg::None,
    "  -H, --no-hydrogens  \tRemove or do not add hydrogens." },
  { KeepHydrogens, 0, "", "keep-hydrogens", Arg::None,
    "  --keep-hydrogens  \tPreserve hydrogens from the input file." },
  { NoZeroOccRestr, 0, "", "no-zero-occ", Arg::None,
    "  --no-zero-occ  \tNo restraints for zero-occupancy atoms." },
  { 0, 0, 0, 0, 0, 0 }
};


// Topology: restraints applied to a model
int count_provenance(const std::vector<Topo::Force>& forces, Topo::Provenance p) {
  return std::count_if(forces.begin(), forces.end(),
                       [&](const Topo::Force& f) { return f.provenance == p; });
}

bool has_anisou(const gemmi::Model& model) {
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& res : chain.residues)
      for (const gemmi::Atom& a : res.atoms)
        if (a.aniso.nonzero())
          return true;
  return false;
}

// for compatibility with makecif, not sure what ccp4_mod_id is really used for
std::string get_ccp4_mod_id(const std::vector<std::string>& mods) {
  for (const std::string& m : mods)
    if (m != "AA-STAND" && !gemmi::starts_with(m, "DEL-OXT") &&
        !gemmi::starts_with(m, "DEL-HN") && m != "DEL-NMH")
      return m;
  return ".";
}

cif::Document make_crd(const gemmi::Structure& st,
                       const gemmi::MonLib& monlib,
                       const Topo& topo) {
  using gemmi::to_str;
  cif::Document crd;
  auto e_id = st.info.find("_entry.id");
  std::string id = cif::quote(e_id != st.info.end() ? e_id->second : st.name);
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
      const Topo::ResInfo::Prev* prev = ri.prev.empty() ? nullptr : &ri.prev[0];
      poly_loop.add_row({ri.res->name,
                         ri.res->seqid.str(),
                         chain_info.entity_id,
                         prev ? prev->link : ".",
                         prev ? prev->get(&ri)->res->seqid.str() : "n/a",
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
  for (const gemmi::Chain& chain : model0.chains) {
    for (const gemmi::Residue& res : chain.residues) {
      std::string auth_seq_id = res.seqid.num.str();
      //std::string ins_code(1, res.icode != ' ' ? res.icode : '?');
      const gemmi::ChemComp& cc = monlib.monomers.at(res.name);
      for (const gemmi::Atom& a : res.atoms) {
        vv.emplace_back("ATOM");
        vv.emplace_back(std::to_string(a.serial));
        vv.emplace_back(a.name);
        vv.emplace_back(1, a.altloc ? a.altloc : '.');
        vv.emplace_back(res.name);
        vv.emplace_back(cif::quote(chain.name));
        vv.emplace_back(auth_seq_id);
        //vv.emplace_back(ins_code);
        vv.emplace_back(to_str(a.pos.x));
        vv.emplace_back(to_str(a.pos.y));
        vv.emplace_back(to_str(a.pos.z));
        vv.emplace_back(to_str(a.occ));
        vv.emplace_back(to_str(a.b_iso));
        vv.emplace_back(a.element.uname());
        vv.emplace_back(1, a.flag ? a.flag : '.'); // calc_flag
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
  }
  return crd;
}

void add_restraints(const Topo::Force force,
                    const Topo& topo, const Restraints& rt,
                    cif::Loop& restr_loop, int (&counters)[5]) {
  //using gemmi::to_str;
  const auto& to_str = gemmi::to_str_prec<3>; // to make comparisons easier
  const auto& to_str3 = gemmi::to_str_prec<3>;
  auto to_str_dot = [&](double x) { return std::isnan(x) ? "." : to_str(x); };
  if (force.rkind == Topo::RKind::Bond) {
    const Topo::Bond& t = topo.bonds[force.index];
    std::string obs = to_str3(t.calculate()) +
                      " # " + t.atoms[0]->name + " " + t.atoms[1]->name;
    restr_loop.add_row({"BOND", std::to_string(++counters[0]),
                        bond_type_to_string(t.restr->type), ".",
                        std::to_string(t.atoms[0]->serial),
                        std::to_string(t.atoms[1]->serial),
                        ".", ".",
                        to_str(t.restr->value), to_str(t.restr->esd),
                        to_str_dot(t.restr->value_nucleus),
                        to_str_dot(t.restr->esd_nucleus),
                        obs});
  } else if (force.rkind == Topo::RKind::Angle) {
    const Topo::Angle& t = topo.angles[force.index];
    std::string obs = to_str3(gemmi::deg(t.calculate()));
    obs += " # " + t.atoms[0]->name + " " +
                   t.atoms[1]->name + " " +
                   t.atoms[2]->name;
    restr_loop.add_row({"ANGL", std::to_string(++counters[1]),
                        ".", ".",
                        std::to_string(t.atoms[0]->serial),
                        std::to_string(t.atoms[1]->serial),
                        std::to_string(t.atoms[2]->serial),
                        ".",
                        to_str(t.restr->value), to_str(t.restr->esd),
                        ".", ".", obs});
  } else if (force.rkind == Topo::RKind::Torsion) {
    const Topo::Torsion& t = topo.torsions[force.index];
    std::string obs = to_str3(gemmi::deg(t.calculate()));
    obs += " # " + t.atoms[0]->name + " " + t.atoms[1]->name +
           " " + t.atoms[2]->name + " " + t.atoms[3]->name;
    restr_loop.add_row({"TORS", std::to_string(++counters[2]),
                        t.restr->label, std::to_string(t.restr->period),
                        std::to_string(t.atoms[0]->serial),
                        std::to_string(t.atoms[1]->serial),
                        std::to_string(t.atoms[2]->serial),
                        std::to_string(t.atoms[3]->serial),
                        to_str(t.restr->value), to_str(t.restr->esd),
                        ".", ".", obs});
  } else if (force.rkind == Topo::RKind::Chirality) {
    const Topo::Chirality& t = topo.chirs[force.index];
    double vol = rt.chiral_abs_volume(*t.restr);
    std::string obs = to_str3(t.calculate()) + " # " + t.atoms[0]->name +
                                                 " " + t.atoms[1]->name +
                                                 " " + t.atoms[2]->name +
                                                 " " + t.atoms[3]->name;
    restr_loop.add_row({"CHIR", std::to_string(++counters[3]),
                        gemmi::chirality_to_string(t.restr->sign), ".",
                        std::to_string(t.atoms[0]->serial),
                        std::to_string(t.atoms[1]->serial),
                        std::to_string(t.atoms[2]->serial),
                        std::to_string(t.atoms[3]->serial),
                        to_str3(vol), "0.020",
                        ".", ".", obs});
  } else if (force.rkind == Topo::RKind::Plane) {
    const Topo::Plane& t = topo.planes[force.index];
    ++counters[4];
    auto coeff = find_best_plane(t.atoms);
    for (const gemmi::Atom* atom : t.atoms) {
      double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
      std::string obs = to_str3(dist) + " # " + atom->name;
      restr_loop.add_row({"PLAN", std::to_string(counters[4]), t.restr->label,
                          ".", std::to_string(atom->serial), ".", ".", ".",
                          to_str(t.restr->esd), ".",
                          ".", ".", obs});
    }
  }
}

cif::Document make_rst(const Topo& topo, const gemmi::MonLib& monlib) {
  using Provenance = Topo::Provenance;
  cif::Document doc;
  doc.blocks.emplace_back("restraints");
  cif::Block& block = doc.blocks[0];
  cif::Loop& restr_loop = block.init_mmcif_loop("_restr.", {
              "record", "number", "label", "period",
              "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
              "value", "dev", "value_nucleus", "dev_nucleus", "val_obs"});
  int counters[5] = {0, 0, 0, 0, 0};
  for (const Topo::ChainInfo& chain_info : topo.chain_infos) {
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      // write link
      for (const Topo::ResInfo::Prev& prev : ri.prev) {
        const gemmi::Residue* prev_res = prev.get(&ri)->res;
        const gemmi::ChemLink* link = monlib.find_link(prev.link);
        if (link && count_provenance(ri.forces, Provenance::PrevLink) > 0) {
          std::string comment = " link " + prev.link + " " +
                                prev_res->seqid.str() + " " +
                                prev_res->name + " - " +
                                ri.res->seqid.str() + " " + ri.res->name;
          restr_loop.add_comment_and_row({comment, "LINK", ".",
                                          cif::quote(prev.link), ".",
                                          ".", ".", ".", ".", ".", ".", ".", ".", "."});
          for (const Topo::Force& force : ri.forces)
            if (force.provenance == Provenance::PrevLink)
              add_restraints(force, topo, link->rt, restr_loop, counters);
        }
      }
      // write monomer
      if (count_provenance(ri.forces, Provenance::Monomer) > 0) {
        std::string res_info = " monomer " + chain_info.name + " " +
                               ri.res->seqid.str() + " " + ri.res->name;
        if (!ri.mods.empty())
          res_info += " modified by " + gemmi::join_str(ri.mods, ", ");

        // need to revisit it later on
        std::string group = cif::quote(ri.chemcomp.group.substr(0, 8));
        if (group == "peptide" || group == "P-peptid" || group == "M-peptid")
          group = "L-peptid";

        restr_loop.add_comment_and_row({res_info, "MONO", ".", group, ".",
                                        ".", ".", ".", ".", ".", ".", ".", ".", "."});
        for (const Topo::Force& force : ri.forces)
          if (force.provenance == Provenance::Monomer)
            add_restraints(force, topo, ri.chemcomp.rt, restr_loop, counters);
      }
    }
  }
  // explicit links
  for (const Topo::ExtraLink& extra_link : topo.extras) {
    const gemmi::ChemLink* chem_link = monlib.find_link(extra_link.link_id);
    assert(chem_link);
    std::string comment = " link " + chem_link->id;
    restr_loop.add_comment_and_row({comment, "LINK", ".",
                                    cif::quote(chem_link->id), ".",
                                    ".", ".", ".", ".", ".", ".", ".", ".", "."});
    for (const Topo::Force& force : extra_link.forces)
      add_restraints(force, topo, chem_link->rt, restr_loop, counters);
  }
  return doc;
}

} // anonymous namespace

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
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  if (p.options[KeepHydrogens] && p.options[NoHydrogens])
    gemmi::fail("cannot use both --no-hydrogens and --keep-hydrogens");
  try {
    gemmi::Structure st = gemmi::read_structure_gz(input,
                                            gemmi::CoorFormat::UnknownAny);
    if (st.input_format == gemmi::CoorFormat::Pdb ||
        st.input_format == gemmi::CoorFormat::ChemComp)
      gemmi::setup_entities(st);
    if (st.models.empty()) {
      fprintf(stderr, "No models found in the input file.\n");
      return 1;
    }
    gemmi::Model& model0 = st.models[0];
    gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir,
                                                model0.get_all_residue_names(),
                                                gemmi::read_cif_gz);

    Topo topo;
    topo.initialize_refmac_topology(st, model0, monlib);

    // add H, sort atoms in residues and assign serial numbers
    int serial = 0;
    for (Topo::ChainInfo& chain_info : topo.chain_infos)
      for (Topo::ResInfo& ri : chain_info.res_infos) {
        const gemmi::ChemComp &cc = ri.chemcomp;
        gemmi::Residue &res = *ri.res;
        if (!p.options[KeepHydrogens]) {
          gemmi::remove_hydrogens(res);
          if (!p.options[NoHydrogens])
            if (cc.name != "HOH") // for compatibility with refmac/makecif
              add_hydrogens_without_positions(cc, res);
        }
        for (gemmi::Atom& atom : res.atoms) {
          auto it = cc.find_atom(atom.name);
          if (it == cc.atoms.end())
            gemmi::fail("No atom " + atom.name + " expected in " + res.name);
          atom.serial = it - cc.atoms.begin(); // temporary, for sorting only
        }
        std::sort(res.atoms.begin(), res.atoms.end(),
                  [](const gemmi::Atom& a, const gemmi::Atom& b) {
                    return a.serial != b.serial ? a.serial < b.serial
                                                : a.altloc < b.altloc;
                  });
        if (1) {  // temporary addition for makecif/refmac compatibility
          if (gemmi::in_vector(std::string("AA-STAND"),  ri.mods) &&
              !ri.res->find_atom("OXT", '*')) {
            gemmi::Atom atom;
            atom.name = "OXT";
            atom.element = gemmi::El::O;
            atom.pos = gemmi::Position(0, 0, 0);
            atom.flag = 'M';
            atom.occ = 0;
            atom.b_iso = 0;
            ri.res->atoms.push_back(atom);
          }
        }
        for (gemmi::Atom& atom : res.atoms)
          atom.serial = ++serial;
      }

    topo.finalize_refmac_topology(monlib);

    if (!p.options[KeepHydrogens] && !p.options[NoHydrogens])
      for (Topo::ChainInfo& chain_info : topo.chain_infos)
        for (Topo::ResInfo& ri : chain_info.res_infos)
          for (gemmi::Atom& atom : ri.res->atoms)
            if (!atom.is_hydrogen()) {
              try {
                place_hydrogens(atom, ri, topo);
              } catch (const std::runtime_error& e) {
                std::string loc = gemmi::atom_str(chain_info.name, *ri.res,
                                                  atom.name, atom.altloc);
                printf("Placing of hydrogen bonded to %s failed:\n  %s\n",
                       loc.c_str(), e.what());
              }
            }

    cif::Document crd = make_crd(st, monlib, topo);
    if (p.options[Verbose])
      printf("Writing coordinates to: %s.crd\n", output.c_str());
    {
      gemmi::Ofstream os(output + ".crd");
      write_cif_to_stream(os.ref(), crd, cif::Style::NoBlankLines);
    }

    if (p.options[NoZeroOccRestr])
      for (gemmi::Chain& chain : model0.chains)
        for (gemmi::Residue& res : chain.residues)
          for (gemmi::Atom& atom : res.atoms)
            if (atom.occ <= 0) {
              if (p.options[Verbose])
                printf("Atom with zero occupancy: %s\n",
                       gemmi::atom_str(chain, res, atom).c_str());
              atom.name += '?';  // hide the atom by mangling the name
            }

    cif::Document rst = make_rst(topo, monlib);
    if (p.options[Verbose])
      printf("Writing restraints to: %s.rst\n", output.c_str());
    {
      gemmi::Ofstream os(output + ".rst");
      write_cif_to_stream(os.ref(), rst, cif::Style::NoBlankLines);
    }
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
