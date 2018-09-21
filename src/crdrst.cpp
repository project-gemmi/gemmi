// Copyright 2017 Global Phasing Ltd.

#include <stdio.h>
#include <cstdlib> // for getenv
#include <cctype>  // for tolower
#include <numeric> // for accumulate
#include <set>
#include <stdexcept>
#include "gemmi/gzread.hpp"
#include "gemmi/chemcomp.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/entstr.hpp"    // for entity_type_to_string
#include "gemmi/to_mmcif.hpp"  // for write_struct_conn
#include "gemmi/sprintf.hpp"   // for to_str, to_str_prec
#include "gemmi/calculate.hpp" // for calculate_angle, find_best_plane
#include "gemmi/polyheur.hpp"  // for are_connected, remove_hydrogens

#define GEMMI_PROG crdrst
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Restraints;

enum OptionIndex { Verbose=3, Monomers, NoHydrogens, NoZeroOccRestr };

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
  { NoHydrogens, 0, "H", "no-hydrogens", Arg::None,
    "  -H, --no-hydrogens  \tRemove or do not add hydrogens." },
  { NoZeroOccRestr, 0, "", "no-zero-occ", Arg::None,
    "  --no-zero-occ  \tNo restraints for zero-occupancy atoms." },
  { 0, 0, 0, 0, 0, 0 }
};


struct MonLib {
  cif::Document mon_lib_list;
  std::map<std::string, gemmi::ChemComp> monomers;
  std::map<std::string, gemmi::ChemLink> links;
  std::map<std::string, gemmi::ChemMod> modifications;

  const gemmi::ChemLink* find_link(const std::string& link_id) const {
    auto link = links.find(link_id);
    return link != links.end() ? &link->second : nullptr;
  }
  const gemmi::ChemMod* find_mod(const std::string& name) const {
    auto modif = modifications.find(name);
    return modif != modifications.end() ? &modif->second : nullptr;
  }
  const gemmi::ChemLink* match_link(const gemmi::ChemLink& link) const {
    for (auto& ml : links)
      if (link.matches(ml.second))
        return &ml.second;
    return nullptr;
  }

};

inline MonLib read_monomers(std::string monomer_dir,
                            const std::set<std::string>& resnames) {
  MonLib monlib;
  assert(!monomer_dir.empty());
  if (monomer_dir.back() != '/' && monomer_dir.back() != '\\')
    monomer_dir += '/';
  monlib.mon_lib_list = gemmi::read_cif_gz(monomer_dir +
                                           "list/mon_lib_list.cif");
  std::string error;
  for (const std::string& name : resnames) {
    std::string path = monomer_dir;
    path += std::tolower(name[0]);
    path += '/';
    path += name + ".cif";
    try {
      cif::Document doc = gemmi::read_cif_gz(path);
      auto cc = gemmi::make_chemcomp_from_cif(name, doc);
      monlib.monomers.emplace(name, cc);
    } catch(std::runtime_error& err) {
      error += "The monomer " + name + " could not be read.\n";
    }
  }
  if (!error.empty())
    gemmi::fail(error + "Please create definitions for missing monomers.");
  monlib.links = gemmi::read_chemlinks(monlib.mon_lib_list);
  monlib.modifications = gemmi::read_chemmods(monlib.mon_lib_list);
  return monlib;
}

struct Linkage {
  struct ResInfo {
    const gemmi::Residue* res;
    std::string prev_link;
    int prev_idx;
    std::vector<std::string> mods;
    gemmi::ChemComp chemcomp;

    ResInfo(const gemmi::Residue* r) : res(r), prev_idx(0) {}
    const gemmi::Residue* prev_res() const {
      return prev_idx != 0 ? (this + prev_idx)->res : nullptr;
    }
  };

  struct ChainInfo {
    std::string name;
    std::string entity_id;
    bool polymer;
    gemmi::PolymerType polymer_type;
    std::vector<ResInfo> residues;
  };

  struct ExtraLink {
    const gemmi::Residue* res1;
    const gemmi::Residue* res2;
    char alt1 = '\0';
    char alt2 = '\0';
    gemmi::ChemLink link;
  };

  std::vector<ChainInfo> chains;
  std::vector<ExtraLink> extra;
};

static const char* get_front_modificat(const gemmi::Residue& res,
                                       gemmi::PolymerType ptype) {
  (void) res;
  if (is_polypeptide(ptype))
    return "NH3";
  if (is_polynucleotide(ptype))
    return "5*END";
  return nullptr;
}

static const char* get_back_modificat(const gemmi::Residue& res,
                                      gemmi::PolymerType ptype) {
  if (is_polypeptide(ptype))
    return res.find_atom("OXT") ? "COO" : "TERMINUS";
  if (is_polynucleotide(ptype))
    return "TERMINUS";
  return nullptr;
}

static Linkage::ChainInfo initialize_chain_info(const gemmi::SubChain& subchain,
                                                const gemmi::Entity* ent) {
  Linkage::ChainInfo lc;
  lc.residues.reserve(subchain.size());
  lc.name = subchain.name();
  lc.entity_id = ent->name;
  lc.polymer = ent && ent->entity_type == gemmi::EntityType::Polymer;
  lc.polymer_type = ent->polymer_type;
  for (const gemmi::Residue& res : subchain) {
    lc.residues.emplace_back(&res);
  }
  return lc;
}

static void setup_polymer_links(Linkage::ChainInfo& ci) {
  if (!ci.polymer || ci.residues.empty())
    return;
  for (auto ri = ci.residues.begin() + 1; ri != ci.residues.end(); ++ri) {
    // For now we ignore microheterogeneity.
    ri->prev_idx = -1;
    const gemmi::Residue* prev_res = (ri + ri->prev_idx)->res;
    if (!prev_res) {
      ri->prev_link = ".";
    } else if (!are_connected(*prev_res, *ri->res, ci.polymer_type)) {
      ri->prev_link = "gap";
    } else if (is_polypeptide(ci.polymer_type)) {
      if (ri->chemcomp.group == "P-peptide")
        ri->prev_link = "P";  // PCIS, PTRANS
      else if (ri->chemcomp.group == "M-peptide")
        ri->prev_link = "NM"; // NMCIS, NMTRANS
      ri->prev_link += prev_res->is_cis ? "CIS" : "TRANS";
    } else if (is_polynucleotide(ci.polymer_type)) {
      ri->prev_link = "p";
    } else {
      ri->prev_link = "?";
    }
  }
}

static bool has_anisou(const gemmi::Model& model) {
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& res : chain.residues)
      for (const gemmi::Atom& a : res.atoms)
        if (a.has_anisou())
          return true;
  return false;
}

// for compatibility with makecif, not sure what ccp4_mod_id is really used for
static std::string get_ccp4_mod_id(const std::vector<std::string>& mods) {
  for (const std::string& m : mods)
    if (m != "AA-STAND" && !gemmi::starts_with(m, "DEL-OXT") &&
        !gemmi::starts_with(m, "DEL-HN") && m != "DEL-NMH")
      return m;
  return ".";
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
  for (const gemmi::Entity& ent : st.entities)
    entity_loop.add_row({ent.name, entity_type_to_string(ent.entity_type)});
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
      std::string prev = res_info.prev_idx ? res_info.prev_res()->seqid.str()
                                           : "n/a";
      std::string mod = get_ccp4_mod_id(res_info.mods);
      poly_loop.add_row({res_info.res->name, res_info.res->seqid.str(),
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
  const gemmi::Model& model0 = st.models.at(0);
  items.emplace_back(cif::CommentArg{"#################\n"
                                     "## STRUCT_ASYM ##\n"
                                     "#################"});
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const gemmi::Chain& chain : model0.chains)
    for (gemmi::SubChain sub : const_cast<gemmi::Chain&>(chain).subchains())
      if (sub.labelled()) {
        const gemmi::Entity* ent = st.get_entity_of(sub);
        asym_loop.add_row({sub.name(), (ent ? ent->name : "?")});
      }
  const auto& connections = model0.connections;
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
      gemmi::ChemComp& cc = monlib.monomers.at(res.name);
      for (const gemmi::Atom& a : res.atoms) {
        vv.emplace_back("ATOM");
        vv.emplace_back(std::to_string(a.serial));
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
        vv.emplace_back(1, a.flag ? a.flag : '.'); // calc_flag
        vv.emplace_back(1, '.'); // label_seg_id
        vv.emplace_back(a.name); // again
        vv.emplace_back(cc.get_atom(a.name).chem_type); // label_chem_id
        if (write_anisou) {
          if (a.has_anisou()) {
            for (float u : {a.u11, a.u22, a.u33, a.u12, a.u13, a.u23})
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

static std::string chirality_to_string(Restraints::Chirality::Type ctype) {
  switch (ctype) {
    case Restraints::Chirality::Positive: return "positive";
    case Restraints::Chirality::Negative: return "negative";
    case Restraints::Chirality::Both: return "both";
    default: return "???";
  }
}

static int add_restraints(const Restraints& rt,
                          const gemmi::Residue& res, const gemmi::Residue* res2,
                          cif::Loop& restr_loop, int (&counters)[5],
                          char altloc='*') {
  //using gemmi::to_str;
  const auto& to_str = gemmi::to_str_prec<3>; // to make comparisons easier
  const auto& to_str3 = gemmi::to_str_prec<3>;

  int init_count = std::accumulate(counters, counters + 5, 0);

  std::string altlocs;
  if (altloc == '*') {
    // find all distinct altlocs
    for (const gemmi::Atom& atom : res.atoms)
      if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
        altlocs += atom.altloc;
    if (res2)
      for (const gemmi::Atom& atom : res2->atoms)
        if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
          altlocs += atom.altloc;
  }
  if (altlocs.empty())
    altlocs += altloc;

  for (const Restraints::Bond& bond : rt.bonds)
    for (char alt : altlocs)
      if (const gemmi::Atom* at1 = bond.id1.get_from(res, res2, alt))
        if (const gemmi::Atom* at2 = bond.id2.get_from(res, res2, alt)) {
          std::string obs = to_str3(at1->pos.dist(at2->pos));
          obs += " # " + at1->name + " " + at2->name;
          restr_loop.add_row({"BOND", std::to_string(++counters[0]),
                              bond_type_to_string(bond.type), ".",
                              std::to_string(at1->serial),
                              std::to_string(at2->serial),
                              ".", ".",
                              to_str(bond.value), to_str(bond.esd), obs});
          if (!at1->altloc && !at2->altloc)
            break;
        }
  for (const Restraints::Angle& angle : rt.angles)
    for (char alt : altlocs)
      if (const gemmi::Atom* at1 = angle.id1.get_from(res, res2, alt))
        if (const gemmi::Atom* at2 = angle.id2.get_from(res, res2, alt))
          if (const gemmi::Atom* at3 = angle.id3.get_from(res, res2, alt)) {
            double a = gemmi::calculate_angle(at1->pos, at2->pos, at3->pos);
            std::string obs = to_str3(gemmi::deg(a));
            obs += " # " + at1->name + " " + at2->name + " " + at3->name;
            restr_loop.add_row({"ANGL", std::to_string(++counters[1]),
                                ".", ".",
                                std::to_string(at1->serial),
                                std::to_string(at2->serial),
                                std::to_string(at3->serial),
                                ".",
                                to_str(angle.value), to_str(angle.esd), obs});
            if (!at1->altloc && !at2->altloc && !at3->altloc)
              break;
          }
  for (const Restraints::Torsion& tor : rt.torsions)
    for (char alt : altlocs)
      if (const gemmi::Atom* at1 = tor.id1.get_from(res, res2, alt))
        if (const gemmi::Atom* at2 = tor.id2.get_from(res, res2, alt))
          if (const gemmi::Atom* at3 = tor.id3.get_from(res, res2, alt))
            if (const gemmi::Atom* at4 = tor.id4.get_from(res, res2, alt)) {
              double d = gemmi::calculate_dihedral(at1->pos, at2->pos,
                                                   at3->pos, at4->pos);
              std::string obs = to_str3(gemmi::deg(d));
              obs += " # " + at1->name + " " + at2->name +
                     " " + at3->name + " " + at4->name;
              restr_loop.add_row({"TORS", std::to_string(++counters[2]),
                                  tor.label, std::to_string(tor.period),
                                  std::to_string(at1->serial),
                                  std::to_string(at2->serial),
                                  std::to_string(at3->serial),
                                  std::to_string(at4->serial),
                                  to_str(tor.value), to_str(tor.esd), obs});
              if (!at1->altloc && !at2->altloc && !at3->altloc && !at4->altloc)
                break;
        }
  for (const Restraints::Chirality& chir : rt.chirs)
    for (char alt : altlocs)
      if (const gemmi::Atom* at1 = chir.id_ctr.get_from(res, res2, alt))
        if (const gemmi::Atom* at2 = chir.id1.get_from(res, res2, alt))
          if (const gemmi::Atom* at3 = chir.id2.get_from(res, res2, alt))
            if (const gemmi::Atom* at4 = chir.id3.get_from(res, res2, alt)) {
              double vol = rt.chiral_abs_volume(chir);
              double obs_vol = gemmi::calculate_chiral_volume(
                                    at1->pos, at2->pos, at3->pos, at4->pos);
              std::string obs = to_str3(obs_vol)
                                + " # " + at1->name + " " + at2->name
                                + " " + at3->name + " " + at4->name;
              restr_loop.add_row({"CHIR", std::to_string(++counters[3]),
                                  chirality_to_string(chir.chir), ".",
                                  std::to_string(at1->serial),
                                  std::to_string(at2->serial),
                                  std::to_string(at3->serial),
                                  std::to_string(at4->serial),
                                  to_str3(vol), "0.020", obs});
              if (!at1->altloc && !at2->altloc && !at3->altloc && !at4->altloc)
                break;
            }
  for (const Restraints::Plane& plane : rt.planes)
    for (char alt : altlocs) {
      std::vector<const gemmi::Atom*> atoms;
      for (const Restraints::AtomId& id : plane.ids)
        if (const gemmi::Atom* atom = id.get_from(res, res2, alt))
          atoms.push_back(atom);
      if (atoms.size() < 4)
        continue;
      ++counters[4];
      auto coeff = find_best_plane(atoms);
      for (const gemmi::Atom* atom : atoms) {
        double dist = coeff[0] * atom->pos.x + coeff[1] * atom->pos.y +
                      coeff[2] * atom->pos.z + coeff[3];
        std::string obs = to_str3(dist) + " # " + atom->name;
        restr_loop.add_row({"PLAN", std::to_string(counters[4]), plane.label,
                            ".", std::to_string(atom->serial), ".", ".", ".",
                            to_str(plane.esd), ".", obs});
      }
      if (std::all_of(atoms.begin(), atoms.end(),
                      [](const gemmi::Atom* a) { return !a->altloc; }))
        break;
    }
  return std::accumulate(counters, counters + 5, 0) - init_count;
}

static void make_unique_link_name(std::string& name, const MonLib& monlib) {
  size_t orig_len = name.size();
  for (int n = 1; monlib.find_link(name) != nullptr; ++n) {
    name.resize(orig_len);
    name += std::to_string(n);
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
      if (const gemmi::Residue* prev = ri.prev_res()) {
        const gemmi::ChemLink* link = monlib.find_link(ri.prev_link);
        if (link && !link->rt.empty()) {
          std::string comment = "# link " + ri.prev_link + " " +
                                 prev->seqid.str() + " " + prev->name + " - " +
                                 ri.res->seqid.str() + " " + ri.res->name;
          restr_loop.add_row({comment + "\nLINK", ".", cif::quote(ri.prev_link),
                              ".", ".", ".", ".", ".", ".", ".", "."});
          int n = add_restraints(link->rt, *prev, ri.res, restr_loop, counters);
          if (n == 0)
            restr_loop.pop_row();  // remove the LINK line
        }
      }
      // write monomer
      if (!ri.chemcomp.rt.empty()) {
        // comments are added relying on how cif writing works
        std::string res_info = "# monomer " + chain_info.name + " " +
                               ri.res->seqid.str() + " " + ri.res->name;

        // need to revisit it later on
        std::string group = cif::quote(ri.chemcomp.group.substr(0, 8));
        if (group == "peptide" || group == "P-peptid" || group == "M-peptid")
          group = "L-peptid";

        restr_loop.add_row({res_info + "\nMONO", ".", group,
                            ".", ".", ".", ".", ".", ".", ".", "."});
        int n = add_restraints(ri.chemcomp.rt, *ri.res, nullptr,
                               restr_loop, counters);
        if (n == 0)
          restr_loop.pop_row();  // remove the MONO line
      }
    }
  }
  // explicit links
  for (const Linkage::ExtraLink& link : linkage.extra) {
    const gemmi::ChemLink* chem_link = monlib.match_link(link.link);
    if (!chem_link) {
      std::string link_name = link.res1->name + "-" + link.res2->name;
      make_unique_link_name(link_name, monlib);
      gemmi::ChemLink& v = monlib.links[link_name];
      v = link.link;
      v.id = link_name;
      chem_link = &v;
    }
    std::string comment = "# link " + chem_link->id;
    restr_loop.add_row({comment + "\nLINK", ".", cif::quote(chem_link->id),
                        ".", ".", ".", ".", ".", ".", ".", "."});
    char altloc = link.alt1 ? link.alt1 : (link.alt2 ? link.alt2 : '*');
    if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
      printf("Warning: LINK between different conformers %c and %c.",
             link.alt1, link.alt2);
    add_restraints(chem_link->rt, *link.res1, link.res2, restr_loop, counters,
                   altloc);
  }
  return doc;
}

static
gemmi::ChemLink connection_to_chemlink(const gemmi::Connection& conn,
                                       const gemmi::Residue& res1,
                                       const gemmi::Residue& res2) {
  gemmi::ChemLink link;
  link.id = res1.name + "-" + res2.name;
  link.comp[0] = res1.name;
  link.comp[1] = res2.name;
  gemmi::Restraints::Bond bond;
  bond.id1 = Restraints::AtomId{1, conn.atom[0].atom_name};
  bond.id2 = Restraints::AtomId{2, conn.atom[1].atom_name};
  bond.type = gemmi::Restraints::Bond::Unspec;
  bond.aromatic = false;
  bond.value = conn.reported_distance;
  bond.esd = 0.02;
  link.rt.bonds.push_back(bond);
  return link;
}

static void add_hydrogens(Linkage::ResInfo& ri) {
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
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  try {
    gemmi::Structure st = gemmi::read_structure_gz(input);
    if (st.input_format == gemmi::CoorFormat::Pdb)
      gemmi::setup_entities(st);
    if (st.models.empty())
      return 1;
    gemmi::Model& model0 = st.models[0];
    if (p.options[NoHydrogens])
      gemmi::remove_hydrogens(model0);
    std::set<std::string> resnames;
    for (const gemmi::Chain& chain : model0.chains)
      for (const gemmi::Residue& res : chain.residues)
        resnames.insert(res.name);

    MonLib monlib = read_monomers(monomer_dir, resnames);

    // sort atoms in residues and assign serial numbers
    int serial = 0;
    for (gemmi::Chain& chain : model0.chains)
      for (gemmi::Residue& res : chain.residues) {
        const gemmi::ChemComp &cc = monlib.monomers.at(res.name);
        for (gemmi::Atom& atom : res.atoms) {
          auto it = cc.find_atom(atom.name);
          if (it == cc.atoms.end())
            gemmi::fail("No atom " + atom.name + " expected in " + res.name);
          atom.serial = it - cc.atoms.begin();
        }
        std::sort(res.atoms.begin(), res.atoms.end(),
                  [](const gemmi::Atom& a, const gemmi::Atom& b) {
                    return a.serial != b.serial ? a.serial < b.serial
                                                : a.altloc < b.altloc;
                  });
        for (gemmi::Atom& atom : res.atoms)
          atom.serial = ++serial;
      }

    Linkage linkage;
    // initialize chains and residues
    for (const gemmi::Chain& chain : model0.chains)
      for (gemmi::SubChain sub : const_cast<gemmi::Chain&>(chain).subchains()) {
        assert(sub.labelled());
        const gemmi::Entity* ent = st.get_entity_of(sub);
        linkage.chains.push_back(initialize_chain_info(sub, ent));
      }
    for (Linkage::ChainInfo& ci : linkage.chains) {
      // copy monomer description
      for (Linkage::ResInfo& ri : ci.residues)
        ri.chemcomp = monlib.monomers.at(ri.res->name);
      // setup polymer links
      setup_polymer_links(ci);
      // add built-in modifications
      if (ci.polymer && !ci.residues.empty()) {
        // we try to get exactly the same numbers that makecif produces
        for (Linkage::ResInfo& ri : ci.residues)
          if (ci.polymer_type == gemmi::PolymerType::PeptideL)
            ri.mods.emplace_back("AA-STAND");
        Linkage::ResInfo& front = ci.residues.front();
        if (const char* mod = get_front_modificat(*front.res, ci.polymer_type))
          front.mods.emplace_back(mod);
        Linkage::ResInfo& back = ci.residues.back();
        if (const char* mod = get_back_modificat(*back.res, ci.polymer_type))
          back.mods.emplace_back(mod);
      }
      // add modifications from standard links
      for (Linkage::ResInfo& ri : ci.residues)
        if (const gemmi::ChemLink* link = monlib.find_link(ri.prev_link)) {
          if (!link->mod[0].empty())
            (&ri + ri.prev_idx)->mods.push_back(link->mod[0]);
          if (!link->mod[1].empty())
            ri.mods.push_back(link->mod[1]);
        }
    }
    // add extra links
    for (const gemmi::Connection& conn : model0.connections) {
      Linkage::ExtraLink extra;
      extra.res1 = model0.find_cra(conn.atom[0]).residue;
      extra.res2 = model0.find_cra(conn.atom[1]).residue;
      extra.alt1 = conn.atom[0].altloc;
      extra.alt2 = conn.atom[1].altloc;
      if (extra.res1 && extra.res2) {
        extra.link = connection_to_chemlink(conn, *extra.res1, *extra.res2);
        linkage.extra.push_back(extra);
      }
    }
    // TODO: automatically determine other links

    for (Linkage::ChainInfo& chain_info : linkage.chains)
      for (Linkage::ResInfo& ri : chain_info.residues) {
        // apply modifications
        for (const std::string& modif : ri.mods) {
          if (const gemmi::ChemMod* chem_mod = monlib.find_mod(modif))
            try {
              chem_mod->apply_to(ri.chemcomp);
            } catch(std::runtime_error& e) {
              printf("Failed to apply modification %s to %s: %s\n",
                     chem_mod->id.c_str(), ri.res->name.c_str(), e.what());
            }
          else
            printf("Modification not found: %s\n", modif.c_str());
        }

        // add hydrogens
        add_hydrogens(ri);
      }

    cif::Document crd = make_crd(st, monlib, linkage);
    if (p.options[Verbose])
      printf("Writing coordinates to: %s.crd\n", output.c_str());
    write_cif_to_file(crd, output + ".crd", cif::Style::NoBlankLines);

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
