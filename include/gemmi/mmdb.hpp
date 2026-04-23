/// @file
/// @brief Conversion between gemmi::Structure and MMDB (CCP4 Macromolecular Data Bank).
///
/// Copyright 2020-2022 Global Phasing Ltd.
///
/// Provides bidirectional conversion functions between Gemmi's Structure
/// representation and MMDB's Manager data structures. This enables
/// interoperability with CCP4's legacy molecular biology library.

#ifndef GEMMI_MMDB_HPP_
#define GEMMI_MMDB_HPP_

#include <cstdlib>           // for atoi
#include <cstring>           // for memcpy
#include "model.hpp"
#include "util.hpp"          // for rtrim_str
#include "polyheur.hpp"      // for assign_subchains
#include <mmdb2/mmdb_manager.h>

namespace gemmi {

/// @brief Copies a Gemmi transformation matrix and vector to MMDB format.
/// @param tr Gemmi transformation object (3x3 rotation + 3D translation).
/// @param mat Output MMDB 3x3 matrix (rotation component).
/// @param vec Output MMDB 3D vector (translation component).
inline void copy_transform_to_mmdb(const Transform& tr,
                                   mmdb::mat33& mat, mmdb::vect3& vec) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      mat[i][j] = tr.mat[i][j];
    vec[i] = tr.vec.at(i);
  }
}

/// @brief Safely copies a string into a fixed-size MMDB character array.
/// @tparam N Size of the destination character array.
/// @param dest Fixed-size MMDB character array (typically a typedef like mmdb::ChainID).
/// @param src String to copy.
/// @throws Throws if src is too long for the destination array.
template<int N>
void strcpy_to_mmdb(char (&dest)[N], const std::string& src) {
  if (src.size() + 1 >= N)
    fail("This string is too long: " + src);
  std::memcpy(dest, src.c_str(), src.size() + 1);
}

/// @brief Sets MMDB sequence number and insertion code from Gemmi SeqId.
/// @param seqnum Pointer to MMDB sequence number (output).
/// @param icode MMDB insertion code array (output).
/// @param seqid Gemmi sequence ID containing number and optional insertion code.
inline void set_seqid_in_mmdb(int* seqnum, mmdb::InsCode& icode, SeqId seqid) {
  *seqnum = *seqid.num;
  icode[0] = icode[1] = '\0';
  if (seqid.has_icode())
    icode[0] = seqid.icode;
}

/// @brief Converts MMDB sequence number and insertion code to Gemmi SeqId.
/// @param seqnum MMDB sequence number.
/// @param inscode MMDB insertion code array.
/// @return Gemmi SeqId with number and optional insertion code.
inline SeqId seqid_from_mmdb(int seqnum, const mmdb::InsCode& inscode) {
  return SeqId(seqnum, inscode[0] ? inscode[0] : ' ');
}

/// @brief Converts Gemmi CisPep to MMDB CisPep structure.
/// @param g Gemmi CisPep record to convert.
/// @param ser_num Serial number for the MMDB record.
/// @param model_num Model number for the MMDB record.
/// @return Newly allocated MMDB CisPep (caller owns memory).
inline mmdb::CisPep* cispep_to_mmdb(const CisPep& g, int ser_num, int model_num) {
  mmdb::CisPep* m = new mmdb::CisPep();
  m->serNum = ser_num;
  strcpy_to_mmdb(m->pep1, g.partner_c.res_id.name);
  strcpy_to_mmdb(m->chainID1, g.partner_c.chain_name);
  set_seqid_in_mmdb(&m->seqNum1, m->icode1, g.partner_c.res_id.seqid);
  strcpy_to_mmdb(m->pep2, g.partner_n.res_id.name);
  strcpy_to_mmdb(m->chainID2, g.partner_n.chain_name);
  set_seqid_in_mmdb(&m->seqNum2, m->icode2, g.partner_n.res_id.seqid);
  m->modNum = model_num;
  m->measure = g.reported_angle;
  return m;
}

/// @brief Converts MMDB CisPep structure to Gemmi CisPep.
/// @param m MMDB CisPep record to convert.
/// @param model_num Model number for the result.
/// @return Gemmi CisPep with data from MMDB record.
inline CisPep cispep_from_mmdb(const mmdb::CisPep& m, int model_num) {
  CisPep g;
  g.partner_c.res_id.name = m.pep1;
  g.partner_c.res_id.seqid = seqid_from_mmdb(m.seqNum1, m.icode1);
  g.partner_c.chain_name = m.chainID1;
  g.partner_n.res_id.name = m.pep2;
  g.partner_n.res_id.seqid = seqid_from_mmdb(m.seqNum2, m.icode2);
  g.partner_n.chain_name = m.chainID2;
  g.model_num = model_num;
  g.reported_angle = m.measure;
  return g;
}

/// @brief Transfers Connection records from Gemmi Structure to MMDB Manager.
/// @param st Gemmi Structure containing Connection records.
/// @param mol MMDB Manager to populate with Link records.
/// Converts Gemmi Connection objects to MMDB Link format and adds them
/// to all models in the MMDB Manager. Links with undefined sequence IDs
/// are skipped.
inline void transfer_links_to_mmdb(const Structure& st,  mmdb::Manager* mol) {
  // based on code provided by Paul Emsley
  for (const Connection& con : st.connections) {
    if (!con.partner1.res_id.seqid.num || !con.partner2.res_id.seqid.num)
      continue;
    mmdb::Link link{};
    // partner1
    strcpy_to_mmdb(link.atName1, con.partner1.atom_name);
    link.aloc1[0] = con.partner1.altloc;
    set_seqid_in_mmdb(&link.seqNum1, link.insCode1, con.partner1.res_id.seqid);
    strcpy_to_mmdb(link.resName1, con.partner1.res_id.name);
    strcpy_to_mmdb(link.chainID1, con.partner1.chain_name);
    // partner2
    strcpy_to_mmdb(link.atName2, con.partner2.atom_name);
    link.aloc2[0] = con.partner2.altloc;
    set_seqid_in_mmdb(&link.seqNum2, link.insCode2, con.partner2.res_id.seqid);
    strcpy_to_mmdb(link.resName2, con.partner2.res_id.name);
    strcpy_to_mmdb(link.chainID2, con.partner2.chain_name);
    if (con.reported_distance > 0)
      link.dist = con.reported_distance;
    if (con.asu == Asu::Different) {
      link.s2 = con.reported_sym[0];
      link.i2 = link.i1 + con.reported_sym[1];
      link.j2 = link.j1 + con.reported_sym[2];
      link.k2 = link.k1 + con.reported_sym[3];
    }
    // add links to models
    for (int imod = 1; imod <= mol->GetNumberOfModels(); imod++)
      if (mmdb::Model* model_p = mol->GetModel(imod))
        model_p->AddLink(new mmdb::Link(link));
  }
}

/// @brief Transfers Link records from MMDB to Gemmi Structure.
/// @param mmdb_links MMDB LinkContainer to read from.
/// @param st Gemmi Structure to populate with Connection records (output).
/// Converts MMDB Link objects to Gemmi Connection format.
/// Note: LinkR (restraint) entries are not transferred.
inline void transfer_links_from_mmdb(mmdb::LinkContainer& mmdb_links, Structure& st) {
  for (int i = 0; i < mmdb_links.Length(); ++i) {
    mmdb::Link& link = *static_cast<mmdb::Link*>(mmdb_links.GetContainerClass(i));
    Connection con;
    // partner1
    con.partner1.atom_name = link.atName1;
    con.partner1.altloc = link.aloc1[0];
    con.partner1.res_id.seqid = seqid_from_mmdb(link.seqNum1, link.insCode1);
    con.partner1.res_id.name = link.resName1;
    con.partner1.chain_name = link.chainID1;
    // partner2
    con.partner2.atom_name = link.atName2;
    con.partner2.altloc = link.aloc2[0];
    con.partner2.res_id.seqid = seqid_from_mmdb(link.seqNum2, link.insCode2);
    con.partner2.res_id.name = link.resName2;
    con.partner2.chain_name = link.chainID2;
    con.reported_distance = link.dist;
    if (link.s1 == link.s2 && link.i1 == link.i2 && link.j1 == link.j2 && link.k1 == link.k2) {
      con.asu = Asu::Same;
    } else {
      con.asu = Asu::Different;
      if (link.s1 == 1)
        con.reported_sym[0] = link.s2;
      con.reported_sym[1] = link.i2 - link.i1;
      con.reported_sym[2] = link.j2 - link.j1;
      con.reported_sym[3] = link.k2 - link.k1;
    }
    // for LinkR we'd also have:
    // con.link_id = link.linkRID;
    st.connections.push_back(con);
  }
}

/// @brief Converts Gemmi sequence vector to MMDB SeqRes format.
/// @param sequence Gemmi polymer sequence (residue names).
/// @param seqres MMDB SeqRes structure to populate (output).
/// Handles microheterogeneity by extracting the first monomer from each
/// sequence entry. Allocates and manages MMDB internal memory.
inline void set_mmdb_seqres(const std::vector<std::string>& sequence,
                            mmdb::SeqRes& seqres) {
  // Free existing sequence data
  delete[] seqres.resName;
  if (seqres.resName)
    seqres.resName = nullptr;

  if (sequence.empty()) {
    seqres.numRes = 0;
    return;
  }

  // Allocate new array
  seqres.numRes = static_cast<int>(sequence.size());
  seqres.resName = new mmdb::ResName[seqres.numRes];

  // Copy sequence data, taking first monomer from each item
  for (int i = 0; i < seqres.numRes; ++i) {
    // Handle microheterogeneity: use gemmi's utility function
    std::string mon = Entity::first_mon(sequence[i]);

    // Copy residue name (max 3 characters for most cases, but ResName allows up to 19)
    strcpy_to_mmdb(seqres.resName[i], mon);
  }
}

/// @brief Checks if Structure has polymer entities with defined sequences.
/// @param st Structure to check.
/// @return True if at least one polymer entity with sequence data exists.
inline bool has_sequences(const Structure& st) {
  // If there are no entities with sequences, try to infer sequences from chains
  for (const Entity& entity : st.entities) {
    if (entity.entity_type == EntityType::Polymer && !entity.full_sequence.empty())
      return true;
  }
  return false;
}

/// @brief Transfers polymer sequence information from Gemmi to MMDB Manager.
/// @param st Gemmi Structure containing entity sequence data.
/// @param manager MMDB Manager to populate (output).
/// Copies full_sequence data from Gemmi entities to corresponding MMDB
/// chain SEQRES records. Does nothing if Structure has no sequences.
inline void transfer_seqres_to_mmdb(const Structure& st, mmdb::Manager* manager) {
  if (!has_sequences(st))
    return;

  // Create a map of chain ID -> entity for faster lookup
  std::map<std::string, const Entity*> chain_to_entity;
  for (const Model& model : st.models)
    for (const Chain& ch : model.chains) {
      const Entity* entity = st.get_entity_of(ch.get_polymer());
      chain_to_entity[ch.name] = entity;
    }

  for (int imodel = 1; imodel <= manager->GetNumberOfModels(); ++imodel) {
    if (mmdb::Model* model = manager->GetModel(imodel)) {
      for (int ichain = 0; ichain < model->GetNumberOfChains(); ++ichain) {
        if (mmdb::Chain* chain = model->GetChain(ichain)) {
          std::string chain_id = chain->GetChainID();
          // MMDB uses empty string (not space) for blank chain IDs

          // Find matching entity
          auto it = chain_to_entity.find(chain_id);
          if (it != chain_to_entity.end()) {
            if (const Entity* entity = it->second)
              set_mmdb_seqres(entity->full_sequence, chain->seqRes);
            // Set chain association to ensure consistency
            chain->seqRes.SetChain(chain);
          }
        }
      }
    }
  }
}

/// @brief Copies Gemmi Structure data into MMDB Manager.
/// @param st Gemmi Structure to convert.
/// @param manager MMDB Manager to populate (output).
/// Transfers all structural information including atoms, residues, chains,
/// crystallographic metadata (cell, spacegroup, NCS operators), SEQRES records,
/// and Connection information. Handles TER records to mark polymer termination.
inline void copy_to_mmdb(const Structure& st, mmdb::Manager* manager) {
  for (const std::string& s : st.raw_remarks) {
    std::string line = rtrim_str(s);
    manager->PutPDBString(line.c_str());
  }

  for (int imodel = 0; imodel < (int) st.models.size(); ++imodel) {
    const Model& model = st.models[imodel];
    mmdb::Model* model2 = mmdb::newModel();
    manager->AddModel(model2);

    for (const Chain& chain : model.chains) {
      mmdb::ChainID chain_id = {};
      strcpy_to_mmdb(chain_id, chain.name);
      mmdb::Chain* chain2 = model2->CreateChain(chain_id);
      for (const Residue& res : chain.residues) {
        mmdb::ResName res_name = {};
        strcpy_to_mmdb(res_name, res.name);
        int seqnum = 0;
        mmdb::InsCode icode = {};
        set_seqid_in_mmdb(&seqnum, icode, res.seqid);
        mmdb::Residue* res2 = chain2->GetResidueCreate(res_name, seqnum,
                                                       icode, true);
        for (const Atom& atom : res.atoms) {
          mmdb::Atom* atom2 = mmdb::newAtom();
          mmdb::AltLoc altloc = {};
          if (atom.altloc != '\0')
            altloc[0] = atom.altloc;
          std::string padded_name = atom.padded_name();
          // padded_name() is padding from the left; MMDB from both sides
          if (padded_name.size() < 4)
            padded_name.resize(4, ' ');
          mmdb::AtomName atom_name = {};
          strcpy_to_mmdb(atom_name, padded_name);
          mmdb::SegID seg_id = {};
          strcpy_to_mmdb(seg_id, res.segment);
          mmdb::Element element = {};
          strcpy_to_mmdb(element, atom.element.uname());
          atom2->SetAtomName(0, atom.serial, atom_name,
                             altloc, seg_id, element);
          atom2->Het = res.het_flag == 'H';
          atom2->SetCharge(atom.charge);
          atom2->SetCoordinates(atom.pos.x, atom.pos.y, atom.pos.z,
                                atom.occ, atom.b_iso);
          if (atom.aniso.nonzero()) {
            atom2->u11 = atom.aniso.u11;
            atom2->u22 = atom.aniso.u22;
            atom2->u33 = atom.aniso.u33;
            atom2->u12 = atom.aniso.u12;
            atom2->u13 = atom.aniso.u13;
            atom2->u23 = atom.aniso.u23;
            atom2->WhatIsSet |= mmdb::ASET_Anis_tFSigma;
          }
          res2->AddAtom(atom2);
        }
        // TER
        if (res.entity_type == EntityType::Polymer &&
            (&res == &chain.residues.back() ||
             (&res + 1)->entity_type != EntityType::Polymer)) {
          mmdb::Atom* atom2 = mmdb::newAtom();
          atom2->MakeTer();
          atom2->serNum = res.atoms.back().serial + 1;
          res2->AddAtom(atom2);
        }
      }
    }
    manager->PutCell(st.cell.a, st.cell.b, st.cell.c,
                     st.cell.alpha, st.cell.beta, st.cell.gamma, 1);
    char spacegroup[64] = {};
    strcpy_to_mmdb(spacegroup, st.spacegroup_hm);
    manager->SetSpaceGroup(spacegroup);
    mmdb::Cryst* cryst = manager->GetCrystData();
    auto z = st.info.find("_cell.Z_PDB");
    if (z != st.info.end() && !z->second.empty()) {
      cryst->Z = std::atoi(z->second.c_str());
      cryst->WhatIsSet |= mmdb::CSET_ZValue;
    }
    if (st.has_origx && !st.origx.is_identity()) {
      copy_transform_to_mmdb(st.origx, cryst->o, cryst->t);
      cryst->WhatIsSet |= mmdb::CSET_OrigMatrix;
    }
    if (st.cell.explicit_matrices) {
      copy_transform_to_mmdb(st.cell.frac, cryst->s, cryst->u);
      cryst->WhatIsSet |= mmdb::CSET_ScaleMatrix;
    }
    if (!st.ncs.empty()) {
      mmdb::mat33 m;
      mmdb::vect3 v;
      if (st.info.find("_struct_ncs_oper.id") != st.info.end()) {
        copy_transform_to_mmdb(Transform{}, m, v);
        cryst->AddNCSMatrix(m, v, 1);
      }
      for (const NcsOp& op : st.ncs) {
        copy_transform_to_mmdb(op.tr, m, v);
        cryst->AddNCSMatrix(m, v, op.given);
      }
    }
  }
  if (!st.cispeps.empty() && !st.models.empty()) {
    int ser_num = 0;
    for (const CisPep& cispep : st.cispeps) {
      // In the PDB, CISPEP records have modNum=0 if there is only one model
      int modnum = st.models.size() > 1 ? cispep.model_num : 0;
      int model_no = std::max(1, modnum);  // GetModel() takes 1-based index
      if (mmdb::Model* m_model = manager->GetModel(model_no))
        m_model->AddCisPep(cispep_to_mmdb(cispep, ++ser_num, modnum));
    }
  }
  transfer_seqres_to_mmdb(st, manager);
  transfer_links_to_mmdb(st, manager);
}


/// @brief Converts MMDB Atom to Gemmi Atom.
/// @param m_atom MMDB Atom to convert (non-const; MMDB API limitation).
/// @return Gemmi Atom with data from MMDB record.
/// Extracts coordinates, occupancy, B-factors, anisotropic temperature factors
/// (if present), element, charge, and atom name/serial information.
inline Atom copy_atom_from_mmdb(mmdb::Atom& m_atom) {
  Atom atom;
  atom.name = m_atom.label_atom_id;
  atom.altloc = m_atom.altLoc[0];
  atom.charge = (signed char) m_atom.charge;
  atom.element = Element(m_atom.element);
  atom.serial = m_atom.serNum;
  atom.pos = Position(m_atom.x, m_atom.y, m_atom.z);
  atom.occ = (float) m_atom.occupancy;
  atom.b_iso = (float) m_atom.tempFactor;
  if (m_atom.WhatIsSet & mmdb::ASET_Anis_tFSigma) {
    atom.aniso.u11 = (float) m_atom.u11;
    atom.aniso.u22 = (float) m_atom.u22;
    atom.aniso.u33 = (float) m_atom.u33;
    atom.aniso.u12 = (float) m_atom.u12;
    atom.aniso.u13 = (float) m_atom.u13;
    atom.aniso.u23 = (float) m_atom.u23;
  }
  return atom;
}

/// @brief Converts MMDB Residue to Gemmi Residue.
/// @param m_res MMDB Residue to convert.
/// @return Gemmi Residue with atoms and metadata.
/// Extracts all atoms, identifies polymer termination via TER pseudo-atoms,
/// and sets het_flag and segment information from first atom.
inline Residue copy_residue_from_mmdb(mmdb::Residue& m_res) {
  Residue res;
  res.name = m_res.name;
  res.seqid = seqid_from_mmdb(m_res.seqNum, m_res.insCode);
  int n = m_res.GetNumberOfAtoms();
  res.atoms.reserve(n);
  bool first = true;
  for (int i = 0; i < n; ++i)
    if (mmdb::Atom* m_atom = m_res.GetAtom(i)) {
      if (m_atom->isTer()) {
        res.entity_type = EntityType::Polymer;
      } else {
        res.atoms.push_back(copy_atom_from_mmdb(*m_atom));
        if (first) {
          res.het_flag = m_atom->Het ? 'H' : 'A';
          res.segment = m_atom->segID;
          first = false;
        }
      }
    }
  return res;
}

/// @brief Converts MMDB SeqRes to Gemmi sequence vector.
/// @param seqres MMDB SeqRes structure to read.
/// @return Vector of residue names (sequence).
/// Extracts sequence data from MMDB SeqRes, handling empty entries and
/// trimming whitespace from each residue name.
inline std::vector<std::string> get_gemmi_sequence(const mmdb::SeqRes& seqres) {
  std::vector<std::string> sequence;
  if (seqres.numRes > 0 && seqres.resName) {
    sequence.reserve(seqres.numRes);
    for (int i = 0; i < seqres.numRes; ++i) {
      const char* s = seqres.resName[i];
      std::string resname(s, rtrim_cstr(s));
      if (!resname.empty())
        sequence.push_back(resname);
    }
  }
  return sequence;
}

/// @brief Converts MMDB Chain to Gemmi Chain and updates Structure entities.
/// @param st Gemmi Structure to update with entity information (output).
/// @param m_chain MMDB Chain to convert.
/// @return Gemmi Chain with residues and sequence data.
/// Extracts chain residues and SEQRES information, creating or updating
/// corresponding polymer entities in the Structure.
inline Chain copy_chain_from_mmdb(Structure& st, mmdb::Chain& m_chain) {
  Chain chain(m_chain.GetChainID());
  int n = m_chain.GetNumberOfResidues();
  chain.residues.reserve(n);
  for (int i = 0; i < n; ++i)
    if (mmdb::Residue* m_res = m_chain.GetResidue(i))
      chain.residues.push_back(copy_residue_from_mmdb(*m_res));
  // in MMDB we may have pseudo-atom TER that marks polymer end
  for (auto i = chain.residues.begin(); i != chain.residues.end(); ++i)
    if (i->entity_type == EntityType::Polymer) {  // residue before TER
      for (auto j = chain.residues.begin(); j != i; ++j)
        j->entity_type = EntityType::Polymer;
      for (auto j = i + 1; j != chain.residues.end(); ++j)
        j->entity_type = j->is_water() ? EntityType::Water
                                       : EntityType::NonPolymer;
      break;
    }
  // Get sequence from SEQRES if available
  std::vector<std::string> sequence = get_gemmi_sequence(m_chain.seqRes);
  if (!sequence.empty()) {
    Entity& ent = impl::find_or_add(st.entities, chain.name);
    ent.entity_type = EntityType::Polymer;
    ent.full_sequence = sequence;
  }
  return chain;
}

/// @brief Converts MMDB Model to Gemmi Model and updates Structure entities.
/// @param st Gemmi Structure to update with entity information (output).
/// @param m_model MMDB Model to convert.
/// @return Gemmi Model with chains.
/// Extracts all chains and their sequence information, updating entity
/// records and deduplicating entity data.
inline Model copy_model_from_mmdb(Structure& st, mmdb::Model& m_model) {
  Model model(m_model.GetSerNum());
  int n = m_model.GetNumberOfChains();
  model.chains.reserve(n);
  for (int i = 0; i < n; ++i)
    if (mmdb::Chain* m_chain = m_model.GetChain(i))
      model.chains.push_back(copy_chain_from_mmdb(st, *m_chain));
  ensure_entities(st);
  for (Model& m : st.models)
    for (Chain& chain : m.chains) {
      Entity& ent = impl::find_or_add(st.entities, chain.name);
      ent.subchains.push_back(chain.get_polymer().subchain_id());
      printf("Dodajemy Entity %s - %s\n",
          chain.name.c_str(), chain.get_polymer().subchain_id().c_str());
    }
  deduplicate_entities(st);
  return model;
}

/// @brief Converts MMDB Manager into Gemmi Structure.
/// @param manager MMDB Manager to convert.
/// @return Gemmi Structure with all models, atoms, and metadata.
/// Transfers crystallographic cell parameters, spacegroup information,
/// all models and chains, CISPEP records, and inter-chain connections.
inline Structure copy_from_mmdb(mmdb::Manager* manager) {
  Structure st;
  const mmdb::Cryst& cryst = *manager->GetCrystData();
  st.cell.set(cryst.a, cryst.b, cryst.c, cryst.alpha, cryst.beta, cryst.gamma);
  if (cryst.WhatIsSet & mmdb::CSET_ZValue)
    st.info.emplace("_cell.Z_PDB", std::to_string(cryst.Z));
  st.spacegroup_hm = cryst.spaceGroup;
  int n = manager->GetNumberOfModels();
  st.models.reserve(n);
  for (int i = 1; i <= n; ++i)
    if (mmdb::Model* m_model = manager->GetModel(i)) {
      st.models.push_back(copy_model_from_mmdb(st, *m_model));
      int model_num = st.models.back().num;
      for (int j = 1; j <= m_model->GetNumberOfCisPeps(); ++j)
        if (const mmdb::CisPep* m_cispep = m_model->GetCisPep(j))
          st.cispeps.push_back(cispep_from_mmdb(*m_cispep, model_num));
    }
  if (n > 0) {
    mmdb::Model* mmdb_model = manager->GetModel(1);
    transfer_links_from_mmdb(*mmdb_model->GetLinks(), st);
  }
  st.input_format = CoorFormat::Pdb;
  return st;
}

} // namespace gemmi
#endif


/*
// Example 1.
// Read a coordinate file using gemmi and write it to pdb using mmdb.

#include <gemmi/mmread_gz.hpp>
#include <gemmi/mmdb.hpp>

// two arguments expected: input and output paths.
int main (int argc, char** argv) {
  if (argc != 3)
    return 1;

  mmdb::InitMatType();
  mmdb::Manager* manager = new mmdb::Manager();
  try {
    gemmi::Structure st = gemmi::read_structure_gz(argv[1]);
    st.merge_chain_parts();
    gemmi::copy_to_mmdb(st, manager);
  } catch(std::runtime_error& e) {
    printf("File reading failed: %s\n", e.what());
    return 1;
  }

  mmdb::ERROR_CODE rc = manager->WritePDBASCII(argv[2]);
  if (rc)
    printf(" ***** ERROR #%i WRITE:\n\n %s\n\n",
           rc, mmdb::GetErrorDescription(rc));
  return 0;
}

// Example 2.
// Read a coordinate file using mmdb and write it to PDB using gemmi.
#include <fstream>
#include <gemmi/mmdb.hpp>
#include <gemmi/to_pdb.hpp>

// two arguments expected: input and output paths.
int main (int argc, char** argv) {
  if (argc != 3)
    return 1;
  mmdb::InitMatType();
  mmdb::Manager manager;
  manager.SetFlag(mmdb::MMDBF_PrintCIFWarnings      |
                  mmdb::MMDBF_FixSpaceGroup         |
                  mmdb::MMDBF_IgnoreHash            |
                  mmdb::MMDBF_IgnoreNonCoorPDBErrors|
                  mmdb::MMDBF_DoNotProcessSpaceGroup);
  mmdb::ERROR_CODE rc = manager.ReadCoorFile(argv[1]);
  if (rc != mmdb::Error_NoError) {
    printf(" ***** ERROR reading the file\n");
    return 1;
  }
  gemmi::Structure st = gemmi::copy_from_mmdb(&manager);
  std::ofstream os(argv[2]);
  gemmi::write_pdb(st, os);
  return 0;
}
*/
