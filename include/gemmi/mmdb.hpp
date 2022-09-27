// Copyright 2020-2022 Global Phasing Ltd.
//
// Converts between gemmi::Structure and mmdb::Manager.

#ifndef GEMMI_MMDB_HPP_
#define GEMMI_MMDB_HPP_

#include <cstdlib>           // for atoi
#include <cstring>           // for memcpy
#include <gemmi/model.hpp>
#include <gemmi/util.hpp>    // for rtrim_str
#include <mmdb2/mmdb_manager.h>

namespace gemmi {

inline void copy_transform_to_mmdb(const Transform& tr,
                                   mmdb::mat33& mat, mmdb::vect3& vec) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      mat[i][j] = tr.mat[i][j];
    vec[i] = tr.vec.at(i);
  }
}

template<int N>
void strcpy_to_mmdb(char (&dest)[N], const std::string& src) {
  if (src.size() >= N)
    fail("This string is too long: " + src);
  std::memcpy(dest, src.c_str(), src.size() + 1);
}

inline void set_seqid_in_mmdb(int* seqnum, mmdb::InsCode& icode, SeqId seqid) {
  *seqnum = *seqid.num;
  icode[0] = seqid.icode;
  icode[1] = '\0';
}

inline mmdb::Manager* copy_to_mmdb(const Structure& st, mmdb::Manager* manager) {
  for (const std::string& s : st.raw_remarks) {
    std::string line = rtrim_str(s);
    manager->PutPDBString(line.c_str());
  }

  for (int imodel = 0; imodel < (int) st.models.size(); ++imodel) {
    const Model& model = st.models[imodel];
    mmdb::PModel model2 = mmdb::newModel();
    model2->SetMMDBManager(manager, imodel);

    for (const Chain& chain : model.chains) {
      mmdb::PChain chain2 = model2->CreateChain(chain.name.c_str());
      for (const Residue& res : chain.residues) {
        const char icode[2] = {res.seqid.icode, '\0'};
        mmdb::PResidue res2 = chain2->GetResidueCreate(res.name.c_str(),
                                                       *res.seqid.num,
                                                       icode,
                                                       true);
        if (res.is_cis) {
          if (const Residue* next = chain.next_residue(res)) {
            mmdb::CisPep* cispep = new mmdb::CisPep();
            cispep->serNum = model2->GetCisPeps()->Length() + 1;
            strcpy_to_mmdb(cispep->chainID1, chain.name);
            strcpy_to_mmdb(cispep->pep1, res.name);
            set_seqid_in_mmdb(&cispep->seqNum1, cispep->icode1, res.seqid);
            strcpy_to_mmdb(cispep->chainID2, chain.name);
            strcpy_to_mmdb(cispep->pep2, next->name);
            set_seqid_in_mmdb(&cispep->seqNum2, cispep->icode2, next->seqid);
            cispep->modNum = imodel;
            model2->AddCisPep(cispep);
          }
        }
        for (const Atom& atom : res.atoms) {
          mmdb::PAtom atom2 = mmdb::newAtom();
          const char altloc[2] = {atom.altloc, '\0'};
          atom2->SetAtomName(0, atom.serial, atom.padded_name().c_str(),
                             altloc, res.segment.c_str(), atom.element.uname());
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
          mmdb::PAtom atom2 = mmdb::newAtom();
          atom2->MakeTer();
          atom2->serNum = res.atoms.back().serial + 1;
          res2->AddAtom(atom2);
        }
      }
    }
    int n_model = manager->AddModel(model2);
    manager->GetModel(n_model)->CopyCisPeps(model2);
    manager->PutCell(st.cell.a, st.cell.b, st.cell.c,
                     st.cell.alpha, st.cell.beta, st.cell.gamma, 1);
    manager->SetSpaceGroup(st.spacegroup_hm.c_str());
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
  return manager;
}


// don't use const for mmdb objects - MMDB doesn't have const method
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

inline Residue copy_residue_from_mmdb(mmdb::Residue& m_res) {
  Residue res;
  res.name = m_res.name;
  res.seqid.num = m_res.seqNum;
  if (m_res.insCode[0])
    res.seqid.icode = m_res.insCode[0];
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

inline Chain copy_chain_from_mmdb(mmdb::Chain& m_chain) {
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
  return chain;
}

inline Model copy_model_from_mmdb(mmdb::Model& m_model) {
  Model model(std::to_string(m_model.GetSerNum()));
  int n = m_model.GetNumberOfChains();
  model.chains.reserve(n);
  for (int i = 0; i < n; ++i)
    if (mmdb::Chain* m_chain = m_model.GetChain(i))
      model.chains.push_back(copy_chain_from_mmdb(*m_chain));
  return model;
}

inline Structure copy_from_mmdb(mmdb::Manager* manager) {
  Structure st;
  const mmdb::Cryst& cryst = *manager->GetCrystData();
  st.cell.set(cryst.a, cryst.b, cryst.c, cryst.alpha, cryst.beta, cryst.gamma);
  st.spacegroup_hm = cryst.spaceGroup;
  int n = manager->GetNumberOfModels();
  st.models.reserve(n);
  for (int i = 0; i < n; ++i)
    if (mmdb::Model* m_model = manager->GetModel(i+1))
      st.models.push_back(copy_model_from_mmdb(*m_model));
  return st;
}

} // namespace gemmi
#endif


/*
// Example 1.
// Read a coordinate file using gemmi and write it to pdb using mmdb.

#include <gemmi/mmread.hpp>
#include <gemmi/mmdb.hpp>

// two arguments expected: input and output paths.
int main (int argc, char** argv) {
  if (argc != 3)
    return 1;

  mmdb::InitMatType();
  mmdb::Manager* manager = new mmdb::Manager();
  try {
    gemmi::Structure st = gemmi::read_structure_file(argv[1]);
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
#define GEMMI_WRITE_IMPLEMENTATION
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
