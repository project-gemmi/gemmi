// Copyright 2026 Global Phasing Ltd.

#include "gemmi/ccp4ener.hpp"
#include "gemmi/acedrg_tables.hpp"
#include "gemmi/ace_graph.hpp"
#include "gemmi/util.hpp"
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace gemmi {
namespace {

struct Ccp4AtomInfo {
  Element el = El::X;
  std::string chem_type;
  std::string ccp4_type;
  int bonding_idx = 0;
  std::map<std::string, int> ring_rep;
  std::vector<int> conn_atoms;
  std::vector<int> conn_atoms_no_metal;
  std::vector<int> conn_h_atoms;
  float par_charge = 0.0f;
  int formal_charge = 0;
};

int ccp4_material_type(Element el) {
  switch (el.elem) {
    case El::H: case El::D:
      return 1;
    case El::C: case El::N: case El::O: case El::P: case El::S: case El::Se:
      return 2;
    case El::Li: case El::Na: case El::K: case El::Rb: case El::Cs: case El::Fr:
      return 3;
    case El::Be: case El::Mg: case El::Ca: case El::Sr: case El::Ba: case El::Ra:
      return 4;
    case El::Sc: case El::Y: case El::Ti: case El::Zr: case El::Hf: case El::Rf:
    case El::V: case El::Nb: case El::Ta: case El::Db:
    case El::Cr: case El::Mo: case El::W: case El::Sg:
    case El::Mn: case El::Tc: case El::Re: case El::Bh:
    case El::Fe: case El::Ru: case El::Os: case El::Hs:
    case El::Co: case El::Rh: case El::Ir: case El::Mt:
    case El::Ni: case El::Pd: case El::Pt: case El::Ds:
    case El::Cu: case El::Ag: case El::Au: case El::Rg:
    case El::Zn: case El::Cd: case El::Hg: case El::Cn:
      return 5;
    case El::Al: case El::Ga: case El::In: case El::Tl: case El::Sn:
    case El::Pb: case El::Bi:
      return 6;
    case El::B: case El::Si: case El::Ge: case El::As: case El::Sb:
    case El::Te: case El::Po:
      return 7;
    case El::F: case El::Cl: case El::Br: case El::I: case El::At:
      return 8;
    case El::La: case El::Ce: case El::Pr: case El::Nd: case El::Pm:
    case El::Sm: case El::Eu: case El::Gd: case El::Tb: case El::Dy:
    case El::Ho: case El::Er: case El::Tm: case El::Yb: case El::Lu:
    case El::Ac: case El::Th: case El::Pa: case El::U: case El::Np:
    case El::Pu: case El::Am: case El::Cm: case El::Bk: case El::Cf:
    case El::Es: case El::Fm: case El::Md: case El::No: case El::Lr:
      return 9;
    case El::He: case El::Ne: case El::Ar: case El::Kr: case El::Xe: case El::Rn:
      return 10;
    default:
      return 0;
  }
}

void set_hydro_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  atom.ccp4_type = "H";
  if (atom.conn_atoms.size() == 1) {
    int nb = atom.conn_atoms[0];
    if (atoms[nb].chem_type == "S")
      atom.ccp4_type = "HSH1";
  }
}

void set_org_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  int r5 = 0;
  int r6 = 0;
  for (const auto& item : atom.ring_rep) {
    if (item.second == 5)
      r5 += 1;
    if (item.second == 6)
      r6 += 1;
  }
  const size_t nconn = atom.conn_atoms_no_metal.size();
  const size_t nh = atom.conn_h_atoms.size();
  if (atom.chem_type == "C") {
    if (atom.bonding_idx == 2) {
      if (r5 && r6) atom.ccp4_type = "CR56";
      else if (r5 == 2) atom.ccp4_type = "CR55";
      else if (r6 == 2) atom.ccp4_type = "CR66";
      else if (r5 == 1) atom.ccp4_type = (nh == 1 ? "CR15" : nh == 0 ? "CR5" : "C");
      else if (r6 == 1) atom.ccp4_type = (nh == 1 ? "CR16" : nh == 0 ? "CR6" : "C");
      else if (nh == 1) atom.ccp4_type = "C1";
      else if (nh == 2) atom.ccp4_type = "C2";
      else if (nh == 0) atom.ccp4_type = "C";
    } else if (atom.bonding_idx == 3) {
      if (nh == 0) atom.ccp4_type = "CT";
      else if (nh == 1) atom.ccp4_type = "CH1";
      else if (nh == 2) {
        bool cage_linker = false;
        int heavy_conn = 0;
        for (int nb : atom.conn_atoms_no_metal) {
          if (atoms[nb].el == El::H)
            continue;
          ++heavy_conn;
          if (atoms[nb].el != El::C)
            continue;
          int heavy_deg = 0;
          bool has_b = false;
          for (int nb2 : atoms[nb].conn_atoms_no_metal) {
            if (atoms[nb2].el == El::H)
              continue;
            ++heavy_deg;
            if (atoms[nb2].el == El::B)
              has_b = true;
          }
          if (heavy_deg >= 5 && has_b)
            cage_linker = true;
        }
        atom.ccp4_type = (cage_linker && heavy_conn == 2) ? "CH3" : "CH2";
      } else if (nh == 3) {
        atom.ccp4_type = "CH3";
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "CSP";
    }
  } else if (atom.chem_type == "N") {
    if (atom.bonding_idx == 2) {
      if (nconn == 3) {
        if (nh == 1) atom.ccp4_type = "NH1";
        else if (nh == 2) atom.ccp4_type = "NH2";
        else if (nh == 0) atom.ccp4_type = "NH0";
        else atom.ccp4_type = "N";
      } else if (nconn == 2) {
        if (nh == 1) atom.ccp4_type = "N21";
        else if (nh == 0) atom.ccp4_type = "N20";
        else atom.ccp4_type = "N";
      }
    } else if (atom.bonding_idx == 3) {
      if (nconn == 4) {
        if (nh == 1) atom.ccp4_type = "NT1";
        else if (nh == 2) atom.ccp4_type = "NT2";
        else if (nh == 3) atom.ccp4_type = "NT3";
        else if (nh == 4) atom.ccp4_type = "NT4";
        else if (nh == 0) atom.ccp4_type = "NT";
        else atom.ccp4_type = "N";
      } else if (nconn == 3) {
        if (nh == 1) atom.ccp4_type = "N31";
        else if (nh == 2) atom.ccp4_type = "N32";
        else if (nh == 3) atom.ccp4_type = "N33";
        else if (nh == 0) atom.ccp4_type = "N30";
        else atom.ccp4_type = "N3";
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "NSP";
    }
  } else if (atom.chem_type == "P") {
    atom.ccp4_type = (nconn == 4 ? "P" : "P1");
  } else if (atom.chem_type == "O") {
    auto has_nb_type = [&](const char* t) {
      return std::any_of(atom.conn_atoms.begin(), atom.conn_atoms.end(),
                         [&](int nb) { return atoms[nb].chem_type == t; });
    };
    bool lP = has_nb_type("P"), lS = has_nb_type("S"), lB = has_nb_type("B");
    bool has_par_charge = std::fabs(atom.par_charge) > 1e-6f;
    bool has_negative_charge = atom.formal_charge < 0;
    auto oc_type = [&]() -> const char* {
      return lP ? "OP" : lS ? "OS" : lB ? "OB" : "OC";
    };
    if (atom.bonding_idx == 2) {
      if (has_par_charge && atom.par_charge < 0) {
        atom.ccp4_type = oc_type();
      } else if (nconn == 2) {
        if (nh == 1) atom.ccp4_type = "OH1";
        else if (nh == 2) atom.ccp4_type = "OH2";
        else atom.ccp4_type = "O";
      } else if (has_negative_charge) {
        atom.ccp4_type = oc_type();
      } else {
        atom.ccp4_type = "O";
      }
    } else if (atom.bonding_idx == 3) {
      bool lC = has_nb_type("C");
      if (lC && nh == 1 && nconn == 2) atom.ccp4_type = "OH1";
      else if (nh == 2) atom.ccp4_type = "OH2";
      else if (nconn == 2) {
        if (has_par_charge) atom.ccp4_type = "OC2";
        else if (nh == 1) atom.ccp4_type = "OH1";
        else atom.ccp4_type = "O2";
      } else if (nconn == 1 && has_negative_charge) {
        atom.ccp4_type = oc_type();
      }
    } else if (nconn == 1) {
      atom.ccp4_type = has_negative_charge ? oc_type() : "O";
    } else {
      atom.ccp4_type = "O";
    }
  } else if (atom.chem_type == "S") {
    if (nconn == 3 || nconn == 4) atom.ccp4_type = (nh == 0 ? "S3" : "SH1");
    else if (nconn == 2) atom.ccp4_type = (nh == 0 ? "S2" : "SH1");
    else if (nconn == 1) atom.ccp4_type = "S1";
    else atom.ccp4_type = (nh == 1 ? "SH1" : "S");
  } else if (atom.chem_type == "Se") {
    atom.ccp4_type = "SE";
  } else {
    atom.ccp4_type = atom.chem_type;
  }
}

void set_one_ccp4_type(std::vector<Ccp4AtomInfo>& atoms, size_t idx) {
  int ntype = ccp4_material_type(atoms[idx].el);
  switch (ntype) {
    case 1: set_hydro_ccp4_type(atoms, idx); break;
    case 2: set_org_ccp4_type(atoms, idx); break;
    default: atoms[idx].ccp4_type = atoms[idx].chem_type; break;
  }
  atoms[idx].ccp4_type = to_upper(atoms[idx].ccp4_type);
}

std::vector<Ccp4AtomInfo> build_ccp4_atoms(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int> >& neighbors) {
  std::vector<Ccp4AtomInfo> atoms;
  atoms.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    Ccp4AtomInfo info;
    info.el = cc.atoms[i].el;
    info.chem_type = cc.atoms[i].el.name();
    info.ccp4_type = info.chem_type;
    info.bonding_idx = atom_info[i].bonding_idx;
    info.ring_rep = atom_info[i].ring_rep;
    info.conn_atoms = neighbors[i];
    for (int nb : neighbors[i]) {
      if (cc.atoms[nb].is_hydrogen())
        info.conn_h_atoms.push_back(nb);
      if (!cc.atoms[nb].el.is_metal())
        info.conn_atoms_no_metal.push_back(nb);
    }
    if (info.el == El::C && info.bonding_idx == 2) {
      size_t total_conn = info.conn_atoms.size();
      if (total_conn >= 4 && total_conn > info.conn_atoms_no_metal.size())
        info.bonding_idx = 3;
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(atom_info[i].charge));
    atoms.emplace_back(std::move(info));
  }
  return atoms;
}

void assign_all_ccp4_types(std::vector<Ccp4AtomInfo>& atoms) {
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);
}

}  // namespace

std::vector<std::string> AcedrgTables::compute_ccp4_types(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int> >& neighbors) const {
  std::vector<Ccp4AtomInfo> atoms = build_ccp4_atoms(cc, atom_info, neighbors);
  assign_all_ccp4_types(atoms);
  std::vector<std::string> out;
  out.reserve(atoms.size());
  for (size_t i = 0; i < atoms.size(); ++i)
    out.push_back(atoms[i].ccp4_type);
  return out;
}

void assign_chemcomp_ccp4_types(ChemComp& cc) {
  if (cc.atoms.empty())
    return;
  AcedrgTables tables;
  std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);
  AceGraphView graph = make_ace_graph_view(cc);
  std::vector<Ccp4AtomInfo> atoms = build_ccp4_atoms(cc, atom_info, graph.neighbors);
  assign_all_ccp4_types(atoms);
  for (size_t i = 0; i < atoms.size(); ++i)
    cc.atoms[i].chem_type = atoms[i].ccp4_type;
}

void AcedrgTables::assign_ccp4_types(ChemComp& cc) const {
  assign_chemcomp_ccp4_types(cc);
}

}  // namespace gemmi
