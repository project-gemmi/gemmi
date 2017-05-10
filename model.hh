// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HH_
#define GEMMI_MODEL_HH_

#include <cstring>
#include <vector>
#include <string>
#include "elem.hh"
#include "unitcell.hh"

namespace gemmi {
namespace mol {

struct Structure;
struct Model;
struct Chain;
struct Residue;

enum class EntityType { Unknown, Polymer, NonPolymer, Water };

struct Atom {
  std::string name;
  char altloc;
  signed char charge;  // [-8, +8]
  Element element = El::X;
  Position pos;
  float occ;
  float b_iso;
  float u11=0, u22=0, u33=0, u12=0, u13=0, u23=0;
  Residue* parent = nullptr;
};

struct Residue {
  enum { UnknownId=-1000 };
  int seq_id = UnknownId;
  int auth_seq_id = UnknownId;
  char ins_code = '\0';
  std::string name;
  std::vector<Atom> atoms;
  Chain* parent = nullptr;

  Residue(int id, int auth_id, std::string rname) noexcept
    : seq_id(id), auth_seq_id(auth_id), name(rname) {}
  int seq_id_for_pdb() const {
    return auth_seq_id != UnknownId ? auth_seq_id : seq_id;
  }
  bool has_standard_pdb_name() const;
};

struct Chain {
  std::string name;
  std::string auth_name; // not guaranteed to be the same for the whole chain?
  EntityType entity_type = EntityType::Unknown;
  std::vector<Residue> residues;
  Model* parent = nullptr;
  explicit Chain(std::string cname) noexcept : name(cname) {}
};

struct Model {
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  Structure* parent = nullptr;
  explicit Model(std::string mname) noexcept : name(mname) {}
};

struct NcsOp {
  bool given;
  // may be changed to tranformation matrix, or quaternion+vector
  Matrix33 rot;
  Position tran;
};

struct Structure {
  UnitCell cell;
  std::string sg_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;

  // Minimal metadata with keys being mmcif tags: _entry.id, _exptl.method, ...
  std::map<std::string, std::string> info;

  const char* get_info(const std::string& tag, const char* def=nullptr) const {
    auto it = info.find(tag);
    return it != info.end() ? it->second.c_str() : def;
  }
};

inline bool Residue::has_standard_pdb_name() const {
#define SR(s) int(#s[0] << 16 | #s[1] << 8 | #s[2])
  const int standard_aa[26] = {
    SR(ALA), SR(ARG), SR(ASN), SR(ASP), SR(ASX), SR(CYS), SR(GLN), SR(GLU),
    SR(GLX), SR(GLY), SR(HIS), SR(ILE), SR(LEU), SR(LYS), SR(MET), SR(PHE),
    SR(PRO), SR(SER), SR(THR), SR(TRP), SR(TYR), SR(UNK), SR(VAL), SR(SEC),
    SR(PYL), 0
  };
  if (name.size() == 3) {
    int n = name[0] << 16 | name[1] << 8 | name[2];
    for (int i = 0; i < 26; i++)
      if (standard_aa[i] == n)
        return true;
    //return std::find(standard_aa, standard_aa + 26, name) != standard_aa + 26;
  } else if (name.size() == 1) {
    return std::strchr("ACGITU", name[0]) != nullptr;
  } else if (name.size() == 2) {
    return (name[0] == '+' || name[0] == 'D') && std::strchr("ACGITU", name[1]) != nullptr;
  }
  return false;
#undef SR
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
