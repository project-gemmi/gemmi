// Copyright 2017 Global Phasing Ltd.
//
// Data structures to keep macromolecular structure model.

#ifndef GEMMI_MODEL_HH_
#define GEMMI_MODEL_HH_

#include <cstring>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include "elem.hh"
#include "unitcell.hh"

namespace gemmi {
namespace mol {

namespace internal {

template<typename T>
T* find_or_null(std::vector<T>& vec, const std::string& name) {
  auto it = std::find_if(vec.begin(), vec.end(), [&name](const T& m) {
      return m.name == name;
  });
  return it != vec.end() ? &*it : nullptr;
}

template<typename T>
T* find_or_add(std::vector<T>& vec, const std::string& name) {
  T* ret = find_or_null(vec, name);
  if (ret)
    return ret;
  vec.emplace_back(name);
  return &vec.back();
}

} // namespace internal


struct Structure;
struct Model;
struct Chain;
struct Residue;

enum class EntityType { Unknown, Polymer, NonPolymer, Water };

struct SequenceItem {
  int num;
  std::string mon;
  SequenceItem(std::string m) noexcept : num(-1), mon(m) {}
  SequenceItem(int n, std::string m) noexcept : num(n), mon(m) {}
};

using Sequence = std::vector<SequenceItem>;

struct Entity {
  std::string id;  // it does not need to be number according to mmCIF spec
  EntityType type;
  Sequence sequence;

  explicit Entity(const std::string& id_, EntityType t) : id(id_), type(t) {}

  std::string type_as_string() {
    switch (type) {
      case EntityType::Polymer: return "polymer";
      case EntityType::NonPolymer: return "non-polymer";
      case EntityType::Water: return "water";
      default /*EntityType::Unknown*/: return "?";
    }
  }
};

EntityType entity_type_from_string(const std::string& t) {
  if (t == "polymer")     return EntityType::Polymer;
  if (t == "non-polymer") return EntityType::NonPolymer;
  if (t == "water")       return EntityType::Water;
  return EntityType::Unknown;
}

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

  Residue(int id, int auth_id, char ins, std::string rname) noexcept
    : seq_id(id), auth_seq_id(auth_id), ins_code(ins), name(rname) {}
  int seq_id_for_pdb() const {
    return auth_seq_id != UnknownId ? auth_seq_id : seq_id;
  }
  bool has_standard_pdb_name() const;
  std::vector<Atom>& children() { return atoms; }
  const std::vector<Atom>& children() const { return atoms; }
};

struct Chain {
  std::string name;
  std::string auth_name; // not guaranteed to be the same for the whole chain?
  std::vector<Residue> residues;
  Entity *entity = nullptr;
  Model* parent = nullptr;

  explicit Chain(std::string cname) noexcept : name(cname) {}
  Residue* find_residue(int seq_id, int auth_seq_id, char icode,
                        const std::string& chem);
  Residue* find_or_add_residue(int seq_id, int auth_seq_id, char icode,
                               const std::string& name);
  std::vector<Residue>& children() { return residues; }
  const std::vector<Residue>& children() const { return residues; }
};

struct Model {
  std::string name;  // actually an integer number
  std::vector<Chain> chains;
  Structure* parent = nullptr;
  explicit Model(std::string mname) noexcept : name(mname) {}

  Chain* find_chain(const std::string& chain_name) {
    return internal::find_or_null(chains, chain_name);
  }
  Chain* find_or_add_chain(const std::string& chain_name) {
    return internal::find_or_add(chains, chain_name);
  }
  std::vector<Chain>& children() { return chains; }
  const std::vector<Chain>& children() const { return chains; }
};

struct NcsOp {
  bool given;
  // may be changed to tranformation matrix, or quaternion+vector
  Matrix33 rot;
  Position tran;
};

struct Structure {
  std::string name;
  UnitCell cell;
  std::string sg_hm;
  std::vector<Model> models;
  std::vector<NcsOp> ncs;
  std::vector<std::unique_ptr<Entity>> entities;

  // Minimal metadata with keys being mmcif tags: _entry.id, _exptl.method, ...
  std::map<std::string, std::string> info;

  const char* get_info(const std::string& tag, const char* def=nullptr) const {
    auto it = info.find(tag);
    return it != info.end() ? it->second.c_str() : def;
  }
  Model* find_model(const std::string& model_name) {
    return internal::find_or_null(models, model_name);
  }
  Model* find_or_add_model(const std::string& model_name) {
    return internal::find_or_add(models, model_name);
  }

  Entity* find_entity(const std::string& ent_id) {
    for (auto& ent : entities)
      if (ent->id == ent_id)
        return ent.get();
    return nullptr;
  }
  Entity* find_or_add_entity(const std::string& ent_id) {
    Entity* ent = find_entity(ent_id);
    if (!ent) {
      ent = new Entity(ent_id, EntityType::Unknown);
      entities.emplace_back(ent);
    }
    return ent;
  }

  std::vector<Model>& children() { return models; }
  const std::vector<Model>& children() const { return models; }
  void finish();
};


inline Residue* Chain::find_residue(int seq_id, int auth_seq_id,
                                    char icode, const std::string& chem) {
  auto it = std::find_if(residues.begin(), residues.end(),
                         [&](const Residue& r) {
      return r.seq_id == seq_id && r.name == chem &&
             (r.seq_id != Residue::UnknownId ||
              (r.auth_seq_id == auth_seq_id && r.ins_code == icode));
  });
  return it != residues.end() ? &*it : nullptr;
}

inline Residue* Chain::find_or_add_residue(int seq_id, int auth_seq_id,
                                         char icode, const std::string& chem) {
  Residue* r = find_residue(seq_id, auth_seq_id, icode, chem);
  if (r)
    return r;
  residues.emplace_back(seq_id, auth_seq_id, icode, chem);
  return &residues.back();
}

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


template<class T> void add_backlinks(T& obj) {
  for (auto& child : obj.children()) {
    child.parent = &obj;
    add_backlinks(child);
  }
}
template<> void add_backlinks(Atom&) {}


template<class T> size_t count_atom_sites(const T& obj) {
  size_t sum = 0;
  for (const auto& child : obj.children())
    sum += count_atom_sites(child);
  return sum;
}
template<> size_t count_atom_sites(const Residue& res) {
  return res.atoms.size();
}


template<class T> double count_occupancies(const T& obj) {
  double sum = 0;
  for (const auto& child : obj.children())
    sum += count_occupancies(child);
  return sum;
}
template<> double count_occupancies(const Atom& atom) { return atom.occ; }

inline void Structure::finish() {
  add_backlinks(*this);
  // if "entities" were not specifed, deduce them based on sequence
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
