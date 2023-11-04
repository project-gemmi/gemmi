// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include <map>           // for multimap
#include <ostream>       // for ostream
#include <memory>        // for unique_ptr
#include <unordered_map> // for unordered_map
#include "chemcomp.hpp"  // for ChemComp
#include "monlib.hpp"    // for MonLib
#include "model.hpp"     // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral

namespace gemmi {

enum class HydrogenChange {
  NoChange, Shift, Remove, ReAdd, ReAddButWater, ReAddKnown
};

struct GEMMI_DLL Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  struct Bond {
    const Restraints::Bond* restr;
    std::array<Atom*, 2> atoms;
    Asu asu;
    double calculate() const {
      return asu != Asu::Different ? atoms[0]->pos.dist(atoms[1]->pos) : NAN;
    }
    double calculate_z_(double d) const { return std::abs(d - restr->value) / restr->esd; }
    double calculate_z() const { return calculate_z_(calculate()); }
  };
  struct Angle {
    const Restraints::Angle* restr;
    std::array<Atom*, 3> atoms;
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    double calculate_z() const { return angle_z(calculate(), *restr); }
  };
  struct Torsion {
    const Restraints::Torsion* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z() const {
      return angle_z(calculate(), *restr, 360. / std::max(1, restr->period));
    }
  };
  struct Chirality {
    const Restraints::Chirality* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z(double ideal_abs_vol, double esd) const {
      double calc = calculate();
      if (restr->sign == ChiralityType::Negative ||
          (restr->sign == ChiralityType::Both && calc < 0))
        ideal_abs_vol *= -1;
      return std::abs(calc - ideal_abs_vol) / esd;
    }
    bool check() const { return !restr->is_wrong(calculate()); }
  };
  struct Plane {
    const Restraints::Plane* restr;
    std::vector<Atom*> atoms;
    bool has(const Atom* atom) const {
      return in_vector(const_cast<Atom*>(atom), atoms);
    }
  };

  enum class RKind { Bond, Angle, Torsion, Chirality, Plane };
  struct Rule {
    RKind rkind;
    size_t index; // index in the respective vector (bonds, ...) in Topo
  };

  struct Link {
    std::string link_id;
    Residue* res1 = nullptr;
    Residue* res2 = nullptr;
    std::vector<Rule> link_rules;
    char alt1 = '\0';
    char alt2 = '\0';
    Asu asu = Asu::Any;  // used only in Links in ChainInfo::extras
    bool is_cis = false;  // helper field for CISPEP record generation

    // aliasing1/2 points to vector element in ChemComp::aliases.
    // The pointers should stay valid even if a ChemComp is moved.
    const ChemComp::Aliasing* aliasing1 = nullptr;
    const ChemComp::Aliasing* aliasing2 = nullptr;

    // only for polymer links, res1 and res2 must be in the same vector (Chain)
    std::ptrdiff_t res_distance() const { return res1 - res2; }
  };

  struct Mod {
    std::string id;         // id of ChemMod from the dictionary (MonLib)
    ChemComp::Group alias;  // alias to be used when applying the modification
    char altloc;            // \0 = all conformers

    bool operator==(const Mod& o) const {
      return id == o.id && alias == o.alias && altloc == o.altloc;
    }
  };

  struct FinalChemComp {
    char altloc;  // Restraints apply to this conformer
    const ChemComp* cc;
  };

  struct ResInfo {
    Residue* res;
    // in case of microheterogeneity we may have 2+ previous residues
    std::vector<Link> prev;
    std::vector<Mod> mods;
    // Pointer to ChemComp in MonLib::monomers.
    const ChemComp* orig_chemcomp = nullptr;
    // Pointer to restraints with modifications applied (if any).
    std::vector<FinalChemComp> chemcomps;
    std::vector<Rule> monomer_rules;

    ResInfo(Residue* r) : res(r) {}
    void add_mod(const std::string& m, const ChemComp::Aliasing* aliasing, char altloc) {
      if (!m.empty()) {
        auto alias_group = aliasing ? aliasing->group : ChemComp::Group::Null;
        Mod mod{m, alias_group, altloc};
        if (!in_vector(mod, mods))
          mods.push_back(mod);
      }
    }

    const ChemComp& get_final_chemcomp(char altloc) const {
      if (chemcomps.size() == 1)
        return *chemcomps[0].cc;
      assert(!chemcomps.empty());
      for (const FinalChemComp& it : chemcomps)
        if (it.altloc == altloc)
          return *it.cc;
      return *chemcomps[0].cc;
    }
  };

  // corresponds to a sub-chain
  struct ChainInfo {
    const Chain& chain_ref;
    std::string subchain_name;
    std::string entity_id;
    bool polymer;
    PolymerType polymer_type;
    std::vector<ResInfo> res_infos;

    ChainInfo(ResidueSpan& subchain, const Chain& chain, const Entity* ent);
    using iterator = std::vector<ResInfo>::iterator;
    iterator group_end(iterator b) const {
      auto e = b + 1;
      while (e != res_infos.end() && e->res->group_key() == b->res->group_key())
        ++e;
      return e;
    }
  };

  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  std::ostream* warnings = nullptr;
  bool only_bonds = false;  // an internal flag for apply_restraints()
  std::vector<ChainInfo> chain_infos;
  std::vector<Link> extras;

  // Restraints applied to Model
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  std::multimap<const Atom*, Bond*> bond_index;       // indexes both atoms
  std::multimap<const Atom*, Angle*> angle_index;     // only middle atom
  std::multimap<const Atom*, Torsion*> torsion_index; // two middle atoms
  std::multimap<const Atom*, Plane*> plane_index;     // all atoms

  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chain_infos)
      for (ResInfo& ri : ci.res_infos)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  Bond* first_bond_in_link(const Link& link) {
    for (const Rule& rule : link.link_rules)
      if (rule.rkind == RKind::Bond)
        return &bonds[rule.index];
    return nullptr;
  }

  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    auto range = bond_index.equal_range(a);
    for (auto i = range.first; i != range.second; ++i) {
      const Bond* bond = i->second;
      if ((bond->atoms[0] == b && bond->atoms[1] == a) ||
          (bond->atoms[1] == b && bond->atoms[0] == a))
        return bond->restr;
    }
    return nullptr;
  }

  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    auto range = angle_index.equal_range(b);
    for (auto i = range.first; i != range.second; ++i) {
      const Angle* ang = i->second;
      if ((ang->atoms[0] == a && ang->atoms[2] == c) ||
          (ang->atoms[0] == c && ang->atoms[2] == a))
        return ang->restr;
    }
    return nullptr;
  }

  const Chirality* get_chirality(const Atom* ctr) const {
    for (const Chirality& chir : chirs)
      if (chir.atoms[0] == ctr)
        return &chir;
    return nullptr;
  }

  double ideal_chiral_abs_volume(const Chirality &ch) const;

  std::vector<Rule> apply_restraints(const Restraints& rt,
                                     Residue& res, Residue* res2, Asu asu,
                                     char altloc1, char altloc2, bool require_alt);
  void apply_restraints_from_link(Link& link, const MonLib& monlib);

  // Structure is non-const b/c connections may have link_id assigned.
  // Model is non-const b/c we store non-const pointers to residues in Topo.
  // Because of the pointers, don't add or remove residues after this step.
  // Monlib may get modified by addition of extra links from the model.
  void initialize_refmac_topology(Structure& st, Model& model0,
                                  MonLib& monlib, bool ignore_unknown_links=false);

  // This step stores pointers to gemmi::Atom's from model0,
  // so after this step don't add or remove atoms.
  // monlib is needed only for links.
  void apply_all_restraints(const MonLib& monlib);

  // prepare bond_index, angle_index, torsion_index, plane_index
  void create_indices();

  Link* find_polymer_link(const AtomAddress& a1, const AtomAddress& a2) {
    for (ChainInfo& ci : chain_infos)
      if (a1.chain_name == ci.chain_ref.name && a2.chain_name == ci.chain_ref.name) {
        for (ResInfo& ri : ci.res_infos)
          for (Link& link : ri.prev) {
            assert(link.res1 && link.res2);
            if ((a1.res_id.matches_noseg(*link.res1) &&
                 a2.res_id.matches_noseg(*link.res2) &&
                 a1.altloc == link.alt1 && a2.altloc == link.alt2) ||
                (a2.res_id.matches_noseg(*link.res1) &&
                 a1.res_id.matches_noseg(*link.res2) &&
                 a1.altloc == link.alt2 && a2.altloc == link.alt1))
              return &link;
          }
      }
    return nullptr;
  }

  void set_cispeps_in_structure(Structure& st);

  GEMMI_COLD void err(const std::string& msg) const {
    if (warnings == nullptr)
      fail(msg);
    *warnings << "Warning: " << msg << std::endl;
  }

private:
  // storage for link restraints modified by aliases
  std::vector<std::unique_ptr<Restraints>> rt_storage;
  // cache for ChemComps after applying modifications
  std::unordered_map<std::string, std::unique_ptr<ChemComp>> cc_cache;
  // storage for ad-hoc ChemComps (placeholders for those missing in MonLib)
  std::vector<std::unique_ptr<ChemComp>> cc_storage;

  void setup_connection(Connection& conn, Model& model0, MonLib& monlib,
                        bool ignore_unknown_links);
};

GEMMI_DLL std::unique_ptr<Topo>
prepare_topology(Structure& st, MonLib& monlib, size_t model_index,
                 HydrogenChange h_change, bool reorder,
                 std::ostream* warnings=nullptr, bool ignore_unknown_links=false,
                 bool use_cispeps=false);


GEMMI_DLL std::unique_ptr<ChemComp> make_chemcomp_with_restraints(const Residue& res);

GEMMI_DLL std::vector<AtomAddress> find_missing_atoms(const Topo& topo,
                                                      bool including_hydrogen=false);

} // namespace gemmi
#endif
