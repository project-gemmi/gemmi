// Copyright 2018 Global Phasing Ltd.
//
// Topo(logy) - restraints (from a monomer library) applied to a model.

#ifndef GEMMI_TOPO_HPP_
#define GEMMI_TOPO_HPP_

#include "chemcomp.hpp" // for ChemComp
#include "monlib.hpp" // for MonLib
#include "model.hpp" // for Residue, Atom
#include "calculate.hpp" // for calculate_angle, calculate_dihedral
#include "polyheur.hpp" // for are_connected

namespace gemmi {

struct Topo {
  // We have internal pointers in this class (pointers setup in
  // apply_restraints() that point to ResInfo::chemcomp.rt),
  // disable copying this class.
  Topo() = default;
  Topo(Topo const&) = delete;
  Topo& operator=(Topo const&) = delete;

  struct Bond {
    const Restraints::Bond* restr;
    std::array<Atom*, 2> atoms;
    double calculate() const { return atoms[0]->pos.dist(atoms[1]->pos); }
    double calculate_z() const {
      return std::abs(calculate() - restr->value) / restr->esd;
    }
  };
  struct Angle {
    const Restraints::Angle* restr;
    std::array<Atom*, 3> atoms;
    double calculate() const {
      return calculate_angle(atoms[0]->pos, atoms[1]->pos, atoms[2]->pos);
    }
    double calculate_z() const {
      return angle_abs_diff(deg(calculate()), restr->value) / restr->esd;
    }
  };
  struct Torsion {
    const Restraints::Torsion* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_dihedral(atoms[0]->pos, atoms[1]->pos,
                                atoms[2]->pos, atoms[3]->pos);
    }
    double calculate_z() const {
      return angle_abs_diff(deg(calculate()), restr->value) / restr->esd;
    }
  };
  struct Chirality {
    const Restraints::Chirality* restr;
    std::array<Atom*, 4> atoms;
    double calculate() const {
      return calculate_chiral_volume(atoms[0]->pos, atoms[1]->pos,
                                     atoms[2]->pos, atoms[3]->pos);
    }
    // positive value for preserved chirality
    double check() const {
      double value = calculate();
      if ((restr->chir == Restraints::Chirality::Type::Positive && value < 0) ||
          (restr->chir == Restraints::Chirality::Type::Negative && value > 0))
        return -std::abs(value);
      return std::abs(value);
    }
  };
  struct Plane {
    const Restraints::Plane* restr;
    std::vector<Atom*> atoms;
  };

  enum class Provenance { None, PrevLink, Monomer, NextLink, ExtraLink };
  enum class RKind { Bond, Angle, Torsion, Chirality, Plane };
  struct Force {
    Provenance provenance;
    RKind rkind;
    size_t index; // index in the respective vector (bonds, ...) in Topo
  };

  struct ResInfo {
    Residue* res;
    std::string prev_link;
    int prev_idx;
    std::vector<std::string> mods;
    ChemComp chemcomp;
    std::vector<Force> forces;

    ResInfo(Residue* r) : res(r), prev_idx(0) {}
    Residue* prev_res() {
      return prev_idx != 0 ? (this + prev_idx)->res : nullptr;
    }
    const Residue* prev_res() const {
      return const_cast<ResInfo*>(this)->prev_res();
    }
  };

  struct ChainInfo {
    std::string name;
    std::string entity_id;
    bool polymer;
    PolymerType polymer_type;
    std::vector<ResInfo> residues;

    void initialize(SubChain& subchain, const Entity* ent);
    void setup_polymer_links();
    void add_refmac_builtin_modifications();
  };

  struct ExtraLink {
    Residue* res1;
    Residue* res2;
    char alt1 = '\0';
    char alt2 = '\0';
    ChemLink link;
    std::vector<Force> forces;
  };

  template<typename T>
  static int has_atom(const Atom* a, const T& t) {
    for (int i = 0; (size_t) i != t.atoms.size(); ++i)
      if (t.atoms[i] == a)
        return i;
    return -1;
  }

  std::vector<ChainInfo> chains;
  std::vector<ExtraLink> extras;

  // Restraints applied to Model
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  ResInfo* find_resinfo(const Residue* res) {
    for (ChainInfo& ci : chains)
      for (ResInfo& ri : ci.residues)
        if (ri.res == res)
          return &ri;
    return nullptr;
  }

  const Restraints::Bond* take_bond(const Atom* a, const Atom* b) const {
    for (const Bond& bond : bonds)
      if ((bond.atoms[0] == a && bond.atoms[1] == b) ||
          (bond.atoms[0] == b && bond.atoms[1] == a))
        return bond.restr;
    return nullptr;
  }

  const Restraints::Angle* take_angle(const Atom* a,
                                      const Atom* b,
                                      const Atom* c) const {
    for (const Angle& ang : angles)
      if (ang.atoms[1] == b && ((ang.atoms[0] == a && ang.atoms[2] == c) ||
                                (ang.atoms[0] == c && ang.atoms[1] == a)))
        return ang.restr;
    return nullptr;
  }

  std::vector<Force> apply_restraints(const Restraints& rt,
                                      Residue& res, Residue* res2,
                                      char altloc='*') {
    std::string altlocs;
    if (altloc == '*') {
      // find all distinct altlocs
      for (const Atom& atom : res.atoms)
        if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
          altlocs += atom.altloc;
      if (res2)
        for (const Atom& atom : res2->atoms)
          if (atom.altloc && altlocs.find(atom.altloc) == std::string::npos)
            altlocs += atom.altloc;
    }
    if (altlocs.empty())
      altlocs += altloc;

    std::vector<Force> forces;
    Provenance pro = Provenance::None;
    for (const Restraints::Bond& bond : rt.bonds)
      for (char alt : altlocs)
        if (Atom* at1 = bond.id1.get_from(res, res2, alt))
          if (Atom* at2 = bond.id2.get_from(res, res2, alt)) {
            forces.push_back({pro, RKind::Bond, bonds.size()});
            bonds.push_back({&bond, {{at1, at2}}});
            if (!at1->altloc && !at2->altloc)
              break;
          }
    for (const Restraints::Angle& angle : rt.angles)
      for (char alt : altlocs)
        if (Atom* at1 = angle.id1.get_from(res, res2, alt))
          if (Atom* at2 = angle.id2.get_from(res, res2, alt))
            if (Atom* at3 = angle.id3.get_from(res, res2, alt)) {
              forces.push_back({pro, RKind::Angle, angles.size()});
              angles.push_back({&angle, {{at1, at2, at3}}});
              if (!at1->altloc && !at2->altloc && !at3->altloc)
                break;
            }
    for (const Restraints::Torsion& tor : rt.torsions)
      for (char alt : altlocs)
        if (Atom* at1 = tor.id1.get_from(res, res2, alt))
          if (Atom* at2 = tor.id2.get_from(res, res2, alt))
            if (Atom* at3 = tor.id3.get_from(res, res2, alt))
              if (Atom* at4 = tor.id4.get_from(res, res2, alt)) {
                forces.push_back({pro, RKind::Torsion, torsions.size()});
                torsions.push_back({&tor, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
          }
    for (const Restraints::Chirality& chir : rt.chirs)
      for (char alt : altlocs)
        if (Atom* at1 = chir.id_ctr.get_from(res, res2, alt))
          if (Atom* at2 = chir.id1.get_from(res, res2, alt))
            if (Atom* at3 = chir.id2.get_from(res, res2, alt))
              if (Atom* at4 = chir.id3.get_from(res, res2, alt)) {
                forces.push_back({pro, RKind::Chirality, chirs.size()});
                chirs.push_back({&chir, {{at1, at2, at3, at4}}});
                if (!at1->altloc && !at2->altloc &&
                    !at3->altloc && !at4->altloc)
                  break;
              }
    for (const Restraints::Plane& plane : rt.planes)
      for (char alt : altlocs) {
        std::vector<Atom*> atoms;
        for (const Restraints::AtomId& id : plane.ids)
          if (Atom* atom = id.get_from(res, res2, alt))
            atoms.push_back(atom);
        if (atoms.size() >= 4) {
          forces.push_back({pro, RKind::Plane, planes.size()});
          planes.push_back({&plane, atoms});
        }
        if (std::all_of(atoms.begin(), atoms.end(),
                        [](Atom* a) { return !a->altloc; }))
          break;
      }
    return forces;
  }

  void apply_internal_restraints_to_residue(ResInfo& ri) {
    auto forces = apply_restraints(ri.chemcomp.rt, *ri.res, nullptr);
    for (const auto& f : forces)
      ri.forces.push_back({Provenance::Monomer, f.rkind, f.index});
  }

  void apply_restraints_to_residue(ResInfo& ri, const MonLib& monlib) {
    if (Residue* prev = ri.prev_res()) {
      if (const ChemLink* link = monlib.find_link(ri.prev_link)) {
        auto forces = apply_restraints(link->rt, *prev, ri.res);
        for (const auto& f : forces) {
          ri.forces.push_back({Provenance::PrevLink, f.rkind, f.index});
          if (ri.prev_idx != 0)
            (&ri + ri.prev_idx)->forces.push_back({Provenance::NextLink,
                                                   f.rkind, f.index});
        }
      }
    }
    apply_internal_restraints_to_residue(ri);
  }

  void apply_restraints_to_extra_link(ExtraLink& link, const MonLib& monlib) {
    if (const ChemLink* cl = monlib.match_link(link.link)) {
      if (link.alt1 && link.alt2 && link.alt1 != link.alt2)
        printf("Warning: LINK between different conformers %c and %c.",
               link.alt1, link.alt2);
      ResInfo* ri1 = find_resinfo(link.res1);
      ResInfo* ri2 = find_resinfo(link.res2);
      auto forces = apply_restraints(cl->rt, *link.res1, link.res2,
                                     link.alt1);
      for (Force& f : forces) {
        f.provenance = Provenance::ExtraLink;
        link.forces.push_back(f);
        if (ri1)
          ri1->forces.push_back(f);
        if (ri2)
          ri2->forces.push_back(f);
      }
    }
  }

  // Model is non-const b/c we store non-const pointers to residues in Topo.
  void prepare_refmac_topology(Model& model0,
                               const std::vector<Entity>& entities,
                               MonLib& monlib);
};

inline void Topo::ChainInfo::initialize(SubChain& subchain, const Entity* ent) {
  residues.reserve(subchain.size());
  name = subchain.name();
  entity_id = ent->name;
  polymer = ent && ent->entity_type == EntityType::Polymer;
  polymer_type = ent->polymer_type;
  for (Residue& res : subchain)
    residues.emplace_back(&res);
}

inline void Topo::ChainInfo::setup_polymer_links() {
  if (!polymer || residues.empty())
    return;
  for (auto ri = residues.begin() + 1; ri != residues.end(); ++ri) {
    // For now we ignore microheterogeneity.
    ri->prev_idx = -1;
    const Residue* prev_res = (ri + ri->prev_idx)->res;
    if (!prev_res) {
      ri->prev_link = ".";
    } else if (!are_connected(*prev_res, *ri->res, polymer_type)) {
      ri->prev_link = "gap";
    } else if (is_polypeptide(polymer_type)) {
      if (ri->chemcomp.group == "P-peptide")
        ri->prev_link = "P";  // PCIS, PTRANS
      else if (ri->chemcomp.group == "M-peptide")
        ri->prev_link = "NM"; // NMCIS, NMTRANS
      ri->prev_link += prev_res->is_cis ? "CIS" : "TRANS";
    } else if (is_polynucleotide(polymer_type)) {
      ri->prev_link = "p";
    } else {
      ri->prev_link = "?";
    }
  }
}

inline void Topo::ChainInfo::add_refmac_builtin_modifications() {
  if (polymer && !residues.empty()) {
    // we try to get exactly the same numbers that makecif produces
    for (Topo::ResInfo& ri : residues)
      if (polymer_type == PolymerType::PeptideL)
        ri.mods.emplace_back("AA-STAND");
    Topo::ResInfo& front = residues.front();
    Topo::ResInfo& back = residues.back();
    if (is_polypeptide(polymer_type)) {
      front.mods.emplace_back("NH3");
      back.mods.emplace_back(back.res->find_atom("OXT") ? "COO" : "TERMINUS");
    } else if (is_polynucleotide(polymer_type)) {
      front.mods.emplace_back("5*END");
      back.mods.emplace_back("TERMINUS");
    }
  }
}

inline ChemLink connection_to_chemlink(const Connection& conn,
                                       const Residue& res1,
                                       const Residue& res2) {
  ChemLink link;
  link.id = res1.name + "-" + res2.name;
  link.comp[0] = res1.name;
  link.comp[1] = res2.name;
  Restraints::Bond bond;
  bond.id1 = Restraints::AtomId{1, conn.atom[0].atom_name};
  bond.id2 = Restraints::AtomId{2, conn.atom[1].atom_name};
  bond.type = Restraints::BondType::Unspec;
  bond.aromatic = false;
  bond.value = conn.reported_distance;
  bond.esd = 0.02;
  link.rt.bonds.push_back(bond);
  return link;
}


// Model is non-const b/c we store non-const pointers to residues in Topo.
inline void Topo::prepare_refmac_topology(Model& model0,
                                          const std::vector<Entity>& entities,
                                          MonLib& monlib) {
  // initialize chains and residues
  for (Chain& chain : model0.chains)
    for (SubChain sub : chain.subchains()) {
      assert(sub.labelled());
      const Entity* ent = get_entity_of(sub, entities);
      chains.emplace_back();
      chains.back().initialize(sub, ent);
    }
  for (ChainInfo& ci : chains) {
    // copy monomer description
    for (ResInfo& ri : ci.residues)
      ri.chemcomp = monlib.monomers.at(ri.res->name);

    ci.setup_polymer_links();

    ci.add_refmac_builtin_modifications();

    // add modifications from standard links
    for (ResInfo& ri : ci.residues)
      if (const ChemLink* link = monlib.find_link(ri.prev_link)) {
        if (!link->mod[0].empty())
          (&ri + ri.prev_idx)->mods.push_back(link->mod[0]);
        if (!link->mod[1].empty())
          ri.mods.push_back(link->mod[1]);
      }
  }
  // add extra links
  for (const Connection& conn : model0.connections) {
    ExtraLink extra;
    extra.res1 = model0.find_cra(conn.atom[0]).residue;
    extra.res2 = model0.find_cra(conn.atom[1]).residue;
    if (extra.res1 && extra.res2) {
      extra.link = connection_to_chemlink(conn, *extra.res1, *extra.res2);
      if (!monlib.match_link(extra.link)) {
        monlib.ensure_unique_link_name(extra.link.id);
        monlib.links.emplace(extra.link.id, extra.link);
      }
      extra.alt1 = conn.atom[0].altloc;
      extra.alt2 = conn.atom[1].altloc;
      extras.push_back(extra);
    }
  }
  // TODO: automatically determine other links

  for (ChainInfo& chain_info : chains)
    for (ResInfo& ri : chain_info.residues) {
      // apply modifications
      for (const std::string& modif : ri.mods) {
        if (const ChemMod* chem_mod = monlib.find_mod(modif))
          try {
            chem_mod->apply_to(ri.chemcomp);
          } catch(std::runtime_error& e) {
            printf("Failed to apply modification %s to %s: %s\n",
                   chem_mod->id.c_str(), ri.res->name.c_str(), e.what());
          }
        else
          printf("Modification not found: %s\n", modif.c_str());
      }
    }

  // finalize
  for (ChainInfo& chain_info : chains)
    for (ResInfo& ri : chain_info.residues)
      apply_restraints_to_residue(ri, monlib);
  for (ExtraLink& link : extras)
    apply_restraints_to_extra_link(link, monlib);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
