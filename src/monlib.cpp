// Copyright 2018-2023 Global Phasing Ltd.

#include <gemmi/monlib.hpp>
#include <gemmi/calculate.hpp>  // for calculate_chiral_volume
#include <gemmi/modify.hpp>     // for rename_atom_names

namespace gemmi {

int ChemLink::calculate_score(const Residue& res1, const Residue* res2,
                              char alt, char alt2,
                              const ChemComp::Aliasing* aliasing1,
                              const ChemComp::Aliasing* aliasing2) const {
  int link_score = side1.specificity() + side2.specificity();

  auto get_from = [&](const Restraints::AtomId& atomid) {
    auto aliasing = (atomid.comp != 2 || res2 == nullptr) ? aliasing1 : aliasing2;
    if (aliasing)
      if (const std::string* real_id = aliasing->name_from_alias(atomid.atom))
        return Restraints::AtomId{atomid.comp, *real_id}.get_from(res1, res2, alt, alt2);
    return atomid.get_from(res1, res2, alt, alt2);
  };

  // check chirality
  for (const Restraints::Chirality& chirality : rt.chirs)
    if (chirality.sign != ChiralityType::Both) {
      const Atom* a1 = get_from(chirality.id_ctr);
      const Atom* a2 = get_from(chirality.id1);
      const Atom* a3 = get_from(chirality.id2);
      const Atom* a4 = get_from(chirality.id3);
      if (a1 && a2 && a3 && a4) {
        double vol = calculate_chiral_volume(a1->pos, a2->pos,
                                             a3->pos, a4->pos);
        if (chirality.is_wrong(vol))
          link_score -= 10;
      }
    }
  // check fixed torsion angle (_chem_link_tor.period == 0)
  for (const Restraints::Torsion& tor : rt.torsions)
    if (tor.period == 0) {
      const Atom* a1 = get_from(tor.id1);
      const Atom* a2 = get_from(tor.id2);
      const Atom* a3 = get_from(tor.id3);
      const Atom* a4 = get_from(tor.id4);
      double z = 10.;
      if (a1 && a2 && a3 && a4)
        z = angle_z(calculate_dihedral(a1->pos, a2->pos, a3->pos, a4->pos),
                    tor);
      link_score -= (int) z;
    }
  return link_score;
}

namespace {

Restraints read_link_restraints(const cif::Block& block_) {
  auto read_aid = [](cif::Table::Row& row, int n) {
    return Restraints::AtomId{cif::as_int(row[n]), row.str(n+1)};
  };
  Restraints rt;
  cif::Block& block = const_cast<cif::Block&>(block_);
  for (auto row : block.find("_chem_link_bond.",
                             {"atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "type",
                              "value_dist", "value_dist_esd"}))
    rt.bonds.push_back({read_aid(row, 0), read_aid(row, 2),
                        bond_type_from_string(row[4]), false,
                        cif::as_number(row[5]), cif::as_number(row[6]),
                        NAN, NAN});
  for (auto row : block.find("_chem_link_angle.",
                             {"atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "value_angle", "value_angle_esd"}))
    rt.angles.push_back({read_aid(row, 0), read_aid(row, 2), read_aid(row, 4),
                         cif::as_number(row[6]), cif::as_number(row[7])});
  for (auto row : block.find("_chem_link_tor.",
                             {"id",
                              "atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "atom_4_comp_id", "atom_id_4",
                              "value_angle", "value_angle_esd",
                              "period"}))
    rt.torsions.push_back({row.str(0),
                           read_aid(row, 1), read_aid(row, 3),
                           read_aid(row, 5), read_aid(row, 7),
                           cif::as_number(row[9]), cif::as_number(row[10]),
                           cif::as_int(row[11])});
  for (auto row : block.find("_chem_link_chir.",
                             {"atom_centre_comp_id", "atom_id_centre",
                              "atom_1_comp_id", "atom_id_1",
                              "atom_2_comp_id", "atom_id_2",
                              "atom_3_comp_id", "atom_id_3",
                              "volume_sign"}))
    if (row[4][0] != 'c') // ignore crossN
      rt.chirs.push_back({read_aid(row, 0), read_aid(row, 2),
                          read_aid(row, 4), read_aid(row, 6),
                          chirality_from_string(row[8])});
  for (auto row : block.find("_chem_link_plane.",
                             {"plane_id", "atom_comp_id", "atom_id",
                              "dist_esd"})) {
    Restraints::Plane& plane = rt.get_or_add_plane(row.str(0));
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[3]);
    plane.ids.push_back(read_aid(row, 1));
  }
  return rt;
}

void insert_chemlinks_into(const cif::Document& doc, std::map<std::string,ChemLink>& links) {
  const cif::Block* list_block = doc.find_block("link_list");
  auto use_chem_link = [](const cif::Block& block, ChemLink& link) {
    for (auto row : const_cast<cif::Block&>(block).find("_chem_link.",
                                   {"id", "?name",
                                    "?comp_id_1", "?mod_id_1", "?group_comp_1",
                                    "?comp_id_2", "?mod_id_2", "?group_comp_2"})) {
      if (row.str(0) == link.id) {
        if (row.has2(1))
          link.name = row.str(1);
        if (row.has2(2))
          link.side1.comp = row.str(2);
        if (row.has2(3))
          link.side1.mod = row.str(3);
        if (row.has2(4))
          link.side1.group = ChemComp::read_group(row[4]);
        if (row.has2(5))
          link.side2.comp = row.str(5);
        if (row.has2(6))
          link.side2.mod = row.str(6);
        link.side2.group = ChemComp::read_group(row[7]);
        break;
      }
    }
  };
  for (const cif::Block& block : doc.blocks)
    if (starts_with(block.name, "link_") && block.name != "link_list") {
      ChemLink link;
      link.id = block.name.substr(5);
      if (list_block)
        use_chem_link(*list_block, link);
      use_chem_link(block, link);
      link.rt = read_link_restraints(block);
      link.block = block;
      links.emplace(link.id, link);
    }
}

// Helper function. str is one of "add", "delete", "change".
int chem_mod_type(const std::string& str) {
  char c = str[0] | 0x20;
  if (c != 'a' && c != 'd' && c != 'c')
    fail("Unexpected value of _chem_mod_*.function: " + str);
  return c;
}

Restraints read_restraint_modifications(const cif::Block& block_) {
  Restraints rt;
  cif::Block& block = const_cast<cif::Block&>(block_);
  for (auto row : block.find("_chem_mod_bond.",
                             {"function", "atom_id_1", "atom_id_2",
                              "new_type",
                              "new_value_dist", "new_value_dist_esd",
                              "?new_value_dist_nucleus", "?new_value_dist_nucleus_esd"}))
    rt.bonds.push_back({Restraints::AtomId{chem_mod_type(row[0]), row.str(1)},
                        Restraints::AtomId{1, row.str(2)},
                        bond_type_from_string(row[3]), false,
                        cif::as_number(row[4]), cif::as_number(row[5]),
                        row.has(6) ? cif::as_number(row[6]) : NAN,
                        row.has(7) ? cif::as_number(row[7]) : NAN});
  for (auto row : block.find("_chem_mod_angle.",
                             {"function", "atom_id_1",
                              "atom_id_2", "atom_id_3",
                              "new_value_angle", "new_value_angle_esd"}))
    rt.angles.push_back({{chem_mod_type(row[0]), row.str(1)},
                         {1, row.str(2)}, {1, row.str(3)},
                         cif::as_number(row[4]), cif::as_number(row[5])});
  for (auto row : block.find("_chem_mod_tor.",
                              {"function", "?id", "atom_id_1",
                               "atom_id_2", "atom_id_3", "atom_id_4",
                               "new_value_angle", "new_value_angle_esd",
                               "?new_period"}))
    rt.torsions.push_back({row.has(1) ? row.str(1) : "",
                           {chem_mod_type(row[0]), row.str(2)},
                           {1, row.str(3)}, {1, row.str(4)}, {1, row.str(5)},
                           cif::as_number(row[6]), cif::as_number(row[7]),
                           row.has(8) ? cif::as_int(row[8], 0) : -1});
  for (auto row : block.find("_chem_mod_chir.",
                             {"function", "atom_id_centre", "atom_id_1",
                              "atom_id_2", "atom_id_3",
                              "new_volume_sign"}))
    rt.chirs.push_back({{1, row.str(1)}, {chem_mod_type(row[0]), row.str(2)},
                        {1, row.str(3)}, {1, row.str(4)},
                        chirality_from_string(row[5])});
  for (auto row : block.find("_chem_mod_plane_atom.",
                             {"function", "plane_id", "atom_id" ,
                              "new_dist_esd"})) {
    Restraints::Plane& plane = rt.get_or_add_plane(row.str(1));
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[3]);
    plane.ids.push_back({chem_mod_type(row[0]), row.str(2)});
  }
  return rt;
}

void insert_chemmods_into(const cif::Document& doc, std::map<std::string, ChemMod>& mods) {
  const cif::Block* list_block = doc.find_block("mod_list");
  auto use_chem_mod = [](const cif::Block& block, ChemMod& mod) {
    for (auto row : const_cast<cif::Block&>(block).find("_chem_mod.",
                                    {"id", "?name", "?comp_id", "?group_id"}))
      if (row.str(0) == mod.id) {
        if (row.has2(1))
          mod.name = row.str(1);
        if (row.has2(2))
          mod.comp_id = row.str(2);
        if (row.has2(3))
          mod.group_id = row.str(3);
        break;
      }
  };
  for (const cif::Block& block : doc.blocks)
    if (starts_with(block.name, "mod_") && block.name != "mod_list") {
      ChemMod mod;
      mod.id = block.name.substr(4);
      if (list_block)
        use_chem_mod(*list_block, mod);
      use_chem_mod(block, mod);
      for (auto ra : const_cast<cif::Block&>(block).find("_chem_mod_atom.",
                                {"function", "atom_id", "new_atom_id",
                                 "new_type_symbol", "new_type_energy",
                                 "?new_charge", "?new_partial_charge"}))
        mod.atom_mods.push_back({chem_mod_type(ra[0]), ra.str(1), ra.str(2),
                                 Element(ra.str(3)),
                                 (float) cif::as_number(ra.one_of(5, 6)),
                                 ra.str(4)});
      mod.rt = read_restraint_modifications(block);
      mod.block = block;
      mods.emplace(mod.id, mod);
    }
}

} // anonymous namespace

void ChemMod::apply_to(ChemComp& chemcomp, ChemComp::Group alias_group) const {
  auto real = [&chemcomp, alias_group](const std::string& atom_id) -> const std::string& {
    if (alias_group != ChemComp::Group::Null) {
      const ChemComp::Aliasing& aliasing = chemcomp.get_aliasing(alias_group);
      if (const std::string* real_id = aliasing.name_from_alias(atom_id))
        return *real_id;
    }
    return atom_id;
  };
  cat_to(chemcomp.name, '+', id);
  // _chem_mod_atom
  for (const AtomMod& mod : atom_mods) {
    if (mod.func == 'a') {
      if (mod.new_id.empty())
        fail("New atom id is not given");
      if (!chemcomp.has_atom(real(mod.new_id)))
        chemcomp.atoms.push_back({mod.new_id, "", mod.el,
                                  std::isnan(mod.charge) ? mod.charge : 0,
                                  mod.chem_type, Position()});
      continue;
    }
    const std::string& old = real(mod.old_id);
    auto it = chemcomp.find_atom(old);
    switch (mod.func) {
      case 'd':
        if (it != chemcomp.atoms.end()) {
          chemcomp.atoms.erase(it);
          // delete restraints containing mod.old_id
          vector_remove_if(chemcomp.rt.bonds, [&](const Restraints::Bond& b) {
              return b.id1 == old || b.id2 == old;
          });
          vector_remove_if(chemcomp.rt.angles, [&](const Restraints::Angle& a) {
              return a.id1 == old || a.id2 == old || a.id3 == old;
          });
          vector_remove_if(chemcomp.rt.torsions,
              [&](const Restraints::Torsion& t) {
                return t.id1 == old || t.id2 == old || t.id3 == old ||
                       t.id4 == old;
          });
          vector_remove_if(chemcomp.rt.chirs,
              [&](const Restraints::Chirality& c) {
                return c.id_ctr == old || c.id1 == old || c.id2 == old ||
                       c.id3 == old;
          });
          for (Restraints::Plane& plane : chemcomp.rt.planes)
            vector_remove_if(plane.ids, [&](const Restraints::AtomId& a) {
                return a.atom == old;
            });
        }
        break;
      case 'c':
        if (it != chemcomp.atoms.end()) {
          // the modification shouldn't change the atom name, so we don't do:
          // if (!mod.new_id.empty())
          //   it->id = mod.new_id;
          if (mod.el != El::X)
            it->el = mod.el;
          if (!std::isnan(mod.charge))
            it->charge = mod.charge;
          if (!mod.chem_type.empty())
            it->chem_type = mod.chem_type;
        }
        break;
    }
  }

  // _chem_mod_bond
  for (const Restraints::Bond& mod : rt.bonds) {
    auto it = chemcomp.rt.find_bond(real(mod.id1.atom), real(mod.id2.atom));
    switch (mod.id1.comp) {
      case 'a':
        if (it == chemcomp.rt.bonds.end()) {
          chemcomp.rt.bonds.push_back(mod);
          // id1.comp was temporarily set to 'a', set it back to 1
          chemcomp.rt.bonds.back().id1.comp = 1;
        }
        break;
      case 'd':
        if (it != chemcomp.rt.bonds.end())
          chemcomp.rt.bonds.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.bonds.end()) {
          if (mod.type != BondType::Unspec)
            it->type = mod.type;
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
          if (!std::isnan(mod.value_nucleus))
            it->value_nucleus = mod.value_nucleus;
          if (!std::isnan(mod.esd_nucleus))
            it->esd_nucleus = mod.esd_nucleus;
        }
        break;
    }
  }

  // _chem_mod_angle
  for (const Restraints::Angle& mod : rt.angles) {
    auto it = chemcomp.rt.find_angle(real(mod.id1.atom),
                                     real(mod.id2.atom),
                                     real(mod.id3.atom));
    switch (mod.id1.comp) {
      case 'a':
        if (it == chemcomp.rt.angles.end()) {
          chemcomp.rt.angles.push_back(mod);
          chemcomp.rt.angles.back().id1.comp = 1;
        }
        break;
      case 'd':
        if (it != chemcomp.rt.angles.end())
          chemcomp.rt.angles.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.angles.end()) {
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
        }
        break;
    }
  }

  // _chem_mod_tor
  for (const Restraints::Torsion& mod : rt.torsions) {
    auto it = chemcomp.rt.find_torsion(real(mod.id1.atom), real(mod.id2.atom),
                                       real(mod.id3.atom), real(mod.id4.atom));
    switch (mod.id1.comp) {
      case 'a':
        if (it == chemcomp.rt.torsions.end()) {
          chemcomp.rt.torsions.push_back(mod);
          chemcomp.rt.torsions.back().id1.comp = 1;
        }
        break;
      case 'd':
        if (it != chemcomp.rt.torsions.end())
          chemcomp.rt.torsions.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.torsions.end()) {
          if (!mod.label.empty())
            it->label = mod.label;
          if (!std::isnan(mod.value))
            it->value = mod.value;
          if (!std::isnan(mod.esd))
            it->esd = mod.esd;
          if (mod.period != -1)
            it->period = mod.period;
        }
        break;
    }
  }

  // _chem_mod_chir
  for (const Restraints::Chirality& mod : rt.chirs) {
    auto it = chemcomp.rt.find_chir(real(mod.id_ctr.atom), real(mod.id1.atom),
                                    real(mod.id2.atom), real(mod.id3.atom));
    switch (mod.id1.comp) {
      case 'a':
        if (it == chemcomp.rt.chirs.end()) {
          chemcomp.rt.chirs.push_back(mod);
          chemcomp.rt.chirs.back().id1.comp = 1;
        }
        break;
      case 'd':
        if (it != chemcomp.rt.chirs.end())
          chemcomp.rt.chirs.erase(it);
        break;
      case 'c':
        if (it != chemcomp.rt.chirs.end())
          it->sign = mod.sign;
        break;
    }
  }

  // _chem_mod_plane_atom
  for (const Restraints::Plane& mod : rt.planes)
    for (const Restraints::AtomId& atom_id : mod.ids) {
      const std::string& real_id = real(atom_id.atom);
      if (atom_id.comp == 'a') {
        Restraints::Plane& plane = chemcomp.rt.get_or_add_plane(mod.label);
        if (plane.esd == 0.0 && !std::isnan(mod.esd))
          plane.esd = mod.esd;
        auto it = std::find(plane.ids.begin(), plane.ids.end(), real_id);
        if (it == plane.ids.end())
          plane.ids.push_back({1, real_id});
      } else if (atom_id.comp == 'd') {
        auto it = chemcomp.rt.get_plane(mod.label);
        if (it != chemcomp.rt.planes.end()) {
          auto item = std::find(it->ids.begin(), it->ids.end(), real_id);
          if (item != it->ids.end())
            it->ids.erase(item);
        }
      }
    }
}

void MonLib::read_monomer_doc(const cif::Document& doc) {
  // ChemComp
  if (const cif::Block* block = doc.find_block("comp_list"))
    for (auto row : const_cast<cif::Block*>(block)->find("_chem_comp.", {"id", "group"}))
      cc_groups.emplace(row.str(0), ChemComp::read_group(row.str(1)));
  for (const cif::Block& block : doc.blocks)
    add_monomer_if_present(block);
  // ChemLink
  insert_chemlinks_into(doc, links);
  // ChemMod
  insert_chemmods_into(doc, modifications);
}

namespace {

const std::string* find_chemtype(const std::map<std::string, ChemComp>& monomers,
                                 const const_CRA& cra) {
  auto cc = monomers.find(cra.residue->name);
  if (cc != monomers.end()) {
    auto cc_atom = cc->second.find_atom(cra.atom->name);
    if (cc_atom != cc->second.atoms.end())
      return &cc_atom->chem_type;
  }
  return nullptr;
}

// Searches data from _lib_atom in ener_lib.cif.
double find_radius(const EnerLib& ener_lib, const std::string& chemtype,
                   EnerLib::RadiusType type) {
  double r = 0;
  auto it = ener_lib.atoms.find(chemtype);
  if (it != ener_lib.atoms.end()) {
    if (type == EnerLib::RadiusType::Ion)
      r = it->second.ion_radius;
    else if (type == EnerLib::RadiusType::Vdw)
      r = it->second.vdw_radius;
    else if (type == EnerLib::RadiusType::Vdwh)
      r = it->second.vdwh_radius;
  }
  return std::isnan(r) ? 0 : r;
}

} // anonymous namespace

double MonLib::find_ideal_distance(const const_CRA& cra1, const const_CRA& cra2) const {
  std::string types[2] = {cra1.atom->element.uname(),
                          cra2.atom->element.uname()};
  if (const std::string* tmp = find_chemtype(monomers, cra1))
    types[0] = *tmp;
  if (const std::string* tmp = find_chemtype(monomers, cra2))
    types[1] = *tmp;
  // if either one is metal, use ion radius + ion radius
  if (cra1.atom->element.is_metal() != cra2.atom->element.is_metal()) {
    // return if ion radius is found for both of them
    double r1 = find_radius(ener_lib, types[0], EnerLib::RadiusType::Ion);
    double r2 = find_radius(ener_lib, types[1], EnerLib::RadiusType::Ion);
    if (r1 > 0 && r2 > 0)
      return r1 + r2;
  }
  // otherwise, look for defined distance or use average of distances
  double r[2] = {0, 0};
  for (int j = 0; j < 2; ++j) {
    auto range = ener_lib.bonds.equal_range(types[j]);
    for (auto i = range.first; i != range.second; ++i)
      if (i->second.atom_type_2 == types[1-j] && !std::isnan(i->second.length))
        return i->second.length;
      else if (i->second.atom_type_2.empty() && i->second.type == BondType::Single)
        r[j] = i->second.length / 2;
    if ((r[j] == 0 || std::isnan(r[j])) && range.first != range.second)
      r[j] = range.first->second.length / 2;
    if (r[j] == 0 || std::isnan(r[j]))
      r[j] = (j==0 ? cra1 : cra2).atom->element.covalent_r();
  }
  return r[0] + r[1];
}

// Returns a multi-line message that can be shown to the user.
// When we have a logging mechanism, the return type will be void.
std::string MonLib::update_old_atom_names(Structure& st) const {
  std::string msg;
  for (const auto& it : monomers)
    // monomers should have only monomers needed for this structure.
    // Few of them (usually none or one) have old names defined.
    if (it.second.has_old_names()) {
      const std::string& resname = it.first;
      const ChemComp& cc = it.second;
      int old_vs_new = 0;
      for (const Model& model : st.models)
        for (const Chain& chain : model.chains)
          for (const Residue& res : chain.residues)
            if (res.name == resname) {
              for (const Atom& atom : res.atoms) {
                if (cc.find_atom(atom.name) == cc.atoms.end())
                  ++old_vs_new;
                if (cc.find_atom_by_old_name(atom.name) == cc.atoms.end())
                  --old_vs_new;
              }
            }
      if (old_vs_new > 0) {
        cat_to(msg, "Updating atom names in ", resname, ':');
        std::map<std::string, std::string> mapping;
        for (const ChemComp::Atom& a : cc.atoms)
          if (!a.old_id.empty() && a.old_id != a.id) {
            mapping.emplace(a.old_id, a.id);
            cat_to(msg, ' ', a.old_id, "->", a.id);
          }
        msg += '\n';
        rename_atom_names(st, resname, mapping);
      }
    }
  return msg;
}

} // namespace gemmi
