// Copyright 2022 Global Phasing Ltd.
//
// Create cif::Block with monomer library _chem_comp* categories
// from struct ChemComp.

#ifndef GEMMI_TO_CHEMCOMP_HPP_
#define GEMMI_TO_CHEMCOMP_HPP_

#include "chemcomp.hpp"  // for ChemComp
#include "sprintf.hpp"   // for to_str

namespace gemmi {

inline void add_chemcomp_to_block(const ChemComp& cc, cif::Block& block) {
  {
    std::vector<std::string> tags =
        {"comp_id", "atom_id", "type_symbol", "type_energy", "charge"};
    if (cc.has_coordinates)
      for (char c = 'x'; c <= 'z'; ++c)
        tags.emplace_back(1, c);
    cif::Table tab = block.find_or_add("_chem_comp_atom.", tags);
    if (!tab.loop_item)
      tab.convert_pair_to_loop();
    size_t pos = tab.length();
    cif::Loop& loop = tab.loop_item->loop;
    loop.values.resize(loop.values.size() + loop.width() * cc.atoms.size(), ".");
    for (const ChemComp::Atom& a : cc.atoms) {
      cif::Table::Row row = tab[pos++];
      row[0] = cc.name;
      row[1] = a.id;
      row[2] = a.el.name();
      row[3] = cif::quote(a.chem_type);
      row[4] = std::to_string(iround(a.charge));
      if (cc.has_coordinates) {
        row[5] = to_str(a.xyz.x);
        row[6] = to_str(a.xyz.y);
        row[7] = to_str(a.xyz.z);
      }
    }
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_bond.",
        {"comp_id", "atom_id_1", "atom_id_2", "type", "aromatic",
         "value_dist", "value_dist_esd", "value_dist_nucleus", "value_dist_nucleus_esd"});
    for (const Restraints::Bond& a : cc.rt.bonds)
      tab.append_row({cc.name, a.id1.atom, a.id2.atom, bond_type_to_string(a.type),
                      std::string(1, a.aromatic ? 'y' : 'n'),
                      to_str(a.value), to_str(a.esd),
                      to_str(a.value_nucleus), to_str(a.esd_nucleus)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_angle.",
        {"comp_id", "atom_id_1", "atom_id_2", "atom_id_3",
         "value_angle", "value_angle_esd"});
    for (const Restraints::Angle& a : cc.rt.angles)
      tab.append_row({cc.name, a.id1.atom, a.id2.atom, a.id3.atom,
                      to_str(a.value), to_str(a.esd)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_tor.",
        {"comp_id", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
         "value_angle", "value_angle_esd", "period"});
    for (const Restraints::Torsion& a : cc.rt.torsions)
      tab.append_row({cc.name, a.label, a.id1.atom, a.id2.atom, a.id3.atom, a.id4.atom,
                      to_str(a.value), to_str(a.esd), std::to_string(a.period)});
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_chir.",
        {"comp_id", "id", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3",
         "volume_sign"});
    for (const Restraints::Chirality& a : cc.rt.chirs) {
      std::string label = "chir_" + std::to_string(tab.length() + 1);
      tab.append_row({cc.name, label, a.id_ctr.atom, a.id1.atom, a.id2.atom, a.id3.atom,
                      chirality_to_string(a.sign)});
    }
  }
  {
    cif::Table tab = block.find_or_add("_chem_comp_plane_atom.",
        {"comp_id", "plane_id", "atom_id", "dist_esd"});
    for (const Restraints::Plane& p : cc.rt.planes)
      for (const Restraints::AtomId& atom_id : p.ids)
        tab.append_row({cc.name, p.label, atom_id.atom, to_str(p.esd)});
  }
}

} // namespace gemmi
#endif
