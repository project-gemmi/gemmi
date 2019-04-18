// Copyright 2019 Global Phasing Ltd.
//
// Searching for links based on the _chem_link table from monomer dictionary.

#ifndef GEMMI_LINKHUNT_HPP_
#define GEMMI_LINKHUNT_HPP_

#include <map>
#include <unordered_map>
#include "monlib.hpp"
#include "subcells.hpp"

namespace gemmi {

struct LinkHunt {
  struct Match {
    const ChemLink* chem_link;
    CRA cra1;
    CRA cra2;
    bool same_asu;
    float bond_length;
  };

  double global_max_dist = 2.34; // ZN-CYS
  std::multimap<std::string, const ChemLink*> links;
  std::unordered_map<std::string, ChemLink::Group> res_group;
  std::map<std::string, double> max_dist_per_atom;

  void index_chem_links(const MonLib& monlib) {
    static const std::vector<std::string> blacklist = {
      "TRANS", "PTRANS", "NMTRANS", "CIS", "PCIS", "NMCIS", "p"
    };
    for (const auto& iter : monlib.links) {
      const ChemLink& link = iter.second;
      if (link.rt.bonds.empty())
        continue;
      if (link.rt.bonds.size() > 1)
        fprintf(stderr, "Note: considering only the first bond in %s\n",
                link.id.c_str());
      if (link.side1.comp.empty() && link.side2.comp.empty())
        if (link.side1.group == ChemLink::Group::Null ||
            link.side2.group == ChemLink::Group::Null ||
            in_vector(link.id, blacklist))
          continue;
      const Restraints::Bond& bond = link.rt.bonds[0];
      if (bond.value > global_max_dist)
        global_max_dist = bond.value;
      for (const std::string& atom_name : {bond.id1.atom, bond.id1.atom}) {
        auto r = max_dist_per_atom.emplace(atom_name, bond.value);
        if (!r.second && r.first->second < bond.value)
          r.first->second = bond.value;
      }
      links.emplace(bond.lexicographic_str(), &link);
    }
    for (const auto& ri : monlib.residue_infos) {
      ChemLink::Group group = ChemLink::Group::Null;
      if (ri.second.is_amino_acid())
        group = ChemLink::Group::Peptide;
      else if (ri.second.is_nucleic_acid())
        group = ChemLink::Group::DnaRna;
      else if (ri.second.kind == ResidueInfo::PYR)
        group = ChemLink::Group::Pyranose;
      res_group.emplace(ri.first, group);
    }
  }

  bool match_link_side(const ChemLink::Side& side,
                       const std::string& resname) const {
    if (!side.comp.empty())
      return side.comp == resname;
    if (side.group == ChemLink::Group::Null)
      return false;
    auto iter = res_group.find(resname);
    return iter != res_group.end() && iter->second == side.group;
  }

  std::vector<Match> find_possible_links(Structure& st, double tolerance) {
    std::vector<Match> results;
    Model& model = st.models.at(0);
    SubCells sc(model, st.cell, std::max(5.0, global_max_dist * tolerance));
    sc.populate(model);
    for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
      Chain& chain = model.chains[n_ch];
      for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
        Residue& res = chain.residues[n_res];
        for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
          Atom& atom = res.atoms[n_atom];
          auto max_dist = max_dist_per_atom.find(atom.name);
          if (max_dist == max_dist_per_atom.end())
            continue;
          sc.for_each(atom.pos, atom.altloc, (float) max_dist->second,
                      [&](SubCells::Mark& m, float dist_sq) {
              // do not consider connections inside a residue
              if (m.image_idx == 0 && m.chain_idx == n_ch &&
                  m.residue_idx == n_res)
                return;
              // avoid reporting connections twice (A-B and B-A)
              if (m.chain_idx < n_ch || (m.chain_idx == n_ch &&
                    (m.residue_idx < n_res || (m.residue_idx == n_res &&
                                               m.atom_idx < n_atom))))
                return;
              // atom can be linked with its image, but if the image
              // is too close the atom is likely on special position.
              if (m.chain_idx == n_ch && m.residue_idx == n_res &&
                  m.atom_idx == n_atom && dist_sq < sq(0.8f))
                return;
              CRA cra = m.to_cra(model);
              auto range = links.equal_range(Restraints::lexicographic_str(
                                                  atom.name, cra.atom->name));
              for (auto iter = range.first; iter != range.second; ++iter) {
                const ChemLink& link = *iter->second;
                const Restraints::Bond& bond = link.rt.bonds[0];
                if (dist_sq > sq(bond.value * tolerance))
                  continue;
                if (bond.id1.atom == atom.name &&
                    match_link_side(link.side1, res.name) &&
                    match_link_side(link.side2, cra.residue->name)) {
                  results.push_back({&link, {&chain, &res, &atom}, cra,
                                     !m.image_idx, std::sqrt(dist_sq)});
                  break;
                }
                if (bond.id2.atom == atom.name &&
                    match_link_side(link.side2, res.name) &&
                    match_link_side(link.side1, cra.residue->name)) {
                  results.push_back({&link, cra, {&chain, &res, &atom},
                                     !m.image_idx, std::sqrt(dist_sq)});
                  break;
                }
              }
          });
        }
      }
    }
    return results;
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
