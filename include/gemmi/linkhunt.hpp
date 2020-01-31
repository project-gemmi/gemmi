// Copyright 2019 Global Phasing Ltd.
//
// Searching for links based on the _chem_link table from monomer dictionary.

#ifndef GEMMI_LINKHUNT_HPP_
#define GEMMI_LINKHUNT_HPP_

#include <map>
#include <unordered_map>
#include "calculate.hpp"  // for calculate_chiral_volume
#include "elem.hpp"
#include "model.hpp"
#include "monlib.hpp"
#include "subcells.hpp"

namespace gemmi {

struct LinkHunt {
  struct Match {
    const ChemLink* chem_link = nullptr;
    int chem_link_count = 0;
    int score = -1000;
    CRA cra1;
    CRA cra2;
    bool same_image;
    float bond_length = 0.f;
    Connection* conn = nullptr;
  };

  double global_max_dist = 2.34; // ZN-CYS
  std::multimap<std::string, const ChemLink*> links;
  std::unordered_map<std::string, ChemLink::Group> res_group;

  void index_chem_links(const MonLib& monlib) {
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
            link.id == "SS")
          continue;
      const Restraints::Bond& bond = link.rt.bonds[0];
      if (bond.value > global_max_dist)
        global_max_dist = bond.value;
      links.emplace(bond.lexicographic_str(), &link);
    }
    for (const auto& ri : monlib.residue_infos)
      res_group.emplace(ri.first, ChemLink::group_from_residue_info(ri.second));
  }

  bool match_link_side(const ChemLink::Side& side,
                       const std::string& resname) const {
    if (!side.comp.empty())
      return side.comp == resname;
    if (side.group == ChemLink::Group::Null)
      return false;
    auto iter = res_group.find(resname);
    return iter != res_group.end() && side.matches_group(iter->second);
  }

  std::vector<Match> find_possible_links(Structure& st,
                                         double bond_margin,
                                         double radius_margin,
                                         bool skip_intra_residue_links=true) {
    std::vector<Match> results;
    Model& model = st.models.at(0);
    double search_radius = std::max(global_max_dist * bond_margin,
                                    /*max r1+r2 ~=*/3.0 * radius_margin);
    SubCells sc(model, st.cell, std::max(5.0, search_radius));
    sc.populate();
    for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
      Chain& chain = model.chains[n_ch];
      for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
        Residue& res = chain.residues[n_res];
        for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
          Atom& atom = res.atoms[n_atom];
          sc.for_each(atom.pos, atom.altloc, (float) search_radius,
                      [&](SubCells::Mark& m, float dist_sq) {
              // do not consider connections inside a residue
              if (skip_intra_residue_links && m.image_idx == 0 &&
                  m.chain_idx == n_ch && m.residue_idx == n_res)
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

              Match match;
              if (bond_margin > 0) {
                // search for a match in chem_links
                auto range = links.equal_range(Restraints::lexicographic_str(
                                                    atom.name, cra.atom->name));
                for (auto iter = range.first; iter != range.second; ++iter) {
                  const ChemLink& link = *iter->second;
                  const Restraints::Bond& bond = link.rt.bonds[0];
                  if (dist_sq > sq(bond.value * bond_margin))
                    continue;
                  bool order1;
                  if (bond.id1.atom == atom.name &&
                      match_link_side(link.side1, res.name) &&
                      match_link_side(link.side2, cra.residue->name))
                    order1 = true;
                  else if (bond.id2.atom == atom.name &&
                      match_link_side(link.side2, res.name) &&
                      match_link_side(link.side1, cra.residue->name))
                    order1 = false;
                  else
                    continue;
                  int link_score = link.side1.specificity() +
                                   link.side2.specificity();
                  // check chirality
                  Residue& res1 = order1 ? res : *cra.residue;
                  Residue* res2 = order1 ? cra.residue : &res;
                  char alt = atom.altloc ? atom.altloc : cra.atom->altloc;
                  for (const Restraints::Chirality& chirality : link.rt.chirs)
                    if (chirality.sign != ChiralityType::Both) {
                      Atom* at1 = chirality.id_ctr.get_from(res1, res2, alt);
                      Atom* at2 = chirality.id1.get_from(res1, res2, alt);
                      Atom* at3 = chirality.id2.get_from(res1, res2, alt);
                      Atom* at4 = chirality.id3.get_from(res1, res2, alt);
                      if (at1 && at2 && at3 && at4) {
                        double vol = calculate_chiral_volume(at1->pos, at2->pos,
                                                             at3->pos, at4->pos);
                        if (chirality.is_wrong(vol))
                          link_score -= 10;
                      }
                    }
                  // check fixed torsion angle (_chem_link_tor.period == 0)
                  for (const Restraints::Torsion& tor : link.rt.torsions)
                    if (tor.period == 0) {
                      Atom* at1 = tor.id1.get_from(res1, res2, alt);
                      Atom* at2 = tor.id2.get_from(res1, res2, alt);
                      Atom* at3 = tor.id3.get_from(res1, res2, alt);
                      Atom* at4 = tor.id4.get_from(res1, res2, alt);
                      double z = 10.;
                      if (at1 && at2 && at3 && at4)
                        z = angle_z(calculate_dihedral(at1->pos, at2->pos,
                                                       at3->pos, at4->pos),
                                    tor);
                      link_score -= (int) z;
                    }
                  match.chem_link_count++;
                  if (link_score < match.score)
                    continue;
                  match.chem_link = &link;
                  match.score = link_score;
                  if (order1) {
                    match.cra1 = {&chain, &res, &atom};
                    match.cra2 = cra;
                  } else {
                    match.cra1 = cra;
                    match.cra2 = {&chain, &res, &atom};
                  }
                }
              }

              // potential other links according to covalent radii
              if (!match.chem_link) {
                float r1 = atom.element.covalent_r();
                float r2 = cra.atom->element.covalent_r();
                if (dist_sq > sq((r1 + r2) * radius_margin))
                  return;
                match.cra1 = {&chain, &res, &atom};
                match.cra2 = cra;
              }

              match.same_image = !m.image_idx;
              match.bond_length = std::sqrt(dist_sq);
              results.push_back(match);
          });
        }
      }
    }
    for (Match& match : results)
      match.conn = st.find_connection_by_cra(match.cra1, match.cra2);
    return results;
  }
};

} // namespace gemmi
#endif
