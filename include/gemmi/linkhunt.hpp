// Copyright 2019 Global Phasing Ltd.
//
// Searching for links based on the _chem_link table from monomer dictionary.

#ifndef GEMMI_LINKHUNT_HPP_
#define GEMMI_LINKHUNT_HPP_

#include <map>
#include <unordered_map>
#include "elem.hpp"
#include "model.hpp"
#include "monlib.hpp"
#include "neighbor.hpp"
#include "contact.hpp"

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
  const MonLib* monlib_ptr = nullptr;
  std::multimap<std::string, const ChemLink*> links;

  void index_chem_links(const MonLib& monlib, bool use_alias=true) {
    std::map<ChemComp::Group, std::map<std::string, std::vector<std::string>>> aliases;
    if (use_alias)
      for (const auto& iter : monlib.monomers)
        for (const ChemComp::Aliasing& a : iter.second.aliases)
          for (const std::pair<std::string, std::string>& r : a.related) {
            const ChemComp::Group& gr = ChemComp::is_nucleotide_group(a.group) ? ChemComp::Group::DnaRna : a.group;
            aliases[gr][r.second].push_back(r.first);
          }

    for (const auto& iter : monlib.links) {
      const ChemLink& link = iter.second;
      if (link.rt.bonds.empty())
        continue;
      if (link.rt.bonds.size() > 1)
        fprintf(stderr, "Note: considering only the first bond in %s\n",
                link.id.c_str());
      if (link.side1.comp.empty() && link.side2.comp.empty())
        if (link.side1.group == ChemComp::Group::Null ||
            link.side2.group == ChemComp::Group::Null ||
            link.id == "SS")
          continue;
      const Restraints::Bond& bond = link.rt.bonds[0];
      if (bond.value > global_max_dist)
        global_max_dist = bond.value;
      links.emplace(bond.lexicographic_str(), &link);

      if (!use_alias || (!link.side1.comp.empty() && !link.side2.comp.empty()))
        continue;
      std::vector<std::string> *names1 = nullptr, *names2 = nullptr;
      if (link.side1.comp.empty()) {
        auto i = aliases.find(link.side1.group);
        if (i != aliases.end()) {
          auto j = i->second.find(bond.id1.atom);
          if (j != i->second.end())
            names1 = &j->second;
        }
      }
      if (link.side2.comp.empty()) {
        auto i = aliases.find(link.side2.group);
        if (i != aliases.end()) {
          auto j = i->second.find(bond.id2.atom);
          if (j != i->second.end())
            names2 = &j->second;
        }
      }
     if (names1 && names2)
       for (const std::string& n1 : *names1)
         for (const std::string& n2 : *names2)
           links.emplace(Restraints::lexicographic_str(n1, n2), &link);
     else if (names1 || names2) {
       const std::string& n1 = names1 ? bond.id2.atom : bond.id1.atom;
       for (const std::string& n2 : (names1 ? *names1 : *names2))
         links.emplace(Restraints::lexicographic_str(n1, n2), &link);
     }
    }
    monlib_ptr = &monlib;
  }

  std::vector<Match> find_possible_links(Structure& st,
                                         double bond_margin,
                                         double radius_margin,
                                         ContactSearch::Ignore ignore) {
    std::vector<Match> results;
    Model& model = st.first_model();
    double search_radius = std::max(global_max_dist * bond_margin,
                                    /*max r1+r2 ~=*/3.0 * radius_margin);
    NeighborSearch ns(model, st.cell, std::max(5.0, search_radius));
    ns.populate();

    ContactSearch contacts((float) search_radius);
    contacts.ignore = ignore;
    contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                      int image_idx, float dist_sq) {
        Match match;

        // search for a match in chem_links
        if (bond_margin > 0) {
          auto range = links.equal_range(Restraints::lexicographic_str(
                                            cra1.atom->name, cra2.atom->name));
          // similar to MonLib::match_link()
          for (auto iter = range.first; iter != range.second; ++iter) {
            const ChemLink& link = *iter->second;
            const Restraints::Bond& bond = link.rt.bonds[0];
            if (dist_sq > sq(bond.value * bond_margin))
              continue;
            const ChemComp::Aliasing* aliasing1 = nullptr;
            const ChemComp::Aliasing* aliasing2 = nullptr;
            bool order1;
            if (monlib_ptr->link_side_matches_residue(link.side1, cra1.residue->name, &aliasing1) &&
                monlib_ptr->link_side_matches_residue(link.side2, cra2.residue->name, &aliasing2) &&
                atom_match_with_alias(bond.id1.atom, cra1.atom->name, aliasing1))
              order1 = true;
            else if (monlib_ptr->link_side_matches_residue(link.side2, cra1.residue->name, &aliasing1) &&
                     monlib_ptr->link_side_matches_residue(link.side1, cra2.residue->name, &aliasing2) &&
                     atom_match_with_alias(bond.id2.atom, cra1.atom->name, aliasing1))
              order1 = false;
            else
              continue;
            int link_score = link.calculate_score(
                    order1 ? *cra1.residue : *cra2.residue,
                    order1 ? cra2.residue : cra1.residue,
                    cra1.atom->altloc_or(cra2.atom->altloc),
                    order1 ? aliasing1 : aliasing2,
                    order1 ? aliasing2 : aliasing1);
            match.chem_link_count++;
            if (link_score > match.score) {
              match.chem_link = &link;
              match.score = link_score;
              if (order1) {
                match.cra1 = cra1;
                match.cra2 = cra2;
              } else {
                match.cra1 = cra2;
                match.cra2 = cra1;
              }
            }
          }
        }

        // potential other links according to covalent radii
        if (!match.chem_link) {
          float r1 = cra1.atom->element.covalent_r();
          float r2 = cra2.atom->element.covalent_r();
          if (dist_sq > sq((r1 + r2) * radius_margin))
            return;
          match.cra1 = cra1;
          match.cra2 = cra2;
        }

        // finalize
        match.same_image = !image_idx;
        match.bond_length = std::sqrt(dist_sq);
        results.push_back(match);
    });

    // add references to st.connections
    for (Match& match : results)
      match.conn = st.find_connection_by_cra(match.cra1, match.cra2);

    return results;
  }
};

} // namespace gemmi
#endif
