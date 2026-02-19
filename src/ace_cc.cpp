// Copyright 2025 Global Phasing Ltd.
//
// Restraint generation helpers for prepare_chemcomp():
// chemical adjustments, torsions, chirality, planes, H naming.

#include "gemmi/ace_cc.hpp"
#include "gemmi/acedrg_tables.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/ace_graph.hpp"
#include "gemmi/resinfo.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <set>
#include <climits>

namespace gemmi {

namespace {

template <typename Range>
int count_missing_values(const Range& range) {
  int missing = 0;
  for (const auto& item : range)
    if (std::isnan(item.value)) missing++;
  return missing;
}

void remove_atom_by_id(ChemComp& cc, const std::string& atom_id) {
  auto is_id = [&](const Restraints::AtomId& id) { return id.atom == atom_id; };
  vector_remove_if(cc.rt.bonds, [&](const Restraints::Bond& b) {
    return is_id(b.id1) || is_id(b.id2);
  });
  vector_remove_if(cc.rt.angles, [&](const Restraints::Angle& a) {
    return is_id(a.id1) || is_id(a.id2) || is_id(a.id3);
  });
  vector_remove_if(cc.rt.torsions, [&](const Restraints::Torsion& t) {
    return is_id(t.id1) || is_id(t.id2) || is_id(t.id3) || is_id(t.id4);
  });
  vector_remove_if(cc.rt.chirs, [&](const Restraints::Chirality& c) {
    return is_id(c.id_ctr) || is_id(c.id1) || is_id(c.id2) || is_id(c.id3);
  });
  for (auto it = cc.rt.planes.begin(); it != cc.rt.planes.end(); ) {
    auto& ids = it->ids;
    vector_remove_if(ids, [&](const Restraints::AtomId& id) { return is_id(id); });
    if (ids.empty())
      it = cc.rt.planes.erase(it);
    else
      ++it;
  }
  vector_remove_if(cc.atoms, [&](const ChemComp::Atom& a) { return a.id == atom_id; });
}

BondType get_bond_type(const ChemComp& cc, const std::string& a, const std::string& b) {
  for (const auto& bond : cc.rt.bonds)
    if ((bond.id1.atom == a && bond.id2.atom == b) ||
        (bond.id2.atom == a && bond.id1.atom == b))
      return bond.type;
  return BondType::Unspec;
}

bool atom_has_unsaturated_bond(const ChemComp& cc, const std::string& atom_id) {
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == atom_id || bond.id2.atom == atom_id) {
      if (bond.type == BondType::Double ||
          bond.type == BondType::Aromatic ||
          bond.type == BondType::Deloc)
        return true;
    }
  }
  return false;
}

bool has_carboxylic_acid_neighbor(const ChemComp& cc,
    const std::map<std::string, std::vector<std::string>>& neighbors,
    const std::string& alpha_carbon_id) {
  auto it = neighbors.find(alpha_carbon_id);
  if (it == neighbors.end())
    return false;
  for (const std::string& c_neighbor : it->second) {
    int cn_idx = cc.find_atom_index(c_neighbor);
    if (cn_idx < 0 || cc.atoms[cn_idx].el != El::C)
      continue;
    int o_count = 0;
    bool has_double_o = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom != c_neighbor && bond.id2.atom != c_neighbor)
        continue;
      std::string other = (bond.id1.atom == c_neighbor) ? bond.id2.atom : bond.id1.atom;
      int other_idx = cc.find_atom_index(other);
      if (other_idx >= 0 && cc.atoms[other_idx].el == El::O) {
        ++o_count;
        if (bond.type == BondType::Double || bond.type == BondType::Deloc)
          has_double_o = true;
      }
    }
    if (o_count >= 2 && has_double_o)
      return true;
  }
  return false;
}

bool is_halogen(Element el) {
  return el == El::F || el == El::Cl || el == El::Br || el == El::I || el == El::At;
}

void adjust_terminal_carboxylate(ChemComp& cc) {
  if (cc.find_atom("N") == cc.atoms.end())
    return;
  auto oxt_it = cc.find_atom("OXT");
  auto hxt_it = cc.find_atom("HXT");
  if (oxt_it == cc.atoms.end() || hxt_it == cc.atoms.end())
    return;
  if (oxt_it->el != El::O || hxt_it->el != El::H)
    return;

  std::string oxt_neighbor;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == "OXT" && bond.id2.atom != "HXT")
      oxt_neighbor = bond.id2.atom;
    else if (bond.id2.atom == "OXT" && bond.id1.atom != "HXT")
      oxt_neighbor = bond.id1.atom;
  }
  if (oxt_neighbor.empty())
    return;

  auto neighbor_it = cc.find_atom(oxt_neighbor);
  if (neighbor_it == cc.atoms.end() || neighbor_it->el != El::C)
    return;

  int o_count = 0;
  bool has_double_o = false;
  bool has_n_neighbor = false;
  for (const auto& bond : cc.rt.bonds) {
    std::string other;
    if (bond.id1.atom == oxt_neighbor)
      other = bond.id2.atom;
    else if (bond.id2.atom == oxt_neighbor)
      other = bond.id1.atom;
    else
      continue;
    auto other_it = cc.find_atom(other);
    if (other_it != cc.atoms.end()) {
      if (other_it->el == El::O) {
        ++o_count;
        if (bond.type == BondType::Double || bond.type == BondType::Deloc)
          has_double_o = true;
      } else if (other_it->el == El::N) {
        has_n_neighbor = true;
      }
    }
  }
  if (o_count < 2 || !has_double_o)
    return;
  if (has_n_neighbor)
    return;

  oxt_it->charge = -1.0f;
  remove_atom_by_id(cc, "HXT");
}

bool add_n_terminal_h3(ChemComp& cc) {
  auto n_it = cc.find_atom("N");
  if (n_it == cc.atoms.end())
    return false;
  if (std::fabs(n_it->charge) > 0.5f)
    return false;
  if (cc.find_atom("H3") != cc.atoms.end())
    return false;
  if (cc.find_atom("H") == cc.atoms.end() || cc.find_atom("H2") == cc.atoms.end())
    return false;

  bool n_has_h = false;
  bool n_has_h2 = false;
  bool n_has_h3 = false;
  std::string ca_atom;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == "N" || bond.id2.atom == "N") {
      const std::string& other = (bond.id1.atom == "N") ? bond.id2.atom : bond.id1.atom;
      if (other == "H")
        n_has_h = true;
      else if (other == "H2")
        n_has_h2 = true;
      else if (other == "H3")
        n_has_h3 = true;
      else if (other[0] != 'H')
        ca_atom = other;
    }
  }
  if (!n_has_h || !n_has_h2 || n_has_h3)
    return false;

  if (!ca_atom.empty()) {
    if (atom_has_unsaturated_bond(cc, ca_atom))
      return false;
    auto neighbors = make_neighbor_names(cc);
    if (!has_carboxylic_acid_neighbor(cc, neighbors, ca_atom))
      return false;
  }

  cc.atoms.push_back(ChemComp::Atom{"H3", "", El::H, 0.0f, "H", "", Position()});
  cc.rt.bonds.push_back({{1, "N"}, {1, "H3"}, BondType::Single, false,
                        NAN, NAN, NAN, NAN});

  auto add_angle_if_present = [&](const std::string& a1, const std::string& a2,
                                  const std::string& a3) {
    if (cc.find_atom(a1) == cc.atoms.end() ||
        cc.find_atom(a2) == cc.atoms.end() ||
        cc.find_atom(a3) == cc.atoms.end())
      return;
    cc.rt.angles.push_back({{1, a1}, {1, a2}, {1, a3}, NAN, NAN});
  };

  add_angle_if_present("CA", "N", "H3");
  add_angle_if_present("H", "N", "H3");
  add_angle_if_present("H2", "N", "H3");
  return true;
}

Restraints::Angle* find_angle(ChemComp& cc, const std::string& center,
                              const std::string& a1, const std::string& a3) {
  for (auto& angle : cc.rt.angles) {
    if (angle.id2.atom != center)
      continue;
    if ((angle.id1.atom == a1 && angle.id3.atom == a3) ||
        (angle.id1.atom == a3 && angle.id3.atom == a1))
      return &angle;
  }
  return nullptr;
}

void sync_n_terminal_h3_angles(ChemComp& cc) {
  if (cc.find_atom("H3") == cc.atoms.end() || cc.find_atom("N") == cc.atoms.end())
    return;

  auto copy_angle = [&](const std::string& a1, const std::string& a3,
                        const std::string& src1, const std::string& src3) {
    Restraints::Angle* target = find_angle(cc, "N", a1, a3);
    if (!target)
      return;
    Restraints::Angle* source = find_angle(cc, "N", src1, src3);
    if (!source || std::isnan(source->value))
      return;
    target->value = source->value;
    target->esd = source->esd;
  };

  if (cc.find_atom("CA") != cc.atoms.end()) {
    copy_angle("CA", "H3", "CA", "H");
    copy_angle("CA", "H3", "CA", "H2");
  }
  copy_angle("H", "H3", "H", "H2");
  copy_angle("H2", "H3", "H", "H2");
}

// Deprotonate O-H groups bonded to a qualifying central atom (P or S).
// For phosphate: central P with 4 neighbors, >=3 O, no direct H.
// For sulfate: central S with 4 neighbors, all 4 O.
void adjust_oxoacid_group(ChemComp& cc, Element central_el,
                          int min_o_count, bool reject_h_on_center) {
  auto neighbors = make_neighbor_names(cc);

  std::set<std::string> qualifying_centers;
  for (const auto& atom : cc.atoms) {
    if (atom.el != central_el)
      continue;
    const auto& nb = neighbors[atom.id];
    if (nb.size() != 4)
      continue;
    int o_count = 0;
    bool has_h = false;
    for (const std::string& nid : nb) {
      int idx = cc.find_atom_index(nid);
      if (idx >= 0) {
        if (cc.atoms[idx].el == El::O)
          ++o_count;
        else if (cc.atoms[idx].el == El::H)
          has_h = true;
      }
    }
    if (o_count >= min_o_count && !(reject_h_on_center && has_h))
      qualifying_centers.insert(atom.id);
  }

  std::vector<std::string> h_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    bool bonded_to_center = false;
    std::string h_neighbor;
    for (const std::string& nid : nb) {
      if (qualifying_centers.count(nid))
        bonded_to_center = true;
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::H)
        h_neighbor = nid;
    }
    if (bonded_to_center && !h_neighbor.empty()) {
      atom.charge = -1.0f;
      h_to_remove.push_back(h_neighbor);
    }
  }

  for (const std::string& h_id : h_to_remove)
    remove_atom_by_id(cc, h_id);
}

void adjust_hexafluorophosphate(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);

  for (auto& atom : cc.atoms) {
    if (atom.el != El::P)
      continue;
    int f_count = 0;
    bool has_h = false;
    for (const std::string& nb : neighbors[atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0)
        continue;
      Element el = cc.atoms[idx].el;
      if (el == El::F)
        ++f_count;
      else if (el == El::H)
        has_h = true;
    }
    if (f_count != 6 || has_h)
      continue;

    std::string new_h = "H";
    if (cc.find_atom(new_h) != cc.atoms.end()) {
      for (int i = 1; i < 100; ++i) {
        new_h = cat("H", i);
        if (cc.find_atom(new_h) == cc.atoms.end())
          break;
      }
    }
    if (cc.find_atom(new_h) != cc.atoms.end())
      continue;

    cc.atoms.push_back({new_h, "", El::H, 0.0f, "H", "", Position()});
    cc.rt.bonds.push_back({{1, atom.id}, {1, new_h}, BondType::Single, false,
                          NAN, NAN, NAN, NAN});
    for (const std::string& nb : neighbors[atom.id]) {
      cc.rt.angles.push_back({{1, nb}, {1, atom.id}, {1, new_h}, NAN, NAN});
    }
  }
}

void adjust_carboxy_asp(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);

  std::set<std::string> matched_atoms;
  for (const auto& atom : cc.atoms) {
    if (atom.el != El::C)
      continue;
    std::string carbonyl_o;
    std::vector<std::string> single_o;
    std::vector<std::string> alpha_c;
    for (const std::string& nid : neighbors[atom.id]) {
      int idx = cc.find_atom_index(nid);
      if (idx < 0)
        continue;
      Element el = cc.atoms[idx].el;
      if (el == El::O) {
        if (get_bond_type(cc, atom.id, nid) == BondType::Double)
          carbonyl_o = nid;
        else
          single_o.push_back(nid);
      } else if (el == El::C) {
        alpha_c.push_back(nid);
      }
    }
    if (carbonyl_o.empty() || single_o.empty() || alpha_c.empty())
      continue;
    for (const std::string& ac : alpha_c) {
      bool has_substituent = false;
      for (const std::string& ac_nb : neighbors[ac]) {
        if (ac_nb == atom.id)
          continue;
        int idx = cc.find_atom_index(ac_nb);
        if (idx >= 0 && cc.atoms[idx].el != El::H) {
          has_substituent = true;
          break;
        }
      }
      if (!has_substituent)
        continue;
      bool is_conjugated = false;
      for (const auto& bond : cc.rt.bonds) {
        std::string other;
        if (bond.id1.atom == ac)
          other = bond.id2.atom;
        else if (bond.id2.atom == ac)
          other = bond.id1.atom;
        else
          continue;
        if (bond.type == BondType::Double && !bond.aromatic) {
          int idx = cc.find_atom_index(other);
          if (idx >= 0 && cc.atoms[idx].el == El::C) {
            if (!atoms_in_same_ring_by_alt_path(ac, other, neighbors)) {
              is_conjugated = true;
              break;
            }
          }
        }
      }
      if (!is_conjugated) {
        for (const auto& so : single_o)
          matched_atoms.insert(so);
      }
      for (const std::string& ac_nb : neighbors[ac]) {
        if (ac_nb != atom.id)
          matched_atoms.insert(ac_nb);
      }
    }
  }

  std::vector<std::string> h_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O || matched_atoms.find(atom.id) == matched_atoms.end())
      continue;
    std::string h_neighbor;
    int h_count = 0;
    for (const std::string& nid : neighbors[atom.id]) {
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::H) {
        h_neighbor = nid;
        ++h_count;
      }
    }
    if (h_count == 1) {
      atom.charge = -1.0f;
      h_to_remove.push_back(h_neighbor);
    }
  }
  for (const std::string& h_id : h_to_remove)
    remove_atom_by_id(cc, h_id);
}

std::string acedrg_h_name(std::set<std::string>& used_names, const std::string& n_id,
                          const std::vector<std::string>& h_on_n) {
  std::string root;
  for (char c : n_id)
    if (std::isalpha(static_cast<unsigned char>(c)))
      root += c;
  std::string h_root = (root.size() == 2) ? "H" + root.substr(1) : "H";

  if (h_root.size() >= 2 && used_names.find(h_root) == used_names.end()) {
    used_names.insert(h_root);
    return h_root;
  }

  int idx_max = 0;
  for (const std::string& h : h_on_n) {
    size_t d = 0;
    for (; d < h.size(); ++d)
      if (std::isdigit(static_cast<unsigned char>(h[d])))
        break;
    if (d < h.size()) {
      const std::string suffix = h.substr(d);
      if (std::all_of(suffix.begin(), suffix.end(),
                      [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
        int val = std::stoi(suffix);
        if (val > idx_max)
          idx_max = val;
      }
    }
  }

  int start = (idx_max > 0) ? idx_max + 1 : 2;
  std::string cand = cat(h_root, start);
  if (used_names.find(cand) == used_names.end()) {
    used_names.insert(cand);
    return cand;
  }

  for (int n = 2; n < 10000; ++n) {
    cand = cat(h_root, n);
    if (used_names.find(cand) == used_names.end()) {
      used_names.insert(cand);
      return cand;
    }
  }
  return {};
}

void adjust_guanidinium_group(ChemComp& cc, std::set<std::string>& used_names) {

  auto neighbors = make_neighbor_names(cc);

  std::map<std::string, bool> has_unsat_bond;
  for (const auto& bond : cc.rt.bonds) {
    bool unsat = bond.aromatic ||
                 bond.type == BondType::Double ||
                 bond.type == BondType::Triple ||
                 bond.type == BondType::Aromatic ||
                 bond.type == BondType::Deloc;
    if (unsat) {
      has_unsat_bond[bond.id1.atom] = true;
      has_unsat_bond[bond.id2.atom] = true;
    }
  }

  auto is_carbonyl_like = [&](const std::string& c_id) {
    int c_idx = cc.find_atom_index(c_id);
    if (c_idx < 0 || cc.atoms[c_idx].el != El::C)
      return false;
    bool has_double_hetero = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom != c_id && bond.id2.atom != c_id)
        continue;
      if (bond.aromatic || bond.type == BondType::Aromatic ||
          bond.type == BondType::Deloc || bond.type == BondType::Triple)
        return false;
      if (bond.type == BondType::Double) {
        const std::string& other = (bond.id1.atom == c_id) ? bond.id2.atom
                                                           : bond.id1.atom;
        int idx = cc.find_atom_index(other);
        if (idx < 0)
          continue;
        Element el = cc.atoms[idx].el;
        if (el == El::O || el == El::S || el == El::Se)
          has_double_hetero = true;
        else
          return false;
      }
    }
    return has_double_hetero;
  };

  for (auto& atom : cc.atoms) {
    if (atom.el != El::C)
      continue;
    const auto& nb = neighbors[atom.id];
    std::vector<std::string> n_neighbors;
    for (const std::string& nid : nb) {
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::N)
        n_neighbors.push_back(nid);
    }
    if (n_neighbors.size() != 3)
      continue;

    for (const std::string& n_id : n_neighbors) {
      bool has_double = false;
      for (const auto& bond : cc.rt.bonds) {
        if ((bond.id1.atom == atom.id && bond.id2.atom == n_id) ||
            (bond.id2.atom == atom.id && bond.id1.atom == n_id)) {
          if (bond.type == BondType::Double)
            has_double = true;
          break;
        }
      }
      if (!has_double)
        continue;

      bool other_n_has_unsat_sub = false;
      for (const std::string& other_n_id : n_neighbors) {
        if (other_n_id == n_id)
          continue;
        const auto& other_nb = neighbors[other_n_id];
        for (const std::string& onb : other_nb) {
          if (onb == atom.id)
            continue;
          int onb_idx = cc.find_atom_index(onb);
          if (onb_idx < 0)
            continue;
          if (cc.atoms[onb_idx].el == El::H)
            continue;
          if (cc.atoms[onb_idx].el != El::C &&
              cc.atoms[onb_idx].el != El::Si &&
              cc.atoms[onb_idx].el != El::Ge) {
            other_n_has_unsat_sub = true;
            break;
          }
          if (has_unsat_bond[onb] && !is_carbonyl_like(onb)) {
            other_n_has_unsat_sub = true;
            break;
          }
        }
        if (other_n_has_unsat_sub)
          break;
      }
      if (other_n_has_unsat_sub)
        continue;

      const auto& n_nb = neighbors[n_id];
      std::vector<std::string> h_neighbors;
      bool has_substituent = false;
      for (const std::string& nid : n_nb) {
        int nid_idx = cc.find_atom_index(nid);
        if (nid_idx < 0)
          continue;
        if (cc.atoms[nid_idx].el == El::H)
          h_neighbors.push_back(nid);
        else if (nid != atom.id)
          has_substituent = true;
      }

      if (has_substituent)
        continue;

      if (h_neighbors.size() == 1) {
        int n_idx = cc.find_atom_index(n_id);
        if (n_idx < 0)
          continue;

        cc.atoms[n_idx].charge = 1.0f;

        std::string new_h_id = acedrg_h_name(used_names, n_id, h_neighbors);
        if (new_h_id.empty())
          continue;

        cc.atoms.push_back(ChemComp::Atom{new_h_id, "", El::H, 0.0f, "H", "", Position()});
        cc.rt.bonds.push_back({{1, n_id}, {1, new_h_id}, BondType::Single, false,
                              NAN, NAN, NAN, NAN});
        // Add angles for the new hydrogen with all existing neighbors
        cc.rt.angles.push_back({{1, atom.id}, {1, n_id}, {1, new_h_id}, NAN, NAN});
        cc.rt.angles.push_back({{1, h_neighbors[0]}, {1, n_id}, {1, new_h_id}, NAN, NAN});
        neighbors[n_id].push_back(new_h_id);
        neighbors[new_h_id].push_back(n_id);
      }
    }
  }
}

void adjust_amino_ter_amine(ChemComp& cc, std::set<std::string>& used_names) {
  auto neighbors = make_neighbor_names(cc);

  for (auto& n1 : cc.atoms) {
    if (n1.el != El::N)
      continue;
    if (std::fabs(n1.charge) > 0.5f)
      continue;

    std::vector<std::string> h_ids;
    std::string c1_id;
    for (const std::string& nb : neighbors[n1.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0)
        continue;
      if (cc.atoms[idx].el == El::H)
        h_ids.push_back(nb);
      else if (cc.atoms[idx].el == El::C)
        c1_id = nb;
    }
    if (h_ids.size() != 2 || c1_id.empty())
      continue;

    bool c1_has_double = false;
    for (const auto& bond : cc.rt.bonds) {
      if ((bond.id1.atom == c1_id || bond.id2.atom == c1_id) &&
          bond.type == BondType::Double) {
        c1_has_double = true;
        break;
      }
    }
    if (c1_has_double)
      continue;

    for (const std::string& c2_id : neighbors[c1_id]) {
      if (c2_id == n1.id)
        continue;
      int c2_idx = cc.find_atom_index(c2_id);
      if (c2_idx < 0 || cc.atoms[c2_idx].el != El::C)
        continue;

      bool has_carbonyl = false, has_amide_n = false;
      for (const std::string& c2_nb : neighbors[c2_id]) {
        int nb_idx = cc.find_atom_index(c2_nb);
        if (nb_idx < 0)
          continue;
        Element el = cc.atoms[nb_idx].el;
        if (el == El::O && get_bond_type(cc, c2_id, c2_nb) == BondType::Double)
          has_carbonyl = true;
        if (el == El::N && c2_nb != n1.id) {
          for (const std::string& amide_n_nb : neighbors[c2_nb]) {
            if (amide_n_nb == c2_id)
              continue;
            int an_idx = cc.find_atom_index(amide_n_nb);
            if (an_idx >= 0 && cc.atoms[an_idx].el != El::H) {
              has_amide_n = true;
              break;
            }
          }
        }
      }

      if (has_carbonyl && has_amide_n) {
        n1.charge = 1.0f;
        std::string new_h = acedrg_h_name(used_names, n1.id, h_ids);
        if (new_h.empty())
          break;
        cc.atoms.push_back({new_h, "", El::H, 0.0f, "H", "", Position()});
        cc.rt.bonds.push_back({{1, n1.id}, {1, new_h}, BondType::Single, false,
                              NAN, NAN, NAN, NAN});
        // Add angles for the new hydrogen with all existing neighbors
        for (const std::string& h_id : h_ids)
          cc.rt.angles.push_back({{1, h_id}, {1, n1.id}, {1, new_h}, NAN, NAN});
        cc.rt.angles.push_back({{1, c1_id}, {1, n1.id}, {1, new_h}, NAN, NAN});
        break;
      }
    }
  }
}

void adjust_terminal_amine(ChemComp& cc, std::set<std::string>& used_names) {
  auto neighbors = make_neighbor_names(cc);

  for (auto& n_atom : cc.atoms) {
    if (n_atom.el != El::N)
      continue;
    if (std::fabs(n_atom.charge) > 0.5f)
      continue;

    int n_carbon = 0, n_hydrogen = 0, n_other = 0;
    std::string alpha_carbon_id;
    std::vector<std::string> h_ids;
    for (const std::string& nb : neighbors[n_atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0)
        continue;
      Element el = cc.atoms[idx].el;
      if (el == El::C) {
        n_carbon++;
        alpha_carbon_id = nb;
      } else if (el == El::H) {
        n_hydrogen++;
        h_ids.push_back(nb);
      } else {
        n_other++;
      }
    }

    if (n_carbon != 1 || n_other != 0)
      continue;
    if (n_hydrogen != 2)
      continue;

    if (atom_has_unsaturated_bond(cc, alpha_carbon_id))
      continue;
    if (atom_has_unsaturated_bond(cc, n_atom.id))
      continue;
    if (!has_carboxylic_acid_neighbor(cc, neighbors, alpha_carbon_id))
      continue;

    n_atom.charge = 1.0f;

    std::string new_h = acedrg_h_name(used_names, n_atom.id, h_ids);
    if (new_h.empty())
      continue;

    cc.atoms.push_back({new_h, "", El::H, 0.0f, "H", "", Position()});
    cc.rt.bonds.push_back({{1, n_atom.id}, {1, new_h}, BondType::Single, false,
                          NAN, NAN, NAN, NAN});

    for (const std::string& h_id : h_ids) {
      cc.rt.angles.push_back({{1, h_id}, {1, n_atom.id}, {1, new_h}, NAN, NAN});
    }
    cc.rt.angles.push_back({{1, alpha_carbon_id}, {1, n_atom.id}, {1, new_h}, NAN, NAN});
  }
}

void add_angles_from_bonds_if_missing(ChemComp& cc) {
  if (!cc.rt.angles.empty())
    return;

  auto atom_index = cc.make_atom_index();

  std::vector<std::vector<size_t>> neighbors(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    neighbors[it1->second].push_back(it2->second);
    neighbors[it2->second].push_back(it1->second);
  }

  std::set<std::tuple<std::string, std::string, std::string>> seen;
  for (size_t center = 0; center < neighbors.size(); ++center) {
    if (cc.atoms[center].el.is_metal())
      continue;
    auto& nbs = neighbors[center];
    if (nbs.size() < 2)
      continue;
    for (size_t i = 0; i + 1 < nbs.size(); ++i) {
      for (size_t j = i + 1; j < nbs.size(); ++j) {
        // AceDRG doesn't generate angles where both flanks are metals
        if (cc.atoms[nbs[i]].el.is_metal() && cc.atoms[nbs[j]].el.is_metal())
          continue;
        const std::string& a1 = cc.atoms[nbs[i]].id;
        const std::string& a3 = cc.atoms[nbs[j]].id;
        std::string first = a1;
        std::string third = a3;
        if (third < first)
          std::swap(first, third);
        auto key = std::make_tuple(cc.atoms[center].id, first, third);
        if (!seen.insert(key).second)
          continue;
        cc.rt.angles.push_back({{1, a1}, {1, cc.atoms[center].id},
                                {1, a3}, NAN, NAN});
      }
    }
  }
}

int chirality_priority(Element el) {
  if (el == El::O) return 0;
  if (el == El::N) return 1;
  if (el == El::S) return 2;
  if (el == El::Br || el == El::I || el == El::At) return 3;
  if (el == El::P) return 4;
  if (el == El::F || el == El::Cl) return 5;
  if (el == El::C) return 6;
  if (el == El::H) return 7;
  return 8;
}

int rdkit_twice_bond_type(BondType bt) {
  switch (bt) {
    case BondType::Single:
      return 2;
    case BondType::Double:
      return 4;
    case BondType::Triple:
      return 6;
    case BondType::Aromatic:
    case BondType::Deloc:
      return 3;
    default:
      return 0;
  }
}

int rdkit_cip_bond_count(const ChemComp& cc, const AceBondAdjacency& adj,
                         const AceBondNeighbor& nb) {
  if (nb.type == BondType::Double && cc.atoms[nb.idx].el == El::P) {
    size_t deg = adj[nb.idx].size();
    if (deg == 3 || deg == 4)
      return 1;
  }
  return rdkit_twice_bond_type(nb.type);
}

struct CipSortableRef {
  std::vector<int>* entry = nullptr;
  size_t atom_idx = 0;
  unsigned curr_rank = 0;
};

bool cip_ref_less(const CipSortableRef& a, const CipSortableRef& b) {
  return *a.entry < *b.entry;
}

void find_cip_tied_segments(std::vector<CipSortableRef>& sorted,
                            std::vector<std::pair<size_t, size_t>>& segments,
                            unsigned& num_ranks) {
  segments.clear();
  num_ranks = static_cast<unsigned>(sorted.size());
  if (sorted.empty())
    return;
  CipSortableRef* current = &sorted[0];
  unsigned running_rank = 0;
  current->curr_rank = running_rank;
  bool in_equal = false;
  for (size_t i = 1; i < sorted.size(); ++i) {
    CipSortableRef& next = sorted[i];
    if (*current->entry == *next.entry) {
      next.curr_rank = running_rank;
      --num_ranks;
      if (!in_equal) {
        in_equal = true;
        segments.push_back(std::make_pair(i - 1, size_t(0)));
      }
    } else {
      ++running_rank;
      next.curr_rank = running_rank;
      current = &next;
      if (in_equal) {
        segments.back().second = i;
        in_equal = false;
      }
    }
  }
  if (in_equal)
    segments.back().second = sorted.size() - 1;
}

std::vector<unsigned> compute_rdkit_legacy_cip_ranks(
    const ChemComp& cc, const AceBondAdjacency& adj) {
  size_t n = cc.atoms.size();
  std::vector<unsigned> ranks(n, 0);
  if (n == 0)
    return ranks;

  std::vector<std::vector<int>> cip_entries(n);
  std::vector<CipSortableRef> sorted;
  sorted.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    cip_entries[i].reserve(16);
    int atomic_num = cc.atoms[i].el.atomic_number() % 128;
    int mass = 512;
    unsigned long invariant = static_cast<unsigned long>(atomic_num);
    invariant = (invariant << 10) | static_cast<unsigned long>(mass);
    invariant = (invariant << 10);  // atom map number is absent in CCD input
    cip_entries[i].push_back(static_cast<int>(invariant));
    sorted.push_back({&cip_entries[i], i, 0});
  }

  std::sort(sorted.begin(), sorted.end(), cip_ref_less);
  std::vector<std::pair<size_t, size_t>> needs_sorting;
  unsigned num_ranks = 0;
  find_cip_tied_segments(sorted, needs_sorting, num_ranks);
  for (size_t i = 0; i < sorted.size(); ++i)
    ranks[sorted[i].atom_idx] = sorted[i].curr_rank;

  for (size_t i = 0; i < n; ++i) {
    cip_entries[i][0] = cc.atoms[i].el.atomic_number();
    cip_entries[i].push_back(static_cast<int>(ranks[i]));
  }

  const size_t cip_rank_index = 2;  // seedWithInvars == false in RDKit
  std::vector<std::vector<std::pair<int, size_t>>> weighted_neighbors(n);
  for (size_t i = 0; i < n; ++i) {
    weighted_neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i]) {
      if (cc.atoms[i].el.is_metal() || cc.atoms[nb.idx].el.is_metal())
        continue;
      weighted_neighbors[i].push_back(
          std::make_pair(rdkit_cip_bond_count(cc, adj, nb), nb.idx));
    }
  }

  unsigned max_iterations = static_cast<unsigned>(n / 2 + 1);
  unsigned iter = 0;
  int last_num_ranks = -1;
  while (!needs_sorting.empty() && iter < max_iterations &&
         (last_num_ranks < 0 ||
          static_cast<unsigned>(last_num_ranks) < num_ranks)) {
    for (size_t i = 0; i < n; ++i) {
      std::vector<std::pair<int, size_t>>& neighbors = weighted_neighbors[i];
      if (neighbors.size() > 1) {
        std::sort(neighbors.begin(), neighbors.end(),
                  [&](const std::pair<int, size_t>& a,
                      const std::pair<int, size_t>& b) {
                    return ranks[a.second] > ranks[b.second];
                  });
      }
      std::vector<int>& cip = cip_entries[i];
      for (const std::pair<int, size_t>& nb : neighbors)
        if (nb.first > 0)
          cip.insert(cip.end(), nb.first, static_cast<int>(ranks[nb.second] + 1));
      // AceDRG uses RDKit's getTotalNumHs(false). With explicit-H atom graphs
      // (the path used here), that contributes zero.
    }

    last_num_ranks = static_cast<int>(num_ranks);
    for (const std::pair<size_t, size_t>& seg : needs_sorting)
      std::sort(sorted.begin() + seg.first, sorted.begin() + seg.second + 1, cip_ref_less);
    find_cip_tied_segments(sorted, needs_sorting, num_ranks);
    for (size_t i = 0; i < sorted.size(); ++i)
      ranks[sorted[i].atom_idx] = sorted[i].curr_rank;

    if (static_cast<unsigned>(last_num_ranks) != num_ranks) {
      for (size_t i = 0; i < n; ++i) {
        cip_entries[i].resize(cip_rank_index + 1);
        cip_entries[i][cip_rank_index] = static_cast<int>(ranks[i]);
      }
    }
    ++iter;
  }
  return ranks;
}

void sort_neighbors_by_rdkit_cip_rank(std::vector<size_t>& neighbors,
                                      const std::vector<unsigned>& cip_ranks) {
  std::stable_sort(neighbors.begin(), neighbors.end(),
                   [&](size_t a, size_t b) { return cip_ranks[a] > cip_ranks[b]; });
}

bool is_oxygen_column(Element el) {
  return el == El::O || el == El::S || el == El::Se ||
         el == El::Te || el == El::Po;
}

// Find a pair of neighbors (one from each side of a bond) that share a ring.
// Returns {SIZE_MAX, SIZE_MAX} if no such pair exists.
std::pair<size_t, size_t> find_ring_sharing_pair(
    const AceBondAdjacency& adj, const std::vector<CodAtomInfo>& atom_info,
    size_t side1, size_t side2) {
  for (const auto& nb1 : adj[side1]) {
    if (nb1.idx == side2) continue;
    for (const auto& nb2 : adj[side2]) {
      if (nb2.idx == side1 || nb2.idx == nb1.idx) continue;
      if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings))
        return {nb1.idx, nb2.idx};
    }
  }
  return {SIZE_MAX, SIZE_MAX};
}

enum class RingParity { Even, Odd, NoFlip };
enum class RingFlip { Even, Odd };

std::map<std::pair<size_t, size_t>, RingParity>
build_ring_bond_parity(const AceBondAdjacency& adj,
                       const std::vector<CodAtomInfo>& atom_info) {
  std::map<std::pair<size_t, size_t>, RingParity> bond_ring_parity;
  int max_ring_id = 0;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      max_ring_id = std::max(max_ring_id, rid + 1);
  std::vector<std::vector<size_t>> rings(max_ring_id);
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      rings[rid].push_back(i);
  for (auto& ring : rings) {
    if (ring.empty())
      continue;
    std::set<size_t> ring_set(ring.begin(), ring.end());
    size_t start = *std::min_element(ring.begin(), ring.end());
    std::vector<size_t> traversal = {start};
    size_t prev = SIZE_MAX;
    size_t cur = start;
    while (traversal.size() < ring.size()) {
      bool found = false;
      for (const auto& nb : adj[cur]) {
        if (nb.idx != prev && ring_set.count(nb.idx) &&
            std::find(traversal.begin(), traversal.end(), nb.idx) == traversal.end()) {
          traversal.push_back(nb.idx);
          prev = cur;
          cur = nb.idx;
          found = true;
          break;
        }
      }
      if (!found)
        break;
    }
    for (size_t i = 0; i + 1 < traversal.size(); ++i) {
      auto key = std::minmax(traversal[i], traversal[i + 1]);
      bond_ring_parity[key] = (i % 2 == 0) ? RingParity::Even : RingParity::Odd;
    }
    if (traversal.size() == ring.size()) {
      auto ckey = std::minmax(traversal.back(), traversal.front());
      if (bond_ring_parity.find(ckey) == bond_ring_parity.end())
        bond_ring_parity[ckey] = RingParity::NoFlip;
    }
  }
  return bond_ring_parity;
}

std::map<std::pair<size_t, size_t>, RingFlip>
build_ring_bond_flip(const ChemComp& cc,
                     const AceBondAdjacency& adj,
                     const std::vector<CodAtomInfo>& atom_info) {
  struct RingWalk {
    std::string rep;
    std::vector<size_t> seq;
  };
  std::map<int, std::vector<size_t>> ring_atoms;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      ring_atoms[rid].push_back(i);

  std::vector<RingWalk> walks;
  walks.reserve(ring_atoms.size());
  for (const auto& kv : ring_atoms) {
    const std::vector<size_t>& r_atoms = kv.second;
    if (r_atoms.size() < 3)
      continue;
    std::set<size_t> rset(r_atoms.begin(), r_atoms.end());
    std::vector<size_t> linked;
    size_t start = *std::min_element(r_atoms.begin(), r_atoms.end());
    linked.push_back(start);
    size_t cur = start;
    int guard = 1;
    while (linked.size() < r_atoms.size() && guard < (int)r_atoms.size()) {
      bool advanced = false;
      for (const auto& nb : adj[cur]) {
        if (rset.count(nb.idx) == 0)
          continue;
        if (std::find(linked.begin(), linked.end(), nb.idx) != linked.end())
          continue;
        linked.push_back(nb.idx);
        cur = nb.idx;
        advanced = true;
        break;
      }
      if (!advanced)
        break;
      ++guard;
    }
    if (linked.size() != r_atoms.size())
      continue;
    std::vector<std::string> names;
    names.reserve(linked.size());
    for (size_t idx : linked)
      names.push_back(cc.atoms[idx].id);
    std::sort(names.begin(), names.end());
    std::string rep;
    for (const auto& n : names)
      rep += n;
    walks.push_back({std::move(rep), std::move(linked)});
  }
  std::sort(walks.begin(), walks.end(),
            [](const RingWalk& a, const RingWalk& b) { return a.rep < b.rep; });

  std::map<std::pair<size_t, size_t>, RingFlip> flips;
  for (const RingWalk& rw : walks) {
    for (size_t i = 0; i + 1 < rw.seq.size(); ++i) {
      auto key = std::minmax(rw.seq[i], rw.seq[i + 1]);
      if (flips.find(key) != flips.end())
        continue;
      flips[key] = (i % 2 == 0) ? RingFlip::Even : RingFlip::Odd;
    }
  }
  return flips;
}

struct SugarRingInfo {
  std::map<int, std::set<size_t>> ring_sets;
  std::set<std::pair<size_t, size_t>> ring_bonds;
  std::map<int, std::vector<size_t>> ring_seq;
};

SugarRingInfo detect_sugar_rings(const ChemComp& cc, const AceBondAdjacency& adj,
                                 const std::vector<CodAtomInfo>& atom_info) {
  SugarRingInfo out;

  auto co_token = [&](size_t idx) -> std::string {
    int nC = 0, nO = 0;
    for (const auto& nb : adj[idx]) {
      if (cc.atoms[nb.idx].el == El::C)
        ++nC;
      else if (cc.atoms[nb.idx].el == El::O)
        ++nO;
    }
    if (nC == 1 && nO == 2) return "CCO2";
    if (nC == 2 && nO == 1) return "CC2O";
    if (nC == 1 && nO == 1) return "CCO";
    if (nC == 2) return "CC2";
    return "";
  };
  auto ring_shape = [&](const std::set<size_t>& rset, size_t o_idx, size_t start,
                        std::vector<size_t>* out_order = nullptr) -> std::string {
    std::vector<size_t> order;
    order.reserve(rset.size());
    order.push_back(o_idx);
    order.push_back(start);
    size_t prev = o_idx;
    size_t curr = start;
    while (order.size() < rset.size()) {
      size_t next = SIZE_MAX;
      for (const auto& nb : adj[curr]) {
        if (nb.idx != prev && rset.count(nb.idx)) {
          next = nb.idx;
          break;
        }
      }
      if (next == SIZE_MAX)
        return "";
      order.push_back(next);
      prev = curr;
      curr = next;
    }
    bool closes = false;
    for (const auto& nb : adj[order.back()]) {
      if (nb.idx == o_idx) {
        closes = true;
        break;
      }
    }
    if (!closes)
      return "";
    if (out_order)
      *out_order = order;
    int o_c = 0;
    for (const auto& nb : adj[o_idx])
      if (cc.atoms[nb.idx].el == El::C)
        ++o_c;
    if (o_c != 2)
      return "";
    std::string s = "(OC2)";
    for (size_t i = 1; i < order.size(); ++i) {
      std::string t = co_token(order[i]);
      if (t.empty())
        return "";
      s += "(" + t + ")";
    }
    return s;
  };
  static const char* k_ace_sugar_shapes[] = {
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CC2O)(CC2O)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)(CCO)",
    "(OC2)(CCO2)(CC2O)(CC2O)(CC2O)",
    "(OC2)(CCO2)(CC2)(CC2O)(CC2O)(CC2O)"
  };
  static const std::set<std::string> allowed_shapes(
      std::begin(k_ace_sugar_shapes), std::end(k_ace_sugar_shapes));
  std::map<int, std::vector<size_t>> ring_atoms;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      ring_atoms[rid].push_back(i);
  for (auto& kv : ring_atoms) {
    auto& atoms = kv.second;
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    int n = (int)atoms.size();
    if (n != 5 && n != 6)
      continue;
    std::set<size_t> rset(atoms.begin(), atoms.end());
    size_t oxy_idx = SIZE_MAX;
    int oxy = 0, carb = 0, sp3 = 0;
    bool ok = true;
    for (size_t idx : rset) {
      if (atom_info[idx].bonding_idx == 3)
        ++sp3;
      if (atom_info[idx].bonding_idx != 3) {
        ok = false;
        break;
      }
      if (cc.atoms[idx].el == El::O) {
        ++oxy;
        oxy_idx = idx;
      } else if (cc.atoms[idx].el == El::C) {
        ++carb;
      } else {
        ok = false;
        break;
      }
      int ring_deg = 0;
      for (const auto& nb : adj[idx])
        if (rset.count(nb.idx))
          ++ring_deg;
      if (ring_deg != 2) {
        ok = false;
        break;
      }
    }
    if (!ok || sp3 != n || oxy != 1 || carb != n - 1 || oxy_idx == SIZE_MAX)
      continue;
    std::vector<size_t> o_ring_nbs;
    for (const auto& nb : adj[oxy_idx])
      if (rset.count(nb.idx))
        o_ring_nbs.push_back(nb.idx);
    if (o_ring_nbs.size() != 2)
      continue;
    std::vector<size_t> seq1, seq2;
    std::string shape1 = ring_shape(rset, oxy_idx, o_ring_nbs[0], &seq1);
    std::string shape2 = ring_shape(rset, oxy_idx, o_ring_nbs[1], &seq2);
    bool ok1 = allowed_shapes.count(shape1) != 0;
    bool ok2 = allowed_shapes.count(shape2) != 0;
    if (!ok1 && !ok2)
      continue;
    out.ring_seq[kv.first] = ok1 ? seq1 : seq2;
    out.ring_sets.emplace(kv.first, std::move(rset));
  }
  for (const auto& kv : out.ring_sets) {
    const auto& rset = kv.second;
    for (size_t idx : rset)
      for (const auto& nb : adj[idx])
        if (rset.count(nb.idx) && idx < nb.idx)
          out.ring_bonds.insert(std::make_pair(idx, nb.idx));
  }
  return out;
}

struct ChiralCenterInfo {
  std::set<size_t> stereo_chiral_centers;
  std::map<size_t, std::map<size_t, std::vector<size_t>>> chir_mut_table;
};

bool is_sp1_like(const CodAtomInfo& ai);
bool is_sp2_like(const CodAtomInfo& ai);
bool is_sp3_like(const CodAtomInfo& ai);

ChiralCenterInfo detect_chiral_centers_and_mut_table(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const std::map<std::string, std::string>& atom_stereo) {
  ChiralCenterInfo out;
  // Prefer explicitly listed/generated chirality rows when available.
  // AceDRG torsion code uses inChirals[0] and its mutTable directly.
  if (!cc.rt.chirs.empty()) {
    auto atom_index = cc.make_atom_index();
    for (size_t i = 0; i < cc.atoms.size(); ++i) {
      auto st_it = atom_stereo.find(cc.atoms[i].id);
      if (st_it == atom_stereo.end() || st_it->second.empty())
        continue;
      char s = lower(st_it->second[0]);
      if ((s == 'r' || s == 's') && cc.atoms[i].el == El::C)
        out.stereo_chiral_centers.insert(i);
    }
    for (const auto& chir : cc.rt.chirs) {
      auto c_it = atom_index.find(chir.id_ctr.atom);
      auto a1_it = atom_index.find(chir.id1.atom);
      auto a2_it = atom_index.find(chir.id2.atom);
      auto a3_it = atom_index.find(chir.id3.atom);
      if (c_it == atom_index.end() || a1_it == atom_index.end() ||
          a2_it == atom_index.end() || a3_it == atom_index.end())
        continue;
      size_t center = c_it->second;
      if (out.chir_mut_table.count(center) != 0)
        continue;  // AceDRG uses the first chiral record for a center.
      size_t a1 = a1_it->second;
      size_t a2 = a2_it->second;
      size_t a3 = a3_it->second;
      size_t missing = SIZE_MAX;
      for (const auto& nb : adj[center])
        if (nb.idx != a1 && nb.idx != a2 && nb.idx != a3) {
          missing = nb.idx;
          break;
        }
      bool negative = (chir.sign == ChiralityType::Negative);
      auto& mt = out.chir_mut_table[center];
      if (!negative) {
        mt[a1] = {a3, a2};
        mt[a2] = {a1, a3};
        mt[a3] = {a2, a1};
      } else {
        mt[a1] = {a2, a3};
        mt[a2] = {a3, a1};
        mt[a3] = {a1, a2};
      }
      if (missing != SIZE_MAX) {
        mt[a1].push_back(missing);
        mt[a2].push_back(missing);
        mt[a3].push_back(missing);
        if (!negative)
          mt[missing] = {a1, a2, a3};
        else
          mt[missing] = {a3, a2, a1};
      }
    }
    return out;
  }

  std::set<size_t> stereo_negative_centers;
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    if (atom_info[i].hybrid != Hybridization::SP3)
      continue;
    std::vector<size_t> non_h_nbs;
    for (const auto& nb : adj[i])
      if (!cc.atoms[nb.idx].is_hydrogen())
        non_h_nbs.push_back(nb.idx);
    std::vector<size_t> chiral_legs = non_h_nbs;
    if (cc.atoms[i].el == El::N && non_h_nbs.size() == 2) {
      bool n31_like = true;
      for (size_t nb : non_h_nbs)
        if (cc.atoms[nb].el != El::C || atom_info[nb].hybrid != Hybridization::SP3) {
          n31_like = false;
          break;
        }
      bool has_sp_neighbor = false;
      for (size_t nb : non_h_nbs)
        if (cc.atoms[nb].el == El::S || cc.atoms[nb].el == El::P) {
          has_sp_neighbor = true;
          break;
        }
      if (n31_like || has_sp_neighbor) {
        for (const auto& nb : adj[i])
          if (cc.atoms[nb.idx].is_hydrogen()) {
            chiral_legs.push_back(nb.idx);
            break;
          }
      }
    }
    auto shared_ring_count = [&](size_t a, size_t b) {
      int n = 0;
      for (int ra : atom_info[a].in_rings)
        for (int rb : atom_info[b].in_rings)
          if (ra == rb)
            ++n;
      return n;
    };
    bool has_halogen_nb = false;
    for (size_t nb : non_h_nbs) {
      Element e = cc.atoms[nb].el;
      if (is_halogen(e)) {
        has_halogen_nb = true;
        break;
      }
    }
    bool use_chiral_priority_sort = !(cc.atoms[i].el == El::C && has_halogen_nb);
    // Save bond-table order before sorting; chir_mut_table uses bond-table order for N.
    std::vector<size_t> bt_chiral_legs = chiral_legs;
    auto bond_to_center_type = [&](size_t nb_idx) {
      for (const auto& nb : adj[i])
        if (nb.idx == nb_idx)
          return nb.type;
      return BondType::Unspec;
    };
    if (use_chiral_priority_sort) {
      if (cc.atoms[i].el == El::P) {
        auto p_nb_rank = [&](size_t nb_idx) {
          if (cc.atoms[nb_idx].el != El::O)
            return 3;
          bool bridged_to_p = false;
          for (const auto& nb2 : adj[nb_idx]) {
            if (nb2.idx == i || cc.atoms[nb2.idx].is_hydrogen())
              continue;
            if (cc.atoms[nb2.idx].el == El::P) {
              bridged_to_p = true;
              break;
            }
          }
          if (bridged_to_p)
            return 0;
          BondType bt = bond_to_center_type(nb_idx);
          if (bt == BondType::Double || bt == BondType::Deloc)
            return 2;
          return 1;
        };
        std::stable_sort(chiral_legs.begin(), chiral_legs.end(),
                         [&](size_t a, size_t b) {
                           int ra = p_nb_rank(a);
                           int rb = p_nb_rank(b);
                           if (ra != rb)
                             return ra < rb;
                           int pa = chirality_priority(cc.atoms[a].el);
                           int pb = chirality_priority(cc.atoms[b].el);
                           return pa < pb;
                         });
      } else {
        auto is_pi_like_nb = [&](size_t nb_idx) {
          for (const auto& nb : adj[nb_idx])
            if (nb.idx != i &&
                (nb.type == BondType::Double ||
                 nb.type == BondType::Deloc ||
                 nb.type == BondType::Aromatic))
              return true;
          return false;
        };
        bool center_has_pi_nonh = false;
        if (cc.atoms[i].el == El::C) {
          for (size_t nb_idx : non_h_nbs)
            if (cc.atoms[nb_idx].el == El::C && is_pi_like_nb(nb_idx)) {
              center_has_pi_nonh = true;
              break;
            }
        }
        std::stable_sort(chiral_legs.begin(), chiral_legs.end(),
                         [&](size_t a, size_t b) {
                           if (cc.atoms[i].el == El::S) {
                             auto s_nb_rank = [&](size_t nb_idx) {
                               if (cc.atoms[nb_idx].el != El::O)
                                 return 2;
                               BondType bt = bond_to_center_type(nb_idx);
                               return (bt == BondType::Double || bt == BondType::Deloc) ? 0 : 1;
                             };
                             int ra = s_nb_rank(a);
                             int rb = s_nb_rank(b);
                             if (ra != rb)
                               return ra < rb;
                           }
                           int pa = chirality_priority(cc.atoms[a].el);
                           int pb = chirality_priority(cc.atoms[b].el);
                           if (pa != pb)
                             return pa < pb;
                           if (cc.atoms[i].el == El::C &&
                               cc.atoms[a].el == El::C &&
                               cc.atoms[b].el == El::C &&
                               center_has_pi_nonh)
                             return cc.atoms[a].id > cc.atoms[b].id;
                           if (cc.atoms[i].el == El::N) {
                             int sa = shared_ring_count(i, a);
                             int sb = shared_ring_count(i, b);
                             if (sa != sb)
                               return sa > sb;
                             return cc.atoms[a].id < cc.atoms[b].id;
                           }
                           return false;
                         });
      }
    }
    if (cc.atoms[i].el == El::N && non_h_nbs.size() == 2 && chiral_legs.size() >= 3) {
      size_t s_idx = SIZE_MAX, o_idx = SIZE_MAX, h_idx = SIZE_MAX;
      for (size_t a : chiral_legs) {
        if (cc.atoms[a].el == El::S)
          s_idx = a;
        else if (cc.atoms[a].el == El::O)
          o_idx = a;
        else if (cc.atoms[a].is_hydrogen() && h_idx == SIZE_MAX)
          h_idx = a;
      }
      if (s_idx != SIZE_MAX && o_idx != SIZE_MAX && h_idx != SIZE_MAX)
        chiral_legs = {s_idx, o_idx, h_idx};
    }
    if (chiral_legs.size() < 3)
      continue;
    bool is_stereo_carbon = false;
    bool is_stereo_r = false;
    if (cc.atoms[i].el == El::C) {
      auto st_it = atom_stereo.find(cc.atoms[i].id);
      if (st_it != atom_stereo.end() && !st_it->second.empty()) {
        char s = lower(st_it->second[0]);
        is_stereo_carbon = (s == 'r' || s == 's');
        is_stereo_r = (s == 'r');
      }
    }
    if (is_stereo_carbon) {
      if (is_stereo_r)
        stereo_negative_centers.insert(i);
      out.stereo_chiral_centers.insert(i);
      std::stable_sort(chiral_legs.begin(), chiral_legs.end(),
                       [&](size_t a, size_t b) {
                         return chirality_priority(cc.atoms[a].el) <
                                chirality_priority(cc.atoms[b].el);
                       });
    }
    if (cc.atoms[i].el == El::C && chiral_legs.size() >= 4) {
      std::vector<size_t> halogens;
      bool has_oxygen = false;
      for (size_t nb : chiral_legs) {
        Element e = cc.atoms[nb].el;
        if (is_halogen(e))
          halogens.push_back(nb);
        if (e == El::O)
          has_oxygen = true;
      }
      if (has_oxygen && halogens.size() >= 3) {
        chiral_legs.clear();
        chiral_legs.push_back(halogens[0]);
        chiral_legs.push_back(halogens[1]);
        chiral_legs.push_back(halogens[2]);
      }
    }
    // For N: COD CIF uses bond-table order; for S/P/C: sorted order matches.
    const std::vector<size_t>& legs_for_mt =
        (cc.atoms[i].el == El::N) ? bt_chiral_legs : chiral_legs;
    size_t a1 = legs_for_mt[0], a2 = legs_for_mt[1], a3 = legs_for_mt[2];
    size_t missing = SIZE_MAX;
    for (const auto& nb : adj[i])
      if (nb.idx != a1 && nb.idx != a2 && nb.idx != a3) {
        missing = nb.idx;
        break;
      }
    auto& mt = out.chir_mut_table[i];
    bool negative_chiral = stereo_negative_centers.count(i) != 0;
    if (!negative_chiral) {
      mt[a1] = {a3, a2};
      mt[a2] = {a1, a3};
      mt[a3] = {a2, a1};
    } else {
      mt[a1] = {a2, a3};
      mt[a2] = {a3, a1};
      mt[a3] = {a1, a2};
    }
    if (missing != SIZE_MAX) {
      mt[a1].push_back(missing);
      mt[a2].push_back(missing);
      mt[a3].push_back(missing);
      if (!negative_chiral)
        mt[missing] = {a1, a2, a3};
      else
        mt[missing] = {a3, a2, a1};
    }
  }
  return out;
}

std::vector<bool> build_aromatic_like_mask(
    const ChemComp& cc, const std::vector<CodAtomInfo>& atom_info,
    const std::map<std::string, size_t>& atom_index) {
  std::vector<bool> aromatic_like(cc.atoms.size(), false);
  for (size_t i = 0; i < atom_info.size(); ++i)
    aromatic_like[i] = atom_info[i].is_aromatic;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.type == BondType::Aromatic || bond.type == BondType::Deloc || bond.aromatic) {
      auto it1 = atom_index.find(bond.id1.atom);
      auto it2 = atom_index.find(bond.id2.atom);
      if (it1 != atom_index.end())
        aromatic_like[it1->second] = true;
      if (it2 != atom_index.end())
        aromatic_like[it2->second] = true;
    }
  }
  return aromatic_like;
}

bool is_sp1_like(const CodAtomInfo& ai) {
  return ai.bonding_idx == 1;
}

bool is_sp2_like(const CodAtomInfo& ai) {
  return ai.bonding_idx == 2;
}

bool is_sp3_like(const CodAtomInfo& ai) {
  return ai.bonding_idx == 3;
}

enum class TvMode { Default, SP3SP3, SP2SP3_SP3, SP3_OXY };

int max_tv_len_for_center(const CodAtomInfo& ai) {
  return is_sp2_like(ai) ? 2 : 3;
}

std::vector<size_t> build_tv_neighbor_order(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info, size_t ctr) {
  (void)cc;
  (void)atom_info;
  std::vector<size_t> nb_order;
  nb_order.reserve(adj[ctr].size());
  for (const auto& nb : adj[ctr])
    nb_order.push_back(nb.idx);
  return nb_order;
}

void append_chiral_cluster_like_acedrg(std::vector<size_t>& out,
                                       const std::vector<size_t>& mut_list) {
  // AceDRG buildChiralCluster2/getMutList behavior:
  // - if output is empty: append mut_list as-is.
  // - if output has one seeded atom and it exists in mut_list:
  //   append mut_list rotated after that atom (excluding it).
  // - otherwise append mut_list as-is.
  if (mut_list.empty())
    return;
  if (out.size() == 1) {
    auto it = std::find(mut_list.begin(), mut_list.end(), out[0]);
    if (it != mut_list.end()) {
      size_t pos = (size_t)std::distance(mut_list.begin(), it);
      for (size_t i = 1; i < mut_list.size(); ++i) {
        size_t j = (pos + i) % mut_list.size();
        if (std::find(out.begin(), out.end(), mut_list[j]) == out.end())
          out.push_back(mut_list[j]);
      }
      return;
    }
  }
  for (size_t cand : mut_list)
    if (std::find(out.begin(), out.end(), cand) == out.end())
      out.push_back(cand);
}

// Build the torsion-value neighbor ordering (tV) list for a center with
// respect to the bonded atom on the other side.
std::vector<size_t> build_tv_list_for_center(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const std::set<size_t>& stereo_chiral_centers,
    const std::map<size_t, std::map<size_t, std::vector<size_t>>>& chir_mut_table,
    size_t ctr, size_t other, TvMode mode, size_t forced_rs = SIZE_MAX) {
  std::vector<size_t> tv;
  size_t rs = forced_rs;
  if (rs == SIZE_MAX)
    rs = find_ring_sharing_pair(adj, atom_info, ctr, other).first;
  if (rs != SIZE_MAX)
    tv.push_back(rs);

  int max_len = max_tv_len_for_center(atom_info[ctr]);
  bool used_chiral = false;
  bool stereo_sp3 = (atom_info[ctr].hybrid == Hybridization::SP3 &&
                     stereo_chiral_centers.count(ctr) != 0);
  if (atom_info[ctr].hybrid == Hybridization::SP3) {
    bool use_chiral_mut_table = !stereo_sp3 && cc.atoms[ctr].el != El::P;
    auto mit = chir_mut_table.find(ctr);
    if (use_chiral_mut_table && mit != chir_mut_table.end()) {
      auto mt_it = mit->second.find(other);
      if (mt_it != mit->second.end()) {
        used_chiral = true;
        std::vector<size_t> mut_filtered;
        mut_filtered.reserve(mt_it->second.size());
        for (size_t cand : mt_it->second)
          if (cand != other && !cc.atoms[cand].el.is_metal())
            mut_filtered.push_back(cand);
        append_chiral_cluster_like_acedrg(tv, mut_filtered);
      }
    }
  }

  if (!used_chiral) {
    bool want_h_first = (mode == TvMode::SP2SP3_SP3 &&
                         atom_info[ctr].hybrid == Hybridization::SP3 &&
                         stereo_chiral_centers.count(ctr) == 0 &&
                         !is_oxygen_column(cc.atoms[other].el) &&
                         rs == SIZE_MAX);
    bool want_non_h_first = (mode == TvMode::SP3_OXY &&
                             atom_info[ctr].hybrid == Hybridization::SP3 &&
                             rs == SIZE_MAX);
    std::vector<size_t> nb_order =
        build_tv_neighbor_order(cc, adj, atom_info, ctr);
    size_t first_non_h = SIZE_MAX;
    if (want_non_h_first && (int)tv.size() < max_len) {
      for (size_t nb_idx : nb_order) {
        if (nb_idx == other || nb_idx == rs)
          continue;
        if (!cc.atoms[nb_idx].is_hydrogen()) {
          first_non_h = nb_idx;
          tv.push_back(first_non_h);
          break;
        }
      }
      if (first_non_h == SIZE_MAX && !adj[ctr].empty()) {
        first_non_h = adj[ctr][0].idx;
        tv.push_back(first_non_h);
      }
    }
    size_t first_h = SIZE_MAX;
    if (want_h_first && (int)tv.size() < max_len) {
      for (size_t nb_idx : nb_order) {
        if (nb_idx == other || nb_idx == rs)
          continue;
        if (cc.atoms[nb_idx].is_hydrogen()) {
          first_h = nb_idx;
          tv.push_back(first_h);
          break;
        }
      }
    }

    for (size_t nb_idx : nb_order) {
      if ((int)tv.size() >= max_len)
        break;
      if (nb_idx == other || nb_idx == rs ||
          nb_idx == first_h || nb_idx == first_non_h)
        continue;
      tv.push_back(nb_idx);
    }
  }

  if ((int)tv.size() > max_len)
    tv.resize(max_len);
  return tv;
}

std::vector<size_t> build_tv_list_sp3sp3_like_acedrg(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const std::set<size_t>& stereo_chiral_centers,
    const std::map<size_t, std::map<size_t, std::vector<size_t>>>& chir_mut_table,
    size_t center, size_t other, size_t rs = SIZE_MAX) {
  std::vector<size_t> tv;
  if (rs != SIZE_MAX)
    tv.push_back(rs);

  bool used_chiral = false;
  bool stereo_sp3 = (atom_info[center].hybrid == Hybridization::SP3 &&
                     stereo_chiral_centers.count(center) != 0);
  if (atom_info[center].hybrid == Hybridization::SP3 && !stereo_sp3) {
    auto mit = chir_mut_table.find(center);
    if (mit != chir_mut_table.end()) {
      auto mt_it = mit->second.find(other);
      if (mt_it != mit->second.end()) {
        used_chiral = true;
        std::vector<size_t> mut_filtered;
        mut_filtered.reserve(mt_it->second.size());
        for (size_t cand : mt_it->second)
          if (cand != other)
            mut_filtered.push_back(cand);
        append_chiral_cluster_like_acedrg(tv, mut_filtered);
      }
    }
  }

  if (!used_chiral) {
    for (const auto& nb : adj[center]) {
      if (nb.idx != other &&
          std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
          !cc.atoms[nb.idx].el.is_metal() &&
          !cc.atoms[nb.idx].is_hydrogen()) {
        tv.push_back(nb.idx);
        break;
      }
    }
  }

  for (const auto& nb : adj[center]) {
    if (nb.idx == other ||
        cc.atoms[nb.idx].el.is_metal() ||
        std::find(tv.begin(), tv.end(), nb.idx) != tv.end())
      continue;
    tv.push_back(nb.idx);
  }
  if (tv.size() > 3)
    tv.resize(3);
  return tv;
}

std::pair<std::vector<size_t>, std::vector<size_t>> build_tv_lists_sp2sp2_like_acedrg(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center1, size_t center2) {
  std::vector<size_t> tv1, tv2;
  size_t s1 = SIZE_MAX, s2 = SIZE_MAX;
  for (const auto& nb1 : adj[center1]) {
    if (nb1.idx == center2)
      continue;
    for (const auto& nb2 : adj[center2]) {
      if (nb2.idx == center1 || nb1.idx == nb2.idx)
        continue;
      if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
        s1 = nb1.idx;
        s2 = nb2.idx;
        tv1.push_back(s1);
        tv2.push_back(s2);
        break;
      }
    }
    if (s1 != SIZE_MAX && s2 != SIZE_MAX)
      break;
  }

  for (const auto& nb : adj[center1])
    if (nb.idx != center2 && nb.idx != s1)
      tv1.push_back(nb.idx);
  for (const auto& nb : adj[center2])
    if (nb.idx != center1 && nb.idx != s2)
      tv2.push_back(nb.idx);

  if (s1 == SIZE_MAX && s2 == SIZE_MAX && tv1.size() == 2 && tv2.size() == 2) {
    auto h_only_excluding = [&](size_t idx, size_t excl) {
      for (const auto& nb : adj[idx])
        if (nb.idx != excl && !cc.atoms[nb.idx].is_hydrogen())
          return false;
      return true;
    };
    auto degree = [&](size_t idx) { return adj[idx].size(); };

    bool h0 = h_only_excluding(tv1[0], center1);
    bool h1 = h_only_excluding(tv1[1], center1);
    if (h0 && !h1)
      std::swap(tv1[0], tv1[1]);
    else if (degree(tv1[0]) < degree(tv1[1]))
      std::swap(tv1[0], tv1[1]);

    h0 = h_only_excluding(tv2[0], center2);
    h1 = h_only_excluding(tv2[1], center2);
    if (h0) {
      if (!h1 || degree(tv2[0]) < degree(tv2[1]))
        std::swap(tv2[0], tv2[1]);
    } else if (degree(tv2[0]) < degree(tv2[1]) && !h1) {
      std::swap(tv2[0], tv2[1]);
    }
  }
  return {std::move(tv1), std::move(tv2)};
}


static void emit_one_torsion(
    ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const AcedrgTables& tables,
    const std::set<size_t>& stereo_chiral_centers,
    const std::map<size_t, std::map<size_t, std::vector<size_t>>>& chir_mut_table,
    const std::map<std::pair<size_t, size_t>, RingParity>& bond_ring_parity,
    const std::map<std::pair<size_t, size_t>, RingFlip>& bond_ring_flip,
    bool peptide_mode,
    size_t center2, size_t center3,
    bool bond_aromatic, int ring_size,
    size_t a1_idx, size_t a4_idx,
    std::vector<Restraints::Torsion>* out_torsions = nullptr) {
  const CodAtomInfo& info2 = atom_info[center2];
  const CodAtomInfo& info3 = atom_info[center3];
  bool sp3_2 = is_sp3_like(info2);
  bool sp3_3 = is_sp3_like(info3);
  bool sp2_2 = is_sp2_like(info2);
  bool sp2_3 = is_sp2_like(info3);
  if (a1_idx == a4_idx || a1_idx == center3 || a4_idx == center2)
    return;
  if (cc.atoms[a1_idx].el.is_metal() || cc.atoms[a4_idx].el.is_metal())
    return;
  if (cc.atoms[a1_idx].id == cc.atoms[a4_idx].id)
    return;
  const ChemComp::Atom* a1 = &cc.atoms[a1_idx];
  const ChemComp::Atom* a4 = &cc.atoms[a4_idx];
  TorsionEntry tors_entry;
  double value = 180.0;
  double esd = 10.0;
  int period = 3;
  bool lookup_found = false;
  if (peptide_mode)
    lookup_found = tables.lookup_pep_tors(a1->id, cc.atoms[center2].id,
                                          cc.atoms[center3].id, a4->id,
                                          tors_entry);
  if (lookup_found) {
    value = tors_entry.value;
    esd = 10.0;
    period = tors_entry.period;
  }

  if (!lookup_found && bond_aromatic && ring_size > 0 && info2.is_aromatic && info3.is_aromatic &&
             !(sp2_2 && sp2_3)) {
    auto shares_ring_across = [&](size_t terminal_idx, size_t opp_center) {
      for (const auto& nb : adj[opp_center])
        if (nb.idx != center2 && nb.idx != center3 &&
            share_ring_ids(atom_info[terminal_idx].in_rings, atom_info[nb.idx].in_rings))
          return true;
      return false;
    };
    bool a1_ring = shares_ring_across(a1_idx, center3);
    bool a4_ring = shares_ring_across(a4_idx, center2);
    value = (a1_ring != a4_ring) ? 180.0 : 0.0;
    esd = 0.0;
    period = 1;
  } else if (!lookup_found && ring_size > 0 && sp3_2 && sp3_3) {
    // SP3-SP3 ring bond: use even/odd/no-flip matrix based on ring traversal parity
    auto parity_key = std::minmax(center2, center3);
    auto pit = bond_ring_parity.find(parity_key);
    size_t side1 = center2;
    size_t side2 = center3;
    std::pair<size_t, size_t> rs_pair =
        find_ring_sharing_pair(adj, atom_info, side1, side2);
    size_t rs1 = rs_pair.first;
    size_t rs2 = rs_pair.second;
    size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
    size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
    std::vector<size_t> tv1 = build_tv_list_sp3sp3_like_acedrg(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side1, side2, rs1);
    std::vector<size_t> tv2 = build_tv_list_sp3sp3_like_acedrg(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side2, side1, rs2);
    if (const char* env = std::getenv("GEMMI_DBG_SP3")) {
      std::string v(env);
      if (v == "1" || v.find(cc.name) != std::string::npos) {
        std::fprintf(stderr, "[sp3 %s ring=%d rs=%s/%s] %s-%s term=%s/%s tv1=",
                     cc.name.c_str(), ring_size, cc.atoms[center2].id.c_str(),
                     rs1 != SIZE_MAX ? cc.atoms[rs1].id.c_str() : "-",
                     rs2 != SIZE_MAX ? cc.atoms[rs2].id.c_str() : "-",
                     cc.atoms[center3].id.c_str(), cc.atoms[term1].id.c_str(),
                     cc.atoms[term2].id.c_str());
        for (size_t i = 0; i < tv1.size(); ++i)
          std::fprintf(stderr, "%s%s", i ? "," : "", cc.atoms[tv1[i]].id.c_str());
        std::fprintf(stderr, " tv2=");
        for (size_t i = 0; i < tv2.size(); ++i)
          std::fprintf(stderr, "%s%s", i ? "," : "", cc.atoms[tv2[i]].id.c_str());
        std::fprintf(stderr, "\n");
      }
    }
    int i_pos = -1, j_pos = -1;
    for (int i = 0; i < (int)tv1.size(); ++i)
      if (tv1[i] == term1) { i_pos = i; break; }
    for (int j = 0; j < (int)tv2.size(); ++j)
      if (tv2[j] == term2) { j_pos = j; break; }
    if (i_pos >= 0 && j_pos >= 0) {
      static const double even_m[3][3] = {
        {60,180,-60}, {-60,60,180}, {180,-60,60}};
      static const double odd_m[3][3] = {
        {-60,60,180}, {180,-60,60}, {60,180,-60}};
      static const double noflip_m[3][3] = {
        {180,-60,60}, {60,180,-60}, {-60,60,180}};
      auto fit = bond_ring_flip.find(parity_key);
      const auto& m = (fit != bond_ring_flip.end() && fit->second == RingFlip::Even) ? even_m
                     : (fit != bond_ring_flip.end() && fit->second == RingFlip::Odd) ? odd_m
                     : (pit != bond_ring_parity.end() && pit->second == RingParity::Even) ? even_m
                     : (pit != bond_ring_parity.end() && pit->second == RingParity::Odd) ? odd_m
                     : noflip_m;
      value = m[i_pos][j_pos];
    }
  } else if (!lookup_found && sp3_2 && sp3_3 && ring_size == 0) {
    // Non-ring SP3-SP3.
    size_t side1 = center2;
    size_t side2 = center3;
    std::pair<size_t, size_t> rs_pair =
        find_ring_sharing_pair(adj, atom_info, side1, side2);
    size_t rs1 = rs_pair.first;
    size_t rs2 = rs_pair.second;
    size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
    size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
    std::vector<size_t> tv1 = build_tv_list_sp3sp3_like_acedrg(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side1, side2, rs1);
    std::vector<size_t> tv2 = build_tv_list_sp3sp3_like_acedrg(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side2, side1, rs2);
    if (const char* env = std::getenv("GEMMI_DBG_SP3")) {
      std::string v(env);
      if (v == "1" || v.find(cc.name) != std::string::npos) {
        std::fprintf(stderr, "[sp3 %s ring=%d rs=%s/%s] %s-%s term=%s/%s tv1=",
                     cc.name.c_str(), ring_size, cc.atoms[center2].id.c_str(),
                     rs1 != SIZE_MAX ? cc.atoms[rs1].id.c_str() : "-",
                     rs2 != SIZE_MAX ? cc.atoms[rs2].id.c_str() : "-",
                     cc.atoms[center3].id.c_str(), cc.atoms[term1].id.c_str(),
                     cc.atoms[term2].id.c_str());
        for (size_t i = 0; i < tv1.size(); ++i)
          std::fprintf(stderr, "%s%s", i ? "," : "", cc.atoms[tv1[i]].id.c_str());
        std::fprintf(stderr, " tv2=");
        for (size_t i = 0; i < tv2.size(); ++i)
          std::fprintf(stderr, "%s%s", i ? "," : "", cc.atoms[tv2[i]].id.c_str());
        std::fprintf(stderr, "\n");
      }
    }
    int i_pos = -1, j_pos = -1;
    for (int i = 0; i < (int)tv1.size(); ++i)
      if (tv1[i] == term1) { i_pos = i; break; }
    for (int j = 0; j < (int)tv2.size(); ++j)
      if (tv2[j] == term2) { j_pos = j; break; }
    if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
      static const double noflip_m[3][3] = {
        {180,-60,60}, {60,180,-60}, {-60,60,180}};
      value = noflip_m[i_pos][j_pos];
    }
  } else if (!lookup_found && ((sp2_2 && sp3_3) || (sp3_2 && sp2_3))) {
    // SP2-SP3:
    // - regular path: tS3 matrix for "two SP2-side neighbours in same ring",
    //   otherwise default matrix (AceDRG SetOneSP2SP3Bond).
    // - Oxy-column SP2 center (O/S/Se/Te/Po): use SetOneSP3OxyColumnBond
    //   style period-3 progression.
    size_t sp2_center = sp2_2 ? center2 : center3;
    size_t sp3_center = sp2_2 ? center3 : center2;
    auto is_nb = [&](size_t ctr, size_t idx) {
      for (const auto& nb : adj[ctr])
        if (nb.idx == idx)
          return true;
      return false;
    };
    size_t sp2_term = SIZE_MAX;
    size_t sp3_term = SIZE_MAX;
    if (is_nb(sp2_center, a1_idx) && is_nb(sp3_center, a4_idx)) {
      sp2_term = a1_idx;
      sp3_term = a4_idx;
    } else if (is_nb(sp2_center, a4_idx) && is_nb(sp3_center, a1_idx)) {
      sp2_term = a4_idx;
      sp3_term = a1_idx;
    } else {
      // Fallback for malformed inputs: keep legacy orientation.
      sp2_term = sp2_2 ? a1_idx : a4_idx;
      sp3_term = sp2_2 ? a4_idx : a1_idx;
    }
    bool dbg_sp23 = false;
    if (const char* env = std::getenv("GEMMI_DBG_SP23")) {
      std::string v(env);
      dbg_sp23 = (v == "1" || v.find(cc.name) != std::string::npos);
    }
    std::pair<size_t, size_t> rs_pair =
        find_ring_sharing_pair(adj, atom_info, sp2_center, sp3_center);
    size_t sp2_rs = rs_pair.first;
    size_t sp3_rs = rs_pair.second;
    bool oxy_col_sp2 = is_oxygen_column(cc.atoms[sp2_center].el);
    bool oxy_col_sp3 = is_oxygen_column(cc.atoms[sp3_center].el);

    auto normalize_acedrg_angle = [](double ang) {
      if (ang > 360.0)
        ang -= 360.0;
      if (ang > 180.0)
        ang -= 360.0;
      return ang;
    };

    if (oxy_col_sp2) {
      // Direct AceDRG SetOneSP3OxyColumnBond(tIdx1=SP3, tIdx2=SP2-oxy).
      std::vector<size_t> tV1, tV2;
      if (sp3_rs != SIZE_MAX && sp2_rs != SIZE_MAX)
        tV1.push_back(sp3_rs);
      if (tV1.empty()) {
        size_t iS = SIZE_MAX;
        for (const auto& nb : adj[sp3_center]) {
          if (nb.idx != sp2_center && !cc.atoms[nb.idx].is_hydrogen()) {
            iS = nb.idx;
            break;
          }
        }
        if (iS == SIZE_MAX && !adj[sp3_center].empty())
          iS = adj[sp3_center][0].idx;
        if (iS != SIZE_MAX)
          tV1.push_back(iS);
      }
      bool used_chiral = false;
      bool stereo_sp3 = (atom_info[sp3_center].hybrid == Hybridization::SP3 &&
                         stereo_chiral_centers.count(sp3_center) != 0);
      if (atom_info[sp3_center].hybrid == Hybridization::SP3 &&
          !stereo_sp3 && cc.atoms[sp3_center].el != El::P) {
        auto mit = chir_mut_table.find(sp3_center);
        if (mit != chir_mut_table.end()) {
          auto mt_it = mit->second.find(sp2_center);
          if (mt_it != mit->second.end()) {
            used_chiral = true;
            std::vector<size_t> mut_filtered;
            mut_filtered.reserve(mt_it->second.size());
            for (size_t cand : mt_it->second)
              if (cand != sp2_center)
                mut_filtered.push_back(cand);
            append_chiral_cluster_like_acedrg(tV1, mut_filtered);
          }
        }
      }
      (void)used_chiral;
      for (const auto& nb : adj[sp3_center]) {
        if (nb.idx == sp2_center ||
            std::find(tV1.begin(), tV1.end(), nb.idx) != tV1.end())
          continue;
        tV1.push_back(nb.idx);
      }
      for (const auto& nb : adj[sp2_center]) {
        if (nb.idx == sp3_center ||
            std::find(tV2.begin(), tV2.end(), nb.idx) != tV2.end())
          continue;
        tV2.push_back(nb.idx);
      }
      int i_pos = -1, j_pos = -1;
      for (int i = 0; i < (int)tV1.size(); ++i)
        if (tV1[i] == sp3_term) { i_pos = i; break; }
      for (int j = 0; j < (int)tV2.size(); ++j)
        if (tV2[j] == sp2_term) { j_pos = j; break; }
      if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
        auto fkey = std::minmax(sp2_center, sp3_center);
        auto fit = bond_ring_flip.find(fkey);
        double ini_value = 180.0;
        if (sp2_rs != SIZE_MAX) {
          if (fit != bond_ring_flip.end())
            ini_value = (fit->second == RingFlip::Even) ? 60.0 : -60.0;
          else
            ini_value = 60.0;
        } else if (oxy_col_sp3) {
          ini_value = 90.0;
        }
        value = normalize_acedrg_angle(ini_value + (i_pos + j_pos) * 120.0);
        esd = 20.0;
        period = 3;
        if (dbg_sp23) {
          std::fprintf(stderr, "[sp23 %s oxy] %s-%s-%s-%s i=%d j=%d ini=%.1f val=%.1f tv_sp3=",
                       cc.name.c_str(), cc.atoms[sp2_term].id.c_str(),
                       cc.atoms[sp2_center].id.c_str(), cc.atoms[sp3_center].id.c_str(),
                       cc.atoms[sp3_term].id.c_str(), i_pos, j_pos, ini_value, value);
          for (size_t ti = 0; ti < tV1.size(); ++ti)
            std::fprintf(stderr, "%s%s", ti ? "," : "", cc.atoms[tV1[ti]].id.c_str());
          std::fprintf(stderr, " tv_sp2=");
          for (size_t tj = 0; tj < tV2.size(); ++tj)
            std::fprintf(stderr, "%s%s", tj ? "," : "", cc.atoms[tV2[tj]].id.c_str());
          std::fprintf(stderr, "\n");
        }
      }
    } else {
      std::vector<size_t> tv_sp2;
      if (sp2_rs != SIZE_MAX && sp3_rs != SIZE_MAX)
        tv_sp2.push_back(sp2_rs);
      // Keep SP2-side ordering deterministic with heavy neighbors first.
      for (const auto& nb : adj[sp2_center]) {
        if (nb.idx == sp3_center || nb.idx == sp2_rs)
          continue;
        if (cc.atoms[nb.idx].el.is_metal())
          continue;
        if (!cc.atoms[nb.idx].is_hydrogen())
          tv_sp2.push_back(nb.idx);
      }
      for (const auto& nb : adj[sp2_center]) {
        if (nb.idx == sp3_center || nb.idx == sp2_rs)
          continue;
        if (cc.atoms[nb.idx].el.is_metal())
          continue;
        if (cc.atoms[nb.idx].is_hydrogen())
          tv_sp2.push_back(nb.idx);
      }

      std::vector<size_t> tv_sp3;
      if (sp2_rs != SIZE_MAX && sp3_rs != SIZE_MAX)
        tv_sp3.push_back(sp3_rs);
      bool used_chiral = false;
      bool stereo_sp3 = (atom_info[sp3_center].hybrid == Hybridization::SP3 &&
                         stereo_chiral_centers.count(sp3_center) != 0);
      if (atom_info[sp3_center].hybrid == Hybridization::SP3 && !stereo_sp3) {
        auto mit = chir_mut_table.find(sp3_center);
        if (mit != chir_mut_table.end()) {
          auto mt_it = mit->second.find(sp2_center);
          if (mt_it != mit->second.end()) {
            used_chiral = true;
            std::vector<size_t> mut_filtered;
            mut_filtered.reserve(mt_it->second.size());
            for (size_t cand : mt_it->second)
              if (cand != sp2_center)
                mut_filtered.push_back(cand);
            append_chiral_cluster_like_acedrg(tv_sp3, mut_filtered);
          }
        }
      }
      if (!used_chiral) {
        int h_count = 0;
        for (const auto& nb : adj[sp3_center])
          if (nb.idx != sp2_center &&
              !cc.atoms[nb.idx].el.is_metal() &&
              cc.atoms[nb.idx].is_hydrogen())
            ++h_count;
        if (!cc.atoms[sp2_center].is_hydrogen() && h_count == 2) {
          for (const auto& nb : adj[sp3_center]) {
            if (nb.idx == sp2_center)
              continue;
            if (cc.atoms[nb.idx].el.is_metal())
              continue;
            if (cc.atoms[nb.idx].is_hydrogen()) {
              tv_sp3.push_back(nb.idx);
              break;
            }
          }
          for (const auto& nb : adj[sp3_center]) {
            if (nb.idx == sp2_center)
              continue;
            if (cc.atoms[nb.idx].el.is_metal())
              continue;
            if (!cc.atoms[nb.idx].is_hydrogen() &&
                std::find(tv_sp3.begin(), tv_sp3.end(), nb.idx) == tv_sp3.end()) {
              tv_sp3.push_back(nb.idx);
              break;
            }
          }
        }
      }
      for (const auto& nb : adj[sp3_center]) {
        if (nb.idx == sp2_center)
          continue;
        if (cc.atoms[nb.idx].el.is_metal())
          continue;
        if (std::find(tv_sp3.begin(), tv_sp3.end(), nb.idx) != tv_sp3.end())
          continue;
        tv_sp3.push_back(nb.idx);
      }
      if (sp3_rs == SIZE_MAX && cc.atoms[sp3_center].el == El::P)
        std::stable_sort(tv_sp3.begin(), tv_sp3.end());

      static const double ts3_m[2][3] = {{150,-90,30}, {-30,90,-150}};
      static const double ts1_m[2][3] = {{0,120,-120}, {180,-60,60}};
      std::vector<size_t> sp2_side_nbs;
      for (const auto& nb : adj[sp2_center])
        if (nb.idx != sp3_center)
          sp2_side_nbs.push_back(nb.idx);
      bool has_ts3_ring = (sp2_side_nbs.size() == 2 &&
                           share_ring_ids(atom_info[sp2_side_nbs[0]].in_rings,
                                          atom_info[sp2_side_nbs[1]].in_rings));
      bool is_ts3 = (sp2_rs == SIZE_MAX && has_ts3_ring);
      if (is_ts3 && tv_sp3.size() == 2 &&
          !cc.atoms[tv_sp3[0]].is_hydrogen() && cc.atoms[tv_sp3[1]].is_hydrogen()) {
        std::swap(tv_sp3[0], tv_sp3[1]);
      }
      const auto& sp2sp3_m = is_ts3 ? ts3_m : ts1_m;
      int i_pos = -1, j_pos = -1;
      bool normal_shape = tv_sp2.size() <= tv_sp3.size();
      if (normal_shape) {
        for (int i = 0; i < (int)tv_sp2.size(); ++i)
          if (tv_sp2[i] == sp2_term) { i_pos = i; break; }
        for (int j = 0; j < (int)tv_sp3.size(); ++j)
          if (tv_sp3[j] == sp3_term) { j_pos = j; break; }
      } else {
        // AceDRG else-branch emits [tv_sp3, sp2, sp3, tv_sp2] but still uses va[i][j].
        for (int i = 0; i < (int)tv_sp3.size(); ++i)
          if (tv_sp3[i] == sp3_term) { i_pos = i; break; }
        for (int j = 0; j < (int)tv_sp2.size(); ++j)
          if (tv_sp2[j] == sp2_term) { j_pos = j; break; }
      }
      if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 3) {
        value = sp2sp3_m[i_pos][j_pos];
        esd = 20.0;
        period = 6;
        if (dbg_sp23) {
          std::fprintf(stderr, "[sp23 %s %s] %s-%s-%s-%s i=%d j=%d val=%.1f rs=%s/%s tv_sp2=",
                       cc.name.c_str(), is_ts3 ? "ts3" : "ts1",
                       cc.atoms[sp2_term].id.c_str(),
                       cc.atoms[sp2_center].id.c_str(), cc.atoms[sp3_center].id.c_str(),
                       cc.atoms[sp3_term].id.c_str(), i_pos, j_pos, value,
                       sp2_rs == SIZE_MAX ? "-" : cc.atoms[sp2_rs].id.c_str(),
                       sp3_rs == SIZE_MAX ? "-" : cc.atoms[sp3_rs].id.c_str());
          for (size_t ti = 0; ti < tv_sp2.size(); ++ti)
            std::fprintf(stderr, "%s%s", ti ? "," : "", cc.atoms[tv_sp2[ti]].id.c_str());
          std::fprintf(stderr, " tv_sp3=");
          for (size_t tj = 0; tj < tv_sp3.size(); ++tj)
            std::fprintf(stderr, "%s%s", tj ? "," : "", cc.atoms[tv_sp3[tj]].id.c_str());
          std::fprintf(stderr, " normal=%d\n", normal_shape ? 1 : 0);
        }
      }
    }
  } else if (!lookup_found && sp2_2 && sp2_3) {
    // Non-aromatic SP2-SP2: use 2x2 matrix.
    // Ring-sharing pair is selected once globally (AceDRG tS1/tS2 logic).
    bool dbg_sp2 = false;
    if (const char* env = std::getenv("GEMMI_DBG_SP2")) {
      std::string v(env);
      dbg_sp2 = (v == "1" || v.find(cc.name) != std::string::npos);
    }
    size_t side1 = center2;
    size_t side2 = center3;
    // DictCifFile::SetOneSP2SP2Bond() does not exclude nb1==nb2 when
    // detecting ring-sharing terminals. Keep that quirk here.
    size_t rs1 = SIZE_MAX, rs2 = SIZE_MAX;
    for (const auto& nb1 : adj[side1]) {
      if (nb1.idx == side2)
        continue;
      for (const auto& nb2 : adj[side2]) {
        if (nb2.idx == side1)
          continue;
        if (share_ring_ids(atom_info[nb1.idx].in_rings,
                           atom_info[nb2.idx].in_rings)) {
          rs1 = nb1.idx;
          rs2 = nb2.idx;
          break;
        }
      }
      if (rs1 != SIZE_MAX)
        break;
    }
    bool has_ring_sharing = (rs1 != SIZE_MAX && rs2 != SIZE_MAX);
    int side1_nb_count = 0;
    int side2_nb_count = 0;
    for (const auto& nb : adj[side1])
      if (nb.idx != side2)
        ++side1_nb_count;
    for (const auto& nb : adj[side2])
      if (nb.idx != side1)
        ++side2_nb_count;
    static const double ring_m[2][2] = {{0,180},{180,0}};
    static const double noring_m[2][2] = {{180,0},{0,180}};
    // SP2-SP2 neighbor ordering: for non-ring case with 2 neighbors per side,
    // AceDRG reorders by "H-only" check then connectivity tiebreaker.
    auto sp2sp2_tv_pos = [&](size_t center, size_t other, size_t target,
                             bool is_side1) -> int {
      auto conn_degree = [&](size_t idx) {
        return (int)adj[idx].size();
      };
      std::vector<size_t> tv;
      // Ring-sharing atom from the globally selected pair.
      if (has_ring_sharing) {
        size_t rs = is_side1 ? rs1 : rs2;
        if (rs != SIZE_MAX)
          tv.push_back(rs);
      }
      // Remaining in adj order
      for (const auto& nb : adj[center])
        if (nb.idx != other &&
            std::find(tv.begin(), tv.end(), nb.idx) == tv.end())
          tv.push_back(nb.idx);
      // Non-ring reordering: H-only and connectivity swap
      if (tv.empty() || has_ring_sharing) goto find_target;
      if (tv.size() == 2 && side1_nb_count == 2 && side2_nb_count == 2) {
        auto is_h_only = [&](size_t idx) {
          for (const auto& nb : adj[idx])
            if (nb.idx != center &&
                !cc.atoms[nb.idx].el.is_metal() &&
                !cc.atoms[nb.idx].is_hydrogen())
              return false;
          return true;
        };
        bool h0 = is_h_only(tv[0]);
        bool h1 = is_h_only(tv[1]);
        if (is_side1) {
          if (h0 && !h1)
            std::swap(tv[0], tv[1]);
          else if (conn_degree(tv[0]) < conn_degree(tv[1]))
            std::swap(tv[0], tv[1]);
        } else {
          if (h0 && (!h1 || conn_degree(tv[0]) < conn_degree(tv[1])))
            std::swap(tv[0], tv[1]);
          else if (!h0 && conn_degree(tv[0]) < conn_degree(tv[1]) && !h1)
            std::swap(tv[0], tv[1]);
        }
      }
      find_target:
      if (dbg_sp2) {
        std::fprintf(stderr, "[sp2 %s] %s-%s side%d target=%s tv=",
                     cc.name.c_str(), cc.atoms[center].id.c_str(),
                     cc.atoms[other].id.c_str(), is_side1 ? 1 : 2,
                     cc.atoms[target].id.c_str());
        for (size_t ti = 0; ti < tv.size(); ++ti) {
          int deg = (int)adj[tv[ti]].size();
          std::fprintf(stderr, "%s%s(deg=%d)%s", ti ? "," : "",
                       cc.atoms[tv[ti]].id.c_str(), deg,
                       cc.atoms[tv[ti]].is_hydrogen() ? "H" : "");
        }
        std::fprintf(stderr, " rs=%s/%s ring=%d nbs=%d/%d\n",
                     rs1 == SIZE_MAX ? "-" : cc.atoms[rs1].id.c_str(),
                     rs2 == SIZE_MAX ? "-" : cc.atoms[rs2].id.c_str(),
                     has_ring_sharing ? 1 : 0, side1_nb_count, side2_nb_count);
      }
      for (int i = 0; i < (int)tv.size(); ++i)
        if (tv[i] == target)
          return i;
      return -1;
    };
    size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
    size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
    const auto& m22 = has_ring_sharing ? ring_m : noring_m;
    int i_pos = sp2sp2_tv_pos(side1, side2, term1, true);
    int j_pos = sp2sp2_tv_pos(side2, side1, term2, false);
    if (dbg_sp2) {
      std::fprintf(stderr, "[sp2 %s] pick %s-%s-%s-%s i=%d j=%d m=%s\n",
                   cc.name.c_str(),
                   cc.atoms[term1].id.c_str(), cc.atoms[side1].id.c_str(),
                   cc.atoms[side2].id.c_str(), cc.atoms[term2].id.c_str(),
                   i_pos, j_pos, has_ring_sharing ? "ring" : "noring");
    }
    if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 2)
      value = m22[i_pos][j_pos];
    bool end1_ring = share_ring_ids(atom_info[term1].in_rings, atom_info[side1].in_rings);
    bool end2_ring = share_ring_ids(atom_info[term2].in_rings, atom_info[side2].in_rings);
    esd = (end1_ring || end2_ring) ? 20.0 : 5.0;
    period = 2;
  }

  if (bond_aromatic && ring_size > 0 && info2.is_aromatic && info3.is_aromatic &&
      sp2_2 && sp2_3) {
    esd = 0.0;
    period = 1;
  }

  Restraints::Torsion tor{
      "auto",
      {1, a1->id},
      {1, cc.atoms[center2].id},
      {1, cc.atoms[center3].id},
      {1, a4->id},
      value, esd, period};
  if (out_torsions)
    out_torsions->push_back(std::move(tor));
  else
    cc.rt.torsions.push_back(std::move(tor));
}

const Restraints::Torsion* find_generated_torsion(
    const std::vector<Restraints::Torsion>& torsions,
    const std::string& a1, const std::string& a2,
    const std::string& a3, const std::string& a4) {
  for (const auto& tor : torsions)
    if ((tor.id1.atom == a1 && tor.id2.atom == a2 &&
         tor.id3.atom == a3 && tor.id4.atom == a4) ||
        (tor.id1.atom == a4 && tor.id2.atom == a3 &&
         tor.id3.atom == a2 && tor.id4.atom == a1))
      return &tor;
  return nullptr;
}

// AceDRG's Python layer (confirmAAandNames) checks if a compound has standard
// AA backbone topology: atoms N, CA, C, O, OXT with correct types and bonding.
// If yes, it passes -C flag to the C++ binary, enabling setPeptideTorsions().
bool confirm_aa_backbone(const ChemComp& cc,
                         const AceBondAdjacency& adj,
                         const std::map<std::string, size_t>& atom_index) {
  auto find = [&](const std::string& name) -> size_t {
    auto it = atom_index.find(name);
    return it != atom_index.end() ? it->second : SIZE_MAX;
  };
  size_t i_n = find("N"), i_ca = find("CA"), i_c = find("C"),
         i_o = find("O"), i_oxt = find("OXT");
  if (i_n == SIZE_MAX || i_ca == SIZE_MAX || i_c == SIZE_MAX ||
      i_o == SIZE_MAX || i_oxt == SIZE_MAX)
    return false;
  if (cc.atoms[i_n].el != El::N || cc.atoms[i_ca].el != El::C ||
      cc.atoms[i_c].el != El::C || cc.atoms[i_o].el != El::O ||
      cc.atoms[i_oxt].el != El::O)
    return false;
  // CA must have exactly 4 bonds, bonded to N and C, exactly one H named HA
  if (adj[i_ca].size() != 4)
    return false;
  bool ca_has_n = false, ca_has_c = false;
  int ca_h_count = 0;
  bool ca_has_ha = false;
  for (const auto& nb : adj[i_ca]) {
    if (nb.idx == i_n) ca_has_n = true;
    else if (nb.idx == i_c) ca_has_c = true;
    if (cc.atoms[nb.idx].is_hydrogen()) {
      ca_h_count++;
      if (cc.atoms[nb.idx].id == "HA")
        ca_has_ha = true;
    }
  }
  if (!ca_has_n || !ca_has_c || ca_h_count != 1 || !ca_has_ha)
    return false;
  // C must have exactly 3 bonds: CA, O, OXT
  if (adj[i_c].size() != 3)
    return false;
  bool c_has_o = false, c_has_oxt = false;
  for (const auto& nb : adj[i_c]) {
    if (nb.idx == i_o) c_has_o = true;
    else if (nb.idx == i_oxt) c_has_oxt = true;
  }
  if (!c_has_o || !c_has_oxt)
    return false;
  // N must have 3-4 bonds with proper H naming
  size_t n_bonds = adj[i_n].size();
  if (n_bonds < 3 || n_bonds > 4)
    return false;
  std::vector<std::string> n_h_names;
  for (const auto& nb : adj[i_n])
    if (cc.atoms[nb.idx].is_hydrogen())
      n_h_names.push_back(cc.atoms[nb.idx].id);
  if (n_bonds == 3) {
    if (n_h_names.empty() || n_h_names.size() > 3)
      return false;
    bool has_h2 = std::find(n_h_names.begin(), n_h_names.end(), "H2") != n_h_names.end();
    if (!has_h2)
      return false;
    if (n_h_names.size() == 1) {
      if (n_h_names[0] != "H2")
        return false;
    } else {
      bool has_h = std::find(n_h_names.begin(), n_h_names.end(), "H") != n_h_names.end();
      if (!has_h || !has_h2)
        return false;
    }
  } else {  // n_bonds == 4
    if (n_h_names.size() != 3)
      return false;
    bool has_h = std::find(n_h_names.begin(), n_h_names.end(), "H") != n_h_names.end();
    bool has_h2 = std::find(n_h_names.begin(), n_h_names.end(), "H2") != n_h_names.end();
    if (!has_h || !has_h2)
      return false;
    std::string third_h;
    for (const std::string& hname : n_h_names)
      if (hname != "H" && hname != "H2") {
        third_h = hname;
        break;
      }
    if (third_h.empty())
      return false;
    if (third_h != "H3" &&
        !(third_h.size() > 1 && third_h[0] == 'H' &&
          std::isdigit(static_cast<unsigned char>(third_h[1]))))
      return false;
  }
  return true;
}

void apply_peptide_tmpchi2_override(
    Restraints::Torsion& tor,
    const std::map<std::string, size_t>& atom_index,
    const std::vector<CodAtomInfo>& atom_info) {
  std::string key = cat(tor.id1.atom, '_', tor.id2.atom, '_', tor.id3.atom, '_',
                        tor.id4.atom);
  if (key != "CA_CB_CG_CD" && key != "CA_CB_CG_CD1" && key != "CA_CB_CG_CD2")
    return;
  auto it2 = atom_index.find(tor.id2.atom);
  auto it3 = atom_index.find(tor.id3.atom);
  if (it2 == atom_index.end() || it3 == atom_index.end())
    return;
  int b2 = atom_info[it2->second].bonding_idx;
  int b3 = atom_info[it3->second].bonding_idx;
  if ((b2 == 2 && b3 == 3) || (b2 == 3 && b3 == 2)) {
    tor.period = 6;
    tor.value = 90.0;
  }
}

void set_peptide_torsion_idx_from_one_bond_like_acedrg(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const AcedrgTables& tables,
    const std::map<std::string, size_t>& atom_index,
    const std::vector<CodAtomInfo>& atom_info,
    size_t idx1, size_t idx2,
    std::vector<Restraints::Torsion>& chi_tors,
    std::vector<Restraints::Torsion>& hh_tors,
    std::vector<Restraints::Torsion>& cst_tors) {
  int best_priority = 10;
  bool found = false;
  Restraints::Torsion best;
  for (const auto& nb1 : adj[idx1]) {
    size_t a1 = nb1.idx;
    for (const auto& nb2 : adj[idx2]) {
      size_t a4 = nb2.idx;
      if (a1 == idx2 || a4 == idx1)
        continue;
      TorsionEntry entry;
      bool forward = true;
      if (!tables.lookup_pep_tors(cc.atoms[a1].id, cc.atoms[idx1].id,
                                  cc.atoms[idx2].id, cc.atoms[a4].id, entry)) {
        if (!tables.lookup_pep_tors(cc.atoms[a4].id, cc.atoms[idx2].id,
                                    cc.atoms[idx1].id, cc.atoms[a1].id, entry))
          continue;
        forward = false;
      }
      Restraints::Torsion tor;
      if (forward) {
        tor.id1 = {1, cc.atoms[a1].id};
        tor.id2 = {1, cc.atoms[idx1].id};
        tor.id3 = {1, cc.atoms[idx2].id};
        tor.id4 = {1, cc.atoms[a4].id};
      } else {
        tor.id1 = {1, cc.atoms[a4].id};
        tor.id2 = {1, cc.atoms[idx2].id};
        tor.id3 = {1, cc.atoms[idx1].id};
        tor.id4 = {1, cc.atoms[a1].id};
      }
      tor.label = entry.id;
      tor.value = entry.value;
      tor.period = entry.period;
      apply_peptide_tmpchi2_override(tor, atom_index, atom_info);
      if (entry.priority < best_priority) {
        best_priority = entry.priority;
        best = tor;
        found = true;
      }
    }
  }
  if (!found)
    return;
  if (best.label.find("chi") != std::string::npos) {
    chi_tors.push_back(std::move(best));
  } else if (best.label.find("hh") != std::string::npos) {
    hh_tors.push_back(std::move(best));
  } else if (best.label.find("CONST") != std::string::npos) {
    best.esd = 0.0;
    cst_tors.push_back(std::move(best));
  }
}

std::vector<Restraints::Torsion> set_peptide_torsions_like_acedrg(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const AcedrgTables& tables,
    const std::map<std::string, size_t>& atom_index,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<Restraints::Torsion>& all_torsions) {
  std::vector<Restraints::Torsion> chi_tors, hh_tors, cst_tors;
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t b1 = it1->second;
    size_t b2 = it2->second;
    std::vector<size_t> pos;
    if (adj[b1].size() > 1)
      pos.push_back(b1);
    if (adj[b2].size() > 1)
      pos.push_back(b2);
    if (pos.size() != 2)
      continue;
    if (cc.atoms[pos[1]].id < cc.atoms[pos[0]].id)
      std::swap(pos[0], pos[1]);
    set_peptide_torsion_idx_from_one_bond_like_acedrg(
        cc, adj, tables, atom_index, atom_info,
        pos[0], pos[1], chi_tors, hh_tors, cst_tors);
  }

  std::vector<Restraints::Torsion> mini_torsions;
  std::map<std::string, size_t> map_id_ser;
  std::vector<std::string> done_bonds;
  auto torsion_id = [](const Restraints::Torsion& t) {
    return cat(t.id1.atom, '_', t.id2.atom, '_', t.id3.atom, '_', t.id4.atom);
  };
  auto bond_id = [](const Restraints::Torsion& t) {
    return cat(t.id2.atom, '_', t.id3.atom);
  };

  if (!chi_tors.empty()) {
    std::map<std::string, size_t> sorted_ids;
    for (size_t i = 0; i < chi_tors.size(); ++i)
      sorted_ids[chi_tors[i].label] = i;
    for (const auto& kv : sorted_ids) {
      const auto& tor = chi_tors[kv.second];
      map_id_ser[torsion_id(tor)] = mini_torsions.size();
      done_bonds.push_back(bond_id(tor));
      mini_torsions.push_back(tor);
    }
  }

  for (size_t i = 0; i < cst_tors.size(); ++i) {
    Restraints::Torsion tor = cst_tors[i];
    tor.label = cat("CONST_", i + 1);
    tor.esd = 0.0;
    map_id_ser[torsion_id(tor)] = mini_torsions.size();
    done_bonds.push_back(bond_id(tor));
    mini_torsions.push_back(std::move(tor));
  }

  for (size_t i = 0; i < hh_tors.size(); ++i) {
    Restraints::Torsion tor = hh_tors[i];
    tor.label = cat(tor.label, i + 1);
    map_id_ser[torsion_id(tor)] = mini_torsions.size();
    done_bonds.push_back(bond_id(tor));
    mini_torsions.push_back(std::move(tor));
  }

  std::vector<Restraints::Torsion> tmp_tors;
  for (const auto& tor : all_torsions) {
    std::string bid12 = cat(tor.id2.atom, '_', tor.id3.atom);
    std::string bid21 = cat(tor.id3.atom, '_', tor.id2.atom);
    std::string tid12 = torsion_id(tor);
    std::string tid21 = cat(tor.id4.atom, '_', tor.id3.atom, '_', tor.id2.atom, '_',
                            tor.id1.atom);
    auto it12 = map_id_ser.find(tid12);
    if (it12 != map_id_ser.end()) {
      mini_torsions[it12->second].period = tor.period;
      continue;
    }
    auto it21 = map_id_ser.find(tid21);
    if (it21 != map_id_ser.end()) {
      mini_torsions[it21->second].period = tor.period;
      continue;
    }
    if (std::find(done_bonds.begin(), done_bonds.end(), bid12) == done_bonds.end() &&
        std::find(done_bonds.begin(), done_bonds.end(), bid21) == done_bonds.end())
      tmp_tors.push_back(tor);
  }

  std::map<std::string, std::vector<Restraints::Torsion>> map_tors2;
  std::map<std::string, int> map_tor_hi;
  for (const auto& tor : tmp_tors) {
    auto it2 = atom_index.find(tor.id2.atom);
    auto it3 = atom_index.find(tor.id3.atom);
    if (it2 == atom_index.end() || it3 == atom_index.end())
      continue;
    int z2 = Element(cc.atoms[it2->second].el).atomic_number();
    int z3 = Element(cc.atoms[it3->second].el).atomic_number();
    if (z2 <= 0 || z3 <= 0)
      continue;
    int av = z2 + z3;
    std::string bid = cat(tor.id2.atom, '_', tor.id3.atom);
    auto h = map_tor_hi.find(bid);
    if (h == map_tor_hi.end()) {
      map_tor_hi[bid] = av;
      map_tors2[bid].push_back(tor);
    } else if (av > h->second) {
      map_tors2[bid].clear();
      map_tor_hi[bid] = av;
      map_tors2[bid].push_back(tor);
    }
  }
  for (const auto& kv : map_tors2)
    for (const auto& tor : kv.second)
      mini_torsions.push_back(tor);
  return mini_torsions;
}

const Restraints::Torsion* select_one_torsion_from_candidates(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center2, size_t center3,
    const std::vector<Restraints::Torsion>& candidates,
    bool* used_path1 = nullptr) {
  bool dbg_sel = false;
  if (const char* env = std::getenv("GEMMI_DBG_SEL")) {
    std::string v(env);
    dbg_sel = v == "1" || v.find(cc.name) != std::string::npos;
  }
  if (dbg_sel) {
    std::fprintf(stderr, "[tor-sel %s] bond %s-%s cand=%zu\n", cc.name.c_str(),
                 cc.atoms[center2].id.c_str(), cc.atoms[center3].id.c_str(), candidates.size());
  }
  if (used_path1)
    *used_path1 = false;
  if (candidates.empty())
    return nullptr;
  std::vector<size_t> idxR1, idxNonH1, idxH1, idxR2, idxNonH2, idxH2;
  for (const auto& nb : adj[center2]) {
    if (nb.idx == center3 || cc.atoms[nb.idx].el.is_metal())
      continue;
    if (!atom_info[nb.idx].in_rings.empty())
      idxR1.push_back(nb.idx);
    else if (!cc.atoms[nb.idx].is_hydrogen())
      idxNonH1.push_back(nb.idx);
    else
      idxH1.push_back(nb.idx);
  }
  for (const auto& nb : adj[center3]) {
    if (nb.idx == center2 || cc.atoms[nb.idx].el.is_metal())
      continue;
    if (!atom_info[nb.idx].in_rings.empty())
      idxR2.push_back(nb.idx);
    else if (!cc.atoms[nb.idx].is_hydrogen())
      idxNonH2.push_back(nb.idx);
    else
      idxH2.push_back(nb.idx);
  }
  const std::string& a2 = cc.atoms[center2].id;
  const std::string& a3 = cc.atoms[center3].id;
  auto pick = [&](size_t a1, size_t a4) -> const Restraints::Torsion* {
    if (dbg_sel)
      std::fprintf(stderr, "  try %s-%s-%s-%s\n", cc.atoms[a1].id.c_str(), a2.c_str(),
                   a3.c_str(), cc.atoms[a4].id.c_str());
    const Restraints::Torsion* t =
        find_generated_torsion(candidates, cc.atoms[a1].id, a2, a3, cc.atoms[a4].id);
    if (t && t->id1.atom != t->id4.atom) {
      if (dbg_sel)
        std::fprintf(stderr, "    -> pick %s-%s-%s-%s\n", t->id1.atom.c_str(),
                     t->id2.atom.c_str(), t->id3.atom.c_str(), t->id4.atom.c_str());
      return t;
    }
    if (dbg_sel && t)
      std::fprintf(stderr, "    -> reject self %s-%s-%s-%s\n", t->id1.atom.c_str(),
                   t->id2.atom.c_str(), t->id3.atom.c_str(), t->id4.atom.c_str());
    return nullptr;
  };
  if (!idxNonH1.empty() && !idxNonH2.empty()) {
    if (used_path1)
      *used_path1 = true;
    if (auto t = pick(idxNonH1[0], idxNonH2[0])) return t;
    return nullptr;
  } else if (!idxNonH1.empty() && idxNonH2.empty()) {
    if (!idxR2.empty())
      if (auto t = pick(idxNonH1[0], idxR2[0])) return t;
    if (!idxH2.empty())
      if (auto t = pick(idxNonH1[0], idxH2[0])) return t;
    return nullptr;
  } else if (idxNonH1.empty() && !idxNonH2.empty()) {
    if (!idxR1.empty())
      if (auto t = pick(idxR1[0], idxNonH2[0])) return t;
    if (!idxH1.empty())
      if (auto t = pick(idxH1[0], idxNonH2[0])) return t;
    return nullptr;
  } else if (!idxR1.empty() && !idxR2.empty()) {
    // AceDRG: if ring+ring fails (e.g. atoms[0]==atoms[3] in 3-ring),
    // it does not fall through to ring+H due to its else-if structure.
    if (auto t = pick(idxR1[0], idxR2[0])) return t;
    return nullptr;
  } else if (!idxR1.empty() && idxR2.empty()) {
    if (!idxH2.empty())
      if (auto t = pick(idxR1[0], idxH2[0])) return t;
    return nullptr;
  } else if (idxR1.empty() && !idxR2.empty()) {
    if (!idxH1.empty())
      if (auto t = pick(idxH1[0], idxR2[0])) return t;
    return nullptr;
  } else if (!idxH1.empty() && !idxH2.empty()) {
    if (auto t = pick(idxH1[0], idxH2[0])) return t;
    return nullptr;
  }
  // Literal torsion.cpp fallback: both centers connect only to H.
  // Pick first non-bond connection from each side regardless of type.
  if (adj[center2].size() > 1 && adj[center3].size() > 1) {
    size_t c1 = SIZE_MAX, c2 = SIZE_MAX;
    for (const auto& nb : adj[center2]) {
      if (nb.idx != center3) {
        c1 = nb.idx;
        break;
      }
    }
    for (const auto& nb : adj[center3]) {
      if (nb.idx != center2) {
        c2 = nb.idx;
        break;
      }
    }
    if (c1 != SIZE_MAX && c2 != SIZE_MAX)
      if (auto t = pick(c1, c2)) return t;
  }
  return nullptr;
}

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables,
                                        const std::vector<CodAtomInfo>& atom_info,
                                        const std::map<std::string, std::string>& atom_stereo,
                                        const AceGraphView& graph,
                                        const std::map<std::string, Position>* sugar_coord_overrides) {
  if (std::getenv("GEMMI_DBG_SEL"))
    std::fprintf(stderr, "[tor-start %s] existing=%zu\n", cc.name.c_str(), cc.rt.torsions.size());
  if (!cc.rt.torsions.empty())
    return;

  auto& atom_index = graph.atom_index;
  auto& adj = graph.adjacency;
  std::string type_upper = to_upper(cc.type_or_group);
  bool peptide_mode = type_upper.find("PEPTIDE") != std::string::npos;
  bool nucleic_mode = (type_upper.find("DNA") != std::string::npos ||
                       type_upper.find("RNA") != std::string::npos);
  // AceDRG applies pepCorr/naCorr only when a descriptor loop is present.
  if (!cc.has_descriptor) {
    peptide_mode = false;
    nucleic_mode = false;
  }
  // AceDRG's Python layer (confirmAAandNames) verifies standard backbone
  // topology before passing -C flag to enable setPeptideTorsions().
  if (peptide_mode && !confirm_aa_backbone(cc, adj, atom_index))
    peptide_mode = false;
  const ResidueInfo& ri = find_tabulated_residue(cc.name);
  bool standard_aa = ri.is_standard() && ri.kind == ResidueKind::AA;
  bool use_peptide_torsions = standard_aa || peptide_mode;
  std::vector<bool> aromatic_like = build_aromatic_like_mask(cc, atom_info, atom_index);

  // Pre-compute ring parity used in SP3-SP3 torsion matrix selection.
  std::map<std::pair<size_t, size_t>, RingParity> bond_ring_parity =
      build_ring_bond_parity(adj, atom_info);
  std::map<std::pair<size_t, size_t>, RingFlip> bond_ring_flip =
      build_ring_bond_flip(cc, adj, atom_info);

  // AceDRG has a dedicated sugar-ring mode: ring bonds are represented by
  // one nu torsion each (from ring geometry), while non-ring bonds keep the
  // full torsion set.
  // Match AceDRG checkOneRingSugar()/getRStr() shape gating:
  // only specific "(OC2)(...)" ring-shape strings are treated as sugar.
  SugarRingInfo sugar_info = detect_sugar_rings(cc, adj, atom_info);
  auto& sugar_ring_bonds = sugar_info.ring_bonds;
  auto& sugar_ring_seq = sugar_info.ring_seq;
  // AceDRG: isPeptide takes priority over sugar detection (resetSystem2).
  bool has_sugar_ring = !peptide_mode && !sugar_ring_bonds.empty();
  if (std::getenv("GEMMI_DBG_SEL"))
    std::fprintf(stderr, "[tor-start %s] sugar=%d rings=%zu\n", cc.name.c_str(),
                 has_sugar_ring ? 1 : 0, sugar_ring_bonds.size());

  ChiralCenterInfo chiral_info =
      detect_chiral_centers_and_mut_table(cc, adj, atom_info, atom_stereo);
  auto& stereo_chiral_centers = chiral_info.stereo_chiral_centers;
  auto& chir_mut_table = chiral_info.chir_mut_table;

  // TEMPORARY: detect 3-membered rings for AceDRG off-by-one bug simulation.
  // AceDRG's getTorsion() returns seriNum (original global index) but callers
  // use it as a vector index.  After self-referencing torsions (from 3-rings)
  // are removed, seriNum != index, causing the wrong torsion to be picked.
  bool has_3ring = false;
  for (size_t ii = 0; ii < atom_info.size(); ++ii)
    if (atom_info[ii].min_ring_size == 3) { has_3ring = true; break; }

  // Per-bond tracking for the 3-ring bug simulation.
  struct BondBugInfo {
    std::vector<Restraints::Torsion> generated;
    std::vector<size_t> atv1, atv2;    // exact per-bond TV lists used for emission
    bool swap_term_emit;               // whether emit loop order was swapped
    size_t rt_idx;                     // index in cc.rt.torsions (SIZE_MAX if none)
    size_t sel_idx1, sel_idx2;         // min/max numeric indices (AceDRG cascade sides)
  };
  std::vector<BondBugInfo> bug_infos;
  std::vector<Restraints::Torsion> peptide_all_torsions;

  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;
    if (cc.atoms[idx1].el.is_metal() || cc.atoms[idx2].el.is_metal())
      continue;

    // Match AceDRG map-ordered bond traversal (fullAtoms map key order).
    size_t center2 = cc.atoms[idx1].id < cc.atoms[idx2].id ? idx1 : idx2;
    size_t center3 = center2 == idx1 ? idx2 : idx1;
    // AceDRG's setupMiniTorsions sorts center indices numerically (min_max),
    // so selectOneTorFromOneBond's idx1=min, idx2=max.
    size_t sel_center2 = std::min(idx1, idx2);
    size_t sel_center3 = std::max(idx1, idx2);
    // Match AceDRG setTorsionFromOneBond() dispatch orientation.
    bool c2_sp2 = is_sp2_like(atom_info[center2]);
    bool c3_sp2 = is_sp2_like(atom_info[center3]);
    bool c2_sp3 = is_sp3_like(atom_info[center2]);
    bool c3_sp3 = is_sp3_like(atom_info[center3]);
    bool c2_oxy_sp2 = c2_sp2 && is_oxygen_column(cc.atoms[center2].el);
    bool c3_oxy_sp2 = c3_sp2 && is_oxygen_column(cc.atoms[center3].el);
    if ((c2_oxy_sp2 && c3_sp3) || (c3_oxy_sp2 && c2_sp3)) {
      // SetOneSP3OxyColumnBond(tIdx1=sp3, tIdx2=sp2-oxy)
      if (c2_oxy_sp2 && c3_sp3)
        std::swap(center2, center3);
    } else if (c2_sp3 && c3_sp2) {
      // SetOneSP2SP3Bond(tIdx1=sp2, tIdx2=sp3)
      std::swap(center2, center3);
    }
    c2_sp2 = is_sp2_like(atom_info[center2]);
    c3_sp2 = is_sp2_like(atom_info[center3]);
    c2_sp3 = is_sp3_like(atom_info[center2]);
    c3_sp3 = is_sp3_like(atom_info[center3]);

    int ring_size = shared_ring_size_from_ring_ids(atom_info[center2].in_rings,
                                                   atom_info[center2].min_ring_size,
                                                   atom_info[center3].in_rings,
                                                   atom_info[center3].min_ring_size);

    bool bond_aromatic = (bond.type == BondType::Aromatic ||
                          bond.type == BondType::Deloc || bond.aromatic ||
                          (aromatic_like[idx1] && aromatic_like[idx2]));

    // AceDRG generates SP1 torsions internally but setupMiniTorsions() filters
    // them all out (buggy check on atoms[1].bondingIdx effectively removes all
    // SP1-centered torsions since SetOneSP1* always puts SP1 as atoms[1]).
    // Skip SP1 bonds entirely to match AceDRG output.
    int b2 = atom_info[center2].bonding_idx;
    int b3 = atom_info[center3].bonding_idx;
    if (b2 <= 1 || b3 <= 1)
      continue;

    std::vector<size_t> tv1_idx, tv2_idx;
    TvMode mode12 = TvMode::Default;
    TvMode mode21 = TvMode::Default;
    bool sp3_like_12 = is_sp3_like(atom_info[center2]);
    bool sp3_like_21 = is_sp3_like(atom_info[center3]);
    bool sp2_like_12 = is_sp2_like(atom_info[center2]);
    bool sp2_like_21 = is_sp2_like(atom_info[center3]);
    if (sp2_like_12 && sp2_like_21) {
      auto tvs = build_tv_lists_sp2sp2_like_acedrg(cc, adj, atom_info, center2, center3);
      tv1_idx = std::move(tvs.first);
      tv2_idx = std::move(tvs.second);
    } else {
      if (sp3_like_12 && sp3_like_21) {
        mode12 = TvMode::SP3SP3;
        mode21 = TvMode::SP3SP3;
      } else if (sp3_like_12 && sp2_like_21) {
        mode12 = TvMode::SP2SP3_SP3;
      } else if (sp2_like_12 && sp3_like_21) {
        mode21 = TvMode::SP2SP3_SP3;
      }
      tv1_idx = build_tv_list_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          center2, center3, mode12);
      tv2_idx = build_tv_list_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          center3, center2, mode21);
    }
    ring_size = shared_ring_size_from_ring_ids(atom_info[center2].in_rings,
                                               atom_info[center2].min_ring_size,
                                               atom_info[center3].in_rings,
                                               atom_info[center3].min_ring_size);
    if (tv1_idx.empty() || tv2_idx.empty())
      continue;



    std::vector<Restraints::Torsion> generated;
    generated.reserve(tv1_idx.size() * tv2_idx.size());
    bool oxy_c2_sp2 = c2_sp2 && is_oxygen_column(cc.atoms[center2].el);
    bool oxy_c3_sp2 = c3_sp2 && is_oxygen_column(cc.atoms[center3].el);
    bool plain_sp2sp3 = ((c2_sp2 && c3_sp3) || (c2_sp3 && c3_sp2)) &&
                        !oxy_c2_sp2 && !oxy_c3_sp2;
    bool swap_term_emit = plain_sp2sp3 && tv1_idx.size() > tv2_idx.size();
    if (swap_term_emit) {
      for (size_t a1_idx : tv2_idx)
        for (size_t a4_idx : tv1_idx)
          emit_one_torsion(cc, adj, atom_info, tables, stereo_chiral_centers,
                           chir_mut_table, bond_ring_parity, bond_ring_flip, peptide_mode,
                           center2, center3, bond_aromatic, ring_size,
                           a1_idx, a4_idx, &generated);
    } else {
      for (size_t a1_idx : tv1_idx)
        for (size_t a4_idx : tv2_idx)
          emit_one_torsion(cc, adj, atom_info, tables, stereo_chiral_centers,
                           chir_mut_table, bond_ring_parity, bond_ring_flip, peptide_mode,
                           center2, center3, bond_aromatic, ring_size,
                           a1_idx, a4_idx, &generated);
    }
    if (generated.empty())
      continue;

    if (const char* env = std::getenv("GEMMI_DBG_SEL")) {
      std::string v(env);
      if (v == "1" || v.find(cc.name) != std::string::npos) {
        std::fprintf(stderr, "[tor-gen %s] bond %s-%s generated=%zu\n",
                     cc.name.c_str(), cc.atoms[center2].id.c_str(),
                     cc.atoms[center3].id.c_str(), generated.size());
        for (const auto& t : generated)
          std::fprintf(stderr, "  cand %s-%s-%s-%s val=%.2f esd=%.2f p.%d\n",
                       t.id1.atom.c_str(), t.id2.atom.c_str(),
                       t.id3.atom.c_str(), t.id4.atom.c_str(),
                       t.value, t.esd, t.period);
      }
    }

    if (has_sugar_ring) {
      if (has_3ring)
        bug_infos.push_back(BondBugInfo{generated, tv1_idx, tv2_idx,
                                        swap_term_emit, SIZE_MAX, sel_center2, sel_center3});
      cc.rt.torsions.insert(cc.rt.torsions.end(), generated.begin(), generated.end());
    } else if (use_peptide_torsions) {
      peptide_all_torsions.insert(peptide_all_torsions.end(),
                                  generated.begin(), generated.end());
    } else {
      const Restraints::Torsion* chosen = select_one_torsion_from_candidates(
          cc, adj, atom_info, sel_center2, sel_center3, generated);
      size_t rt_idx = SIZE_MAX;
      if (chosen) {
        rt_idx = cc.rt.torsions.size();
        cc.rt.torsions.push_back(*chosen);
      }
      if (has_3ring)
        bug_infos.push_back(BondBugInfo{generated, tv1_idx, tv2_idx,
                                        swap_term_emit, rt_idx, sel_center2, sel_center3});
    }
  }

  if (use_peptide_torsions) {
    cc.rt.torsions = set_peptide_torsions_like_acedrg(
        cc, adj, tables, atom_index, atom_info, peptide_all_torsions);
  }

  // TEMPORARY: apply AceDRG off-by-one bug for 3-ring compounds.
  // AceDRG's getTorsion() returns seriNum (original global index) but callers
  // use it as a vector index.  After self-referencing torsions (from 3-rings)
  // are removed, seriNum != index, causing the wrong torsion to be picked.
  // ALL cascade paths use the same buggy getTorsion.
  //
  // Key: AceDRG classifies atoms as Ring if inRings is non-empty (any ring size).
  // The cascade uses full connAtoms lists (not limited TV lists).
  // idx1 side uses chemType != "H", idx2 side uses chemType.find("H")==npos
  // (asymmetric check — AceDRG bug).
  if (has_3ring && !bug_infos.empty()) {
    bool dbg_tor = false;
    if (const char* env = std::getenv("GEMMI_DBG_TOR")) {
      std::string v(env);
      dbg_tor = v == "1" || v.find(cc.name) != std::string::npos;
    }
    // Step 1: Sort bug_infos by AceDRG's TorsionSetOneBond map key order.
    // Key = IntToStr(min_idx) + "_" + IntToStr(max_idx), string-lexicographic.
    auto make_sel_key = [](size_t idx1, size_t idx2) -> std::string {
      return std::to_string(idx1) + "_" + std::to_string(idx2);
    };
    std::vector<size_t> sorted_order(bug_infos.size());
    std::iota(sorted_order.begin(), sorted_order.end(), 0);
    std::sort(sorted_order.begin(), sorted_order.end(), [&](size_t a, size_t b) {
      return make_sel_key(bug_infos[a].sel_idx1, bug_infos[a].sel_idx2) <
             make_sel_key(bug_infos[b].sel_idx1, bug_infos[b].sel_idx2);
    });
    if (dbg_tor) {
      std::fprintf(stderr, "[tor-bug %s] sorted keys:\n", cc.name.c_str());
      for (size_t bi : sorted_order)
        std::fprintf(stderr, "  key=%s rt=%zu bi=%zu\n",
                     make_sel_key(bug_infos[bi].sel_idx1, bug_infos[bi].sel_idx2).c_str(),
                     bug_infos[bi].rt_idx, bi);
    }

    // Step 2: Build AceDRG cross product in CORRECT Phase 1 / Phase 2 order.
    // AceDRG's setAllTorsions2():
    //   Phase 1: ring-walk bonds (N-1 per ring, allRingsV string-sort order).
    //     Uses SetOneSP3SP3Bond(flip) which has "if (tV1[i]!=tV2[j])" check
    //     → phantom entries (a1==a4) are SKIPPED.
    //   Phase 2: remaining bonds in CIF/allBonds order.
    //     Uses SetOneSP3SP3Bond (no flip) → phantoms ARE included.
    //     Only 3-ring closing bonds produce phantoms (the shared third vertex
    //     appears in both TV lists, giving a1==a4 at position 0).
    //
    // Build ring walks to classify Phase 1 vs Phase 2 closing bonds.
    int max_rid = 0;
    for (size_t ii = 0; ii < atom_info.size(); ++ii)
      for (int rid : atom_info[ii].in_rings)
        max_rid = std::max(max_rid, rid + 1);
    std::vector<std::vector<size_t>> rings_by_id(max_rid);
    for (size_t ii = 0; ii < atom_info.size(); ++ii)
      for (int rid : atom_info[ii].in_rings)
        rings_by_id[rid].push_back(ii);

    struct RingWalkBug {
      std::string rep;                               // allRingsV sort key
      std::vector<std::pair<size_t,size_t>> p1bonds; // Phase 1 walk bonds
      std::pair<size_t,size_t> closing;              // Phase 2 closing bond
      int size;
    };
    std::vector<RingWalkBug> rwalks;
    for (auto& r_atoms : rings_by_id) {
      if (r_atoms.size() < 3) continue;
      std::set<size_t> rset(r_atoms.begin(), r_atoms.end());
      size_t start = *std::min_element(r_atoms.begin(), r_atoms.end());
      std::vector<size_t> tr = {start};
      size_t cur = start;
      while (tr.size() < r_atoms.size()) {
        bool found = false;
        for (const auto& nb : adj[cur])
          if (rset.count(nb.idx) &&
              std::find(tr.begin(), tr.end(), nb.idx) == tr.end()) {
            tr.push_back(nb.idx); cur = nb.idx; found = true; break;
          }
        if (!found) break;
      }
      if (tr.size() != r_atoms.size()) continue;
      std::vector<std::string> names;
      for (size_t idx : tr) names.push_back(cc.atoms[idx].id);
      std::sort(names.begin(), names.end());
      std::string rep;
      for (const auto& n : names) rep += n;
      RingWalkBug rw;
      rw.rep = std::move(rep);
      rw.size = (int)tr.size();
      for (size_t i = 0; i + 1 < tr.size(); ++i)
        rw.p1bonds.push_back(std::minmax(tr[i], tr[i+1]));
      rw.closing = std::minmax(tr.back(), tr.front());
      rwalks.push_back(std::move(rw));
    }
    std::sort(rwalks.begin(), rwalks.end(),
              [](const RingWalkBug& a, const RingWalkBug& b){ return a.rep < b.rep; });

    // Collect Phase 1 SP3-SP3 bond set.
    // AceDRG Phase 1 processes ALL ring walk bonds (SP2-SP2, SP2-SP3, SP3-SP3).
    // Only SP3-SP3 walk bonds use SetOneSP3SP3Bond(flip) which has a skip check
    // for phantoms (a1==a4). SP2-SP3 walk bonds use SetOneSP2SP3Bond which has
    // NO skip check, so they generate phantoms.
    std::set<std::pair<size_t,size_t>> phase1_sp3sp3_set;
    for (const auto& rw : rwalks) {
      for (auto& pb : rw.p1bonds)
        if (is_sp3_like(atom_info[pb.first]) && is_sp3_like(atom_info[pb.second]))
          phase1_sp3sp3_set.insert(pb);
    }

    // Build bug_info key→index map.
    std::map<std::pair<size_t,size_t>, size_t> bkey_to_bi;
    for (size_t bi = 0; bi < bug_infos.size(); ++bi)
      bkey_to_bi[std::minmax(bug_infos[bi].sel_idx1, bug_infos[bi].sel_idx2)] = bi;

    // Ordering: Phase 1 (ALL ring walk) bonds first, then Phase 2 (CIF order).
    std::vector<size_t> ordered_bis;
    {
      std::set<size_t> ordered_set;
      for (const auto& rw : rwalks)
        for (auto& pb : rw.p1bonds) {
          auto it = bkey_to_bi.find(pb);
          if (it != bkey_to_bi.end() && !ordered_set.count(it->second)) {
            ordered_bis.push_back(it->second);
            ordered_set.insert(it->second);
          }
        }
      for (size_t bi = 0; bi < bug_infos.size(); ++bi)
        if (!ordered_set.count(bi)) {
          ordered_bis.push_back(bi);
          ordered_set.insert(bi);
        }
    }

    struct BondCPInfo {
      std::vector<std::pair<size_t, size_t>> entries;
      size_t global_start;
    };
    std::vector<BondCPInfo> bond_cps(bug_infos.size());
    size_t global_pos = 0;
    bool any_phantom = false;
    for (size_t bi : ordered_bis) {
      auto& info = bug_infos[bi];
      auto& bcp = bond_cps[bi];
      bcp.global_start = global_pos;
      auto bkey = std::minmax(info.sel_idx1, info.sel_idx2);
      const std::vector<size_t>* p_outer;
      const std::vector<size_t>* p_inner;
      p_outer = info.swap_term_emit ? &info.atv2 : &info.atv1;
      p_inner = info.swap_term_emit ? &info.atv1 : &info.atv2;
      for (size_t a1 : *p_outer)
        for (size_t a4 : *p_inner) {
          // Only SP3-SP3 ring walk bonds have the skip check in AceDRG
          // (SetOneSP3SP3Bond with flip). SP2-SP3 walk bonds (SetOneSP2SP3Bond)
          // do NOT skip, so they generate phantoms.  Closing bonds and other
          // Phase 2 bonds also do NOT skip.
          if (phase1_sp3sp3_set.count(bkey) > 0 && a1 == a4) continue;
          bcp.entries.push_back({a1, a4});
          if (a1 == a4) {
            any_phantom = true;
            if (dbg_tor)
              std::fprintf(stderr, "  phantom key=%zu_%zu at global=%zu atom=%s\n",
                           bkey.first, bkey.second, global_pos + bcp.entries.size() - 1,
                           cc.atoms[a1].id.c_str());
          }
        }
      global_pos += bcp.entries.size();
    }
    if (dbg_tor) {
      std::fprintf(stderr, "[tor-bug %s] cross order:\n", cc.name.c_str());
      for (size_t bi : ordered_bis) {
        auto key = std::minmax(bug_infos[bi].sel_idx1, bug_infos[bi].sel_idx2);
        std::fprintf(stderr, "  key=%zu_%zu start=%zu n=%zu bi=%zu\n",
                     key.first, key.second, bond_cps[bi].global_start,
                     bond_cps[bi].entries.size(), bi);
        for (size_t ei = 0; ei < bond_cps[bi].entries.size(); ++ei) {
          const auto& e = bond_cps[bi].entries[ei];
          std::fprintf(stderr, "    e%zu=%s-%s\n", ei,
                       cc.atoms[e.first].id.c_str(), cc.atoms[e.second].id.c_str());
        }
      }
    }

    if (any_phantom) {
      // Step 3: Build global phantom mask and new_to_orig mapping.
      std::vector<bool> is_phantom(global_pos, false);
      for (size_t bi = 0; bi < bond_cps.size(); ++bi)
        for (size_t i = 0; i < bond_cps[bi].entries.size(); ++i)
          if (bond_cps[bi].entries[i].first == bond_cps[bi].entries[i].second)
            is_phantom[bond_cps[bi].global_start + i] = true;
      std::vector<size_t> new_to_orig;
      new_to_orig.reserve(global_pos);
      for (size_t i = 0; i < global_pos; ++i)
        if (!is_phantom[i])
          new_to_orig.push_back(i);
      if (dbg_tor)
        std::fprintf(stderr, "[tor-bug %s] global=%zu nonphantom=%zu\n",
                     cc.name.c_str(), global_pos, new_to_orig.size());


      // Step 4: For each bond, simulate AceDRG's cascade to find the pair
      // it would pick, then apply the seriNum shift.
      // AceDRG removes phantoms at the END of selectOneTorFromOneBond.
      // The first bond processed has phantoms still present (seriNum == index).
      // Subsequent bonds see the shortened vector (seriNum != index).
      std::string min_sel_key = make_sel_key(
          bug_infos[sorted_order[0]].sel_idx1,
          bug_infos[sorted_order[0]].sel_idx2);
      bool first_bond_in_map = true;
      for (size_t si = 0; si < sorted_order.size(); ++si) {
        size_t bi = sorted_order[si];
        auto& info = bug_infos[bi];
        if (info.rt_idx == SIZE_MAX)
          continue;  // no torsion selected for this bond

        std::string this_key = make_sel_key(info.sel_idx1, info.sel_idx2);

        // Simulate AceDRG's cascade classification using FULL neighbor lists.
        // AceDRG: inRings.size()!=0 → Ring; not hydrogen → nonH; else → H.
        std::vector<size_t> nonH1, R1, H1, nonH2, R2, H2;
        for (const auto& nb : adj[info.sel_idx1]) {
          if (nb.idx == info.sel_idx2 || cc.atoms[nb.idx].el.is_metal()) continue;
          if (!atom_info[nb.idx].in_rings.empty())
            R1.push_back(nb.idx);
          else if (!cc.atoms[nb.idx].is_hydrogen())
            nonH1.push_back(nb.idx);
          else
            H1.push_back(nb.idx);
        }
        for (const auto& nb : adj[info.sel_idx2]) {
          if (nb.idx == info.sel_idx1 || cc.atoms[nb.idx].el.is_metal()) continue;
          if (!atom_info[nb.idx].in_rings.empty())
            R2.push_back(nb.idx);
          else if (!cc.atoms[nb.idx].is_hydrogen())
            nonH2.push_back(nb.idx);
          else
            H2.push_back(nb.idx);
        }

        // Follow AceDRG's exact cascade structure (if-else-if, NOT fall-through).
        // ALL paths use getTorsion (buggy seriNum return).
        size_t cas_a1 = SIZE_MAX, cas_a4 = SIZE_MAX;
        if (!nonH1.empty() && !nonH2.empty()) {
          // Path 1: nonH+nonH
          cas_a1 = nonH1[0]; cas_a4 = nonH2[0];
        } else if (!nonH1.empty() && nonH2.empty()) {
          // Paths 2-3: nonH1 present, nonH2 empty
          if (!R2.empty()) {
            cas_a1 = nonH1[0]; cas_a4 = R2[0];  // Path 2: nonH+R
          } else if (!H2.empty()) {
            cas_a1 = nonH1[0]; cas_a4 = H2[0];  // Path 3: nonH+H
          }
        } else if (nonH1.empty() && !nonH2.empty()) {
          // Paths 4-5: nonH2 present, nonH1 empty
          if (!R1.empty()) {
            cas_a1 = R1[0]; cas_a4 = nonH2[0];  // Path 4: R+nonH
          } else if (!H1.empty()) {
            cas_a1 = H1[0]; cas_a4 = nonH2[0];  // Path 5: H+nonH
          }
        } else if (!R1.empty() || !R2.empty()) {
          // Paths 6-8: R atoms present, no nonH on either side
          if (!R1.empty() && !R2.empty()) {
            cas_a1 = R1[0]; cas_a4 = R2[0];     // Path 6: R+R
          } else if (!R1.empty() && R2.empty()) {
            if (!H2.empty()) {
              cas_a1 = R1[0]; cas_a4 = H2[0];   // Path 7: R+H
            }
          } else if (!R2.empty() && R1.empty()) {
            if (!H1.empty()) {
              cas_a1 = H1[0]; cas_a4 = R2[0];   // Path 8: H+R
            }
          }
        } else {
          // Path 9: H+H fallback — pick first non-other-center neighbor
          for (const auto& nb : adj[info.sel_idx1]) {
            if (nb.idx != info.sel_idx2) { cas_a1 = nb.idx; break; }
          }
          for (const auto& nb : adj[info.sel_idx2]) {
            if (nb.idx != info.sel_idx1) { cas_a4 = nb.idx; break; }
          }
        }
        if (cas_a1 == SIZE_MAX || cas_a4 == SIZE_MAX)
          continue;

        // Check if this is the very first bond in the FULL map iteration.
        // The first bond's cascade runs with phantoms present (no shift).
        if (first_bond_in_map && this_key == min_sel_key) {
          first_bond_in_map = false;
          if (dbg_tor)
            std::fprintf(stderr, "  key=%s first map bond: no shift\n", this_key.c_str());
          continue;
        }

        // Find cascade pair in cross product (search both orderings since
        // generation order center2/center3 may differ from cascade sel_idx1/sel_idx2).
        auto& bcp = bond_cps[bi];
        if (dbg_tor) {
          std::fprintf(stderr, "  key=%s cas-pair=%s-%s entries=%zu\n",
                       this_key.c_str(), cc.atoms[cas_a1].id.c_str(),
                       cc.atoms[cas_a4].id.c_str(), bcp.entries.size());
        }
        size_t cross_pos = SIZE_MAX;
        for (size_t i = 0; i < bcp.entries.size(); ++i)
          if ((bcp.entries[i].first == cas_a1 && bcp.entries[i].second == cas_a4) ||
              (bcp.entries[i].first == cas_a4 && bcp.entries[i].second == cas_a1)) {
            cross_pos = i;
            break;
          }
        if (cross_pos == SIZE_MAX)
          continue;

        size_t seriNum = bcp.global_start + cross_pos;
        if (seriNum >= new_to_orig.size())
          continue;
        size_t shifted_orig = new_to_orig[seriNum];
        if (dbg_tor) {
          std::fprintf(stderr, "  key=%s cas=%s-%s cross=%zu seri=%zu shifted=%zu\n",
                       this_key.c_str(), cc.atoms[cas_a1].id.c_str(),
                       cc.atoms[cas_a4].id.c_str(), cross_pos, seriNum, shifted_orig);
        }
        if (shifted_orig == bcp.global_start + cross_pos)
          continue;

        // Find the replacement torsion at the shifted position.
        size_t r_a1 = SIZE_MAX, r_a4 = SIZE_MAX, src_bi = SIZE_MAX;
        for (size_t bi2 = 0; bi2 < bug_infos.size(); ++bi2) {
          auto& bcp2 = bond_cps[bi2];
          if (shifted_orig >= bcp2.global_start &&
              shifted_orig < bcp2.global_start + bcp2.entries.size()) {
            size_t cp2 = shifted_orig - bcp2.global_start;
            r_a1 = bcp2.entries[cp2].first;
            r_a4 = bcp2.entries[cp2].second;
            src_bi = bi2;
            break;
          }
        }
        if (r_a1 == SIZE_MAX || src_bi == SIZE_MAX)
          continue;

        const std::string& r_a1_name = cc.atoms[r_a1].id;
        const std::string& r_a4_name = cc.atoms[r_a4].id;
        auto try_apply = [&](const std::vector<Restraints::Torsion>& src) {
          for (const auto& t : src)
            if ((t.id1.atom == r_a1_name && t.id4.atom == r_a4_name) ||
                (t.id1.atom == r_a4_name && t.id4.atom == r_a1_name)) {
              cc.rt.torsions[info.rt_idx] = t;
              return true;
            }
          return false;
        };
        // Prefer the actual shifted source bond; fall back to same-bond lookup.
        bool applied = try_apply(bug_infos[src_bi].generated);
        if (!applied)
          applied = try_apply(info.generated);
        if (dbg_tor) {
          std::fprintf(stderr, "    -> repl=%s-%s src_bi=%zu applied=%d\n",
                       r_a1_name.c_str(), r_a4_name.c_str(), src_bi, applied ? 1 : 0);
        }
      }
    }
  }

  if (has_sugar_ring) {
    std::vector<Restraints::Torsion> nu_torsions;
    std::set<std::pair<size_t, size_t>> selected_bonds;
    auto sugar_pos = [&](size_t idx) -> Position {
      if (sugar_coord_overrides) {
        auto it = sugar_coord_overrides->find(cc.atoms[idx].id);
        if (it != sugar_coord_overrides->end())
          return it->second;
      }
      return cc.atoms[idx].xyz;
    };
    for (const auto& kv : sugar_ring_seq) {
      const std::vector<size_t>& seq = kv.second;
      size_t n = seq.size();
      if (n != 5 && n != 6)
        continue;
      for (size_t i = 0; i < n; ++i) {
        size_t a = seq[(i + n - 1) % n];
        size_t b = seq[i];
        size_t c = seq[(i + 1) % n];
        size_t d = seq[(i + 2) % n];
        selected_bonds.insert(std::minmax(b, c));
        Position pa = sugar_pos(a);
        Position pb = sugar_pos(b);
        Position pc = sugar_pos(c);
        Position pd = sugar_pos(d);
        double dih = deg(calculate_dihedral(pa, pb, pc, pd));
        if (!std::isfinite(dih))
          continue;
        Restraints::Torsion nu{
            "auto",
            {1, cc.atoms[a].id},
            {1, cc.atoms[b].id},
            {1, cc.atoms[c].id},
            {1, cc.atoms[d].id},
            dih, 10.0, 3};
        nu_torsions.push_back(std::move(nu));
      }
    }

    vector_remove_if(cc.rt.torsions, [&](const Restraints::Torsion& tor) {
      auto itb = atom_index.find(tor.id2.atom);
      auto itc = atom_index.find(tor.id3.atom);
      size_t b = itb != atom_index.end() ? itb->second : SIZE_MAX;
      size_t c = itc != atom_index.end() ? itc->second : SIZE_MAX;
      if (b == SIZE_MAX || c == SIZE_MAX)
        return false;
      auto bkey = std::minmax(b, c);
      return selected_bonds.count(bkey) != 0;
    });
    cc.rt.torsions.insert(cc.rt.torsions.end(), nu_torsions.begin(), nu_torsions.end());
  }

  // Replace torsions with nucleic-acid-specific table entries when applicable.
  const bool apply_nucl_tors = false;  // AceDRG does not apply nucl_tors in CCD outputs
  if (apply_nucl_tors && nucleic_mode && !cc.rt.torsions.empty()) {
    std::vector<Restraints::Torsion> replaced;
    replaced.reserve(cc.rt.torsions.size());
    std::unordered_set<std::string> seen_keys;
    for (const auto& tor : cc.rt.torsions) {
      std::vector<TorsionEntry> entries;
      if (tables.lookup_nucl_tors(tor.id1.atom, tor.id2.atom, tor.id3.atom, tor.id4.atom, entries) ||
          tables.lookup_nucl_tors(tor.id4.atom, tor.id3.atom, tor.id2.atom, tor.id1.atom, entries)) {
        for (const auto& e : entries) {
          Restraints::Torsion t = tor;
          t.label = e.id;
          t.value = e.value;
          t.esd = e.sigma;
          t.period = e.period;
          std::string key = cat(t.id1.atom, '|', t.id2.atom, '|', t.id3.atom, '|',
                                t.id4.atom, '|', t.label);
          if (seen_keys.insert(key).second)
            replaced.push_back(std::move(t));
        }
      } else {
        std::string key = cat(tor.id1.atom, '|', tor.id2.atom, '|', tor.id3.atom, '|',
                              tor.id4.atom, '|', tor.label);
        if (seen_keys.insert(key).second)
          replaced.push_back(tor);
      }
    }
    cc.rt.torsions.swap(replaced);
  }

}

void add_chirality_if_missing(
    ChemComp& cc, const std::map<std::string, std::string>& atom_stereo,
    const std::vector<CodAtomInfo>& atom_info,
    const AceGraphView& graph) {
  if (!cc.rt.chirs.empty())
    return;

  auto& adj = graph.adjacency;
  std::vector<unsigned> cip_ranks = compute_rdkit_legacy_cip_ranks(cc, adj);
  std::set<size_t> stereo_centers;

  auto stereo_sign = [&](size_t center) -> ChiralityType {
    if (cc.atoms[center].el != El::C)
      return ChiralityType::Both;
    auto st_it = atom_stereo.find(cc.atoms[center].id);
    if (st_it == atom_stereo.end() || st_it->second.empty())
      return ChiralityType::Both;
    char s = lower(st_it->second[0]);
    if (s == 's')
      return ChiralityType::Positive;
    if (s == 'r')
      return ChiralityType::Negative;
    return ChiralityType::Both;
  };

  for (size_t center = 0; center < cc.atoms.size(); ++center) {
    if (atom_info[center].hybrid != Hybridization::SP3)
      continue;

    std::vector<size_t> non_h = non_hydrogen_neighbors(cc, adj, center);
    std::vector<size_t> h = hydrogen_neighbors(cc, adj, center);
    if (cc.atoms[center].el == El::C) {
      bool has_metal_neighbor = false;
      for (size_t nb : non_h)
        if (cc.atoms[nb].el.is_metal()) {
          has_metal_neighbor = true;
          break;
        }
      if (has_metal_neighbor)
        continue;
    }
    bool n_sp3_2h1_case = false;
    if (cc.atoms[center].el == El::N && non_h.size() == 2 && h.size() == 1) {
      bool n31_like = true;
      for (size_t nb : non_h)
        if (cc.atoms[nb].el != El::C || atom_info[nb].hybrid != Hybridization::SP3) {
          n31_like = false;
          break;
        }
      bool has_sp_neighbor = false;
      for (size_t nb : non_h) {
        if (cc.atoms[nb].el == El::S || cc.atoms[nb].el == El::P ||
            cc.atoms[nb].el == El::O) {
          has_sp_neighbor = true;
          break;
        }
        if (cc.atoms[nb].el == El::N) {
          bool n_has_h = false;
          for (const auto& nb2 : adj[nb])
            if (nb2.idx != center && cc.atoms[nb2.idx].is_hydrogen()) {
              n_has_h = true;
              break;
            }
          if (n_has_h) {
            has_sp_neighbor = true;
            break;
          }
        }
      }
      n_sp3_2h1_case = n31_like || has_sp_neighbor;
    }
    if (non_h.size() < 3 && !n_sp3_2h1_case)
      continue;

    ChiralityType sign = stereo_sign(center);
    bool is_stereo_carbon = (sign != ChiralityType::Both && cc.atoms[center].el == El::C);
    bool is_stereo_s_carbon = (is_stereo_carbon && sign == ChiralityType::Positive);
    if (is_stereo_carbon)
      stereo_centers.insert(center);
    auto bond_to_center_type = [&](size_t nb_idx) {
      for (const auto& nb : adj[center])
        if (nb.idx == nb_idx)
          return nb.type;
      return BondType::Unspec;
    };
    if (cc.atoms[center].el == El::S && non_h.size() == 4) {
      int n_c = 0, n_n = 0, n_o = 0, n_other = 0;
      bool all_single = true;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        if (el == El::C)
          ++n_c;
        else if (el == El::N)
          ++n_n;
        else if (el == El::O)
          ++n_o;
        else
          ++n_other;
        if (bond_to_center_type(idx) != BondType::Single)
          all_single = false;
      }
      // AceDRG does not emit chirality for sulfamates/sulfamides with
      // S(single)-{C,N,O,O} neighbors.
      if (all_single && n_c == 1 && n_n == 1 && n_o == 2 && n_other == 0)
        continue;
    }
    if (cc.atoms[center].el == El::P) {
      int n_c = 0, n_metal = 0, n_other = 0;
      bool all_single = true;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        if (el == El::C)
          ++n_c;
        else if (el.is_metal())
          ++n_metal;
        else
          ++n_other;
        if (bond_to_center_type(idx) != BondType::Single)
          all_single = false;
      }
      // AceDRG omits chirality for trigonal phosphines.
      if (all_single && non_h.size() == 3 && n_c == 3 && n_metal == 0 && n_other == 0)
        continue;
      // AceDRG omits chirality for phosphines bound to one metal and three carbons.
      if (all_single && non_h.size() == 4 && n_c == 3 && n_metal == 1 && n_other == 0)
        continue;
    }
    sort_neighbors_by_rdkit_cip_rank(non_h, cip_ranks);
    if (cc.atoms[center].el == El::C || cc.atoms[center].el == El::N) {
      auto stereo_label_rank = [&](size_t idx) {
        if (cc.atoms[idx].el != El::C)
          return 0;
        auto st_it = atom_stereo.find(cc.atoms[idx].id);
        if (st_it == atom_stereo.end() || st_it->second.empty())
          return 0;
        char s = lower(st_it->second[0]);
        if (s == 'r')
          return 2;
        if (s == 's')
          return 1;
        return 0;
      };
      auto branch_stereo_rank = [&](size_t idx) {
        int rank = stereo_label_rank(idx);
        if (rank == 2)
          return rank;
        for (const auto& nb2 : adj[idx]) {
          if (nb2.idx == center || cc.atoms[nb2.idx].is_hydrogen())
            continue;
          rank = std::max(rank, stereo_label_rank(nb2.idx));
          if (rank == 2)
            break;
        }
        return rank;
      };
      std::stable_sort(non_h.begin(), non_h.end(), [&](size_t a, size_t b) {
        unsigned ra = cip_ranks[a];
        unsigned rb = cip_ranks[b];
        if (ra != rb)
          return ra > rb;
        int sa = branch_stereo_rank(a);
        int sb = branch_stereo_rank(b);
        if (sa != sb)
          return sa > sb;
        return false;
      });
    }
    sort_neighbors_by_rdkit_cip_rank(h, cip_ranks);
    bool assign_noncarbon_sign = false;
    if (!is_stereo_carbon &&
        (cc.atoms[center].el == El::P || cc.atoms[center].el == El::S)) {
      int terminal_oxygen_count = 0;
      for (const auto& nb : adj[center]) {
        if (cc.atoms[nb.idx].el != El::O)
          continue;
        // AceDRG's "one-bond O" check uses total O degree, so O-H does not
        // count as terminal here.
        if (adj[nb.idx].size() == 1)
          ++terminal_oxygen_count;
      }
      assign_noncarbon_sign = terminal_oxygen_count <= 1;
    }
    if (!is_stereo_carbon && !assign_noncarbon_sign)
      sign = ChiralityType::Both;
    bool force_negative_cationic_n = false;
    auto branch_non_h_count = [&](size_t idx) {
      int count = 0;
      for (const auto& nb2 : adj[idx])
        if (nb2.idx != center && !cc.atoms[nb2.idx].is_hydrogen())
          ++count;
      return count;
    };

    std::vector<size_t> chosen;
    if (cc.atoms[center].el == El::C) {
      std::vector<size_t> halogens;
      for (size_t nb : non_h)
        if (is_halogen(cc.atoms[nb].el))
          halogens.push_back(nb);
      if (halogens.size() >= 3) {
        chosen.push_back(halogens[0]);
        chosen.push_back(halogens[1]);
        chosen.push_back(halogens[2]);
      }
    }
    if (chosen.empty() && n_sp3_2h1_case) {
      if (non_h.size() >= 2 && !h.empty())
        chosen = {non_h[0], non_h[1], h[0]};
    } else if (chosen.empty()) {
      for (size_t idx : non_h) {
        if (chosen.size() == 3)
          break;
        chosen.push_back(idx);
      }
      for (size_t idx : h) {
        if (chosen.size() == 3)
          break;
        chosen.push_back(idx);
      }
    }
    if (cc.atoms[center].el == El::As && non_h.size() > 4 && chosen.size() >= 3) {
      std::vector<size_t> oxy;
      for (size_t idx : non_h)
        if (cc.atoms[idx].el == El::O)
          oxy.push_back(idx);
      if (oxy.size() >= 3) {
        sort_neighbors_by_rdkit_cip_rank(oxy, cip_ranks);
        chosen = {oxy[0], oxy[1], oxy[2]};
      }
    }
    if (cc.atoms[center].el == El::B && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> oxy_single;
      for (size_t idx : non_h)
        if (cc.atoms[idx].el == El::O && bond_to_center_type(idx) == BondType::Single)
          oxy_single.push_back(idx);
      if (oxy_single.size() == 4) {
        std::vector<size_t> branched;
        std::vector<size_t> terminal;
        for (size_t idx : oxy_single) {
          bool has_non_h_other = false;
          for (const auto& nb2 : adj[idx])
            if (nb2.idx != center && !cc.atoms[nb2.idx].is_hydrogen()) {
              has_non_h_other = true;
              break;
            }
          (has_non_h_other ? branched : terminal).push_back(idx);
        }
        if (branched.size() >= 2 && !terminal.empty())
          chosen = {branched[0], branched[1], terminal[0]};
      }
    }
    if (cc.atoms[center].el == El::As && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> dbl_o;
      std::vector<size_t> sing_o;
      for (size_t idx : non_h) {
        if (cc.atoms[idx].el != El::O)
          continue;
        BondType bt = bond_to_center_type(idx);
        if (bt == BondType::Double || bt == BondType::Deloc)
          dbl_o.push_back(idx);
        else if (bt == BondType::Single)
          sing_o.push_back(idx);
      }
      if (dbl_o.size() == 1 && sing_o.size() >= 2)
        chosen = {dbl_o[0], sing_o[0], sing_o[1]};
    }
    if (cc.atoms[center].el == El::N &&
        cc.atoms[center].charge > 0 &&
        non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> carbons;
      std::vector<size_t> nitrogens;
      bool all_single = true;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        if (el == El::C)
          carbons.push_back(idx);
        else if (el == El::N)
          nitrogens.push_back(idx);
        if (bond_to_center_type(idx) != BondType::Single)
          all_single = false;
      }
      if (all_single && carbons.size() == 4) {
        std::stable_sort(carbons.begin(), carbons.end(), [&](size_t a, size_t b) {
          return branch_non_h_count(a) > branch_non_h_count(b);
        });
        chosen = {carbons[0], carbons[1], carbons[2]};
      } else if (all_single && carbons.size() == 3 && nitrogens.size() == 1) {
        std::stable_sort(carbons.begin(), carbons.end(), [&](size_t a, size_t b) {
          return branch_non_h_count(a) > branch_non_h_count(b);
        });
        chosen = {nitrogens[0], carbons[0], carbons[1]};
        force_negative_cationic_n = true;
      }
    }
    if (is_stereo_carbon &&
        cc.atoms[center].el == El::C &&
        non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> oxygens;
      std::vector<size_t> carbons;
      for (size_t idx : non_h) {
        if (cc.atoms[idx].el == El::O)
          oxygens.push_back(idx);
        else if (cc.atoms[idx].el == El::C)
          carbons.push_back(idx);
      }
      if (oxygens.size() == 1 && carbons.size() == 3) {
        auto carbon_rank = [&](size_t idx) {
          bool has_triple_other = false;
          bool has_hetero_other = false;
          bool has_pi_other = false;
          int non_h_other = 0;
          int second_shell_non_h = 0;
          for (const auto& nb2 : adj[idx]) {
            if (nb2.idx == center)
              continue;
            if (!cc.atoms[nb2.idx].is_hydrogen()) {
              ++non_h_other;
              if (cc.atoms[nb2.idx].el != El::C)
                has_hetero_other = true;
              for (const auto& nb3 : adj[nb2.idx])
                if (nb3.idx != idx && !cc.atoms[nb3.idx].is_hydrogen())
                  ++second_shell_non_h;
            }
            if (nb2.type == BondType::Double ||
                nb2.type == BondType::Deloc ||
                nb2.type == BondType::Aromatic ||
                nb2.type == BondType::Triple)
              has_pi_other = true;
            if (nb2.type == BondType::Triple)
              has_triple_other = true;
          }
          return std::make_tuple(has_triple_other, has_hetero_other,
                                 has_pi_other, non_h_other,
                                 second_shell_non_h,
                                 cc.atoms[idx].id);
        };
        bool any_triple = false;
        for (size_t c_idx : carbons)
          if (std::get<0>(carbon_rank(c_idx))) {
            any_triple = true;
            break;
          }
        std::stable_sort(carbons.begin(), carbons.end(), [&](size_t a, size_t b) {
          auto ma = carbon_rank(a);
          auto mb = carbon_rank(b);
          if (std::get<1>(ma) != std::get<1>(mb))
            return std::get<1>(ma) > std::get<1>(mb);
          if (any_triple) {
            if (std::get<0>(ma) != std::get<0>(mb))
              return std::get<0>(ma) > std::get<0>(mb);
            if (std::get<2>(ma) != std::get<2>(mb))
              return std::get<2>(ma) > std::get<2>(mb);
            if (std::get<3>(ma) != std::get<3>(mb))
              return std::get<3>(ma) > std::get<3>(mb);
          } else {
            if (std::get<3>(ma) != std::get<3>(mb))
              return std::get<3>(ma) > std::get<3>(mb);
            if (std::get<4>(ma) != std::get<4>(mb))
              return std::get<4>(ma) > std::get<4>(mb);
            if (std::get<2>(ma) != std::get<2>(mb))
              return std::get<2>(ma) > std::get<2>(mb);
          }
          return std::get<5>(ma) < std::get<5>(mb);
        });
        size_t c1 = carbons[0];
        size_t c2 = carbons[1];
        bool oxygen_has_non_h_other = false;
        for (const auto& nb2 : adj[oxygens[0]]) {
          if (nb2.idx == center || cc.atoms[nb2.idx].is_hydrogen())
            continue;
          oxygen_has_non_h_other = true;
          break;
        }
        double vol = calculate_chiral_volume(cc.atoms[center].xyz,
                                             cc.atoms[oxygens[0]].xyz,
                                             cc.atoms[c1].xyz,
                                             cc.atoms[c2].xyz);
        if (std::isfinite(vol) && std::fabs(vol) > 1e-8) {
          ChiralityType vol_sign = (vol > 0.0) ? ChiralityType::Positive
                                               : ChiralityType::Negative;
          if (vol_sign != sign && !oxygen_has_non_h_other)
            std::swap(c1, c2);
        }
        chosen = {oxygens[0], c1, c2};
      }
    }
    if (cc.atoms[center].el == El::C && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> oxygens;
      std::vector<size_t> carbons;
      for (size_t idx : non_h) {
        if (cc.atoms[idx].el == El::O)
          oxygens.push_back(idx);
        else if (cc.atoms[idx].el == El::C)
          carbons.push_back(idx);
      }
      if (oxygens.size() == 2 && carbons.size() == 2) {
        auto carbon_has_oxygen_branch = [&](size_t c_idx) {
          for (const auto& nb2 : adj[c_idx])
            if (nb2.idx != center && cc.atoms[nb2.idx].el == El::O)
              return true;
          return false;
        };
        bool c0_oxy = carbon_has_oxygen_branch(carbons[0]);
        bool c1_oxy = carbon_has_oxygen_branch(carbons[1]);
        if (c0_oxy != c1_oxy)
          chosen = {oxygens[0], oxygens[1], c0_oxy ? carbons[0] : carbons[1]};
      }
    }
    if (cc.atoms[center].el == El::C && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> heavy_halogens;
      std::vector<size_t> phosphorus;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        if (el == El::Br || el == El::I || el == El::At)
          heavy_halogens.push_back(idx);
        else if (cc.atoms[idx].el == El::P)
          phosphorus.push_back(idx);
      }
      if (heavy_halogens.size() == 2 && phosphorus.size() == 2) {
        auto p_bridge_score = [&](size_t p_idx) {
          int bridged_o_to_p = 0;
          int non_h_other = 0;
          for (const auto& nb2 : adj[p_idx]) {
            if (nb2.idx == center || cc.atoms[nb2.idx].is_hydrogen())
              continue;
            ++non_h_other;
            if (cc.atoms[nb2.idx].el != El::O)
              continue;
            for (const auto& nb3 : adj[nb2.idx]) {
              if (nb3.idx == p_idx || cc.atoms[nb3.idx].is_hydrogen())
                continue;
              if (cc.atoms[nb3.idx].el == El::P) {
                ++bridged_o_to_p;
                break;
              }
            }
          }
          return std::make_tuple(bridged_o_to_p, non_h_other, cc.atoms[p_idx].id);
        };
        auto p0 = p_bridge_score(phosphorus[0]);
        auto p1 = p_bridge_score(phosphorus[1]);
        size_t picked_p = phosphorus[0];
        if (std::get<0>(p1) > std::get<0>(p0) ||
            (std::get<0>(p1) == std::get<0>(p0) &&
             std::get<1>(p1) > std::get<1>(p0)))
          picked_p = phosphorus[1];
        chosen = {heavy_halogens[0], heavy_halogens[1], picked_p};
      }
    }
    if (cc.atoms[center].el == El::S && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> dbl_o;
      std::vector<size_t> sing_n;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        BondType bt = bond_to_center_type(idx);
        if (el == El::O && (bt == BondType::Double || bt == BondType::Deloc))
          dbl_o.push_back(idx);
        else if (el == El::N && bt == BondType::Single)
          sing_n.push_back(idx);
      }
      if (dbl_o.size() == 2 && sing_n.size() == 2) {
        std::sort(dbl_o.begin(), dbl_o.end(),
                  [&](size_t a, size_t b) { return cc.atoms[a].id < cc.atoms[b].id; });
        std::sort(sing_n.begin(), sing_n.end(),
                  [&](size_t a, size_t b) { return cc.atoms[a].id < cc.atoms[b].id; });
        chosen = {dbl_o[0], dbl_o[1], sing_n.back()};
      }
    }
    if (cc.atoms[center].el == El::P && non_h.size() == 4 && chosen.size() >= 3) {
      std::vector<size_t> sulfur;
      std::vector<size_t> oxygens;
      for (size_t idx : non_h) {
        if (cc.atoms[idx].el == El::S)
          sulfur.push_back(idx);
        else if (cc.atoms[idx].el == El::O)
          oxygens.push_back(idx);
      }
      if (sulfur.size() == 2 && oxygens.size() == 2) {
        size_t single_o = SIZE_MAX;
        int n_double_o = 0;
        for (size_t o_idx : oxygens) {
          BondType bt = bond_to_center_type(o_idx);
          if (bt == BondType::Single)
            single_o = o_idx;
          if (bt == BondType::Double || bt == BondType::Deloc)
            ++n_double_o;
        }
        if (single_o != SIZE_MAX && n_double_o == 1) {
          std::sort(sulfur.begin(), sulfur.end(), [&](size_t a, size_t b) {
            return cc.atoms[a].id < cc.atoms[b].id;
          });
          chosen = {sulfur[0], sulfur[1], single_o};
          assign_noncarbon_sign = false;
          sign = ChiralityType::Both;
        }
      }
      if (sulfur.size() == 1 && oxygens.size() == 3) {
        std::stable_sort(oxygens.begin(), oxygens.end(), [&](size_t a, size_t b) {
          unsigned ra = cip_ranks[a];
          unsigned rb = cip_ranks[b];
          if (ra != rb)
            return ra > rb;
          return false;
        });
        chosen = {sulfur[0], oxygens[0], oxygens[1]};
      }
    }
    if (chosen.size() < 3)
      continue;
    if (cc.atoms[center].el == El::P) {
      int sulfur_idx = -1;
      int sulfur_count = 0;
      int oxygen_count = 0;
      for (int i = 0; i != 3; ++i) {
        if (cc.atoms[chosen[i]].el == El::S) {
          ++sulfur_count;
          sulfur_idx = i;
        } else if (cc.atoms[chosen[i]].el == El::O) {
          ++oxygen_count;
        }
      }
      if (sulfur_count == 1 && oxygen_count == 2 && sulfur_idx > 0) {
        size_t sulfur = chosen[sulfur_idx];
        chosen.erase(chosen.begin() + sulfur_idx);
        chosen.insert(chosen.begin(), sulfur);
      }
      if (cc.atoms[chosen[0]].el == El::O &&
          cc.atoms[chosen[1]].el == El::O &&
          cc.atoms[chosen[2]].el == El::O) {
        auto branch_rank = [&](size_t o_idx) {
          bool has_p_other = false;
          bool has_metal_other = false;
          bool has_h_other = false;
          bool has_non_metal_other = false;
          for (const auto& nb2 : adj[o_idx]) {
            if (nb2.idx == center || cc.atoms[nb2.idx].is_hydrogen())
              continue;
            if (!cc.atoms[nb2.idx].el.is_metal())
              has_non_metal_other = true;
            if (cc.atoms[nb2.idx].el == El::P)
              has_p_other = true;
            if (cc.atoms[nb2.idx].el.is_metal())
              has_metal_other = true;
          }
          for (const auto& nb2 : adj[o_idx])
            if (nb2.idx != center && cc.atoms[nb2.idx].is_hydrogen())
              has_h_other = true;
          BondType bt = bond_to_center_type(o_idx);
          if (has_p_other)
            return 0;
          if (bt == BondType::Single && has_h_other)
            return 1;
          if (bt == BondType::Single && has_non_metal_other && !has_metal_other)
            return 1;
          if (has_metal_other)
            return 2;
          if (bt == BondType::Double || bt == BondType::Deloc)
            return 3;
          return 1;
        };
        bool any_metal = false;
        for (size_t idx : chosen)
          for (const auto& nb2 : adj[idx])
            if (nb2.idx != center && !cc.atoms[nb2.idx].is_hydrogen() &&
                cc.atoms[nb2.idx].el.is_metal()) {
              any_metal = true;
              break;
            }
        if (any_metal) {
          std::stable_sort(chosen.begin(), chosen.end(), [&](size_t a, size_t b) {
            return branch_rank(a) < branch_rank(b);
          });
        }
      }
    }
    if (assign_noncarbon_sign) {
      bool has_metal_branch = false;
      int protonated_terminal_oxy = 0;
      for (size_t idx : chosen) {
        bool has_non_h_other = false;
        bool has_h_other = false;
        bool has_metal_other = false;
        for (const auto& nb2 : adj[idx]) {
          if (nb2.idx == center)
            continue;
          if (cc.atoms[nb2.idx].is_hydrogen()) {
            has_h_other = true;
          } else {
            has_non_h_other = true;
            if (cc.atoms[nb2.idx].el.is_metal())
              has_metal_other = true;
          }
        }
        if (has_metal_other)
          has_metal_branch = true;
        if (cc.atoms[idx].el == El::O &&
            bond_to_center_type(idx) == BondType::Single &&
            has_h_other && !has_non_h_other)
          ++protonated_terminal_oxy;
      }
      if (has_metal_branch || protonated_terminal_oxy >= 2)
        assign_noncarbon_sign = false;
    }
    bool force_negative_noncarbon = false;
    if (assign_noncarbon_sign && cc.atoms[center].el == El::P) {
      int carbon_count = 0;
      int double_oxygen_count = 0;
      int protonated_single_oxygen_count = 0;
      int other_count = 0;
      for (size_t idx : non_h) {
        Element el = cc.atoms[idx].el;
        BondType bt = bond_to_center_type(idx);
        if (el == El::C) {
          ++carbon_count;
          continue;
        }
        if (el == El::O) {
          bool has_non_h_other = false;
          bool has_h_other = false;
          for (const auto& nb2 : adj[idx]) {
            if (nb2.idx == center)
              continue;
            if (cc.atoms[nb2.idx].is_hydrogen())
              has_h_other = true;
            else
              has_non_h_other = true;
          }
          if (bt == BondType::Double || bt == BondType::Deloc)
            ++double_oxygen_count;
          else if (bt == BondType::Single && has_h_other && !has_non_h_other)
            ++protonated_single_oxygen_count;
          else
            ++other_count;
          continue;
        }
        ++other_count;
      }
      auto st_it = atom_stereo.find(cc.atoms[center].id);
      bool stereo_n = st_it != atom_stereo.end() && !st_it->second.empty() &&
                      lower(st_it->second[0]) == 'n';
      force_negative_noncarbon = stereo_n &&
                                 carbon_count == 2 &&
                                 double_oxygen_count == 1 &&
                                 protonated_single_oxygen_count == 1 &&
                                 other_count == 0;
    }

    if (is_stereo_carbon && sign != ChiralityType::Both) {
      double vol = calculate_chiral_volume(cc.atoms[center].xyz,
                                           cc.atoms[chosen[0]].xyz,
                                           cc.atoms[chosen[1]].xyz,
                                           cc.atoms[chosen[2]].xyz);
      if (std::isfinite(vol) && std::fabs(vol) > 1e-8) {
        sign = (vol > 0.0) ? ChiralityType::Positive : ChiralityType::Negative;
        auto chosen_stereo_label = [&](size_t idx) {
          auto st_it = atom_stereo.find(cc.atoms[idx].id);
          if (st_it == atom_stereo.end() || st_it->second.empty())
            return '\0';
          char s = lower(st_it->second[0]);
          if (s == 'r' || s == 's')
            return s;
          return '\0';
        };
        if (sign == ChiralityType::Positive &&
            !is_stereo_s_carbon &&
            non_h.size() == 4 &&
            cc.atoms[chosen[0]].el == El::C &&
            cc.atoms[chosen[1]].el == El::C &&
            cc.atoms[chosen[2]].el == El::C &&
            chosen_stereo_label(chosen[0]) == '\0' &&
            chosen_stereo_label(chosen[1]) == '\0' &&
            chosen_stereo_label(chosen[2]) == '\0') {
          sign = ChiralityType::Negative;
        }
        if (sign == ChiralityType::Positive &&
            !is_stereo_s_carbon &&
            cc.atoms[chosen[0]].el == El::O &&
            cc.atoms[chosen[1]].el == El::C &&
            cc.atoms[chosen[2]].el == El::C) {
          char s1 = chosen_stereo_label(chosen[1]);
          char s2 = chosen_stereo_label(chosen[2]);
          if ((s1 == 'r' && s2 == 's') || (s1 == 's' && s2 == 'r'))
            sign = ChiralityType::Negative;
        }
        if (sign == ChiralityType::Positive &&
            is_stereo_s_carbon &&
            non_h.size() == 3 &&
            h.size() == 1 &&
            std::fabs(vol) < 0.02 &&
            cc.atoms[chosen[0]].el == El::C &&
            cc.atoms[chosen[1]].el == El::C &&
            cc.atoms[chosen[2]].el == El::C &&
            chosen_stereo_label(chosen[0]) == '\0' &&
            chosen_stereo_label(chosen[1]) == '\0' &&
            chosen_stereo_label(chosen[2]) == '\0') {
          sign = ChiralityType::Negative;
        }
        if (is_stereo_s_carbon && sign == ChiralityType::Positive &&
            std::fabs(vol) < 0.03 && h.size() == 1) {
          bool has_oxygen = false;
          bool has_nitrogen = false;
          for (size_t idx : chosen) {
            if (cc.atoms[idx].el == El::O)
              has_oxygen = true;
            if (cc.atoms[idx].el == El::N)
              has_nitrogen = true;
          }
          if (has_oxygen && !has_nitrogen)
            sign = ChiralityType::Negative;
        }
        if (sign == ChiralityType::Negative &&
            std::fabs(vol) < 0.02 && h.size() == 1) {
          bool all_carbon = true;
          for (size_t idx : chosen)
            if (cc.atoms[idx].el != El::C) {
              all_carbon = false;
              break;
            }
          if (all_carbon && cc.atoms[chosen[1]].id > cc.atoms[chosen[2]].id)
            std::swap(chosen[1], chosen[2]);
        }
      } else {
        sign = ChiralityType::Both;
      }
    } else if (assign_noncarbon_sign) {
      double vol = calculate_chiral_volume(cc.atoms[center].xyz,
                                           cc.atoms[chosen[0]].xyz,
                                           cc.atoms[chosen[1]].xyz,
                                           cc.atoms[chosen[2]].xyz);
      if (std::isfinite(vol) && std::fabs(vol) > 1e-8)
        sign = (vol > 0.0) ? ChiralityType::Positive : ChiralityType::Negative;
      else
        sign = ChiralityType::Both;
      if (force_negative_noncarbon)
        sign = ChiralityType::Negative;
      if (cc.atoms[center].el == El::S &&
          sign == ChiralityType::Positive &&
          cc.atoms[chosen[0]].el == El::O &&
          cc.atoms[chosen[1]].el == El::C &&
          cc.atoms[chosen[2]].el == El::C &&
          (bond_to_center_type(chosen[0]) == BondType::Double ||
           bond_to_center_type(chosen[0]) == BondType::Deloc)) {
        bool carbon_has_pi = false;
        for (size_t c_idx : {chosen[1], chosen[2]}) {
          for (const auto& nb2 : adj[c_idx]) {
            if (nb2.idx == center)
              continue;
            if (nb2.type == BondType::Double ||
                nb2.type == BondType::Deloc ||
                nb2.type == BondType::Aromatic ||
                nb2.type == BondType::Triple) {
              carbon_has_pi = true;
              break;
            }
          }
          if (carbon_has_pi)
            break;
        }
        if (carbon_has_pi)
          sign = ChiralityType::Negative;
      }
    }
    if (force_negative_cationic_n)
      sign = ChiralityType::Negative;

    cc.rt.chirs.push_back({{1, cc.atoms[center].id},
                           {1, cc.atoms[chosen[0]].id},
                           {1, cc.atoms[chosen[1]].id},
                           {1, cc.atoms[chosen[2]].id},
                           sign});
  }

  std::string type = to_upper(cc.type_or_group);
  ChiralityType sign = ChiralityType::Both;
  if (type.find("L-PEPTIDE") != std::string::npos)
    sign = ChiralityType::Positive;
  else if (type.find("D-PEPTIDE") != std::string::npos)
    sign = ChiralityType::Negative;
  if (sign == ChiralityType::Both)
    return;

  auto ca = cc.find_atom("CA");
  if (ca == cc.atoms.end())
    return;
  if (cc.find_atom("N") == cc.atoms.end() ||
      cc.find_atom("C") == cc.atoms.end() ||
      cc.find_atom("CB") == cc.atoms.end())
    return;

  size_t ca_idx = std::distance(cc.atoms.begin(), ca);
  if (stereo_centers.count(ca_idx) == 0)
    cc.rt.chirs.push_back({{1, "CA"}, {1, "N"}, {1, "C"}, {1, "CB"}, sign});
}

void add_planes_if_missing(ChemComp& cc,
                           const std::vector<CodAtomInfo>& atom_info,
                           const AceGraphView& graph) {
  auto& adj = graph.adjacency;

  std::vector<std::set<Restraints::AtomId>> plane_sets;
  plane_sets.reserve(cc.rt.planes.size());
  for (const auto& plane : cc.rt.planes)
    plane_sets.emplace_back(plane.ids.begin(), plane.ids.end());

  auto add_plane = [&](const std::set<size_t>& idxs, size_t min_size = 4,
                       double esd = 0.02) {
    if (idxs.size() < min_size)
      return;
    std::set<Restraints::AtomId> ids;
    for (size_t idx : idxs)
      ids.insert({1, cc.atoms[idx].id});
    for (size_t i = 0; i < plane_sets.size(); ++i)
      if (plane_sets[i] == ids) {
        if (i < cc.rt.planes.size() && cc.rt.planes[i].esd < esd)
          cc.rt.planes[i].esd = esd;
        return;
      }
    Restraints::Plane plane;
    plane.label = cat("plan-", cc.rt.planes.size() + 1);
    plane.esd = esd;
    std::vector<std::string> names;
    names.reserve(ids.size());
    for (const auto& id : ids)
      names.push_back(id.atom);
    std::sort(names.begin(), names.end());
    for (const std::string& name : names)
      plane.ids.push_back({1, name});
    cc.rt.planes.push_back(std::move(plane));
    plane_sets.push_back(std::move(ids));
  };

  std::map<int, std::vector<size_t>> rings;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int ring_id : atom_info[i].in_rings)
      rings[ring_id].push_back(i);

  std::set<size_t> aromatic_ring_atoms;
  for (const auto& ring : rings) {
    const std::vector<size_t>& ring_atoms = ring.second;
    if (ring_atoms.size() < 4)
      continue;
    bool aromatic = true;
    for (size_t idx : ring_atoms) {
      if (!atom_info[idx].is_aromatic) {
        aromatic = false;
        break;
      }
    }
    if (!aromatic)
      continue;
    std::set<size_t> plane_idx;
    for (size_t idx : ring_atoms) {
      plane_idx.insert(idx);
      aromatic_ring_atoms.insert(idx);
      // AceDRG: only add neighbors if the ring atom has != 4 connections
      if (adj[idx].size() != 4) {
        for (const auto& nb : adj[idx])
          if (!cc.atoms[nb.idx].el.is_metal())
            plane_idx.insert(nb.idx);
      }
    }
    add_plane(plane_idx);
  }

  // SP2 planes for atoms not already in aromatic ring planes (AceDRG order: atom serial)
  for (size_t idx = 0; idx < cc.atoms.size(); ++idx) {
    if (atom_info[idx].hybrid != Hybridization::SP2)
      continue;
    if (cc.atoms[idx].is_hydrogen())
      continue;
    if (aromatic_ring_atoms.count(idx))
      continue;
    std::set<size_t> plane_idx;
    plane_idx.insert(idx);
    for (const auto& nb : adj[idx])
      plane_idx.insert(nb.idx);
    add_plane(plane_idx);
  }

  // Metal coordination planes (AceDRG metalMode.py logic):
  // For each metal-nonmetal bond (M-X) where X is SP2 per metalMode rules,
  // create a plane with M + X + X's filtered neighbors.
  // metalMode.py uses its own hybridization:
  //   O: always SP3 (never generates metal planes)
  //   N: SP2 only if has an SP2 non-metal neighbor AND total connections < 4
  //   C: SP2 if non-metal connections < 3, or non-metal == 3 with no metal
  auto metal_mode_is_sp2 = [&](size_t x) -> bool {
    Element el = cc.atoms[x].el;
    if (el == El::O || el == El::P || el == El::S || el == El::Se)
      return false;
    int non_metal = 0, metal = 0;
    for (const auto& nb2 : adj[x]) {
      if (cc.atoms[nb2.idx].el.is_metal())
        metal++;
      else
        non_metal++;
    }
    if (el == El::C || el == El::Si || el == El::Ge) {
      if (non_metal == 3 && metal > 0) return false;  // SP3 in metalMode
      return atom_info[x].bonding_idx == 2;
    }
    if (el == El::N || el == El::As) {
      // metalMode second round: N is SP2 if has SP2 non-metal neighbor
      // and total connections < 4
      if (non_metal + metal >= 4) return false;
      for (const auto& nb2 : adj[x])
        if (!cc.atoms[nb2.idx].el.is_metal() && atom_info[nb2.idx].bonding_idx == 2)
          return true;
      return false;
    }
    return atom_info[x].bonding_idx == 2;
  };
  std::set<size_t> metal_indices;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    if (cc.atoms[i].el.is_metal())
      metal_indices.insert(i);
  for (size_t m : metal_indices) {
    std::set<size_t> m_neighbors;
    for (const auto& nb : adj[m])
      m_neighbors.insert(nb.idx);
    for (const auto& nb : adj[m]) {
      size_t x = nb.idx;
      if (cc.atoms[x].el.is_metal())
        continue;
      if (!metal_mode_is_sp2(x))
        continue;
      std::set<size_t> plane_idx;
      plane_idx.insert(m);
      plane_idx.insert(x);
      for (const auto& nb2 : adj[x]) {
        size_t y = nb2.idx;
        if (y == m)
          continue;
        if (cc.atoms[y].el.is_metal())
          continue;
        if (m_neighbors.count(y))
          continue;
        plane_idx.insert(y);
      }
      add_plane(plane_idx, 3, 0.06);
    }
  }
}

void apply_metal_charge_corrections(ChemComp& cc) {
  if (cc.atoms.empty())
    return;
  auto graph = make_ace_graph_view(cc);
  auto& adj = graph.adjacency;
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    const Element& el = cc.atoms[i].el;
    if (el.is_metal() || el == El::H)
      continue;
    if (!has_metal_and_non_metal_heavy_neighbor(cc, adj, i))
      continue;
    int expected_valence = expected_valence_for_nonmetal(el);
    if (expected_valence == 0)
      continue;
    float sum_bo = sum_non_metal_bond_order(cc, adj, i);
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    cc.atoms[i].charge = static_cast<float>(-rem_v);
  }
}

}  // namespace

void prepare_chemcomp(ChemComp& cc, const AcedrgTables& tables,
                      const std::map<std::string, std::string>& atom_stereo,
                      bool no_angles,
                      const std::map<std::string, Position>* sugar_coord_overrides) {
  if (!no_angles)
    add_angles_from_bonds_if_missing(cc);

  // Collect all original atom names for H naming collision checks.
  std::set<std::string> used_names;
  for (const auto& atom : cc.atoms)
    used_names.insert(atom.id);

  adjust_oxoacid_group(cc, El::P, 3, true);   // phosphate
  adjust_oxoacid_group(cc, El::S, 4, false);  // sulfate
  adjust_hexafluorophosphate(cc);
  adjust_carboxy_asp(cc);
  adjust_terminal_carboxylate(cc);
  adjust_guanidinium_group(cc, used_names);
  adjust_amino_ter_amine(cc, used_names);
  adjust_terminal_amine(cc, used_names);

  bool added_h3 = add_n_terminal_h3(cc);

  apply_metal_charge_corrections(cc);

  int missing_bonds = count_missing_values(cc.rt.bonds);
  int missing_angles = no_angles ? 0 : count_missing_values(cc.rt.angles);
  bool need_fill = (missing_bonds > 0 || missing_angles > 0);

  if (need_fill) {
    tables.fill_restraints(cc);
    if (!no_angles) {
      if (added_h3)
        sync_n_terminal_h3_angles(cc);
    }
    std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);
    AceGraphView graph = make_ace_graph_view(cc);
    add_chirality_if_missing(cc, atom_stereo, atom_info, graph);
    add_torsions_from_bonds_if_missing(cc, tables, atom_info, atom_stereo, graph,
                                       sugar_coord_overrides);
    add_planes_if_missing(cc, atom_info, graph);
  } else {
    if (added_h3 && !no_angles)
      sync_n_terminal_h3_angles(cc);
  }

  tables.assign_ccp4_types(cc);
}

} // namespace gemmi
