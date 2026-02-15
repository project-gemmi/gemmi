// Copyright 2025 Global Phasing Ltd.
//
// Restraint generation helpers for prepare_chemcomp():
// chemical adjustments, torsions, chirality, planes, H naming.

#include "gemmi/ace_cc.hpp"
#include "gemmi/acedrg_tables.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/ace_graph.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
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

int element_priority(Element el) {
  if (el == El::N) return 0;
  if (el == El::C) return 1;
  if (el == El::O) return 2;
  if (el == El::S) return 3;
  if (el == El::P) return 4;
  if (el == El::Se) return 5;
  return 6;
}

int chirality_priority(Element el) {
  if (el == El::O) return 0;
  if (el == El::N) return 1;
  if (el == El::S) return 2;
  if (el == El::P) return 3;
  if (el == El::C) return 4;
  if (el == El::H) return 5;
  return 6;
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

struct SugarRingInfo {
  std::map<int, std::vector<size_t>> ring_orders;
  std::set<std::pair<size_t, size_t>> ring_bonds;
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
  auto trace_ring_order = [&](const std::set<size_t>& rset, size_t o_idx,
                             size_t start) -> std::vector<size_t> {
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
        return {};
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
      return {};
    return order;
  };
  auto ring_shape = [&](const std::vector<size_t>& order) -> std::string {
    int o_c = 0;
    if (order.empty())
      return "";
    size_t o_idx = order.front();
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
      if (atom_info[idx].hybrid == Hybridization::SP3)
        ++sp3;
      if (atom_info[idx].is_aromatic || atom_info[idx].hybrid != Hybridization::SP3) {
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
    std::vector<size_t> order1 = trace_ring_order(rset, oxy_idx, o_ring_nbs[0]);
    std::vector<size_t> order2 = trace_ring_order(rset, oxy_idx, o_ring_nbs[1]);
    std::string shape1 = ring_shape(order1);
    std::string shape2 = ring_shape(order2);
    if (allowed_shapes.count(shape1) == 0 && allowed_shapes.count(shape2) == 0)
      continue;
    if (allowed_shapes.count(shape1) != 0)
      out.ring_orders.emplace(kv.first, order1);
    else
      out.ring_orders.emplace(kv.first, order2);
  }
  for (const auto& kv : out.ring_orders) {
    const auto& rset = kv.second;
    std::set<size_t> rset2(rset.begin(), rset.end());
    for (size_t idx : rset)
      for (const auto& nb : adj[idx])
        if (rset2.count(nb.idx) && idx < nb.idx)
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
    size_t a1 = chiral_legs[0], a2 = chiral_legs[1], a3 = chiral_legs[2];
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
  return ai.bonding_idx > 0 ? ai.bonding_idx == 1
                            : ai.hybrid == Hybridization::SP1;
}

bool is_sp2_like(const CodAtomInfo& ai) {
  return ai.bonding_idx > 0 ? ai.bonding_idx == 2
                            : ai.hybrid == Hybridization::SP2;
}

bool is_sp3_like(const CodAtomInfo& ai) {
  return ai.bonding_idx > 0 ? ai.bonding_idx == 3
                            : ai.hybrid == Hybridization::SP3;
}

enum class TvMode { Default, SP3SP3, SP2SP3_SP3, SP3_OXY };

int max_tv_len_for_center(const CodAtomInfo& ai) {
  return is_sp2_like(ai) ? 2 : 3;
}

std::vector<size_t> build_tv_neighbor_order(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info, size_t ctr) {
  std::vector<size_t> nb_order;
  nb_order.reserve(adj[ctr].size());
  if (is_sp2_like(atom_info[ctr])) {
    for (const auto& nb : adj[ctr])
      if (!cc.atoms[nb.idx].is_hydrogen())
        nb_order.push_back(nb.idx);
    for (const auto& nb : adj[ctr])
      if (cc.atoms[nb.idx].is_hydrogen())
        nb_order.push_back(nb.idx);
  } else {
    for (const auto& nb : adj[ctr])
      nb_order.push_back(nb.idx);
  }
  return nb_order;
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
    bool use_chiral_mut_table = !stereo_sp3;
    auto mit = chir_mut_table.find(ctr);
    if (use_chiral_mut_table && mit != chir_mut_table.end()) {
      auto mt_it = mit->second.find(other);
      if (mt_it != mit->second.end()) {
        used_chiral = true;
        for (size_t cand : mt_it->second) {
          if ((int)tv.size() >= max_len)
            break;
          if (cand == other || cand == rs)
            continue;
          if (std::find(tv.begin(), tv.end(), cand) != tv.end())
            continue;
          tv.push_back(cand);
        }
      }
    }
  }

  if (!used_chiral) {
    bool want_h_first = false;
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

int compute_tv_position_for_center(
    const ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const std::set<size_t>& stereo_chiral_centers,
    const std::map<size_t, std::map<size_t, std::vector<size_t>>>& chir_mut_table,
    size_t center, size_t other_center, size_t target,
    TvMode mode = TvMode::Default, size_t forced_rs = SIZE_MAX) {
  auto tv = build_tv_list_for_center(
      cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
      center, other_center, mode, forced_rs);
  for (int i = 0; i < (int)tv.size(); ++i)
    if (tv[i] == target)
      return i;
  return -1;
}


static void emit_one_torsion(
    ChemComp& cc, const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    const AcedrgTables& tables,
    const std::set<size_t>& stereo_chiral_centers,
    const std::map<size_t, std::map<size_t, std::vector<size_t>>>& chir_mut_table,
    const std::map<std::pair<size_t, size_t>, RingParity>& bond_ring_parity,
    bool peptide_mode,
    size_t center2, size_t center3,
    bool bond_aromatic, int ring_size,
    size_t a1_idx, size_t a4_idx) {
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

  bool phenol_oh =
      ((cc.atoms[center2].id == "CZ" && cc.atoms[center3].id == "OH") ||
       (cc.atoms[center3].id == "CZ" && cc.atoms[center2].id == "OH")) &&
      a4->id == "HH";

  if (phenol_oh) {
    value = 0.0;
    esd = 5.0;
    period = 2;
  } else if (!lookup_found && bond_aromatic && ring_size > 0 && info2.is_aromatic && info3.is_aromatic &&
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
    int i_pos = compute_tv_position_for_center(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side1, side2, term1, TvMode::SP3SP3, rs1);
    int j_pos = compute_tv_position_for_center(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side2, side1, term2, TvMode::SP3SP3, rs2);
    if (pit != bond_ring_parity.end() && i_pos >= 0 && j_pos >= 0) {
      static const double even_m[3][3] = {
        {60,180,-60}, {-60,60,180}, {180,-60,60}};
      static const double odd_m[3][3] = {
        {-60,60,180}, {180,-60,60}, {60,180,-60}};
      static const double noflip_m[3][3] = {
        {180,-60,60}, {60,180,-60}, {-60,60,180}};
      const auto& m = (pit->second == RingParity::Even) ? even_m
                     : (pit->second == RingParity::Odd) ? odd_m
                     : noflip_m;
      value = m[i_pos][j_pos];
    }
  } else if (!lookup_found && sp3_2 && sp3_3 && ring_size == 0) {
    // Non-ring SP3-SP3: always use no-flip matrix
    size_t side1 = center2;
    size_t side2 = center3;
    std::pair<size_t, size_t> rs_pair =
        find_ring_sharing_pair(adj, atom_info, side1, side2);
    size_t rs1 = rs_pair.first;
    size_t rs2 = rs_pair.second;
    size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
    size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
    int i_pos = compute_tv_position_for_center(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side1, side2, term1, TvMode::SP3SP3, rs1);
    int j_pos = compute_tv_position_for_center(
        cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
        side2, side1, term2, TvMode::SP3SP3, rs2);
    if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
      static const double noflip_m[3][3] = {
        {180,-60,60}, {60,180,-60}, {-60,60,180}};
      value = noflip_m[i_pos][j_pos];
    }
    // Peptide-like chi1 around CA-CB...
    size_t ca_idx = SIZE_MAX, cb_idx = SIZE_MAX;
    if (cc.atoms[center2].id == "CA" && cc.atoms[center3].id == "CB") {
      ca_idx = center2; cb_idx = center3;
    } else if (cc.atoms[center3].id == "CA" && cc.atoms[center2].id == "CB") {
      ca_idx = center3; cb_idx = center2;
    }
    if (ca_idx != SIZE_MAX) {
      size_t n_term = (ca_idx == center2) ? a1_idx : a4_idx;
      size_t cg_term = (cb_idx == center2) ? a1_idx : a4_idx;
      if (cc.atoms[n_term].id == "N" && cc.atoms[cg_term].id == "CG") {
        int cb_h = 0;
        for (const auto& nb : adj[cb_idx])
          if (nb.idx != ca_idx && cc.atoms[nb.idx].is_hydrogen())
            ++cb_h;
        if (cb_h >= 2) {
          std::vector<size_t> cg_nonh_other;
          for (const auto& nb : adj[cg_term])
            if (nb.idx != cb_idx && !cc.atoms[nb.idx].is_hydrogen())
              cg_nonh_other.push_back(nb.idx);
          bool only_n_branch = (cg_nonh_other.size() == 1 &&
                                cc.atoms[cg_nonh_other[0]].el == El::N);
          if (only_n_branch)
            value = -60.0;
          else
            value = 180.0;
        }
      }
    }
    // AceDRG: for N(P)2(H)-P bonds...
    auto has_double_oxo = [&](size_t p_idx, size_t term_idx) {
      if (cc.atoms[p_idx].el != El::P || cc.atoms[term_idx].el != El::O)
        return false;
      for (const auto& nb : adj[p_idx])
        if (nb.idx == term_idx)
          return nb.type == BondType::Double || nb.type == BondType::Deloc;
      return false;
    };
    auto n_diphos_h = [&](size_t n_idx) {
      if (cc.atoms[n_idx].el != El::N)
        return false;
      int p_cnt = 0, h_cnt = 0;
      for (const auto& nb : adj[n_idx]) {
        if (cc.atoms[nb.idx].el == El::P) ++p_cnt;
        if (cc.atoms[nb.idx].is_hydrogen()) ++h_cnt;
      }
      return p_cnt == 2 && h_cnt >= 1;
    };
    size_t n_center = SIZE_MAX;
    size_t n_term = SIZE_MAX, p_term = SIZE_MAX;
    if (cc.atoms[center2].el == El::N && cc.atoms[center3].el == El::P) {
      n_center = center2; n_term = a1_idx; p_term = a4_idx;
    } else if (cc.atoms[center3].el == El::N && cc.atoms[center2].el == El::P) {
      n_center = center3; n_term = a4_idx; p_term = a1_idx;
    }
    if (n_center != SIZE_MAX && n_diphos_h(n_center) &&
        cc.atoms[n_term].el == El::P &&
        cc.atoms[p_term].el == El::O) {
      value = 60.0;
    }
    size_t o_center = SIZE_MAX;
    size_t p_center2 = SIZE_MAX;
    size_t o_term = SIZE_MAX;
    size_t p_term2 = SIZE_MAX;
    if (cc.atoms[center2].el == El::O && cc.atoms[center3].el == El::P) {
      o_center = center2; p_center2 = center3; o_term = a1_idx; p_term2 = a4_idx;
    } else if (cc.atoms[center3].el == El::O && cc.atoms[center2].el == El::P) {
      o_center = center3; p_center2 = center2; o_term = a4_idx; p_term2 = a1_idx;
    }
    if (o_center != SIZE_MAX &&
        cc.atoms[o_term].el == El::C &&
        cc.atoms[p_term2].el == El::O &&
        has_double_oxo(p_center2, p_term2)) {
      bool has_n_diphos_h_nb = false;
      for (const auto& nb : adj[p_center2]) {
        if (nb.idx != o_center && n_diphos_h(nb.idx)) {
          has_n_diphos_h_nb = true;
          break;
        }
      }
      if (has_n_diphos_h_nb)
        value = 180.0;
    }
  } else if (!lookup_found && ((sp2_2 && sp3_3) || (sp3_2 && sp2_3))) {
    // SP2-SP3:
    // - regular path: tS3 matrix for "two SP2-side neighbours in same ring",
    //   otherwise default matrix (AceDRG SetOneSP2SP3Bond).
    // - Oxy-column SP2 center (O/S/Se/Te/Po): use SetOneSP3OxyColumnBond
    //   style period-3 progression.
    size_t sp2_center = sp2_2 ? center2 : center3;
    size_t sp3_center = sp2_2 ? center3 : center2;
    size_t sp2_term = sp2_2 ? a1_idx : a4_idx;
    size_t sp3_term = sp2_2 ? a4_idx : a1_idx;
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
      // AceDRG SetOneSP3OxyColumnBond(tIdx1=SP3, tIdx2=SP2-oxy):
      // row index on SP3 side (tV1), column index on SP2 side (tV2),
      // torsion = init + (row+col)*120 with init depending on ring/oxy.
      int i_pos = compute_tv_position_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          sp3_center, sp2_center, sp3_term, TvMode::SP3_OXY, sp3_rs);
      int j_pos = compute_tv_position_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          sp2_center, sp3_center, sp2_term, TvMode::Default, sp2_rs);
      if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
        double ini_value = (sp2_rs != SIZE_MAX) ? 60.0
                       : (oxy_col_sp3 ? 90.0 : 180.0);
        value = normalize_acedrg_angle(ini_value + (i_pos + j_pos) * 120.0);
        esd = 20.0;
        period = 3;
      }
    } else {
      int i_pos = compute_tv_position_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          sp2_center, sp3_center, sp2_term, TvMode::Default, sp2_rs);
      int j_pos = compute_tv_position_for_center(
          cc, adj, atom_info, stereo_chiral_centers, chir_mut_table,
          sp3_center, sp2_center, sp3_term, TvMode::SP2SP3_SP3, sp3_rs);
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
      const auto& sp2sp3_m = is_ts3 ? ts3_m : ts1_m;
      if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 3) {
        value = sp2sp3_m[i_pos][j_pos];
        esd = 20.0;
        period = 6;
      }
    }
  } else if (!lookup_found && sp2_2 && sp2_3) {
    // Non-aromatic SP2-SP2: use 2x2 matrix.
    // Ring-sharing pair is selected once globally (AceDRG tS1/tS2 logic).
    size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                    ? center2 : center3;
    size_t side2 = (side1 == center2) ? center3 : center2;
    auto rs_pair = find_ring_sharing_pair(adj, atom_info, side1, side2);
    size_t rs1 = rs_pair.first;
    size_t rs2 = rs_pair.second;
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
      auto nonmetal_degree = [&](size_t idx) {
        int deg = 0;
        for (const auto& nb : adj[idx])
          if (!cc.atoms[nb.idx].el.is_metal())
            ++deg;
        return deg;
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
          else if (nonmetal_degree(tv[0]) < nonmetal_degree(tv[1]))
            std::swap(tv[0], tv[1]);
        } else {
          if (h0 && (!h1 || nonmetal_degree(tv[0]) < nonmetal_degree(tv[1])))
            std::swap(tv[0], tv[1]);
          else if (!h0 && nonmetal_degree(tv[0]) < nonmetal_degree(tv[1]) && !h1)
            std::swap(tv[0], tv[1]);
        }
      }
      find_target:
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

  if (!lookup_found && sp3_2 && sp3_3) {
    auto is_sulfone_oxo = [&](size_t center, size_t term) {
      if (cc.atoms[center].el != El::S || cc.atoms[term].el != El::O)
        return false;
      int n_oxo = 0;
      bool term_is_oxo = false;
      for (const auto& nb : adj[center]) {
        bool is_oxo = (cc.atoms[nb.idx].el == El::O &&
                       (nb.type == BondType::Double || nb.type == BondType::Deloc));
        if (is_oxo) {
          ++n_oxo;
          if (nb.idx == term)
            term_is_oxo = true;
        }
      }
      return term_is_oxo && n_oxo >= 2;
    };
    auto n_is_substituted = [&](size_t n_idx, size_t s_idx) {
      if (cc.atoms[n_idx].el != El::N)
        return false;
      for (const auto& nb : adj[n_idx])
        if (nb.idx != s_idx && !cc.atoms[nb.idx].is_hydrogen())
          return true;
      return false;
    };
    bool ns_sulfone_oxo = ((cc.atoms[center2].el == El::N &&
                            is_sulfone_oxo(center3, a4_idx) &&
                            n_is_substituted(center2, center3)) ||
                           (cc.atoms[center3].el == El::N &&
                            is_sulfone_oxo(center2, a1_idx) &&
                            n_is_substituted(center3, center2)));
    if (ns_sulfone_oxo) {
      value = 180.0;
      esd = 10.0;
      period = 3;
    }
    auto is_sulfone_single_o = [&](size_t s_idx, size_t term_idx) {
      if (cc.atoms[s_idx].el != El::S || cc.atoms[term_idx].el != El::O)
        return false;
      bool term_single = false;
      int n_oxo = 0;
      for (const auto& nb : adj[s_idx]) {
        bool is_oxo = (cc.atoms[nb.idx].el == El::O &&
                       (nb.type == BondType::Double || nb.type == BondType::Deloc));
        if (is_oxo)
          ++n_oxo;
        if (nb.idx == term_idx)
          term_single = !is_oxo;
      }
      return term_single && n_oxo >= 2;
    };
    size_t n_center = SIZE_MAX, s_center = SIZE_MAX;
    size_t n_term = SIZE_MAX, s_term = SIZE_MAX;
    if (cc.atoms[center2].el == El::N && cc.atoms[center3].el == El::S) {
      n_center = center2; s_center = center3; n_term = a1_idx; s_term = a4_idx;
    } else if (cc.atoms[center3].el == El::N && cc.atoms[center2].el == El::S) {
      n_center = center3; s_center = center2; n_term = a4_idx; s_term = a1_idx;
    }
    if (n_center != SIZE_MAX &&
        is_sulfone_single_o(s_center, s_term) &&
        !cc.atoms[n_term].is_hydrogen()) {
      value = -60.0;
      esd = 10.0;
      period = 3;
    }
  }
  cc.rt.torsions.push_back({"auto",
                            {1, a1->id},
                            {1, cc.atoms[center2].id},
                            {1, cc.atoms[center3].id},
                            {1, a4->id},
                            value, esd, period});
}

double sugar_ring_torsion_value(const ChemComp& cc,
                               size_t a1, size_t a2, size_t a3, size_t a4) {
  double value = deg(calculate_dihedral(cc.atoms[a1].xyz, cc.atoms[a2].xyz,
                                       cc.atoms[a3].xyz, cc.atoms[a4].xyz));
  return std::isfinite(value) ? value : 0.0;
}

void add_sugar_ring_torsions(ChemComp& cc,
                            const std::vector<size_t>& ring_order) {
  if (ring_order.empty())
    return;
  if (ring_order.size() != 5 && ring_order.size() != 6)
    return;

  size_t n = ring_order.size();
  size_t oxygen = SIZE_MAX;
  for (size_t i = 0; i < n; ++i) {
    if (cc.atoms[ring_order[i]].el == El::O) {
      oxygen = i;
      break;
    }
  }
  if (oxygen == SIZE_MAX)
    return;

  std::vector<size_t> seq = ring_order;
  if (oxygen != 0)
    std::rotate(seq.begin(), seq.begin() + static_cast<long>(oxygen), seq.end());

  for (size_t i = 0; i < n; ++i) {
    size_t a1 = seq[(i + n - 1) % n];
    size_t a2 = seq[i];
    size_t a3 = seq[(i + 1) % n];
    size_t a4 = seq[(i + 2) % n];
    double value = sugar_ring_torsion_value(cc, a1, a2, a3, a4);
    cc.rt.torsions.push_back({"auto",
                              {1, cc.atoms[a1].id},
                              {1, cc.atoms[a2].id},
                              {1, cc.atoms[a3].id},
                              {1, cc.atoms[a4].id},
                              value, 10.0, 3});
  }
}

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables,
                                        const std::vector<CodAtomInfo>& atom_info,
                                        const std::map<std::string, std::string>& atom_stereo,
                                        const AceGraphView& graph) {
  if (!cc.rt.torsions.empty())
    return;

  auto& atom_index = graph.atom_index;
  auto& adj = graph.adjacency;
  bool peptide_mode = ChemComp::is_peptide_group(cc.group);
  std::vector<bool> aromatic_like = build_aromatic_like_mask(cc, atom_info, atom_index);

  // Pre-compute ring parity used in SP3-SP3 torsion matrix selection.
  std::map<std::pair<size_t, size_t>, RingParity> bond_ring_parity =
      build_ring_bond_parity(adj, atom_info);

  // AceDRG has a dedicated sugar-ring mode: ring bonds are represented by
  // one nu torsion each (from ring geometry), while non-ring bonds keep the
  // full torsion set.
  // Match AceDRG checkOneRingSugar()/getRStr() shape gating:
  // only specific "(OC2)(...)" ring-shape strings are treated as sugar.
  SugarRingInfo sugar_info = detect_sugar_rings(cc, adj, atom_info);
  auto& sugar_ring_orders = sugar_info.ring_orders;
  auto& sugar_ring_bonds = sugar_info.ring_bonds;
  bool has_sugar_ring = !sugar_ring_bonds.empty();

  ChiralCenterInfo chiral_info =
      detect_chiral_centers_and_mut_table(cc, adj, atom_info, atom_stereo);
  auto& stereo_chiral_centers = chiral_info.stereo_chiral_centers;
  auto& chir_mut_table = chiral_info.chir_mut_table;

  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;
    if (cc.atoms[idx1].el.is_metal() || cc.atoms[idx2].el.is_metal())
      continue;

    int ring_size = shared_ring_size_from_ring_ids(atom_info[idx1].in_rings,
                                                   atom_info[idx1].min_ring_size,
                                                   atom_info[idx2].in_rings,
                                                   atom_info[idx2].min_ring_size);
    size_t center2 = idx1;
    bool c1_carbonyl = is_carbonyl_carbon(cc, adj, idx1);
    bool c2_carbonyl = is_carbonyl_carbon(cc, adj, idx2);
    bool idx1_in_ring = (atom_info[idx1].min_ring_size > 0 || aromatic_like[idx1]);
    bool idx2_in_ring = (atom_info[idx2].min_ring_size > 0 || aromatic_like[idx2]);
    if (c1_carbonyl != c2_carbonyl) {
      center2 = c1_carbonyl ? idx1 : idx2;
    } else if (cc.atoms[idx1].id == "CA" || cc.atoms[idx2].id == "CA") {
      center2 = (cc.atoms[idx1].id == "CA") ? idx1 : idx2;
    } else if (idx1_in_ring && idx2_in_ring) {
      center2 = (cc.atoms[idx1].id < cc.atoms[idx2].id) ? idx1 : idx2;
    } else if (idx1_in_ring != idx2_in_ring) {
      if (cc.atoms[idx1].el == El::O || cc.atoms[idx2].el == El::O)
        center2 = (cc.atoms[idx1].el == El::O) ? idx2 : idx1;
      else
        center2 = idx1_in_ring ? idx2 : idx1;
    } else if (ring_size > 0 &&
               ((cc.atoms[idx1].el == El::C && cc.atoms[idx2].el == El::N) ||
                (cc.atoms[idx2].el == El::C && cc.atoms[idx1].el == El::N))) {
      center2 = (cc.atoms[idx1].el == El::C) ? idx1 : idx2;
    } else if (cc.atoms[idx1].el != cc.atoms[idx2].el) {
      center2 = (element_priority(cc.atoms[idx1].el) <
                 element_priority(cc.atoms[idx2].el)) ? idx1 : idx2;
    } else if (cc.atoms[idx2].id < cc.atoms[idx1].id) {
      center2 = idx2;
    }

    if (ring_size > 0 && idx1_in_ring && idx2_in_ring &&
        cc.atoms[idx1].el == El::C && cc.atoms[idx2].el == El::C)
      center2 = (cc.atoms[idx1].id < cc.atoms[idx2].id) ? idx1 : idx2;
    size_t center3 = (center2 == idx1) ? idx2 : idx1;
    if (ring_size == 0 && idx1_in_ring && !idx2_in_ring &&
        is_sp3_like(atom_info[idx1]) && is_sp3_like(atom_info[idx2])) {
      center2 = idx1;
      center3 = idx2;
    }

    bool bond_aromatic = (bond.type == BondType::Aromatic ||
                          bond.type == BondType::Deloc || bond.aromatic ||
                          (aromatic_like[idx1] && aromatic_like[idx2]));

    std::vector<size_t> tv1_idx, tv2_idx;
    TvMode mode12 = TvMode::Default;
    TvMode mode21 = TvMode::Default;
    bool sp3_like_12 = is_sp3_like(atom_info[center2]);
    bool sp3_like_21 = is_sp3_like(atom_info[center3]);
    bool sp2_like_12 = is_sp2_like(atom_info[center2]);
    bool sp2_like_21 = is_sp2_like(atom_info[center3]);
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
    ring_size = shared_ring_size_from_ring_ids(atom_info[center2].in_rings,
                                               atom_info[center2].min_ring_size,
                                               atom_info[center3].in_rings,
                                               atom_info[center3].min_ring_size);
    if (tv1_idx.empty() || tv2_idx.empty())
      continue;

    if (is_sp1_like(atom_info[center2]) || is_sp1_like(atom_info[center3]))
      continue;

    auto emit_torsion = [&](size_t a1_idx, size_t a4_idx) {
      emit_one_torsion(cc, adj, atom_info, tables, stereo_chiral_centers,
                        chir_mut_table, bond_ring_parity, peptide_mode,
                        center2, center3, bond_aromatic, ring_size,
                        a1_idx, a4_idx);
    };

    auto pick_pair = [&]() -> std::pair<size_t, size_t> {
      size_t side1 = center2;
      size_t side2 = center3;

      std::vector<size_t> r1, n1, h1, r2, n2, h2;
      for (const auto& nb : adj[side1]) {
        if (nb.idx == side2) continue;
        if (!atom_info[nb.idx].in_rings.empty())
          r1.push_back(nb.idx);
        else if (cc.atoms[nb.idx].is_hydrogen())
          h1.push_back(nb.idx);
        else if (cc.atoms[nb.idx].el.is_metal()) {
        } else
          n1.push_back(nb.idx);
      }
      for (const auto& nb : adj[side2]) {
        if (nb.idx == side1) continue;
        if (!atom_info[nb.idx].in_rings.empty())
          r2.push_back(nb.idx);
        else if (cc.atoms[nb.idx].is_hydrogen())
          h2.push_back(nb.idx);
        else if (cc.atoms[nb.idx].el.is_metal()) {
        } else
          n2.push_back(nb.idx);
      }

      // Match AceDRG CodClassify::selectOneTorFromOneBond():
      // classify terminals as ring/non-H/H and pick the first candidate
      // in each class by connection order.
      if (!n1.empty() && !n2.empty())
        return {n1.front(), n2.front()};

      if (!n1.empty() && n2.empty()) {
        if (!r2.empty())
          return {n1.front(), r2.front()};
        if (!h2.empty())
          return {n1.front(), h2.front()};
      }

      if (n1.empty() && !n2.empty()) {
        if (!r1.empty())
          return {r1.front(), n2.front()};
        if (!h1.empty())
          return {h1.front(), n2.front()};
      }

      if (!r1.empty() || !r2.empty()) {
        if (!r1.empty() && !r2.empty())
          return {r1.front(), r2.front()};
        if (!r1.empty() && !h2.empty())
          return {r1.front(), h2.front()};
        if (!h1.empty() && !r2.empty())
          return {h1.front(), r2.front()};
      }

      if (!h1.empty() && !h2.empty())
        return {h1.front(), h2.front()};

      return {tv1_idx.front(), tv2_idx.front()};
    };

    auto chosen = pick_pair();
    emit_torsion(chosen.first, chosen.second);
  }

  if (has_sugar_ring) {
    vector_remove_if(cc.rt.torsions, [&](const Restraints::Torsion& tor) {
      auto b_it = atom_index.find(tor.id2.atom);
      auto c_it = atom_index.find(tor.id3.atom);
      if (b_it == atom_index.end() || c_it == atom_index.end())
        return false;
      auto bkey = std::minmax(b_it->second, c_it->second);
      return sugar_ring_bonds.count(bkey) != 0;
    });
    for (const auto& kv : sugar_ring_orders)
      add_sugar_ring_torsions(cc, kv.second);
  }

}

void add_chirality_if_missing(
    ChemComp& cc, const std::map<std::string, std::string>& atom_stereo,
    const std::vector<CodAtomInfo>& atom_info,
    const AceGraphView& graph) {
  if (!cc.rt.chirs.empty())
    return;

  auto& adj = graph.adjacency;
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
    bool n_sp3_2h1_case = (cc.atoms[center].el == El::N &&
                           non_h.size() == 2 && !h.empty());
    if (non_h.size() < 3 && !n_sp3_2h1_case)
      continue;

    ChiralityType sign = stereo_sign(center);
    bool is_stereo_carbon = (sign != ChiralityType::Both && cc.atoms[center].el == El::C);
    if (is_stereo_carbon)
      stereo_centers.insert(center);
    std::stable_sort(non_h.begin(), non_h.end(), [&](size_t a, size_t b) {
      return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
    });
    std::stable_sort(h.begin(), h.end(), [&](size_t a, size_t b) {
      return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
    });
    if (!is_stereo_carbon)
      sign = ChiralityType::Both;

    std::vector<size_t> chosen;
    if (n_sp3_2h1_case) {
      std::vector<size_t> n_non_h = non_h;
      std::stable_sort(n_non_h.begin(), n_non_h.end(), [&](size_t a, size_t b) {
        return cc.atoms[a].id > cc.atoms[b].id;
      });
      chosen.push_back(n_non_h[0]);
      chosen.push_back(n_non_h[1]);
      chosen.push_back(h[0]);
    } else {
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
    if (chosen.size() < 3)
      continue;

    if (is_stereo_carbon && sign != ChiralityType::Both) {
      double vol = calculate_chiral_volume(cc.atoms[center].xyz,
                                           cc.atoms[chosen[0]].xyz,
                                           cc.atoms[chosen[1]].xyz,
                                           cc.atoms[chosen[2]].xyz);
      if (std::isfinite(vol) && std::fabs(vol) > 1e-8) {
        bool should_be_pos = (sign == ChiralityType::Positive);
        bool is_pos = vol > 0.0;
        if (should_be_pos != is_pos)
          std::swap(chosen[0], chosen[1]);
      }
    }

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

  auto add_plane = [&](const std::set<size_t>& idxs) {
    if (idxs.size() < 4)
      return;
    std::set<Restraints::AtomId> ids;
    for (size_t idx : idxs)
      ids.insert({1, cc.atoms[idx].id});
    if (in_vector(ids, plane_sets))
      return;
    Restraints::Plane plane;
    plane.label = cat("plan-", cc.rt.planes.size() + 1);
    plane.esd = 0.02;
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

  for (const auto& ring : rings) {
    const std::vector<size_t>& ring_atoms = ring.second;
    if (ring_atoms.size() < 5)
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
      for (const auto& nb : adj[idx])
        plane_idx.insert(nb.idx);
    }
    add_plane(plane_idx);
  }

  // SP2 nitrogen planes first, then SP2 carbon planes (matching AceDRG order)
  for (Element sp2_el : {El::N, El::C}) {
    for (size_t idx = 0; idx < cc.atoms.size(); ++idx) {
      if (cc.atoms[idx].el != sp2_el)
        continue;
      if (atom_info[idx].hybrid != Hybridization::SP2)
        continue;
      std::set<size_t> plane_idx;
      plane_idx.insert(idx);
      for (const auto& nb : adj[idx])
        plane_idx.insert(nb.idx);
      add_plane(plane_idx);
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
                      bool only_bonds) {
  if (!only_bonds)
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
  int missing_angles = only_bonds ? 0 : count_missing_values(cc.rt.angles);
  bool need_fill = (missing_bonds > 0 || missing_angles > 0);

  if (need_fill) {
    tables.fill_restraints(cc);
    if (!only_bonds) {
      if (added_h3)
        sync_n_terminal_h3_angles(cc);
      std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);
      AceGraphView graph = make_ace_graph_view(cc);
      add_torsions_from_bonds_if_missing(cc, tables, atom_info, atom_stereo, graph);
      add_chirality_if_missing(cc, atom_stereo, atom_info, graph);
      add_planes_if_missing(cc, atom_info, graph);
    }
  } else {
    if (added_h3 && !only_bonds)
      sync_n_terminal_h3_angles(cc);
  }

  tables.assign_ccp4_types(cc);
}

} // namespace gemmi
