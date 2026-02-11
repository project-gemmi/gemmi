// Copyright 2025 Global Phasing Ltd.
//
// Restraint generation helpers for prepare_chemcomp():
// chemical adjustments, torsions, chirality, planes, H naming.

#include "gemmi/ace_cc.hpp"
#include "gemmi/acedrg_tables.hpp"
#include "gemmi/calculate.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <set>

namespace gemmi {

namespace {

template <typename Range>
int count_missing_values(const Range& range) {
  int missing = 0;
  for (const auto& item : range)
    if (std::isnan(item.value)) missing++;
  return missing;
}

std::map<std::string, std::vector<std::string>> make_neighbor_names(const ChemComp& cc) {
  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }
  return neighbors;
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

// Check if two bonded atoms are in the same ring by looking for an alternative path
bool atoms_in_same_ring(const std::string& atom1, const std::string& atom2,
                        const std::map<std::string, std::vector<std::string>>& neighbors) {
  // BFS from atom1 to atom2, excluding the direct bond between them
  std::set<std::string> visited;
  std::vector<std::string> queue;
  visited.insert(atom1);
  for (const std::string& nb : neighbors.at(atom1)) {
    if (nb != atom2) {
      queue.push_back(nb);
      visited.insert(nb);
    }
  }
  while (!queue.empty()) {
    std::string current = queue.back();
    queue.pop_back();
    if (current == atom2)
      return true;  // Found alternative path -> in same ring
    auto it = neighbors.find(current);
    if (it == neighbors.end())
      continue;
    for (const std::string& nb : it->second) {
      if (visited.find(nb) == visited.end()) {
        visited.insert(nb);
        queue.push_back(nb);
      }
    }
  }
  return false;
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
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom == ca_atom || bond.id2.atom == ca_atom) {
        if (bond.type == BondType::Double ||
            bond.type == BondType::Aromatic ||
            bond.type == BondType::Deloc) {
          return false;
        }
      }
    }
  }

  if (!ca_atom.empty()) {
    bool has_carboxylic_acid = false;
    for (const auto& bond : cc.rt.bonds) {
      std::string other;
      if (bond.id1.atom == ca_atom)
        other = bond.id2.atom;
      else if (bond.id2.atom == ca_atom)
        other = bond.id1.atom;
      else
        continue;

      auto it = cc.find_atom(other);
      if (it == cc.atoms.end() || it->el != El::C)
        continue;

      int o_count = 0;
      bool has_double_o = false;
      for (const auto& b2 : cc.rt.bonds) {
        if (b2.id1.atom != other && b2.id2.atom != other)
          continue;
        std::string o_atom = (b2.id1.atom == other) ? b2.id2.atom : b2.id1.atom;
        auto o_it = cc.find_atom(o_atom);
        if (o_it != cc.atoms.end() && o_it->el == El::O) {
          ++o_count;
          if (b2.type == BondType::Double || b2.type == BondType::Deloc)
            has_double_o = true;
        }
      }
      if (o_count >= 2 && has_double_o) {
        has_carboxylic_acid = true;
        break;
      }
    }
    if (!has_carboxylic_acid)
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

void adjust_phosphate_group(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);

  std::set<std::string> po4_phosphorus;
  for (const auto& atom : cc.atoms) {
    if (atom.el != El::P)
      continue;
    const auto& nb = neighbors[atom.id];
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
    if (nb.size() == 4 && o_count >= 3 && !has_h)
      po4_phosphorus.insert(atom.id);
  }

  std::vector<std::string> h_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    bool bonded_to_po4 = false;
    std::string h_neighbor;
    for (const std::string& nid : nb) {
      if (po4_phosphorus.count(nid))
        bonded_to_po4 = true;
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::H)
        h_neighbor = nid;
    }
    if (bonded_to_po4 && !h_neighbor.empty()) {
      atom.charge = -1.0f;
      h_to_remove.push_back(h_neighbor);
    }
  }

  for (const std::string& h_id : h_to_remove)
    remove_atom_by_id(cc, h_id);
}

void adjust_sulfate_group(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);

  std::set<std::string> sulfate_sulfur;
  for (const auto& atom : cc.atoms) {
    if (atom.el != El::S)
      continue;
    const auto& nb = neighbors[atom.id];
    if (nb.size() != 4)
      continue;
    int o_count = 0;
    for (const std::string& nid : nb) {
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::O)
        ++o_count;
    }
    if (o_count == 4)
      sulfate_sulfur.insert(atom.id);
  }

  std::vector<std::string> h_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    bool bonded_to_sulfate = false;
    std::string h_neighbor;
    for (const std::string& nid : nb) {
      if (sulfate_sulfur.count(nid))
        bonded_to_sulfate = true;
      int idx = cc.find_atom_index(nid);
      if (idx >= 0 && cc.atoms[idx].el == El::H)
        h_neighbor = nid;
    }
    if (bonded_to_sulfate && !h_neighbor.empty()) {
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

  auto get_bond_type = [&](const std::string& a, const std::string& b) {
    for (const auto& bond : cc.rt.bonds)
      if ((bond.id1.atom == a && bond.id2.atom == b) ||
          (bond.id2.atom == a && bond.id1.atom == b))
        return bond.type;
    return BondType::Unspec;
  };

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
        if (get_bond_type(atom.id, nid) == BondType::Double)
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
            if (!atoms_in_same_ring(ac, other, neighbors)) {
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

  auto get_bond_type = [&](const std::string& a, const std::string& b) {
    for (const auto& bond : cc.rt.bonds)
      if ((bond.id1.atom == a && bond.id2.atom == b) ||
          (bond.id2.atom == a && bond.id1.atom == b))
        return bond.type;
    return BondType::Unspec;
  };

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
        if (el == El::O && get_bond_type(c2_id, c2_nb) == BondType::Double)
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
    if (n_hydrogen < 2)
      continue;
    if (n_hydrogen >= 3)
      continue;

    bool c_is_sp2 = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom == alpha_carbon_id || bond.id2.atom == alpha_carbon_id) {
        if (bond.type == BondType::Double ||
            bond.type == BondType::Aromatic ||
            bond.type == BondType::Deloc) {
          c_is_sp2 = true;
          break;
        }
      }
    }
    if (c_is_sp2)
      continue;

    bool n_has_double = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom == n_atom.id || bond.id2.atom == n_atom.id) {
        if (bond.type == BondType::Double ||
            bond.type == BondType::Aromatic ||
            bond.type == BondType::Deloc) {
          n_has_double = true;
          break;
        }
      }
    }
    if (n_has_double)
      continue;

    bool has_carboxylic_acid = false;
    for (const std::string& c_neighbor : neighbors[alpha_carbon_id]) {
      int cn_idx = cc.find_atom_index(c_neighbor);
      if (cn_idx < 0)
        continue;
      if (cc.atoms[cn_idx].el != El::C)
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
      if (o_count >= 2 && has_double_o) {
        has_carboxylic_acid = true;
        break;
      }
    }
    if (!has_carboxylic_acid)
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

struct NeighborBond {
  size_t idx;
  BondType type;
  bool aromatic;
};

std::vector<std::vector<NeighborBond>> build_bond_adjacency(const ChemComp& cc,
                                                            const std::map<std::string, size_t>& atom_index) {
  std::vector<std::vector<NeighborBond>> adj(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;
    bool aromatic = (bond.type == BondType::Aromatic ||
                     bond.type == BondType::Deloc || bond.aromatic);
    adj[idx1].push_back({idx2, bond.type, aromatic});
    adj[idx2].push_back({idx1, bond.type, aromatic});
  }
  return adj;
}

bool is_carbonyl_carbon(size_t idx, const ChemComp& cc,
                        const std::vector<std::vector<NeighborBond>>& adj) {
  if (cc.atoms[idx].el != El::C)
    return false;
  for (const auto& nb : adj[idx]) {
    if (cc.atoms[nb.idx].el == El::O &&
        (nb.type == BondType::Double || nb.type == BondType::Deloc))
      return true;
  }
  return false;
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

bool share_ring(const CodAtomInfo& a, const CodAtomInfo& b) {
  for (int ring_id : a.in_rings)
    for (int other_id : b.in_rings)
      if (ring_id == other_id)
        return true;
  return false;
}

bool is_oxygen_column(Element el) {
  return el == El::O || el == El::S || el == El::Se ||
         el == El::Te || el == El::Po;
}

int shared_ring_size(const CodAtomInfo& a, const CodAtomInfo& b) {
  if (share_ring(a, b)) {
    int size = 999;
    if (a.min_ring_size > 0)
      size = std::min(size, a.min_ring_size);
    if (b.min_ring_size > 0)
      size = std::min(size, b.min_ring_size);
    return size == 999 ? 1 : size;
  }
  return 0;
}

const ChemComp::Atom* pick_torsion_neighbor(
    const ChemComp& cc,
    const std::vector<std::vector<NeighborBond>>& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> candidates;
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx)
      candidates.push_back(nb.idx);
  if (candidates.empty())
    return nullptr;

  // AceDRG priority: non-ring non-H > ring non-H > H
  // Separate non-H candidates into ring vs non-ring, prefer non-ring.
  {
    std::vector<size_t> non_ring_non_h;
    std::vector<size_t> ring_non_h;
    for (size_t idx : candidates) {
      if (cc.atoms[idx].is_hydrogen()) {
        // keep in candidates (will be filtered later)
      } else if (atom_info[idx].min_ring_size > 0) {
        ring_non_h.push_back(idx);
      } else {
        non_ring_non_h.push_back(idx);
      }
    }
    if (!non_ring_non_h.empty())
      candidates.swap(non_ring_non_h);
    else if (!ring_non_h.empty())
      candidates.swap(ring_non_h);
    // else only H atoms remain — keep candidates as-is
  }

  // AceDRG tiebreaker: first in connAtoms = bond-table order = candidates order
  return &cc.atoms[candidates.front()];
}

// AceDRG selectOneTorFromOneBond: classify as Ring / NonH / H,
// priority NonH > Ring > H, tiebreaker = first in connAtoms (bond-table order).
const ChemComp::Atom* pick_aromatic_ring_neighbor(
    const ChemComp& cc,
    const std::vector<std::vector<NeighborBond>>& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> nonring_nonh;
  std::vector<size_t> ring_nb;
  std::vector<size_t> h;
  for (const auto& nb : adj[center_idx]) {
    if (nb.idx == exclude_idx)
      continue;
    if (cc.atoms[nb.idx].is_hydrogen()) {
      h.push_back(nb.idx);
    } else if (atom_info[nb.idx].min_ring_size > 0) {
      ring_nb.push_back(nb.idx);
    } else {
      nonring_nonh.push_back(nb.idx);
    }
  }
  // tiebreaker: bond-table order (= adj iteration order, already preserved)
  if (!nonring_nonh.empty())
    return &cc.atoms[nonring_nonh.front()];
  if (!ring_nb.empty())
    return &cc.atoms[ring_nb.front()];
  if (!h.empty())
    return &cc.atoms[h.front()];
  return nullptr;
}

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables,
                                        const std::vector<CodAtomInfo>& atom_info,
                                        const std::map<std::string, std::string>& atom_stereo) {
  if (!cc.rt.torsions.empty())
    return;

  auto atom_index = cc.make_atom_index();
  auto adj = build_bond_adjacency(cc, atom_index);
  std::map<int, std::vector<size_t>> ring_map;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int rid : atom_info[i].in_rings)
      ring_map[rid].push_back(i);
  std::vector<std::vector<size_t>> sugar_rings;
  std::set<size_t> sugar_ring_atoms;
  for (const auto& r : ring_map) {
    const auto& atoms = r.second;
    if (atoms.size() != 5 && atoms.size() != 6)
      continue;
    int n_o = 0, n_c = 0;
    bool ok = true;
    for (size_t idx : atoms) {
      if (atom_info[idx].hybrid != Hybridization::SP3) {
        ok = false;
        break;
      }
      if (cc.atoms[idx].el == El::O)
        ++n_o;
      else if (cc.atoms[idx].el == El::C)
        ++n_c;
      else {
        ok = false;
        break;
      }
    }
    if (!ok || n_o != 1 || n_c != (int)atoms.size() - 1)
      continue;
    sugar_rings.push_back(atoms);
    for (size_t idx : atoms)
      sugar_ring_atoms.insert(idx);
  }
  bool sugar_mode = !sugar_rings.empty();
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

  // Pre-compute ring bond parity (even/odd/no-flip) by traversing each ring
  // starting from the smallest-index atom, following bond-table (adj) order.
  // AceDRG uses this parity to select the SP3-SP3 torsion value matrix.
  enum class RingParity { Even, Odd, NoFlip };
  std::map<std::pair<size_t,size_t>, RingParity> bond_ring_parity;
  {
    int max_ring_id = 0;
    for (size_t i = 0; i < atom_info.size(); ++i)
      for (int rid : atom_info[i].in_rings)
        max_ring_id = std::max(max_ring_id, rid + 1);
    std::vector<std::vector<size_t>> rings(max_ring_id);
    for (size_t i = 0; i < atom_info.size(); ++i)
      for (int rid : atom_info[i].in_rings)
        rings[rid].push_back(i);
    for (auto& ring : rings) {
      if (ring.empty()) continue;
      std::set<size_t> ring_set(ring.begin(), ring.end());
      size_t start = *std::min_element(ring.begin(), ring.end());
      std::vector<size_t> traversal = {start};
      size_t prev = SIZE_MAX;
      size_t cur = start;
      while (traversal.size() < ring.size()) {
        bool found = false;
        for (const auto& nb : adj[cur]) {
          if (nb.idx != prev && ring_set.count(nb.idx) &&
              std::find(traversal.begin(), traversal.end(), nb.idx)
                == traversal.end()) {
            traversal.push_back(nb.idx);
            prev = cur;
            cur = nb.idx;
            found = true;
            break;
          }
        }
        if (!found) break;
      }
      for (size_t i = 0; i + 1 < traversal.size(); ++i) {
        auto key = std::minmax(traversal[i], traversal[i+1]);
        bond_ring_parity[key] = (i % 2 == 0) ? RingParity::Even
                                              : RingParity::Odd;
      }
      if (traversal.size() == ring.size()) {
        auto ckey = std::minmax(traversal.back(), traversal.front());
        bond_ring_parity[ckey] = RingParity::NoFlip;
      }
    }
  }

  // Detect chiral centers and build mutTable for "both" chirality.
  // AceDRG generates chirality for SP3 atoms with >=3 non-H neighbors.
  // R/S carbon stereocenters get positive/negative (empty mutTable in AceDRG
  // due to string comparison bug). All other eligible atoms get "both"
  // chirality with populated mutTable.
  // Chirality affects torsion neighbor ordering: chiral atoms use plain adj
  // order (or mutTable reorder for "both"), non-chiral use first-non-H/first-H.
  std::set<size_t> chiral_centers;
  std::map<size_t, std::map<size_t, std::vector<size_t>>> chir_mut_table;
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    if (atom_info[i].hybrid != Hybridization::SP3)
      continue;
    std::vector<size_t> non_h_nbs;
    for (const auto& nb : adj[i])
      if (!cc.atoms[nb.idx].is_hydrogen())
        non_h_nbs.push_back(nb.idx);
    std::stable_sort(non_h_nbs.begin(), non_h_nbs.end(),
                     [&](size_t a, size_t b) {
                       int pa = chirality_priority(cc.atoms[a].el);
                       int pb = chirality_priority(cc.atoms[b].el);
                       return pa < pb;
                     });
    if (non_h_nbs.size() < 3)
      continue;
    chiral_centers.insert(i);
    // R/S carbon stereocenters get positive/negative chirality (empty mutTable).
    // All other eligible atoms get "both" chirality with populated mutTable.
    bool is_stereo_carbon = false;
    if (cc.atoms[i].el == El::C) {
      auto st_it = atom_stereo.find(cc.atoms[i].id);
      if (st_it != atom_stereo.end()) {
        char s = lower(st_it->second[0]);
        if (s == 'r' || s == 's')
          is_stereo_carbon = true;
      }
    }
    if (is_stereo_carbon)
      continue;
    // Build mutTable for "both" chirality: legs are first 3 non-H in adj order
    size_t a1 = non_h_nbs[0], a2 = non_h_nbs[1], a3 = non_h_nbs[2];
    size_t missing = SIZE_MAX;
    for (const auto& nb : adj[i])
      if (nb.idx != a1 && nb.idx != a2 && nb.idx != a3) {
        missing = nb.idx;
        break;
      }
    auto& mt = chir_mut_table[i];
    mt[a1] = {a3, a2};
    mt[a2] = {a1, a3};
    mt[a3] = {a2, a1};
    if (missing != SIZE_MAX) {
      mt[a1].push_back(missing);
      mt[a2].push_back(missing);
      mt[a3].push_back(missing);
      mt[missing] = {a1, a2, a3};
    }
  }

  // Torsion-value neighbor ordering mode for non-chiral atoms.
  // SP3SP3: first non-H neighbor gets priority (AceDRG SetOneSP3SP3Bond)
  // SP2SP3_SP3: first H neighbor gets priority (AceDRG SetOneSP2SP3Bond SP3 side)
  enum class TvMode { Default, SP3SP3, SP2SP3_SP3 };

  // Compute the position of a terminal atom in the torsion-value neighbor
  // ordering (tV) for a given bond center. Accounts for ring-sharing priority,
  // "both" chirality mutTable reordering, and non-chiral atom reordering.
  auto compute_tv_position = [&](size_t center, size_t other_center,
                                  size_t target,
                                  TvMode mode = TvMode::Default) -> int {
    size_t ring_sharing = SIZE_MAX;
    for (const auto& nb1 : adj[center]) {
      if (nb1.idx == other_center) continue;
      for (const auto& nb2 : adj[other_center]) {
        if (nb2.idx == center || nb2.idx == nb1.idx) continue;
        if (share_ring(atom_info[nb1.idx], atom_info[nb2.idx])) {
          ring_sharing = nb1.idx;
          goto found_rs;
        }
      }
    }
    found_rs:
    std::vector<size_t> tv;
    if (ring_sharing != SIZE_MAX)
      tv.push_back(ring_sharing);
    auto mt_it = chir_mut_table.find(center);
    if (mt_it != chir_mut_table.end()) {
      auto exc_it = mt_it->second.find(other_center);
      if (exc_it != mt_it->second.end()) {
        for (size_t a : exc_it->second)
          if (std::find(tv.begin(), tv.end(), a) == tv.end())
            tv.push_back(a);
      }
    }
    // Non-chiral atom reordering: AceDRG adds a specific atom type first
    if (chiral_centers.count(center) == 0) {
      if (mode == TvMode::SP3SP3) {
        for (const auto& nb : adj[center])
          if (nb.idx != other_center &&
              std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
              !cc.atoms[nb.idx].is_hydrogen()) {
            tv.push_back(nb.idx);
            break;
          }
      } else if (mode == TvMode::SP2SP3_SP3) {
        for (const auto& nb : adj[center])
          if (nb.idx != other_center &&
              std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
              cc.atoms[nb.idx].is_hydrogen()) {
            tv.push_back(nb.idx);
            break;
          }
      }
    }
    for (const auto& nb : adj[center]) {
      if (nb.idx != other_center &&
          std::find(tv.begin(), tv.end(), nb.idx) == tv.end())
        tv.push_back(nb.idx);
    }
    for (int i = 0; i < (int)tv.size(); ++i)
      if (tv[i] == target)
        return i;
    return -1;
  };
  auto build_tv_list = [&](size_t center, size_t other_center,
                           TvMode mode = TvMode::Default) -> std::vector<size_t> {
    size_t ring_sharing = SIZE_MAX;
    for (const auto& nb1 : adj[center]) {
      if (nb1.idx == other_center) continue;
      for (const auto& nb2 : adj[other_center]) {
        if (nb2.idx == center || nb2.idx == nb1.idx) continue;
        if (share_ring(atom_info[nb1.idx], atom_info[nb2.idx])) {
          ring_sharing = nb1.idx;
          goto found_rs;
        }
      }
    }
    found_rs:
    std::vector<size_t> tv;
    if (ring_sharing != SIZE_MAX)
      tv.push_back(ring_sharing);
    auto mt_it = chir_mut_table.find(center);
    if (mt_it != chir_mut_table.end()) {
      auto exc_it = mt_it->second.find(other_center);
      if (exc_it != mt_it->second.end()) {
        for (size_t a : exc_it->second)
          if (std::find(tv.begin(), tv.end(), a) == tv.end())
            tv.push_back(a);
      }
    }
    if (chiral_centers.count(center) == 0) {
      if (mode == TvMode::SP3SP3) {
        for (const auto& nb : adj[center])
          if (nb.idx != other_center &&
              std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
              !cc.atoms[nb.idx].is_hydrogen()) {
            tv.push_back(nb.idx);
            break;
          }
      } else if (mode == TvMode::SP2SP3_SP3) {
        for (const auto& nb : adj[center])
          if (nb.idx != other_center &&
              std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
              cc.atoms[nb.idx].is_hydrogen()) {
            tv.push_back(nb.idx);
            break;
          }
      }
    }
    for (const auto& nb : adj[center]) {
      if (nb.idx != other_center &&
          std::find(tv.begin(), tv.end(), nb.idx) == tv.end())
        tv.push_back(nb.idx);
    }
    return tv;
  };

  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;

    int ring_size = shared_ring_size(atom_info[idx1], atom_info[idx2]);
    size_t center2 = idx1;
    size_t center3 = idx2;
    bool c1_carbonyl = is_carbonyl_carbon(idx1, cc, adj);
    bool c2_carbonyl = is_carbonyl_carbon(idx2, cc, adj);
    bool idx1_in_ring = (atom_info[idx1].min_ring_size > 0 || aromatic_like[idx1]);
    bool idx2_in_ring = (atom_info[idx2].min_ring_size > 0 || aromatic_like[idx2]);
    if (c1_carbonyl != c2_carbonyl) {
      center2 = c1_carbonyl ? idx1 : idx2;
      center3 = c1_carbonyl ? idx2 : idx1;
    } else if (cc.atoms[idx1].id == "CA" || cc.atoms[idx2].id == "CA") {
      center2 = (cc.atoms[idx1].id == "CA") ? idx1 : idx2;
      center3 = (center2 == idx1) ? idx2 : idx1;
    } else if (idx1_in_ring && idx2_in_ring) {
      center2 = (cc.atoms[idx1].id > cc.atoms[idx2].id) ? idx1 : idx2;
      center3 = (center2 == idx1) ? idx2 : idx1;
    } else if (idx1_in_ring != idx2_in_ring) {
      if (cc.atoms[idx1].el == El::O || cc.atoms[idx2].el == El::O) {
        center2 = (cc.atoms[idx1].el == El::O) ? idx2 : idx1;
      } else {
        center2 = idx1_in_ring ? idx2 : idx1;
      }
      center3 = (center2 == idx1) ? idx2 : idx1;
    } else if (ring_size > 0 &&
               ((cc.atoms[idx1].el == El::C && cc.atoms[idx2].el == El::N) ||
                (cc.atoms[idx2].el == El::C && cc.atoms[idx1].el == El::N))) {
      center2 = (cc.atoms[idx1].el == El::C) ? idx1 : idx2;
      center3 = (center2 == idx1) ? idx2 : idx1;
    } else if (cc.atoms[idx1].el != cc.atoms[idx2].el) {
      center2 = (element_priority(cc.atoms[idx1].el) <
                 element_priority(cc.atoms[idx2].el)) ? idx1 : idx2;
      center3 = (center2 == idx1) ? idx2 : idx1;
    } else if (cc.atoms[idx2].id < cc.atoms[idx1].id) {
      center2 = idx2;
      center3 = idx1;
    }

    if (ring_size > 0 && idx1_in_ring && idx2_in_ring &&
        cc.atoms[idx1].el == El::C && cc.atoms[idx2].el == El::C) {
      center2 = (cc.atoms[idx1].id > cc.atoms[idx2].id) ? idx1 : idx2;
      center3 = (center2 == idx1) ? idx2 : idx1;
    }

    bool bond_aromatic = (bond.type == BondType::Aromatic ||
                          bond.type == BondType::Deloc || bond.aromatic ||
                          (aromatic_like[idx1] && aromatic_like[idx2]));
    const ChemComp::Atom* a1 = nullptr;
    const ChemComp::Atom* a4 = nullptr;
    if (bond_aromatic && idx1_in_ring && idx2_in_ring) {
      a1 = pick_aromatic_ring_neighbor(cc, adj, atom_info, center2, center3);
      a4 = pick_aromatic_ring_neighbor(cc, adj, atom_info, center3, center2);
    } else {
      a1 = pick_torsion_neighbor(cc, adj, atom_info, center2, center3);
      a4 = pick_torsion_neighbor(cc, adj, atom_info, center3, center2);
    }
    ring_size = shared_ring_size(atom_info[center2], atom_info[center3]);
    if (!a1 || !a4)
      continue;

    const CodAtomInfo& h2 = atom_info[center2];
    const CodAtomInfo& h3 = atom_info[center3];
    bool sp3_2 = (h2.hybrid == Hybridization::SP3);
    bool sp3_3 = (h3.hybrid == Hybridization::SP3);
    bool sp2_2 = (h2.hybrid == Hybridization::SP2);
    bool sp2_3 = (h3.hybrid == Hybridization::SP2);
    if (h2.hybrid == Hybridization::SP1 || h3.hybrid == Hybridization::SP1)
      continue;
    if (sugar_mode &&
        sugar_ring_atoms.count(center2) &&
        sugar_ring_atoms.count(center3))
      continue;
    if (sugar_mode && sp3_2 && sp3_3) {
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      std::vector<size_t> tv1 = build_tv_list(side1, side2, TvMode::SP3SP3);
      std::vector<size_t> tv2 = build_tv_list(side2, side1, TvMode::SP3SP3);
      static const double noflip_m[3][3] = {
        {180,-60,60}, {60,180,-60}, {-60,60,180}};
      int n1 = std::min(3, (int)tv1.size());
      int n2 = std::min(3, (int)tv2.size());
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
          cc.rt.torsions.push_back({"auto",
                                    {1, cc.atoms[tv1[i]].id},
                                    {1, cc.atoms[side1].id},
                                    {1, cc.atoms[side2].id},
                                    {1, cc.atoms[tv2[j]].id},
                                    noflip_m[i][j], 10.0, 3});
        }
      }
      continue;
    }

    TorsionEntry tors_entry;
    double value = 180.0;
    double esd = 10.0;
    int period = 3;
    bool lookup_found = tables.lookup_pep_tors(a1->id, cc.atoms[center2].id,
                                               cc.atoms[center3].id, a4->id,
                                               tors_entry);
    if (lookup_found) {
      value = tors_entry.value;
      esd = 10.0;
      period = tors_entry.period;
    }
    auto it_a1 = atom_index.find(a1->id);
    auto it_a4 = atom_index.find(a4->id);
    if (it_a1 == atom_index.end() || it_a4 == atom_index.end())
      continue;
    size_t a1_idx = it_a1->second;
    size_t a4_idx = it_a4->second;

    bool chi2_aromatic =
        (a1->id == "CA" &&
         ((cc.atoms[center2].id == "CB" && cc.atoms[center3].id == "CG") ||
          (cc.atoms[center3].id == "CB" && cc.atoms[center2].id == "CG")) &&
         (a4->id == "CD1" || a4->id == "CD2"));
    bool phenol_oh =
        ((cc.atoms[center2].id == "CZ" && cc.atoms[center3].id == "OH") ||
         (cc.atoms[center3].id == "CZ" && cc.atoms[center2].id == "OH")) &&
        a4->id == "HH";

    if (chi2_aromatic) {
      value = 90.0;
      esd = 20.0;
      period = 6;
    } else if (phenol_oh) {
      value = 0.0;
      esd = 5.0;
      period = 2;
    } else if (bond_aromatic && ring_size > 0 && h2.is_aromatic && h3.is_aromatic) {
      auto shares_ring_across = [&](size_t terminal_idx, size_t opp_center) {
        for (const auto& nb : adj[opp_center])
          if (nb.idx != center2 && nb.idx != center3 &&
              share_ring(atom_info[terminal_idx], atom_info[nb.idx]))
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
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
      size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
      int i_pos = compute_tv_position(side1, side2, term1, TvMode::SP3SP3);
      int j_pos = compute_tv_position(side2, side1, term2, TvMode::SP3SP3);
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
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
      size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
      int i_pos = compute_tv_position(side1, side2, term1, TvMode::SP3SP3);
      int j_pos = compute_tv_position(side2, side1, term2, TvMode::SP3SP3);
      if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
        static const double noflip_m[3][3] = {
          {180,-60,60}, {60,180,-60}, {-60,60,180}};
        value = noflip_m[i_pos][j_pos];
      }
    } else if (!lookup_found && ((sp2_2 && sp3_3) || (sp3_2 && sp2_3))) {
      // SP2-SP3: tS3 uses {150,-90,30; -30,90,-150},
      // all other cases use {0,120,-120; 180,-60,60}
      size_t sp2_center = sp2_2 ? center2 : center3;
      size_t sp3_center = sp2_2 ? center3 : center2;
      size_t sp2_term = sp2_2 ? a1_idx : a4_idx;
      size_t sp3_term = sp2_2 ? a4_idx : a1_idx;
      // Cross-ring sharing: SP2 neighbor shares ring with SP3 neighbor
      size_t sp2_rs = SIZE_MAX;
      for (const auto& nb1 : adj[sp2_center]) {
        if (nb1.idx == sp3_center) continue;
        for (const auto& nb2 : adj[sp3_center]) {
          if (nb2.idx == sp2_center || nb2.idx == nb1.idx) continue;
          if (share_ring(atom_info[nb1.idx], atom_info[nb2.idx])) {
            sp2_rs = nb1.idx;
            goto sp2sp3_found_rs;
          }
        }
      }
      sp2sp3_found_rs:
      // SP2 side tV: [ring_sharing?] + [remaining in adj order]
      {
        std::vector<size_t> sp2_tv;
        if (sp2_rs != SIZE_MAX)
          sp2_tv.push_back(sp2_rs);
        for (const auto& nb : adj[sp2_center])
          if (nb.idx != sp3_center &&
              std::find(sp2_tv.begin(), sp2_tv.end(), nb.idx) == sp2_tv.end())
            sp2_tv.push_back(nb.idx);
        int i_pos = -1;
        for (int k = 0; k < (int)sp2_tv.size(); ++k)
          if (sp2_tv[k] == sp2_term) { i_pos = k; break; }
        // tS3: SP2's non-center neighbors share a ring with each other
        bool is_ts3 = (sp2_rs == SIZE_MAX && sp2_tv.size() == 2 &&
                       share_ring(atom_info[sp2_tv[0]], atom_info[sp2_tv[1]]));
        static const double ts3_m[2][3] = {
          {150,-90,30}, {-30,90,-150}};
        static const double ts1_m[2][3] = {
          {0,120,-120}, {180,-60,60}};
        const auto& sp2sp3_m = is_ts3 ? ts3_m : ts1_m;
        int j_pos = compute_tv_position(sp3_center, sp2_center, sp3_term,
                                        TvMode::SP2SP3_SP3);
        if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 3)
          value = sp2sp3_m[i_pos][j_pos];
      }
      esd = 20.0;
      period = 6;
      int h_count = 0;
      for (const auto& nb : adj[sp3_center])
        if (nb.idx != sp2_center && cc.atoms[nb.idx].is_hydrogen())
          ++h_count;
      if (is_oxygen_column(cc.atoms[sp2_center].el) &&
          h_count >= 3 &&
          cc.atoms[sp3_term].is_hydrogen()) {
        value = -60.0;
        period = 3;
      }
    } else if (!lookup_found && sp2_2 && sp2_3) {
      // Non-aromatic SP2-SP2: use 2x2 matrix
      // Ring-sharing: {0,180; 180,0}. Non-ring: {180,0; 0,180}.
      auto sp2sp2_ring_sharing = [&](size_t c2, size_t c3) {
        for (const auto& nb1 : adj[c2]) {
          if (nb1.idx == c3) continue;
          for (const auto& nb2 : adj[c3]) {
            if (nb2.idx == c2 || nb2.idx == nb1.idx) continue;
            if (share_ring(atom_info[nb1.idx], atom_info[nb2.idx]))
              return true;
          }
        }
        return false;
      };
      bool has_ring_sharing = sp2sp2_ring_sharing(center2, center3);
      int side1_nb_count = 0;
      int side2_nb_count = 0;
      for (const auto& nb : adj[center2])
        if (nb.idx != center3)
          ++side1_nb_count;
      for (const auto& nb : adj[center3])
        if (nb.idx != center2)
          ++side2_nb_count;
      static const double ring_m[2][2] = {{0,180},{180,0}};
      static const double noring_m[2][2] = {{180,0},{0,180}};
      const auto& m22 = has_ring_sharing ? ring_m : noring_m;
      // SP2-SP2 neighbor ordering: for non-ring case with 2 neighbors per side,
      // AceDRG reorders by "H-only" check then connectivity tiebreaker.
      auto sp2sp2_tv_pos = [&](size_t center, size_t other, size_t target,
                               bool is_side1) -> int {
        std::vector<size_t> tv;
        // Ring-sharing atom first
        for (const auto& nb1 : adj[center]) {
          if (nb1.idx == other) continue;
          for (const auto& nb2 : adj[other]) {
            if (nb2.idx == center || nb2.idx == nb1.idx) continue;
            if (share_ring(atom_info[nb1.idx], atom_info[nb2.idx])) {
              tv.push_back(nb1.idx);
              goto done_rs;
            }
          }
        }
        done_rs:
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
              if (nb.idx != center && !cc.atoms[nb.idx].is_hydrogen())
                return false;
            return true;
          };
          bool h0 = is_h_only(tv[0]);
          bool h1 = is_h_only(tv[1]);
          if (is_side1) {
            if (h0 && !h1)
              std::swap(tv[0], tv[1]);
            else if ((int)adj[tv[0]].size() < (int)adj[tv[1]].size())
              std::swap(tv[0], tv[1]);
          } else {
            if (h0 && (!h1 || (int)adj[tv[0]].size() < (int)adj[tv[1]].size()))
              std::swap(tv[0], tv[1]);
            else if (!h0 && (int)adj[tv[0]].size() < (int)adj[tv[1]].size() && !h1)
              std::swap(tv[0], tv[1]);
          }
        }
        find_target:
        for (int i = 0; i < (int)tv.size(); ++i)
          if (tv[i] == target)
            return i;
        return -1;
      };
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
      size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
      int i_pos = sp2sp2_tv_pos(side1, side2, term1, true);
      int j_pos = sp2sp2_tv_pos(side2, side1, term2, false);
      if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 2)
        value = m22[i_pos][j_pos];
      bool end1_ring = share_ring(atom_info[term1], atom_info[side1]);
      bool end2_ring = share_ring(atom_info[term2], atom_info[side2]);
      esd = (end1_ring || end2_ring) ? 20.0 : 5.0;
      period = 2;
    }

    cc.rt.torsions.push_back({"auto",
                              {1, a1->id},
                              {1, cc.atoms[center2].id},
                              {1, cc.atoms[center3].id},
                              {1, a4->id},
                              value, esd, period});
  }

  if (sugar_mode) {
    for (const auto& ring_atoms : sugar_rings) {
      size_t o_idx = SIZE_MAX;
      for (size_t idx : ring_atoms)
        if (cc.atoms[idx].el == El::O) {
          o_idx = idx;
          break;
        }
      if (o_idx == SIZE_MAX)
        continue;
      std::set<size_t> ring_set(ring_atoms.begin(), ring_atoms.end());
      std::vector<size_t> ring_nbs;
      for (const auto& nb : adj[o_idx])
        if (ring_set.count(nb.idx))
          ring_nbs.push_back(nb.idx);
      if (ring_nbs.size() != 2)
        continue;
      size_t next = (cc.atoms[ring_nbs[0]].id < cc.atoms[ring_nbs[1]].id)
                    ? ring_nbs[0] : ring_nbs[1];
      std::vector<size_t> seq;
      seq.reserve(ring_atoms.size());
      seq.push_back(o_idx);
      seq.push_back(next);
      size_t prev = o_idx;
      size_t cur = next;
      while (seq.size() < ring_atoms.size()) {
        size_t chosen = SIZE_MAX;
        for (const auto& nb : adj[cur]) {
          if (!ring_set.count(nb.idx) || nb.idx == prev)
            continue;
          if (std::find(seq.begin(), seq.end(), nb.idx) != seq.end())
            continue;
          if (chosen == SIZE_MAX || cc.atoms[nb.idx].id < cc.atoms[chosen].id)
            chosen = nb.idx;
        }
        if (chosen == SIZE_MAX)
          break;
        seq.push_back(chosen);
        prev = cur;
        cur = chosen;
      }
      if (seq.size() != ring_atoms.size())
        continue;
      auto add_nu = [&](int nu_idx, size_t i1, size_t i2, size_t i3, size_t i4) {
        // AceDRG getTorsion() uses opposite sign convention than calculate_dihedral().
        double value = -deg(calculate_dihedral(cc.atoms[i1].xyz, cc.atoms[i2].xyz,
                                               cc.atoms[i3].xyz, cc.atoms[i4].xyz));
        cc.rt.torsions.push_back({cat("nu", nu_idx),
                                  {1, cc.atoms[i1].id},
                                  {1, cc.atoms[i2].id},
                                  {1, cc.atoms[i3].id},
                                  {1, cc.atoms[i4].id},
                                  value, 10.0, 3});
      };
      if (seq.size() == 6) {
        add_nu(0, seq[5], seq[0], seq[1], seq[2]);
        add_nu(1, seq[0], seq[1], seq[2], seq[3]);
        add_nu(2, seq[1], seq[2], seq[3], seq[4]);
        add_nu(3, seq[2], seq[3], seq[4], seq[5]);
        add_nu(4, seq[3], seq[4], seq[5], seq[0]);
        add_nu(5, seq[4], seq[5], seq[0], seq[1]);
      } else if (seq.size() == 5) {
        add_nu(0, seq[4], seq[0], seq[1], seq[2]);
        add_nu(1, seq[0], seq[1], seq[2], seq[3]);
        add_nu(2, seq[1], seq[2], seq[3], seq[4]);
        add_nu(3, seq[2], seq[3], seq[4], seq[0]);
        add_nu(4, seq[3], seq[4], seq[0], seq[1]);
      }
    }
  }
}

void add_chirality_if_missing(
    ChemComp& cc, const std::map<std::string, std::string>& atom_stereo,
    const std::vector<CodAtomInfo>& atom_info) {
  if (!cc.rt.chirs.empty())
    return;

  auto atom_index = cc.make_atom_index();
  auto adj = build_bond_adjacency(cc, atom_index);
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

    std::vector<size_t> non_h;
    std::vector<size_t> h;
    for (const auto& nb : adj[center]) {
      if (cc.atoms[nb.idx].el == El::H)
        h.push_back(nb.idx);
      else
        non_h.push_back(nb.idx);
    }
    if (non_h.size() < 3)
      continue;

    ChiralityType sign = stereo_sign(center);
    bool is_stereo_carbon = (sign != ChiralityType::Both && cc.atoms[center].el == El::C);
    if (is_stereo_carbon)
      stereo_centers.insert(center);

    if (is_stereo_carbon) {
      std::stable_sort(non_h.begin(), non_h.end(), [&](size_t a, size_t b) {
        return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
      });
      std::stable_sort(h.begin(), h.end(), [&](size_t a, size_t b) {
        return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
      });
    } else {
      std::stable_sort(non_h.begin(), non_h.end(), [&](size_t a, size_t b) {
        return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
      });
      std::stable_sort(h.begin(), h.end(), [&](size_t a, size_t b) {
        return chirality_priority(cc.atoms[a].el) < chirality_priority(cc.atoms[b].el);
      });
      sign = ChiralityType::Both;
    }

    std::vector<size_t> chosen;
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
    if (chosen.size() < 3)
      continue;

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
                           const std::vector<CodAtomInfo>& atom_info) {
  auto atom_index = cc.make_atom_index();
  auto adj = build_bond_adjacency(cc, atom_index);

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

  for (size_t idx = 0; idx < cc.atoms.size(); ++idx) {
    if (cc.atoms[idx].el != El::N)
      continue;
    if (atom_info[idx].hybrid != Hybridization::SP2)
      continue;
    if (atom_info[idx].min_ring_size > 0)
      continue;
    std::set<size_t> plane_idx;
    plane_idx.insert(idx);
    for (const auto& nb : adj[idx])
      plane_idx.insert(nb.idx);
    add_plane(plane_idx);
  }

  for (size_t idx = 0; idx < cc.atoms.size(); ++idx) {
    if (cc.atoms[idx].el != El::C)
      continue;
    if (atom_info[idx].hybrid != Hybridization::SP2)
      continue;
    if (atom_info[idx].min_ring_size > 0)
      continue;
    std::set<size_t> plane_idx;
    plane_idx.insert(idx);
    for (const auto& nb : adj[idx])
      plane_idx.insert(nb.idx);
    add_plane(plane_idx);
  }

  for (size_t idx = 0; idx < cc.atoms.size(); ++idx) {
    if (cc.atoms[idx].el != El::C)
      continue;
    std::vector<std::string> oxy;
    std::vector<std::string> other;
    for (const auto& nb : adj[idx]) {
      if (cc.atoms[nb.idx].el == El::O)
        oxy.push_back(cc.atoms[nb.idx].id);
      else
        other.push_back(cc.atoms[nb.idx].id);
    }
    if (oxy.size() < 2 || other.empty())
      continue;
    std::set<size_t> plane_idx;
    plane_idx.insert(idx);
    plane_idx.insert(atom_index[other.front()]);
    for (const std::string& o : oxy)
      plane_idx.insert(atom_index[o]);
    add_plane(plane_idx);
  }
}

void apply_metal_charge_corrections(ChemComp& cc) {
  if (cc.atoms.empty())
    return;
  auto atom_index = cc.make_atom_index();
  auto adj = build_bond_adjacency(cc, atom_index);
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    const Element& el = cc.atoms[i].el;
    if (el.is_metal() || el == El::H)
      continue;
    // check if atom has any metal neighbor
    bool has_metal = false;
    bool has_non_metal_heavy = false;
    for (const auto& nb : adj[i]) {
      if (cc.atoms[nb.idx].el.is_metal())
        has_metal = true;
      else if (!cc.atoms[nb.idx].is_hydrogen())
        has_non_metal_heavy = true;
    }
    if (!has_metal || !has_non_metal_heavy)
      continue;
    int expected_valence = 0;
    if (el == El::O) expected_valence = 2;
    else if (el == El::N) expected_valence = 3;
    else if (el == El::S) expected_valence = 2;
    else if (el == El::C) expected_valence = 4;
    else if (el == El::P) expected_valence = 3;
    else continue;
    float sum_bo = 0.0f;
    for (const auto& nb : adj[i]) {
      if (!cc.atoms[nb.idx].el.is_metal())
        sum_bo += order_of_bond_type(nb.type);
    }
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

  adjust_phosphate_group(cc);
  adjust_sulfate_group(cc);
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
      add_torsions_from_bonds_if_missing(cc, tables, atom_info, atom_stereo);
      add_chirality_if_missing(cc, atom_stereo, atom_info);
      add_planes_if_missing(cc, atom_info);
    }
  } else {
    if (added_h3 && !only_bonds)
      sync_n_terminal_h3_angles(cc);
  }

  tables.assign_ccp4_types(cc);
}

} // namespace gemmi
