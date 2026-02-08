// Copyright 2025 Global Phasing Ltd.
//
// Restraint generation helpers for prepare_chemcomp():
// chemical adjustments, torsions, chirality, planes, H naming.

#include "gemmi/ace_cc.hpp"
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

std::map<std::string, size_t> make_atom_index(const ChemComp& cc) {
  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;
  return atom_index;
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
        new_h = "H" + std::to_string(i);
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
  std::string cand = h_root + std::to_string(start);
  if (used_names.find(cand) == used_names.end()) {
    used_names.insert(cand);
    return cand;
  }

  for (int n = 2; n < 10000; ++n) {
    cand = h_root + std::to_string(n);
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

  auto atom_index = make_atom_index(cc);

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
    auto& nbs = neighbors[center];
    if (nbs.size() < 2)
      continue;
    for (size_t i = 0; i + 1 < nbs.size(); ++i) {
      for (size_t j = i + 1; j < nbs.size(); ++j) {
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

int torsion_neighbor_priority(Element el) {
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
    const std::vector<bool>& ring_like,
    bool bond_aromatic,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> candidates;
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx)
      candidates.push_back(nb.idx);
  if (candidates.empty())
    return nullptr;

  bool prefer_ring = ring_like[center_idx];
  if (prefer_ring) {
    bool share_other = bond_aromatic;
    if (share_other) {
      std::vector<size_t> non_ring_non_h;
      std::vector<size_t> ring_non_h;
      for (size_t idx : candidates) {
        if (cc.atoms[idx].is_hydrogen())
          continue;
        if (share_ring(atom_info[center_idx], atom_info[idx]) || ring_like[idx])
          ring_non_h.push_back(idx);
        else
          non_ring_non_h.push_back(idx);
      }

      if (non_ring_non_h.size() == 1 && ring_non_h.size() == 1) {
        Element el = cc.atoms[non_ring_non_h.front()].el;
        if (el == El::O || el == El::N || el == El::S || el == El::P) {
          candidates.swap(ring_non_h);
        } else if (cc.atoms[exclude_idx].id < cc.atoms[ring_non_h.front()].id) {
          candidates.swap(non_ring_non_h);
        } else {
          candidates.swap(ring_non_h);
        }
      } else if (!non_ring_non_h.empty()) {
        candidates.swap(non_ring_non_h);
      } else if (!ring_non_h.empty()) {
        candidates.swap(ring_non_h);
      }
    } else {
      std::vector<size_t> ring_candidates;
      for (size_t idx : candidates)
        if (share_ring(atom_info[center_idx], atom_info[idx]))
          ring_candidates.push_back(idx);
      if (!ring_candidates.empty())
        candidates.swap(ring_candidates);
    }
  }

  bool has_non_h = false;
  for (size_t idx : candidates)
    if (!cc.atoms[idx].is_hydrogen())
      has_non_h = true;

  auto better_h = [&](size_t a, size_t b) {
    const std::string& na = cc.atoms[a].id;
    const std::string& nb = cc.atoms[b].id;
    bool a_is_h = (na == "H");
    bool b_is_h = (nb == "H");
    if (a_is_h != b_is_h)
      return a_is_h;
    if (!a_is_h && !b_is_h)
      return na > nb;
    return na < nb;
  };

  size_t best = candidates.front();
  for (size_t idx : candidates) {
    if (has_non_h && cc.atoms[idx].is_hydrogen())
      continue;
    if (has_non_h && cc.atoms[best].is_hydrogen())
      best = idx;
    if (cc.atoms[idx].is_hydrogen()) {
      if (better_h(idx, best))
        best = idx;
      continue;
    }
    if (cc.atoms[best].is_hydrogen()) {
      best = idx;
      continue;
    }
    if (prefer_ring) {
      int ring_idx = share_ring(atom_info[center_idx], atom_info[idx])
                         ? shared_ring_size(atom_info[center_idx], atom_info[idx])
                         : 999;
      int ring_best = share_ring(atom_info[center_idx], atom_info[best])
                          ? shared_ring_size(atom_info[center_idx], atom_info[best])
                          : 999;
      if (ring_idx != ring_best) {
        if (ring_idx < ring_best)
          best = idx;
        continue;
      }
    }
    if (bond_aromatic) {
      bool idx_arom_bond = false;
      for (const auto& nb : adj[center_idx]) {
        if (nb.idx == idx && nb.aromatic) {
          idx_arom_bond = true;
          break;
        }
      }
      bool best_arom_bond = false;
      for (const auto& nb : adj[center_idx]) {
        if (nb.idx == best && nb.aromatic) {
          best_arom_bond = true;
          break;
        }
      }
      if (idx_arom_bond != best_arom_bond) {
        if (!idx_arom_bond)
          best = idx;
        continue;
      }
    }
    int p_idx = torsion_neighbor_priority(cc.atoms[idx].el);
    int p_best = torsion_neighbor_priority(cc.atoms[best].el);
    if (p_idx != p_best) {
      if (p_idx < p_best)
        best = idx;
      continue;
    }
    if (cc.atoms[idx].el == El::C && cc.atoms[best].el == El::C) {
      bool c_idx = is_carbonyl_carbon(idx, cc, adj);
      bool c_best = is_carbonyl_carbon(best, cc, adj);
      if (c_idx != c_best) {
        if (c_idx)
          best = idx;
        continue;
      }
    }
    if (cc.atoms[idx].id < cc.atoms[best].id)
      best = idx;
  }

  return &cc.atoms[best];
}

const ChemComp::Atom* pick_torsion_neighbor_in_ring(
    const ChemComp& cc,
    const std::vector<std::vector<NeighborBond>>& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t other_center_idx,
    size_t exclude_idx) {
  std::vector<size_t> ring_candidates;
  for (const auto& nb : adj[center_idx]) {
    if (nb.idx == exclude_idx)
      continue;
    if (shared_ring_size(atom_info[nb.idx], atom_info[other_center_idx]) > 0)
      ring_candidates.push_back(nb.idx);
  }
  if (ring_candidates.empty())
    return nullptr;

  bool has_non_h = false;
  for (size_t idx : ring_candidates)
    if (!cc.atoms[idx].is_hydrogen())
      has_non_h = true;

  size_t best = ring_candidates.front();
  for (size_t idx : ring_candidates) {
    if (has_non_h && cc.atoms[idx].is_hydrogen())
      continue;
    if (has_non_h && cc.atoms[best].is_hydrogen())
      best = idx;
    if (cc.atoms[idx].is_hydrogen()) {
      if (cc.atoms[idx].id > cc.atoms[best].id)
        best = idx;
      continue;
    }
    if (cc.atoms[best].is_hydrogen()) {
      best = idx;
      continue;
    }
    int ring_idx = shared_ring_size(atom_info[center_idx], atom_info[idx]);
    int ring_best = shared_ring_size(atom_info[center_idx], atom_info[best]);
    if (ring_idx != ring_best) {
      if (ring_idx < ring_best)
        best = idx;
      continue;
    }
    int p_idx = torsion_neighbor_priority(cc.atoms[idx].el);
    int p_best = torsion_neighbor_priority(cc.atoms[best].el);
    if (p_idx != p_best) {
      if (p_idx < p_best)
        best = idx;
      continue;
    }
    if (cc.atoms[idx].id < cc.atoms[best].id)
      best = idx;
  }

  return &cc.atoms[best];
}

const ChemComp::Atom* pick_aromatic_ring_neighbor(
    const ChemComp& cc,
    const std::vector<std::vector<NeighborBond>>& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> nonring_c;
  std::vector<size_t> nonring_other;
  std::vector<size_t> ring_nb;
  std::vector<size_t> h;
  for (const auto& nb : adj[center_idx]) {
    if (nb.idx == exclude_idx)
      continue;
    if (cc.atoms[nb.idx].is_hydrogen()) {
      h.push_back(nb.idx);
    } else if (share_ring(atom_info[center_idx], atom_info[nb.idx])) {
      ring_nb.push_back(nb.idx);
    } else if (cc.atoms[nb.idx].el == El::C) {
      nonring_c.push_back(nb.idx);
    } else {
      nonring_other.push_back(nb.idx);
    }
  }
  auto sort_ids = [&](std::vector<size_t>& v) {
    std::sort(v.begin(), v.end(),
              [&](size_t a, size_t b) { return cc.atoms[a].id < cc.atoms[b].id; });
  };
  sort_ids(nonring_c);
  sort_ids(nonring_other);
  sort_ids(ring_nb);
  sort_ids(h);
  if (cc.atoms[center_idx].id == "CG" &&
      cc.atoms[exclude_idx].id == "CD2") {
    for (size_t idx : ring_nb)
      if (cc.atoms[idx].id == "CD1")
        return &cc.atoms[idx];
  }
  if (!nonring_c.empty())
    return &cc.atoms[nonring_c.front()];
  if (!ring_nb.empty())
    return &cc.atoms[ring_nb.front()];
  if (!nonring_other.empty())
    return &cc.atoms[nonring_other.front()];
  if (!h.empty())
    return &cc.atoms[h.front()];
  return nullptr;
}

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables,
                                        const std::vector<CodAtomInfo>& atom_info) {
  if (!cc.rt.torsions.empty())
    return;

  auto atom_index = make_atom_index(cc);
  auto adj = build_bond_adjacency(cc, atom_index);
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

  for (const auto& bond : cc.rt.bonds) {
    auto it1 = atom_index.find(bond.id1.atom);
    auto it2 = atom_index.find(bond.id2.atom);
    if (it1 == atom_index.end() || it2 == atom_index.end())
      continue;
    size_t idx1 = it1->second;
    size_t idx2 = it2->second;

    bool pep_seeded = false;
    TorsionEntry pep_entry;
    bool pep_forward = true;
    size_t pep_a1 = 0;
    size_t pep_a4 = 0;
    int pep_priority = 9999;
    if (adj[idx1].size() > 1 && adj[idx2].size() > 1) {
      for (const auto& nb1 : adj[idx1]) {
        if (nb1.idx == idx2)
          continue;
        for (const auto& nb2 : adj[idx2]) {
          if (nb2.idx == idx1)
            continue;
          const std::string& a1 = cc.atoms[nb1.idx].id;
          const std::string& a2 = cc.atoms[idx1].id;
          const std::string& a3 = cc.atoms[idx2].id;
          const std::string& a4 = cc.atoms[nb2.idx].id;
          TorsionEntry entry;
          bool forward = false;
          if (tables.lookup_pep_tors(a1, a2, a3, a4, entry)) {
            forward = true;
          } else if (tables.lookup_pep_tors(a4, a3, a2, a1, entry)) {
            forward = false;
          } else {
            continue;
          }
          if (entry.priority < pep_priority) {
            pep_priority = entry.priority;
            pep_entry = entry;
            pep_forward = forward;
            pep_a1 = nb1.idx;
            pep_a4 = nb2.idx;
            pep_seeded = true;
          }
        }
      }
    }

    if (pep_seeded) {
      std::string id_lower = to_lower(pep_entry.id);
      double value = pep_entry.value;
      int period = pep_entry.period;
      double esd = 10.0;
      if (id_lower.find("const") != std::string::npos)
        esd = 0.0;
      else if (id_lower.find("hh") != std::string::npos)
        esd = 5.0;

      const std::string& a1 = cc.atoms[pep_a1].id;
      const std::string& a2 = cc.atoms[idx1].id;
      const std::string& a3 = cc.atoms[idx2].id;
      const std::string& a4 = cc.atoms[pep_a4].id;
      auto is_chi2_combo = [&](const std::string& x1, const std::string& x2,
                               const std::string& x3, const std::string& x4) {
        return x1 == "CA" && x2 == "CB" && x3 == "CG" &&
               (x4 == "CD" || x4 == "CD1" || x4 == "CD2");
      };
      bool chi2_match = is_chi2_combo(a1, a2, a3, a4) ||
                        is_chi2_combo(a4, a3, a2, a1);
      int sp1 = atom_info[idx1].bonding_idx;
      int sp2 = atom_info[idx2].bonding_idx;
      if (chi2_match && ((sp1 == 2 && sp2 == 3) || (sp1 == 3 && sp2 == 2))) {
        period = 6;
        value = 90.0;
      }

      if (pep_forward) {
        cc.rt.torsions.push_back({pep_entry.id,
                                  {1, a1},
                                  {1, a2},
                                  {1, a3},
                                  {1, a4},
                                  value, esd, period});
      } else {
        cc.rt.torsions.push_back({pep_entry.id,
                                  {1, a4},
                                  {1, a3},
                                  {1, a2},
                                  {1, a1},
                                  value, esd, period});
      }
      continue;
    }

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
      a1 = pick_torsion_neighbor(cc, adj, atom_info, aromatic_like,
                                 bond_aromatic, center2, center3);
      a4 = pick_torsion_neighbor(cc, adj, atom_info, aromatic_like,
                                 bond_aromatic, center3, center2);
    }
    ring_size = shared_ring_size(atom_info[center2], atom_info[center3]);
    if (ring_size > 0) {
      const ChemComp::Atom* ring_a1 =
          pick_torsion_neighbor_in_ring(cc, adj, atom_info, center2, center3, center3);
      const ChemComp::Atom* ring_a4 =
          pick_torsion_neighbor_in_ring(cc, adj, atom_info, center3, center2, center2);
      auto it_a1_idx = atom_index.find(a1->id);
      auto it_a4_idx = atom_index.find(a4->id);
      if (ring_a1 && it_a1_idx != atom_index.end() &&
          aromatic_like[it_a1_idx->second])
        a1 = ring_a1;
      if (ring_a4 && it_a4_idx != atom_index.end() &&
          aromatic_like[it_a4_idx->second])
        a4 = ring_a4;
    }
    if (!a1 || !a4)
      continue;

    const CodAtomInfo& h2 = atom_info[center2];
    const CodAtomInfo& h3 = atom_info[center3];
    bool sp3_2 = (h2.hybrid == Hybridization::SP3);
    bool sp3_3 = (h3.hybrid == Hybridization::SP3);
    bool sp2_2 = (h2.hybrid == Hybridization::SP2);
    bool sp2_3 = (h3.hybrid == Hybridization::SP2);

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
    bool a1_arom = (it_a1 != atom_index.end() &&
                    atom_info[it_a1->second].is_aromatic);
    bool a4_arom = (it_a4 != atom_index.end() &&
                    atom_info[it_a4->second].is_aromatic);
    bool phospho_torsion =
        ((cc.atoms[center2].el == El::O && cc.atoms[center3].el == El::P) ||
         (cc.atoms[center2].el == El::P && cc.atoms[center3].el == El::O)) &&
        ((a1->el == El::P && a4->el == El::O) ||
         (a1->el == El::O && a4->el == El::P));
    bool outer_non_ring = false;
    if (it_a1 != atom_index.end() && it_a4 != atom_index.end())
      outer_non_ring = !aromatic_like[it_a1->second] || !aromatic_like[it_a4->second];

    bool chi2_aromatic =
        (a1->id == "CA" &&
         ((cc.atoms[center2].id == "CB" && cc.atoms[center3].id == "CG") ||
          (cc.atoms[center3].id == "CB" && cc.atoms[center2].id == "CG")) &&
         (a4->id == "CD1" || a4->id == "CD2"));
    bool phenol_oh =
        ((cc.atoms[center2].id == "CZ" && cc.atoms[center3].id == "OH") ||
         (cc.atoms[center3].id == "CZ" && cc.atoms[center2].id == "OH")) &&
        a4->id == "HH";

    if (phospho_torsion) {
      value = 60.0;
    } else if (chi2_aromatic) {
      value = 90.0;
      esd = 20.0;
      period = 6;
    } else if (phenol_oh) {
      value = 0.0;
      esd = 5.0;
      period = 2;
    } else if (bond_aromatic && h2.is_aromatic && h3.is_aromatic && outer_non_ring) {
      value = 0.0;
      esd = 0.0;
      period = 1;
    } else if (!lookup_found && ring_size > 0 && h2.is_aromatic && h3.is_aromatic &&
               a1_arom && a4_arom) {
      value = 0.0;
      esd = 0.0;
      period = 2;
    } else if (!lookup_found && ring_size == 5 && sp3_2 && sp3_3) {
      bool ca_n_bond =
          (cc.atoms[center2].id == "CA" && cc.atoms[center3].el == El::N) ||
          (cc.atoms[center3].id == "CA" && cc.atoms[center2].el == El::N);
      if (ca_n_bond) {
        value = 60.0;
      } else if (a1->el == El::N || a4->el == El::N) {
        value = -60.0;
      }
    } else if (!lookup_found && ((sp2_2 && sp3_3) || (sp3_2 && sp2_3))) {
      value = 0.0;
      esd = 20.0;
      period = 6;
    } else if (!lookup_found && sp2_2 && sp2_3) {
      value = 180.0;
      esd = 20.0;
      period = 2;
    }

    cc.rt.torsions.push_back({"auto",
                              {1, a1->id},
                              {1, cc.atoms[center2].id},
                              {1, cc.atoms[center3].id},
                              {1, a4->id},
                              value, esd, period});
  }
}

void add_chirality_if_missing(
    ChemComp& cc, const std::map<std::string, std::string>& atom_stereo,
    const std::vector<CodAtomInfo>& atom_info) {
  if (!cc.rt.chirs.empty())
    return;

  if (!atom_stereo.empty()) {
    auto atom_index = make_atom_index(cc);
    auto adj = build_bond_adjacency(cc, atom_index);

    for (const auto& entry : atom_stereo) {
      char stereo = lower(entry.second[0]);
      if (stereo != 'r' && stereo != 's')
        continue;
      auto it = atom_index.find(entry.first);
      if (it == atom_index.end())
        continue;
      size_t center = it->second;
      if (cc.atoms[center].el != El::C)
        continue;

      auto by_priority = [&](size_t a, size_t b) {
        int pa = chirality_priority(cc.atoms[a].el);
        int pb = chirality_priority(cc.atoms[b].el);
        if (pa != pb)
          return pa < pb;
        bool center_in_ring = (atom_info[center].min_ring_size > 0);
        if (center_in_ring)
          return cc.atoms[a].id > cc.atoms[b].id;
        return cc.atoms[a].id < cc.atoms[b].id;
      };

      std::vector<size_t> non_h;
      std::vector<size_t> h;
      for (const auto& nb : adj[center]) {
        if (cc.atoms[nb.idx].el == El::H)
          h.push_back(nb.idx);
        else
          non_h.push_back(nb.idx);
      }
      std::sort(non_h.begin(), non_h.end(), by_priority);
      std::sort(h.begin(), h.end(), by_priority);

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

      ChiralityType sign = (stereo == 's') ? ChiralityType::Positive
                                           : ChiralityType::Negative;
      cc.rt.chirs.push_back({{1, cc.atoms[center].id},
                             {1, cc.atoms[chosen[0]].id},
                             {1, cc.atoms[chosen[1]].id},
                             {1, cc.atoms[chosen[2]].id},
                             sign});
    }
  }

  if (!cc.rt.chirs.empty())
    return;

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

  cc.rt.chirs.push_back({{1, "CA"}, {1, "N"}, {1, "C"}, {1, "CB"}, sign});
}

void add_planes_if_missing(ChemComp& cc,
                           const std::vector<CodAtomInfo>& atom_info) {
  if (!cc.rt.planes.empty())
    return;

  auto atom_index = make_atom_index(cc);
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
  auto atom_index = make_atom_index(cc);
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
                      const std::map<std::string, std::string>& atom_stereo) {
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
  int missing_angles = count_missing_values(cc.rt.angles);
  bool need_fill = (missing_bonds > 0 || missing_angles > 0);

  if (need_fill) {
    tables.fill_restraints(cc);
    if (added_h3)
      sync_n_terminal_h3_angles(cc);
    std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);
    add_torsions_from_bonds_if_missing(cc, tables, atom_info);
    add_chirality_if_missing(cc, atom_stereo, atom_info);
    add_planes_if_missing(cc, atom_info);
  } else {
    if (added_h3)
      sync_n_terminal_h3_angles(cc);
  }

  tables.assign_ccp4_types(cc);
}

} // namespace gemmi
