// Copyright 2025 Global Phasing Ltd.
//
// AcedrgTables - COD/CSD-based atom classification and restraint value lookup
// Port of AceDRG codClassify system to gemmi.

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <numeric>
#include <sstream>
#include <iostream>
#include <set>
#include "gemmi/acedrg_tables.hpp"
#include "gemmi/elem.hpp"
#include "gemmi/fail.hpp"
#include "gemmi/util.hpp"
#include "gemmi/read_cif.hpp"

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
  cc.rt.bonds.erase(std::remove_if(cc.rt.bonds.begin(), cc.rt.bonds.end(),
                                  [&](const Restraints::Bond& b) {
                                    return is_id(b.id1) || is_id(b.id2);
                                  }),
                    cc.rt.bonds.end());
  cc.rt.angles.erase(std::remove_if(cc.rt.angles.begin(), cc.rt.angles.end(),
                                   [&](const Restraints::Angle& a) {
                                     return is_id(a.id1) || is_id(a.id2) || is_id(a.id3);
                                   }),
                     cc.rt.angles.end());
  cc.rt.torsions.erase(std::remove_if(cc.rt.torsions.begin(), cc.rt.torsions.end(),
                                     [&](const Restraints::Torsion& t) {
                                       return is_id(t.id1) || is_id(t.id2) ||
                                              is_id(t.id3) || is_id(t.id4);
                                     }),
                       cc.rt.torsions.end());
  cc.rt.chirs.erase(std::remove_if(cc.rt.chirs.begin(), cc.rt.chirs.end(),
                                  [&](const Restraints::Chirality& c) {
                                    return is_id(c.id_ctr) || is_id(c.id1) ||
                                           is_id(c.id2) || is_id(c.id3);
                                  }),
                    cc.rt.chirs.end());
  for (auto it = cc.rt.planes.begin(); it != cc.rt.planes.end(); ) {
    auto& ids = it->ids;
    ids.erase(std::remove_if(ids.begin(), ids.end(),
                             [&](const Restraints::AtomId& id) { return is_id(id); }),
              ids.end());
    if (ids.empty())
      it = cc.rt.planes.erase(it);
    else
      ++it;
  }
  cc.atoms.erase(std::remove_if(cc.atoms.begin(), cc.atoms.end(),
                               [&](const ChemComp::Atom& a) { return a.id == atom_id; }),
                 cc.atoms.end());
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
  auto atom_index = make_atom_index(cc);
  auto neighbors = make_neighbor_names(cc);

  std::set<std::string> po4_phosphorus;
  for (const auto& atom : cc.atoms) {
    if (atom.el != El::P)
      continue;
    const auto& nb = neighbors[atom.id];
    int o_count = 0;
    bool has_h = false;
    for (const std::string& nid : nb) {
      auto it = atom_index.find(nid);
      if (it != atom_index.end()) {
        if (cc.atoms[it->second].el == El::O)
          ++o_count;
        else if (cc.atoms[it->second].el == El::H)
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
      auto it = atom_index.find(nid);
      if (it != atom_index.end() && cc.atoms[it->second].el == El::H)
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
  auto atom_index = make_atom_index(cc);
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
      auto it = atom_index.find(nid);
      if (it != atom_index.end() && cc.atoms[it->second].el == El::O)
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
      auto it = atom_index.find(nid);
      if (it != atom_index.end() && cc.atoms[it->second].el == El::H)
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
  auto atom_index = make_atom_index(cc);
  auto neighbors = make_neighbor_names(cc);

  for (auto& atom : cc.atoms) {
    if (atom.el != El::P)
      continue;
    int f_count = 0;
    bool has_h = false;
    for (const std::string& nb : neighbors[atom.id]) {
      auto it = atom_index.find(nb);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
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
  auto atom_index = make_atom_index(cc);
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
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
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
        auto it = atom_index.find(ac_nb);
        if (it != atom_index.end() && cc.atoms[it->second].el != El::H) {
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
          auto it = atom_index.find(other);
          if (it != atom_index.end() && cc.atoms[it->second].el == El::C) {
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
      auto it = atom_index.find(nid);
      if (it != atom_index.end() && cc.atoms[it->second].el == El::H) {
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

  auto atom_index = make_atom_index(cc);
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
    auto itc = atom_index.find(c_id);
    if (itc == atom_index.end() || cc.atoms[itc->second].el != El::C)
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
        auto it = atom_index.find(other);
        if (it == atom_index.end())
          continue;
        Element el = cc.atoms[it->second].el;
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
      auto it = atom_index.find(nid);
      if (it != atom_index.end() && cc.atoms[it->second].el == El::N)
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
          auto it = atom_index.find(onb);
          if (it == atom_index.end())
            continue;
          if (cc.atoms[it->second].el == El::H)
            continue;
          if (cc.atoms[it->second].el != El::C &&
              cc.atoms[it->second].el != El::Si &&
              cc.atoms[it->second].el != El::Ge) {
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
        auto it = atom_index.find(nid);
        if (it == atom_index.end())
          continue;
        if (cc.atoms[it->second].el == El::H)
          h_neighbors.push_back(nid);
        else if (nid != atom.id)
          has_substituent = true;
      }

      if (has_substituent)
        continue;

      if (h_neighbors.size() == 1) {
        auto n_it = atom_index.find(n_id);
        if (n_it == atom_index.end())
          continue;

        cc.atoms[n_it->second].charge = 1.0f;

        std::string new_h_id = acedrg_h_name(used_names, n_id, h_neighbors);
        if (new_h_id.empty())
          continue;

        cc.atoms.push_back(ChemComp::Atom{new_h_id, "", El::H, 0.0f, "H", "", Position()});
        cc.rt.bonds.push_back({{1, n_id}, {1, new_h_id}, BondType::Single, false,
                              NAN, NAN, NAN, NAN});
        atom_index[new_h_id] = cc.atoms.size() - 1;
        neighbors[n_id].push_back(new_h_id);
        neighbors[new_h_id].push_back(n_id);
      }
    }
  }
}

void adjust_amino_ter_amine(ChemComp& cc, std::set<std::string>& used_names) {
  auto atom_index = make_atom_index(cc);
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
      auto it = atom_index.find(nb);
      if (it == atom_index.end())
        continue;
      if (cc.atoms[it->second].el == El::H)
        h_ids.push_back(nb);
      else if (cc.atoms[it->second].el == El::C)
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
      auto it = atom_index.find(c2_id);
      if (it == atom_index.end() || cc.atoms[it->second].el != El::C)
        continue;

      bool has_carbonyl = false, has_amide_n = false;
      for (const std::string& c2_nb : neighbors[c2_id]) {
        auto it2 = atom_index.find(c2_nb);
        if (it2 == atom_index.end())
          continue;
        Element el = cc.atoms[it2->second].el;
        if (el == El::O && get_bond_type(c2_id, c2_nb) == BondType::Double)
          has_carbonyl = true;
        if (el == El::N && c2_nb != n1.id) {
          for (const std::string& amide_n_nb : neighbors[c2_nb]) {
            if (amide_n_nb == c2_id)
              continue;
            auto it3 = atom_index.find(amide_n_nb);
            if (it3 != atom_index.end() && cc.atoms[it3->second].el != El::H) {
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
  auto atom_index = make_atom_index(cc);
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
      auto it = atom_index.find(nb);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
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
      auto it = atom_index.find(c_neighbor);
      if (it == atom_index.end())
        continue;
      if (cc.atoms[it->second].el != El::C)
        continue;
      int o_count = 0;
      bool has_double_o = false;
      for (const auto& bond : cc.rt.bonds) {
        if (bond.id1.atom != c_neighbor && bond.id2.atom != c_neighbor)
          continue;
        std::string other = (bond.id1.atom == c_neighbor) ? bond.id2.atom : bond.id1.atom;
        auto other_it = atom_index.find(other);
        if (other_it != atom_index.end() && cc.atoms[other_it->second].el == El::O) {
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
      std::string id_lower = pep_entry.id;
      for (char& c : id_lower)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
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
      char stereo = static_cast<char>(std::tolower(static_cast<unsigned char>(entry.second[0])));
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

  std::string type = cc.type_or_group;
  for (char& c : type)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
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
    if (std::find(plane_sets.begin(), plane_sets.end(), ids) != plane_sets.end())
      return;
    Restraints::Plane plane;
    plane.label = "plan-" + std::to_string(cc.rt.planes.size() + 1);
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

}  // namespace

// ============================================================================
// Below: AcedrgTables member functions — table loading, atom classification,
// bond/angle search, CCP4 type assignment.
// Above anonymous namespace: restraint generation helpers — torsions,
// chirality, planes, H naming, chemical adjustments, prepare_chemcomp().
// These two parts are candidates for splitting into separate source files.
// ============================================================================

const char* hybridization_to_string(Hybridization h) {
  switch (h) {
    case Hybridization::SP1: return "SP1";
    case Hybridization::SP2: return "SP2";
    case Hybridization::SP3: return "SP3";
    case Hybridization::SPD5: return "SPD5";
    case Hybridization::SPD6: return "SPD6";
    case Hybridization::SPD7: return "SPD7";
    case Hybridization::SPD8: return "SPD8";
    case Hybridization::SP_NON: default: return "SP-NON";
  }
}

Hybridization hybridization_from_string(const std::string& s) {
  if (s == "SP1") return Hybridization::SP1;
  if (s == "SP2") return Hybridization::SP2;
  if (s == "SP3") return Hybridization::SP3;
  if (s == "SPD5") return Hybridization::SPD5;
  if (s == "SPD6") return Hybridization::SPD6;
  if (s == "SPD7") return Hybridization::SPD7;
  if (s == "SPD8") return Hybridization::SPD8;
  return Hybridization::SP_NON;
}


void AcedrgTables::load_tables(const std::string& tables_dir) {
  tables_dir_ = tables_dir;

  // Load hash code mapping
  load_hash_codes(tables_dir + "/allOrgLinkedHashCode.table");

  // Load HRS (summary) tables
  load_bond_hrs(tables_dir + "/allOrgBondsHRS.table");
  //load_angle_hrs(tables_dir + "/allOrgAnglesHRS.table");  // temporarily disabled

  // Load element+hybridization fallback
  load_en_bonds(tables_dir + "/allOrgBondEN.table");

  // Load protonated hydrogen distances
  load_prot_hydr_dists(tables_dir + "/prot_hydr_dists.table");

  // Load metal tables
  load_metal_tables(tables_dir);
  load_covalent_radii(tables_dir + "/radii.table");

  // Load CCP4 energetic library bonds (for AceDRG fallback)
  load_ccp4_bonds(tables_dir + "/ener_lib.cif");

  // Load detailed indexed tables
  load_atom_type_codes(tables_dir + "/allAtomTypesFromMolsCoded.list");
  load_bond_index(tables_dir + "/allOrgBondTables/bond_idx.table");
  load_bond_tables(tables_dir + "/allOrgBondTables");
  //load_angle_index(tables_dir + "/allOrgAngleTables/angle_idx.table");  // temporarily disabled
  //load_angle_tables(tables_dir + "/allOrgAngleTables");  // temporarily disabled
  load_pep_tors(tables_dir + "/pep_tors.table");

  tables_loaded_ = true;
}

void AcedrgTables::load_hash_codes(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash_code, linked_hash;
    std::string footprint;

    if (iss >> hash_code >> footprint >> linked_hash) {
      digit_keys_[hash_code] = footprint;
      linked_hash_[hash_code] = linked_hash;
    }
  }
}

void AcedrgTables::load_bond_hrs(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    fail("Cannot open bond HRS table: ", path);

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash1, hash2;
    std::string hybrid_pair, in_ring;
    double value, sigma;
    int count;

    if (iss >> hash1 >> hash2 >> hybrid_pair >> in_ring >> value >> sigma >> count) {
      BondHRSKey key;
      key.hash1 = std::min(hash1, hash2);
      key.hash2 = std::max(hash1, hash2);
      key.hybrid_pair = hybrid_pair;
      key.in_ring = (in_ring == "Y" || in_ring == "y") ? "Y" : "N";
      ValueStats vs(value, sigma, count);
      bond_hrs_[key] = vs;
      // AceDRG levels 9-11 use the HRS table hierarchy.
      bond_hasp_2d_[key.hash1][key.hash2][key.hybrid_pair][key.in_ring].push_back(vs);
      bond_hasp_1d_[key.hash1][key.hash2][key.hybrid_pair].push_back(vs);
      bond_hasp_0d_[key.hash1][key.hash2].push_back(vs);
    }
  }
}

void AcedrgTables::load_angle_hrs(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    fail("Cannot open angle HRS table: ", path);

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    int hash1, hash2, hash3;
    std::string value_key;
    std::string a1_cod, a2_cod, a3_cod;
    double value1, sigma1;
    int count1;
    double value2, sigma2;
    int count2;

    if (iss >> hash1 >> hash2 >> hash3 >> value_key
        >> a1_cod >> a2_cod >> a3_cod
        >> value1 >> sigma1 >> count1
        >> value2 >> sigma2 >> count2) {
      AngleHRSKey key;
      // For angles, center is always hash2, but we canonicalize flanking atoms
      bool swap_flanks = hash1 > hash3;
      key.hash1 = swap_flanks ? hash3 : hash1;
      key.hash2 = hash2;
      key.hash3 = swap_flanks ? hash1 : hash3;
      if (swap_flanks) {
        size_t colon = value_key.find(':');
        if (colon != std::string::npos) {
          std::string ring_part = value_key.substr(0, colon);
          std::string hybrid_part = value_key.substr(colon + 1);
          std::vector<std::string> parts = split_str(hybrid_part, '_');
          if (parts.size() == 3) {
            std::swap(parts[0], parts[2]);
            hybrid_part = parts[0] + "_" + parts[1] + "_" + parts[2];
          }
          value_key = ring_part + ":" + hybrid_part;
        }
      }
      key.value_key = value_key;
      angle_hrs_[key] = ValueStats(value1, sigma1, count1);
    }
  }
}

void AcedrgTables::load_en_bonds(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string elem1, sp1, elem2, sp2;
    double value, sigma;
    int count;

    if (iss >> elem1 >> sp1 >> elem2 >> sp2 >> value >> sigma >> count) {
      // Canonicalize order
      if (elem1 > elem2 || (elem1 == elem2 && sp1 > sp2)) {
        std::swap(elem1, elem2);
        std::swap(sp1, sp2);
      }
      en_bonds_[elem1][sp1][elem2][sp2].emplace_back(value, sigma, count);
    }
  }
}

void AcedrgTables::load_prot_hydr_dists(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Remove stray commas that may appear in numeric fields
    for (char& c : line)
      if (c == ',') c = ' ';

    std::istringstream iss(line);
    std::string h_elem, heavy_elem, type_key;
    double nucleus_val, nucleus_sigma, v1, s1, electron_val, electron_sigma;

    // Format: H  C  H_sp3_C  1.092  0.010  1.093107  0.003875  0.988486  0.005251  ...
    // Column 4-5: nucleus distance (ener_lib-like)
    // Column 6-7: refined nucleus distance
    // Column 8-9: electron distance (X-ray) - this is what acedrg uses for value_dist
    if (iss >> h_elem >> heavy_elem >> type_key >> nucleus_val >> nucleus_sigma
            >> v1 >> s1 >> electron_val >> electron_sigma) {
      ProtHydrDist& phd = prot_hydr_dists_[type_key];
      phd.electron_val = electron_val;
      phd.electron_sigma = electron_sigma;
      phd.nucleus_val = nucleus_val;
      phd.nucleus_sigma = nucleus_sigma;
    }
  }
}

void AcedrgTables::load_metal_tables(const std::string& dir) {
  // Load allMetalBonds.table
  std::ifstream f(dir + "/allMetalBonds.table");
  if (!f)
    return;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token)
      tokens.push_back(token);

    if (tokens.size() == 10 || tokens.size() == 11) {
      MetalBondEntry entry;
      entry.metal = Element(tokens[0]);
      entry.metal_coord = str_to_int(tokens[1]);
      entry.ligand = Element(tokens[2]);
      entry.ligand_coord = str_to_int(tokens[3]);
      entry.pre_value = std::stod(tokens[4]);
      entry.pre_sigma = std::stod(tokens[5]);
      entry.pre_count = str_to_int(tokens[6]);
      if (tokens.size() == 11) {
        entry.ligand_class = tokens[7];
        entry.value = std::stod(tokens[8]);
        entry.sigma = std::stod(tokens[9]);
        entry.count = str_to_int(tokens[10]);
      } else {
        entry.ligand_class = "NONE";
        entry.value = std::stod(tokens[7]);
        entry.sigma = std::stod(tokens[8]);
        entry.count = str_to_int(tokens[9]);
      }
      metal_bonds_.push_back(entry);
    }
  }

  // Load metal coordination geometry
  std::ifstream f2(dir + "/allMetalDefCoordGeos.table");
  if (f2) {
    while (std::getline(f2, line)) {
      if (line.empty() || line[0] == '#')
        continue;

      std::istringstream iss(line);
      std::string metal_str, geo_str;
      int coord;

      if (iss >> metal_str >> coord >> geo_str) {
        Element metal(metal_str);
        CoordGeometry geo = CoordGeometry::UNKNOWN;
        if (geo_str == "LINEAR") geo = CoordGeometry::LINEAR;
        else if (geo_str == "TRIGONAL_PLANAR") geo = CoordGeometry::TRIGONAL_PLANAR;
        else if (geo_str == "T_SHAPED") geo = CoordGeometry::T_SHAPED;
        else if (geo_str == "TETRAHEDRAL") geo = CoordGeometry::TETRAHEDRAL;
        else if (geo_str == "SQUARE_PLANAR") geo = CoordGeometry::SQUARE_PLANAR;
        else if (geo_str == "TRIGONAL_BIPYRAMIDAL") geo = CoordGeometry::TRIGONAL_BIPYRAMIDAL;
        else if (geo_str == "SQUARE_PYRAMIDAL") geo = CoordGeometry::SQUARE_PYRAMIDAL;
        else if (geo_str == "OCTAHEDRAL") geo = CoordGeometry::OCTAHEDRAL;
        else if (geo_str == "TRIGONAL_PRISM") geo = CoordGeometry::TRIGONAL_PRISM;
        else if (geo_str == "PENTAGONAL_BIPYRAMIDAL") geo = CoordGeometry::PENTAGONAL_BIPYRAMIDAL;
        else if (geo_str == "CAPPED_OCTAHEDRAL") geo = CoordGeometry::CAPPED_OCTAHEDRAL;
        else if (geo_str == "SQUARE_ANTIPRISM") geo = CoordGeometry::SQUARE_ANTIPRISM;
        metal_coord_geo_[metal][coord] = geo;
      }
    }
  }

  // Load metal coordination angles
  std::ifstream f3(dir + "/allMetalCoordGeoAngles.table");
  if (f3) {
    while (std::getline(f3, line)) {
      if (line.empty() || line[0] == '#')
        continue;

      std::istringstream iss(line);
      std::string metal_str, geo_str;
      int coord;
      double angle, sigma;

      if (iss >> metal_str >> coord >> geo_str >> angle >> sigma) {
        MetalAngleEntry entry;
        entry.metal = Element(metal_str);
        entry.coord_number = coord;
        entry.angle = angle;
        entry.sigma = sigma;
        metal_angles_.push_back(entry);
      }
    }
  }
}

void AcedrgTables::load_covalent_radii(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    std::string elem, kind;
    double value = NAN;
    if (!(iss >> elem >> kind >> value))
      continue;
    if (kind != "cova")
      continue;
    std::string key = elem;
    for (char& c : key)
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    covalent_radii_[key] = value;
  }
}

void AcedrgTables::load_atom_type_codes(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: code<whitespace>full_type
    std::istringstream iss(line);
    std::string code, full_type;
    if (iss >> code >> full_type) {
      atom_type_codes_[code] = full_type;
    }
  }
}

void AcedrgTables::load_bond_index(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: ha1 ha2 fileNum
    std::istringstream iss(line);
    int ha1, ha2, file_num;
    if (iss >> ha1 >> ha2 >> file_num) {
      bond_file_index_[ha1][ha2] = file_num;
    }
  }
}

void AcedrgTables::load_bond_tables(const std::string& dir) {
  // Load each bond table file referenced in the index
  std::set<int> loaded_files;

  for (const auto& ha1_pair : bond_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      int file_num = ha2_pair.second;
      if (loaded_files.count(file_num))
        continue;
      loaded_files.insert(file_num);

      std::string path = dir + "/" + std::to_string(file_num) + ".table";
      std::ifstream f(path);
      if (!f)
        continue;

      std::string line;
      while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#')
          continue;

        // Format: ha1 ha2 hybrComb inRing a1NB2 a2NB2 a1NB a2NB atomCode1 atomCode2
        //         value sigma count val2 sig2 cnt2
        std::istringstream iss(line);
        int ha1, ha2;
        std::string hybr_comb, in_ring, a1_nb2, a2_nb2, a1_nb, a2_nb;
        std::string atom_code1, atom_code2;
        double value, sigma;
        int count;
        double value2, sigma2;
        int count2;

        if (!(iss >> ha1 >> ha2 >> hybr_comb >> in_ring
                  >> a1_nb2 >> a2_nb2 >> a1_nb >> a2_nb
                  >> atom_code1 >> atom_code2
                  >> value >> sigma >> count
                  >> value2 >> sigma2 >> count2))
          continue;

        // Get main atom types from codes
        std::string a1_type_m, a2_type_m;
        std::string a1_type_f, a2_type_f;
        auto it1 = atom_type_codes_.find(atom_code1);
        auto it2 = atom_type_codes_.find(atom_code2);
        if (it1 != atom_type_codes_.end()) {
          a1_type_f = it1->second;
          // Extract main type (before '{' if present)
          a1_type_m = a1_type_f;
          size_t brace = a1_type_m.find('{');
          if (brace != std::string::npos)
            a1_type_m = a1_type_m.substr(0, brace);
        }
        if (it2 != atom_type_codes_.end()) {
          a2_type_f = it2->second;
          a2_type_m = a2_type_f;
          size_t brace = a2_type_m.find('{');
          if (brace != std::string::npos)
            a2_type_m = a2_type_m.substr(0, brace);
        }

        // Extract root types (element + ring annotation, e.g., "C[5a]")
        std::string a1_root, a2_root;
        if (!a1_type_m.empty()) {
          size_t paren = a1_type_m.find('(');
          a1_root = (paren != std::string::npos) ? a1_type_m.substr(0, paren) : a1_type_m;
        }
        if (!a2_type_m.empty()) {
          size_t paren = a2_type_m.find('(');
          a2_root = (paren != std::string::npos) ? a2_type_m.substr(0, paren) : a2_type_m;
        }

        ValueStats vs(value, sigma, count);
        ValueStats vs1d(value2, sigma2, count2);

        // Populate 1D structure (full detail)
        bond_idx_1d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb][a1_type_m][a2_type_m].push_back(vs1d);

        // Store full COD-class stats for exact matches (AceDRG uses full codClass first).
        if (!a1_type_f.empty() && !a2_type_f.empty()) {
          bond_idx_full_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                         [a1_nb][a2_nb][a1_type_m][a2_type_m]
                         [a1_type_f][a2_type_f] = vs;
        }

        // Populate 2D structure (no atom types)
        bond_idx_2d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                    [a1_nb][a2_nb].push_back(vs);

        // Populate Nb2D structure (nb2 only, no nb1nb2_sp)
        bond_nb2d_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2].push_back(vs);

        // Populate Nb2DType structure (nb2 + root types, no nb1nb2_sp)
        if (!a1_root.empty() && !a2_root.empty())
          bond_nb2d_type_[ha1][ha2][hybr_comb][in_ring][a1_nb2][a2_nb2]
                         [a1_root][a2_root].push_back(vs1d);

        // Levels 9-11 are populated from allOrgBondsHRS.table in load_bond_hrs().
      }
    }
  }
}

void AcedrgTables::load_angle_index(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return; // Optional file

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;

    // Format: ha1 ha2 ha3 fileNum
    std::istringstream iss(line);
    int ha1, ha2, ha3, file_num;
    if (iss >> ha1 >> ha2 >> ha3 >> file_num) {
      angle_file_index_[ha1][ha2][ha3] = file_num;
    }
  }
}

void AcedrgTables::load_angle_tables(const std::string& dir) {
  // Load each angle table file referenced in the index
  std::set<int> loaded_files;

  for (const auto& ha1_pair : angle_file_index_) {
    for (const auto& ha2_pair : ha1_pair.second) {
      for (const auto& ha3_pair : ha2_pair.second) {
        int file_num = ha3_pair.second;
        if (loaded_files.count(file_num))
          continue;
        loaded_files.insert(file_num);

        std::string path = dir + "/" + std::to_string(file_num) + ".table";
        std::ifstream f(path);
        if (!f)
          continue;

        std::string line;
        while (std::getline(f, line)) {
          if (line.empty() || line[0] == '#')
            continue;

          // 34-column format:
          // 1-3: ha1 ha2 ha3
          // 4: valueKey (ring:hybr_tuple, e.g. "0:SP2_SP2_SP3")
          // 5-7: a1_root a2_root a3_root
          // 8-10: a1_nb2 a2_nb2 a3_nb2
          // 11-13: a1_nb a2_nb a3_nb
          // 14-16: a1_code a2_code a3_code
          // 17-19: AxC value, sigma, count (unused)
          // 20-22: AxM value, sigma, count (1D)
          // 23-25: A_NB value, sigma, count (2D)
          // 26-28: A_NB2 value, sigma, count (3D)
          // 29-31: a1R/a2R/a3R value, sigma, count (4D)
          // 32-34: R3A value, sigma, count (5D)
          std::istringstream iss(line);
          int ha1, ha2, ha3;
          std::string value_key;
          std::string a1_root, a2_root, a3_root;
          std::string a1_nb2, a2_nb2, a3_nb2;
          std::string a1_nb, a2_nb, a3_nb;
          std::string a1_code, a2_code, a3_code;

          if (!(iss >> ha1 >> ha2 >> ha3 >> value_key
                    >> a1_root >> a2_root >> a3_root
                    >> a1_nb2 >> a2_nb2 >> a3_nb2
                    >> a1_nb >> a2_nb >> a3_nb
                    >> a1_code >> a2_code >> a3_code))
            continue;

          // Read 6 sets of value/sigma/count
          double values[6], sigmas[6];
          int counts[6];
          bool ok = true;
          for (int lvl = 0; lvl < 6 && ok; ++lvl) {
            if (!(iss >> values[lvl] >> sigmas[lvl] >> counts[lvl]))
              ok = false;
          }
          if (!ok)
            continue;

          // Get main atom types from codes
          std::string a1_type, a2_type, a3_type;
          auto it1 = atom_type_codes_.find(a1_code);
          auto it2 = atom_type_codes_.find(a2_code);
          auto it3 = atom_type_codes_.find(a3_code);
          if (it1 != atom_type_codes_.end()) {
            a1_type = it1->second;
            size_t brace = a1_type.find('{');
            if (brace != std::string::npos)
              a1_type = a1_type.substr(0, brace);
          }
          if (it2 != atom_type_codes_.end()) {
            a2_type = it2->second;
            size_t brace = a2_type.find('{');
            if (brace != std::string::npos)
              a2_type = a2_type.substr(0, brace);
          }
          if (it3 != atom_type_codes_.end()) {
            a3_type = it3->second;
            size_t brace = a3_type.find('{');
            if (brace != std::string::npos)
              a3_type = a3_type.substr(0, brace);
          }

          // Populate structures at each level with corresponding pre-computed values.
          // AceDRG keeps only the first entry for each key (no aggregation).
          // Level 1D: full detail with atom types
          ValueStats vs1(values[1], sigmas[1], counts[1]);
          auto& angle_1d_vec = angle_idx_1d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb]
                               [a1_type][a2_type][a3_type];
          if (angle_1d_vec.empty())
            angle_1d_vec.push_back(vs1);

          // Level 2D: no atom types
          ValueStats vs2(values[2], sigmas[2], counts[2]);
          auto& angle_2d_vec = angle_idx_2d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2]
                               [a1_nb][a2_nb][a3_nb];
          if (angle_2d_vec.empty())
            angle_2d_vec.push_back(vs2);

          // Level 3D: hash + valueKey + NB2 only
          ValueStats vs3(values[3], sigmas[3], counts[3]);
          auto& angle_3d_vec = angle_idx_3d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root]
                               [a1_nb2][a2_nb2][a3_nb2];
          if (angle_3d_vec.empty())
            angle_3d_vec.push_back(vs3);

          // Level 4D: hash + valueKey + roots
          ValueStats vs4(values[4], sigmas[4], counts[4]);
          auto& angle_4d_vec = angle_idx_4d_[ha1][ha2][ha3][value_key]
                               [a1_root][a2_root][a3_root];
          if (angle_4d_vec.empty())
            angle_4d_vec.push_back(vs4);

          // Level 5D: hash + valueKey only
          ValueStats vs5(values[5], sigmas[5], counts[5]);
          auto& angle_5d_vec = angle_idx_5d_[ha1][ha2][ha3][value_key];
          if (angle_5d_vec.empty())
            angle_5d_vec.push_back(vs5);

          // Level 6D: hash only (leave empty for 34-column data)
        }
      }
    }
  }
}

void AcedrgTables::load_pep_tors(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    return;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    std::string tors_id, label, a1, a2, a3, a4;
    int period = 0;
    int idx = 0;
    double value = 0.0;
    if (!(iss >> tors_id >> label >> a1 >> a2 >> a3 >> a4 >> period >> idx >> value))
      continue;
    TorsionEntry entry;
    entry.value = value;
    entry.period = period;
    entry.priority = idx;
    entry.id = label;
    pep_tors_[a1 + "_" + a2 + "_" + a3 + "_" + a4] = std::move(entry);
  }
}

// ============================================================================
// Implementation - Atom classification
// ============================================================================

std::vector<CodAtomInfo> AcedrgTables::classify_atoms(const ChemComp& cc) const {
  std::vector<CodAtomInfo> atoms(cc.atoms.size());

  // Initialize basic info
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    atoms[i].index = static_cast<int>(i);
    atoms[i].id = cc.atoms[i].id;
    atoms[i].el = cc.atoms[i].el;
    atoms[i].is_metal = is_metal(cc.atoms[i].el);
    atoms[i].charge = cc.atoms[i].charge;
    atoms[i].par_charge = cc.atoms[i].charge;
  }

  // Build adjacency and neighbor lists
  std::vector<std::vector<BondInfo>> adj = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adj);

  // Connectivity counts
  for (size_t i = 0; i < atoms.size(); ++i) {
    atoms[i].connectivity = static_cast<int>(neighbors[i].size());
    atoms[i].conn_atoms_no_metal.clear();
    int metal_conn = 0;
    for (int nb : neighbors[i]) {
      if (atoms[nb].is_metal) {
        ++metal_conn;
      } else {
        atoms[i].conn_atoms_no_metal.push_back(nb);
      }
    }
    atoms[i].metal_connectivity = metal_conn;
  }

  // AceDRG adjusts charges for atoms bonded to metals via valence bookkeeping.
  // Use non-metal bond orders to recompute charge for such atoms so that
  // bonding_idx and aromaticity match AceDRG behavior (e.g., AIV N2).
  for (size_t i = 0; i < atoms.size(); ++i) {
    CodAtomInfo& atom = atoms[i];
    if (atom.is_metal || atom.el == El::H)
      continue;
    bool has_metal_neighbor = atom.connectivity > static_cast<int>(atom.conn_atoms_no_metal.size());
    if (!has_metal_neighbor)
      continue;
    bool has_non_metal_heavy_neighbor = false;
    for (int nb : atom.conn_atoms_no_metal) {
      if (atoms[nb].el != El::H) {
        has_non_metal_heavy_neighbor = true;
        break;
      }
    }
    if (!has_non_metal_heavy_neighbor)
      continue;
    int expected_valence = 0;
    if (atom.el == El::O) expected_valence = 2;
    else if (atom.el == El::N) expected_valence = 3;
    else if (atom.el == El::S) expected_valence = 2;
    else if (atom.el == El::C) expected_valence = 4;
    else if (atom.el == El::P) expected_valence = 3;
    else continue;
    float sum_bo = 0.0f;
    for (const BondInfo& bi : adj[i]) {
      if (!atoms[bi.neighbor_idx].is_metal)
        sum_bo += order_of_bond_type(bi.type);
    }
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    atom.charge = static_cast<float>(-rem_v);
  }

  // Bonding/planarity info is needed for AceDRG ring aromaticity rules.
  set_atoms_bonding_and_chiral_center(atoms, neighbors);

  // Detect rings and populate ring representations
  std::vector<RingInfo> rings;
  detect_rings_acedrg(neighbors, atoms, rings);
  set_ring_aromaticity_from_bonds(adj, atoms, rings);
  set_atoms_ring_rep_s(atoms, rings);

  // Build COD class names (AceDRG style)
  for (size_t i = 0; i < atoms.size(); ++i)
    set_atom_cod_class_name_new2(atoms[i], atoms[i], 2, atoms, neighbors);

  for (size_t i = 0; i < atoms.size(); ++i)
    set_special_3nb_symb2(atoms[i], atoms, neighbors);

  for (size_t i = 0; i < atoms.size(); ++i)
    cod_class_to_atom2(atoms[i].cod_class, atoms[i]);

  // Hybridization and NB1/NB2_SP
  set_atoms_nb1nb2_sp(atoms, neighbors);

  // Ring props from codClass and hash codes
  for (auto& atom : atoms) {
    atom.min_ring_size = get_min_ring2_from_cod_class(atom.cod_class);
    atom.is_aromatic = cod_class_is_aromatic(atom.cod_class);
    atom.hybrid = hybrid_from_bonding_idx(atom.bonding_idx, atom.is_metal,
                                          atom.connectivity);
    compute_hash(atom);
    // Store no-charge variant for COD table lookup consistency.
    // cod_class itself doesn't encode charges (it's based on ring aromaticity and
    // neighbor topology). The charge-sensitive fields are bonding_idx/hybrid/nb1nb2_sp,
    // but those are used separately in lookups. For now, cod_class_no_charge equals
    // cod_class since the COD class string format doesn't include charge info.
    atom.cod_class_no_charge = atom.cod_class;
  }

  // Phase 2: Rebuild codClass with permissive aromaticity for output.
  // AceDRG's reDoAtomCodClassNames() sets isAromatic = isAromaticP after lookups.
  // cod_class_no_charge (used for table lookups) retains the strict version.
  bool has_permissive_diff = false;
  for (const auto& ring : rings)
    if (ring.is_aromatic_permissive != ring.is_aromatic)
      has_permissive_diff = true;
  if (has_permissive_diff) {
    for (auto& ring : rings)
      ring.is_aromatic = ring.is_aromatic_permissive;
    for (auto& atom : atoms)
      atom.ring_rep_s.clear();
    set_atoms_ring_rep_s(atoms, rings);
    for (size_t i = 0; i < atoms.size(); ++i)
      set_atom_cod_class_name_new2(atoms[i], atoms[i], 2, atoms, neighbors);
    for (size_t i = 0; i < atoms.size(); ++i)
      set_special_3nb_symb2(atoms[i], atoms, neighbors);
    // Don't call cod_class_to_atom2 here — it would overwrite cod_main, cod_root,
    // nb_symb, nb2_symb, nb3_symb with permissive-aromaticity values, but these
    // fields must retain strict-aromaticity values for COD table lookups.
  }

  return atoms;
}

// AceDrg hash used in *HRS.table files
void AcedrgTables::compute_hash(CodAtomInfo& atom) const {
  static const int primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229
  };

  // d1: aromaticity (0 or 1)
  int d1 = atom.is_aromatic ? 1 : 0;

  // d2: min ring size mapped to 2-7
  int d2;
  switch (atom.min_ring_size) {
    case 0: d2 = 2; break;
    case 3: d2 = 3; break;
    case 4: d2 = 4; break;
    case 5: d2 = 5; break;
    case 6: d2 = 6; break;
    default: d2 = 7; break;
  }

  // d3: non-metal connectivity + 8
  // AceDRG excludes metal neighbors from connectivity when computing hash.
  int d3 = 8 + static_cast<int>(atom.conn_atoms_no_metal.size());

  // d4: periodic row + 16
  int d4 = 16 + element_row(atom.el);

  // d5: periodic group + 24
  int d5 = 24 + element_group(atom.el);

  // Compute hash as product of primes mod HASH_SIZE
  int64_t prime_product = static_cast<int64_t>(primes[d1]) *
                          primes[d2] * primes[d3] * primes[d4] * primes[d5];
  atom.hashing_value = static_cast<int>(prime_product % HASH_SIZE);

  // If we have hash tables loaded, resolve collisions
  if (!digit_keys_.empty()) {
    std::string footprint = std::to_string(d1) + "_" + std::to_string(d2) +
                            "_" + std::to_string(d3) + "_" + std::to_string(d4) +
                            "_" + std::to_string(d5);

    int pseudo_hash = atom.hashing_value;
    auto it = digit_keys_.find(pseudo_hash);

    if (it != digit_keys_.end()) {
      if (it->second == footprint) {
        atom.hashing_value = pseudo_hash;
      } else {
        // Follow linked hash chain
        bool found = false;
        while (!found) {
          auto link_it = linked_hash_.find(pseudo_hash);
          if (link_it == linked_hash_.end() || link_it->second == -1) {
            break;
          }
          pseudo_hash = link_it->second;
          auto key_it = digit_keys_.find(pseudo_hash);
          if (key_it != digit_keys_.end() && key_it->second == footprint) {
            atom.hashing_value = pseudo_hash;
            found = true;
          }
        }
      }
    }
  }
}

std::vector<std::vector<AcedrgTables::BondInfo>>
AcedrgTables::build_adjacency(const ChemComp& cc) const {
  std::vector<std::vector<BondInfo>> adj(cc.atoms.size());
  for (const auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom);
    int idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 >= 0 && idx2 >= 0) {
      adj[idx1].push_back({idx2, bond.type, bond.aromatic});
      adj[idx2].push_back({idx1, bond.type, bond.aromatic});
    }
  }
  return adj;
}

std::vector<std::vector<int>> AcedrgTables::build_neighbors(
    const std::vector<std::vector<BondInfo>>& adj) const {
  std::vector<std::vector<int>> neighbors(adj.size());
  for (size_t i = 0; i < adj.size(); ++i) {
    neighbors[i].reserve(adj[i].size());
    for (const auto& nb : adj[i])
      neighbors[i].push_back(nb.neighbor_idx);
  }
  return neighbors;
}

void AcedrgTables::set_ring_aromaticity_from_bonds(
    const std::vector<std::vector<BondInfo>>& adj,
    const std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) const {
  // AceDRG has two-phase aromaticity:
  // - Strict (mode 0): only NoMetal pi count → isAromatic (used for COD table lookup)
  // - Permissive (mode 1): NoMetal+All pi counts → isAromaticP (used for output CIF)
  // Ring must be planar (bondingIdx==2, N allowed).
  auto count_non_mc = [&](int idx) -> int {
    int non_mc = 0;
    for (const auto& nb : adj[idx]) {
      const auto& at = atoms[nb.neighbor_idx];
      if (!at.is_metal)
        ++non_mc;
    }
    return non_mc;
  };

  auto is_atom_planar = [&](int idx) -> bool {
    const auto& atom = atoms[idx];
    if (atom.el == El::N)
      return true;
    return atom.bonding_idx == 2;
  };

  auto is_ring_planar = [&](const RingInfo& ring) -> bool {
    for (int idx : ring.atoms)
      if (!is_atom_planar(idx))
        return false;
    return true;
  };

  auto count_atom_pi_no_metal = [&](int idx, int mode) -> double {
    const auto& atom = atoms[idx];
    int non_mc = count_non_mc(idx);
    double aN = 0.0;

    if (atom.bonding_idx == 2) {
      if (atom.charge == 0.0f) {
        if (atom.el == El::C) {
          if (non_mc == 3) {
            bool has_exo = false;
            for (int nb : atom.conn_atoms_no_metal) {
              if (atoms[nb].el == El::O &&
                  atoms[nb].conn_atoms_no_metal.size() == 1 &&
                  atoms[nb].charge == 0.0f) {
                has_exo = true;
              } else if (atoms[nb].el == El::C &&
                         atoms[nb].conn_atoms_no_metal.size() == 3) {
                int h_count = 0;
                for (int nb2 : atoms[nb].conn_atoms_no_metal)
                  if (atoms[nb2].el == El::H)
                    ++h_count;
                if (h_count >= 2)
                  has_exo = true;
              }
            }
            if (!has_exo)
              aN = 1;
          } else if (non_mc == 2) {
            aN = 0.0;
          }
        } else if (atom.el == El::N) {
          if (non_mc == 2)
            aN = 1;
          else if (non_mc == 3)
            aN = 2;
        } else if (atom.el == El::B) {
          if (non_mc == 2)
            aN = 1;
          else if (non_mc == 3)
            aN = 0.0;
        } else if (atom.el == El::O) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.el == El::S) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.el == El::P) {
          if (non_mc == 3)
            aN = 2;
        }
      } else {
        if (atom.el == El::C) {
          if (atom.charge == -1.0f) {
            if (non_mc == 3)
              aN = 2;
            else if (non_mc == 2)
              aN = (mode == 1) ? 1.0 : 2.0;
          } else if (atom.charge == -2.0f) {
            if (non_mc == 2)
              aN = 2;
          }
        } else if (atom.el == El::N) {
          if (atom.charge == -1.0f) {
            if (non_mc == 2)
              aN = 2;
          } else if (atom.charge == 1.0f) {
            if (non_mc == 3)
              aN = (mode == 1) ? 1.0 : 2.0;
          }
        } else if (atom.el == El::O) {
          if (atom.charge == 1.0f && non_mc == 2)
            aN = 1;
        } else if (atom.el == El::B) {
          if (atom.charge == -1.0f && non_mc == 3)
            aN = 1;
        }
      }
    } else if (atom.bonding_idx == 3 &&
               (atom.el == El::N || atom.el == El::B)) {
      if (atom.el == El::N) {
        if (atom.charge == -1.0f) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.charge == 1.0f) {
          if (non_mc == 3)
            aN = 1;
        } else {
          aN = 2;
        }
      } else if (atom.el == El::B) {
        aN = 0.0;
      }
    }

    return aN;
  };

  for (size_t i = 0; i < rings.size(); ++i) {
    RingInfo& ring = rings[i];
    ring.is_aromatic = false;

    if (!is_ring_planar(ring))
      continue;

    // AceDRG uses strict aromaticity (mode 0) for the COD table lookup.
    // Mode 0 differs from mode 1 for charged N+ with 3 connections:
    // mode 0 gives 2 pi electrons, mode 1 gives 1.
    // AceDRG only checks the NoMetal pi count in strict mode.
    double pi1 = 0.0;
    for (int idx : ring.atoms)
      pi1 += count_atom_pi_no_metal(idx, 0);
    if (pi1 > 0.0 && std::fabs(std::fmod(pi1, 4.0) - 2.0) < 0.001)
      ring.is_aromatic = true;
    if (verbose >= 2) {
      std::fprintf(stderr, "    ring %zu (size=%zu): pi1=%.1f aromatic=%d atoms:",
                   i, ring.atoms.size(), pi1, ring.is_aromatic ? 1 : 0);
      for (int idx : ring.atoms)
        std::fprintf(stderr, " %s", atoms[idx].id.c_str());
      std::fprintf(stderr, "\n");
    }
  }

  // AceDRG pyrole rule: if there are exactly 4 five-member rings with 4C+1N
  // and they are planar, mark them aromatic even if pi-count failed.
  std::vector<size_t> pyrole_rings;
  for (size_t i = 0; i < rings.size(); ++i) {
    const RingInfo& ring = rings[i];
    if (ring.atoms.size() != 5)
      continue;
    int num_c = 0;
    int num_n = 0;
    for (int idx : ring.atoms) {
      if (atoms[idx].el == El::C)
        ++num_c;
      else if (atoms[idx].el == El::N)
        ++num_n;
    }
    if (num_c == 4 && num_n == 1)
      pyrole_rings.push_back(i);
  }
  if (pyrole_rings.size() == 4) {
    for (size_t i : pyrole_rings) {
      if (is_ring_planar(rings[i]))
        rings[i].is_aromatic = true;
    }
  }

  // Permissive aromaticity (AceDRG mode 1): used for output codClass names.
  // Differs from strict for charged N+ (mode 1 gives 1 pi, mode 0 gives 2),
  // and also checks the "all" pi count (which always gives 1 for N+).
  auto count_atom_pi_all = [&](int idx) -> double {
    const auto& atom = atoms[idx];
    int non_mc = count_non_mc(idx);
    double aN = 0.0;

    if (atom.bonding_idx == 2) {
      if (atom.charge == 0.0f) {
        if (atom.el == El::C) {
          if (non_mc == 3) {
            bool has_exo = false;
            for (int nb : atom.conn_atoms_no_metal) {
              if (atoms[nb].el == El::O &&
                  atoms[nb].conn_atoms_no_metal.size() == 1 &&
                  atoms[nb].charge == 0.0f) {
                has_exo = true;
              } else if (atoms[nb].el == El::C &&
                         atoms[nb].conn_atoms_no_metal.size() == 3) {
                int h_count = 0;
                for (int nb2 : atoms[nb].conn_atoms_no_metal)
                  if (atoms[nb2].el == El::H)
                    ++h_count;
                if (h_count >= 2)
                  has_exo = true;
              }
            }
            if (!has_exo)
              aN = 1;
          } else if (non_mc == 2) {
            aN = 0.0;
          }
        } else if (atom.el == El::N) {
          if (non_mc == 2)
            aN = 1;
          else if (non_mc == 3)
            aN = 2;
        } else if (atom.el == El::B) {
          if (non_mc == 2)
            aN = 1;
          else if (non_mc == 3)
            aN = 0.0;
        } else if (atom.el == El::O) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.el == El::S) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.el == El::P) {
          if (non_mc == 3)
            aN = 2;
        }
      } else {
        if (atom.el == El::C) {
          if (atom.charge == -1.0f) {
            if (non_mc == 3)
              aN = 2;
            else if (non_mc == 2)
              aN = 1.0;
          }
        } else if (atom.el == El::N) {
          if (atom.charge == -1.0f) {
            if (non_mc == 2)
              aN = 2;
          } else if (atom.charge == 1.0f) {
            if (non_mc == 3)
              aN = 1.0;
          }
        } else if (atom.el == El::O) {
          if (atom.charge == 1.0f && non_mc == 2)
            aN = 1;
        } else if (atom.el == El::B) {
          if (atom.charge == -1.0f && non_mc == 3)
            aN = 1;
        }
      }
    } else if (atom.bonding_idx == 3 &&
               (atom.el == El::N || atom.el == El::B)) {
      if (atom.el == El::N) {
        if (atom.charge == -1.0f) {
          if (non_mc == 2)
            aN = 2;
        } else if (atom.charge == 1.0f) {
          if (non_mc == 3)
            aN = 1;
        } else {
          aN = 2;
        }
      } else if (atom.el == El::B) {
        aN = 0.0;
      }
    }

    return aN;
  };

  for (auto& ring : rings) {
    ring.is_aromatic_permissive = ring.is_aromatic;
    if (ring.is_aromatic)
      continue;
    if (!is_ring_planar(ring))
      continue;
    double pi1 = 0.0;
    double pi2 = 0.0;
    for (int idx : ring.atoms) {
      pi1 += count_atom_pi_no_metal(idx, 1);
      pi2 += count_atom_pi_all(idx);
    }
    if ((pi1 > 0.0 && std::fabs(std::fmod(pi1, 4.0) - 2.0) < 0.001) ||
        (pi2 > 0.0 && std::fabs(std::fmod(pi2, 4.0) - 2.0) < 0.001))
      ring.is_aromatic_permissive = true;
  }
}

void AcedrgTables::detect_rings_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) const {
  rings.clear();
  std::map<std::string, int> ring_index;

  for (size_t i = 0; i < atoms.size(); ++i) {
    if (atoms[i].is_metal || atoms[i].el == El::H)
      continue;
    std::map<int, std::string> seen;
    std::map<int, std::string> path;
    check_one_path_acedrg(neighbors, atoms, rings, ring_index,
                          static_cast<int>(i), static_cast<int>(i), -999,
                          1, 7, seen, path);
  }
}

void AcedrgTables::check_one_path_acedrg(
    const std::vector<std::vector<int>>& neighbors,
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings,
    std::map<std::string, int>& ring_index,
    int ori_idx, int cur_idx, int prev_idx, int cur_lev, int max_ring,
    std::map<int, std::string>& seen_atom_ids,
    std::map<int, std::string>& atom_ids_in_path) const {
  if (cur_lev >= max_ring)
    return;

  bool path_collision = false;
  for (int nb : neighbors[cur_idx]) {
    if (nb != ori_idx && nb != prev_idx &&
        seen_atom_ids.find(nb) != seen_atom_ids.end()) {
      path_collision = true;
      break;
    }
  }

  if (!path_collision) {
    for (int nb : neighbors[cur_idx]) {
      if (nb == ori_idx && nb != prev_idx && cur_lev > 2 && atoms[nb].el != El::H) {
        atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
        std::list<std::string> all_ids;
        std::list<std::string> all_seris;
        std::vector<int> ring_atoms;
        for (const auto& it : atom_ids_in_path) {
          all_seris.push_back(std::to_string(it.first));
          all_ids.push_back(it.second);
          ring_atoms.push_back(it.first);
        }
        all_seris.sort(compare_no_case);
        all_ids.sort(compare_no_case);

        std::string rep;
        for (const auto& id : all_ids)
          rep += id;

        std::string s_rep;
        int nrs = 0;
        for (const auto& seri : all_seris) {
          if (nrs == 0)
            s_rep += seri;
          else
            s_rep += "_" + seri;
          ++nrs;
        }

        atoms[ori_idx].ring_rep[rep] = static_cast<int>(atom_ids_in_path.size());

        if (ring_index.find(s_rep) == ring_index.end()) {
          RingInfo ring;
          ring.atoms = ring_atoms;
          ring.rep = rep;
          ring.s_rep = s_rep;
          ring.is_aromatic = true;
          for (int idx : ring_atoms)
            if (!atoms[idx].is_aromatic)
              ring.is_aromatic = false;

          int idx = static_cast<int>(rings.size());
          ring_index[s_rep] = idx;
          rings.push_back(ring);
          for (int atom_idx : ring_atoms)
            atoms[atom_idx].in_rings.push_back(idx);
        }

        atom_ids_in_path.erase(cur_idx);
        path_collision = true;
        break;
      }
    }
  }

  if (!path_collision) {
    int next_lev = cur_lev + 1;
    seen_atom_ids[cur_idx] = atoms[cur_idx].id;
    if (next_lev < max_ring) {
      if (cur_lev == 1) {
        seen_atom_ids.clear();
        atom_ids_in_path.clear();
        seen_atom_ids[cur_idx] = atoms[cur_idx].id;
        atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
      }
      atom_ids_in_path[cur_idx] = atoms[cur_idx].id;
      for (int nb : neighbors[cur_idx]) {
        if (nb != prev_idx && !atoms[nb].is_metal) {
          check_one_path_acedrg(neighbors, atoms, rings, ring_index,
                                ori_idx, nb, cur_idx, next_lev, max_ring,
                                seen_atom_ids, atom_ids_in_path);
        }
      }
      atom_ids_in_path.erase(cur_idx);
      seen_atom_ids.erase(cur_idx);
    }
    atom_ids_in_path.erase(cur_idx);
    seen_atom_ids.erase(cur_idx);
  }
}

void AcedrgTables::set_atoms_ring_rep_s(
    std::vector<CodAtomInfo>& atoms,
    std::vector<RingInfo>& rings) const {
  for (const auto& ring : rings) {
    std::string size = std::to_string(ring.atoms.size());
    std::string rep_id;
    std::list<std::string> all_seris;
    for (int idx : ring.atoms)
      all_seris.push_back(std::to_string(idx));
    all_seris.sort(compare_no_case);
    int nrs = 0;
    for (const auto& seri : all_seris) {
      if (nrs == 0)
        rep_id += seri;
      else
        rep_id += "_" + seri;
      ++nrs;
    }

    for (int idx : ring.atoms) {
      if (ring.is_aromatic)
        atoms[idx].ring_rep_s[rep_id] = size + "a";
      else
        atoms[idx].ring_rep_s[rep_id] = size;
    }
  }
}

// Build ring size_map for cod_class annotation.
// For fused aromatic rings, collect all sizes (e.g., C[5a,6a] for indole fusion).
static void build_ring_size_map(const std::map<std::string, std::string>& ring_rep_s,
                                std::map<std::string, int>& size_map) {
  for (const auto& it : ring_rep_s) {
    const std::string& ring_str = it.second;
    size_map[ring_str] += 1;
  }
}

void AcedrgTables::set_atom_cod_class_name_new2(
    CodAtomInfo& atom, const CodAtomInfo& ori_atom, int lev,
    const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  if (lev == 1) {
    atom.cod_class.clear();
    atom.cod_class.append(atom.el.name());

    if (!atom.ring_rep_s.empty()) {
      std::map<std::string, int> size_map;
      build_ring_size_map(atom.ring_rep_s, size_map);

      atom.cod_class.append("[");
      int i = 0;
      int j = static_cast<int>(size_map.size());
      for (const auto& it : size_map) {
        std::string size = it.first;
        std::string num = std::to_string(it.second);
        if (it.second >= 3)
          atom.cod_class.append(num + "x" + size);
        else if (it.second == 2)
          atom.cod_class.append(size + "," + size);
        else
          atom.cod_class.append(size);
        if (i != j - 1)
          atom.cod_class.append(",");
        else
          atom.cod_class.append("]");
        ++i;
      }
    }

    std::string t_str;
    std::list<std::string> t_str_list;
    std::map<std::string, int> comps;
    for (int nb : neighbors[atom.index]) {
      if (nb == ori_atom.index || atoms[nb].is_metal)
        continue;
      std::string nb_type = atoms[nb].el.name();
      if (!atoms[nb].ring_rep_s.empty()) {
        std::map<std::string, int> size_map;
        build_ring_size_map(atoms[nb].ring_rep_s, size_map);
        nb_type.append("[");
        int i = 0;
        int j = static_cast<int>(size_map.size());
        for (const auto& it : size_map) {
          std::string size = it.first;
          std::string num = std::to_string(it.second);
          if (it.second >= 3)
            nb_type.append(num + "x" + size);
          else if (it.second == 2)
            nb_type.append(size + "," + size);
          else
            nb_type.append(size);
          if (i != j - 1)
            nb_type.append(",");
          else
            nb_type.append("]");
          ++i;
        }
      }
      comps[nb_type] += 1;
    }

    std::vector<SortMap> sorted;
    for (const auto& it : comps) {
      SortMap sm;
      sm.key = it.first;
      sm.val = it.second;
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), desc_sort_map_key);
    for (const auto& sm : sorted) {
      std::string s1 = sm.key + std::to_string(sm.val);
      std::string s2;
      for (int i = 0; i < sm.val; ++i)
        s2.append(sm.key);
      if (s1.size() < s2.size())
        t_str_list.push_back(s1);
      else
        t_str_list.push_back(s2);
    }
    for (const auto& s : t_str_list)
      t_str.append(s);

    atom.cod_class.append(t_str);
  } else if (lev == 2) {
    atom.cod_class.clear();
    atom.cod_class.append(atom.el.name());

    if (!atom.ring_rep_s.empty()) {
      std::map<std::string, int> size_map;
      build_ring_size_map(atom.ring_rep_s, size_map);

      atom.cod_class.append("[");
      int i = 0;
      int j = static_cast<int>(size_map.size());
      for (const auto& it : size_map) {
        std::string size = it.first;
        std::string num = std::to_string(it.second);
        if (it.second >= 3)
          atom.cod_class.append(num + "x" + size);
        else if (it.second == 2)
          atom.cod_class.append(size + "," + size);
        else
          atom.cod_class.append(size);
        if (i != j - 1)
          atom.cod_class.append(",");
        else
          atom.cod_class.append("]");
        ++i;
      }
    }

    int low_lev = lev - 1;
    std::map<std::string, std::vector<int>> id_map;
    for (int nb : neighbors[atom.index]) {
      if (atoms[nb].is_metal)
        continue;
      CodAtomInfo nb_atom = atoms[nb];
      set_atom_cod_class_name_new2(nb_atom, ori_atom, low_lev, atoms, neighbors);
      auto& entry = id_map[nb_atom.cod_class];
      if (entry.empty()) {
        entry.push_back(1);
        int non_metal = 0;
        for (int nb2 : neighbors[nb])
          if (!atoms[nb2].is_metal)
            ++non_metal;
        entry.push_back(non_metal);
      } else {
        entry[0] += 1;
      }
    }

    std::vector<SortMap2> sorted;
    for (const auto& it : id_map) {
      SortMap2 sm;
      sm.key = it.first;
      sm.val = it.second[0];
      sm.nNB = it.second[1];
      sorted.push_back(sm);
    }
    std::sort(sorted.begin(), sorted.end(), desc_sort_map_key2);
    for (const auto& sm : sorted) {
      if (sm.val == 1)
        atom.cod_class.append("(" + sm.key + ")");
      else
        atom.cod_class.append("(" + sm.key + ")" + std::to_string(sm.val));
    }
  }
}

void AcedrgTables::set_special_3nb_symb2(
    CodAtomInfo& atom, const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  if (atom.ring_rep.empty())
    return;

  std::vector<int> ser_num_nb123;
  std::map<std::string, int> nb3_props;

  for (int nb1 : neighbors[atom.index]) {
    if (atoms[nb1].is_metal)
      continue;
    if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb1) == ser_num_nb123.end())
      ser_num_nb123.push_back(nb1);
    for (int nb2 : neighbors[nb1]) {
      if (atoms[nb2].is_metal)
        continue;
      if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb2) == ser_num_nb123.end() &&
          nb2 != atom.index) {
        ser_num_nb123.push_back(nb2);
      }
    }
  }

  for (int nb1 : neighbors[atom.index]) {
    if (atoms[nb1].ring_rep.empty())
      continue;
    for (int nb2 : neighbors[nb1]) {
      if (atoms[nb2].ring_rep.empty())
        continue;
      for (int nb3 : neighbors[nb2]) {
        if (atoms[nb3].is_metal)
          continue;
        if (std::find(ser_num_nb123.begin(), ser_num_nb123.end(), nb3) == ser_num_nb123.end() &&
            nb3 != atom.index) {
          std::string prop = atoms[nb3].el.name();
          int deg = 0;
          for (int nbx : neighbors[nb3])
            if (!atoms[nbx].is_metal)
              ++deg;
          prop.append("<" + std::to_string(deg) + ">");
          nb3_props[prop] += 1;
          ser_num_nb123.push_back(nb3);
        }
      }
    }
  }

  std::list<std::string> comps;
  for (const auto& it : nb3_props) {
    std::string id = std::to_string(it.second) + "|" + it.first;
    comps.push_back(id);
  }
  comps.sort(compare_no_case2);

  if (!comps.empty()) {
    std::string all3 = "{";
    unsigned i = 0;
    unsigned n = comps.size();
    for (const auto& id : comps) {
      if (i < n - 1)
        all3.append(id + ",");
      else
        all3.append(id);
      ++i;
    }
    all3.append("}");
    atom.cod_class.append(all3);
  }
}

void AcedrgTables::cod_class_to_atom2(const std::string& cod_class,
                                             CodAtomInfo& atom) const {
  std::string t_cod = trim_str(cod_class);
  atom.cod_class = t_cod;
  atom.nb_symb.clear();
  atom.nb2_symb.clear();
  atom.nb3_symb.clear();

  std::vector<std::string> two_parts;
  if (t_cod.find('{') != std::string::npos) {
    two_parts = split_str(t_cod, '{');
    if (two_parts.size() == 2) {
      atom.cod_main = two_parts[0];
      std::vector<std::string> nb3 = split_str(two_parts[1], '}');
      if (!nb3.empty())
        atom.nb3_symb = nb3[0];
    } else {
      atom.cod_main = t_cod;
    }
  } else {
    atom.cod_main = t_cod;
  }

  std::vector<std::string> atm_strs = split_str(atom.cod_main, '(');
  if (!atm_strs.empty()) {
    atom.cod_root = trim_str(atm_strs[0]);
  }

  std::vector<NB1stFam> all_nbs;
  for (size_t i = 1; i < atm_strs.size(); ++i) {
    std::string tS = trim_str(atm_strs[i]);
    std::vector<std::string> nb1 = split_str(tS, ')');
    NB1stFam fam;
    if (nb1.size() > 1) {
      fam.repN = str_to_int(nb1[1]);
      if (fam.repN == 0)
        fam.repN = 1;
    } else {
      fam.repN = 1;
    }

    std::string tS1 = trim_str(nb1[0]);
    get_small_family(tS1, fam);
    all_nbs.push_back(fam);
  }

  for (const auto& fam : all_nbs) {
    for (int j = 0; j < fam.repN; ++j) {
      std::string sN = std::to_string(static_cast<int>(fam.NB2ndList.size()) + 1);
      atom.nb_symb += fam.name + "-" + sN + ":";
      atom.nb2_symb += sN + ":";
    }
  }
}

void AcedrgTables::set_atoms_nb1nb2_sp(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  for (auto& atom : atoms) {
    std::vector<std::string> nb1_nb2_sp_set;
    for (int nb1 : neighbors[atom.index]) {
      if (atoms[nb1].is_metal)
        continue;
      std::string nb1_main = atoms[nb1].cod_root;
      std::vector<int> nb2_sp_set;
      for (int nb2 : neighbors[nb1]) {
        if (atoms[nb2].is_metal)
          continue;
        nb2_sp_set.push_back(atoms[nb2].bonding_idx);
      }
      std::sort(nb2_sp_set.begin(), nb2_sp_set.end(), std::greater<int>());
      std::string nb2_sp_str;
      for (size_t i = 0; i < nb2_sp_set.size(); ++i) {
        nb2_sp_str.append(std::to_string(nb2_sp_set[i]));
        if (i != nb2_sp_set.size() - 1)
          nb2_sp_str.append("_");
      }
      nb1_nb2_sp_set.emplace_back(nb1_main + "-" + nb2_sp_str);
    }
    // Sort alphabetically by the string (same order as AceDRG tables)
    std::sort(nb1_nb2_sp_set.begin(), nb1_nb2_sp_set.end(),
              [](const auto& a, const auto& b) {
                return compare_no_case(a, b);
              });
    // Build nb1nb2_sp in alphabetical order
    atom.nb1nb2_sp.clear();
    for (size_t i = 0; i < nb1_nb2_sp_set.size(); ++i) {
      atom.nb1nb2_sp.append(nb1_nb2_sp_set[i]);
      if (i != nb1_nb2_sp_set.size() - 1)
        atom.nb1nb2_sp.append(":");
    }
  }
}

void AcedrgTables::set_atoms_nb_symb_from_neighbors(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  for (auto& atom : atoms) {
    // Collect neighbor info: cod_root and connectivity
    // The value used in nb_symb/nb2_symb is the neighbor's connectivity
    // (number of bonded atoms), which matches the table format used in
    // acedrg's indexed angle tables (e.g., "C-4:" means carbon with 4 bonds).
    struct NbInfo {
      std::string root;
      int connectivity;
    };
    std::vector<NbInfo> nb_info;
    for (int nb_idx : neighbors[atom.index]) {
      if (atoms[nb_idx].is_metal)
        continue;
      int non_metal_conn = static_cast<int>(atoms[nb_idx].conn_atoms_no_metal.size());
      nb_info.push_back({atoms[nb_idx].cod_root,
                         non_metal_conn});
    }

    // Build "root-connectivity" strings for sorting and lookup
    std::vector<std::string> nb_strs;
    nb_strs.reserve(nb_info.size());
    for (const auto& info : nb_info) {
      nb_strs.push_back(info.root + "-" + std::to_string(info.connectivity));
    }

    // Sort by: 1) string length (longer first), 2) connectivity (higher first)
    // This matches desc_sort_map_key2 sorting used in set_atom_cod_class_name_new2
    std::vector<size_t> indices(nb_info.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&nb_strs, &nb_info](size_t a, size_t b) {
                if (nb_strs[a].length() > nb_strs[b].length())
                  return true;
                if (nb_strs[a].length() < nb_strs[b].length())
                  return false;
                // Same length: sort by connectivity (higher first)
                return nb_info[a].connectivity > nb_info[b].connectivity;
              });

    // Build nb_symb and nb2_symb using sorted order
    atom.nb_symb.clear();
    atom.nb2_symb.clear();
    for (size_t i : indices) {
      atom.nb_symb += nb_strs[i] + ":";
      atom.nb2_symb += std::to_string(nb_info[i].connectivity) + ":";
    }
  }
}

void AcedrgTables::set_atoms_bonding_and_chiral_center(
    std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  std::map<int, std::vector<int>> num_conn_map;

  for (auto& atom : atoms) {
    int t_len = 0;
    int t_m_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
      else
        t_m_len++;
    }
    if (atom.metal_connectivity > 0)
      t_m_len = atom.metal_connectivity;

    if (atom.el == El::C) {
      if (t_len == 3 && t_m_len == 2)
        t_len = 3;
    }

    num_conn_map[atom.index].push_back(t_len);
    num_conn_map[atom.index].push_back(t_m_len);

    if (t_len > 4) {
      atom.bonding_idx = t_len;
    } else if (atom.el == El::C || atom.el == El::Si || atom.el == El::Ge) {
      if (t_len == 4) {
        atom.bonding_idx = 3;
      } else if (t_len == 3) {
        atom.bonding_idx = 2;
      } else if (t_len == 2) {
        if (get_num_oxy_connect(atoms, atom, neighbors) == 1)
          atom.bonding_idx = 1;
        else if (t_m_len == 1 || atom.charge == -1.0f)
          atom.bonding_idx = 2;
        else if (t_m_len == 1 || atom.charge == -2.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 1;
      }
    } else if (atom.el == El::N || atom.el == El::As) {
      if (t_len == 4 || t_len == 3) {
        atom.bonding_idx = 3;
      } else if (t_len == 2) {
        if (atom.charge == 1.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 2;
      } else if (t_len == 1) {
        atom.bonding_idx = 1;
      }
    } else if (atom.el == El::B) {
      if (t_len == 4) {
        atom.bonding_idx = 3;
      } else if (t_len == 3) {
        atom.bonding_idx = 2;
      } else if (t_len == 2) {
        if (atom.charge == 1.0f)
          atom.bonding_idx = 1;
        else
          atom.bonding_idx = 2;
      } else if (t_len == 1) {
        atom.bonding_idx = 1;
      }
    } else if (atom.el == El::O) {
      // Use t_len (non-metal connection count) like AceDRG, which
      // excludes metal bonds from connAtoms before hybridization assignment.
      if (t_len == 2) {
        atom.bonding_idx = 3;
      } else if (t_len == 1) {
        // Find the first non-metal neighbor and check its non-metal connections
        int first_nb = -1;
        for (int nb : neighbors[atom.index])
          if (!atoms[nb].is_metal) { first_nb = nb; break; }
        int nb_non_metal = 0;
        if (first_nb >= 0)
          for (int nnb : neighbors[first_nb])
            if (!atoms[nnb].is_metal)
              nb_non_metal++;
        if (nb_non_metal != 1)
          atom.bonding_idx = 2;
        else
          atom.bonding_idx = 3;
      } else {
        atom.bonding_idx = 3;
      }
    } else if (atom.el == El::P) {
      if (t_len == 4 || t_len == 3 || t_len == 2 || t_len == 5)
        atom.bonding_idx = 3;
    } else if (atom.el == El::S) {
      if (t_len == 2 || t_len == 3 || t_len == 4)
        atom.bonding_idx = 3;
      else if (t_len == 6)
        atom.bonding_idx = 5;
      else if (t_len == 1)
        atom.bonding_idx = 3;
    } else if (atom.el == El::Se) {
      if (t_len == 4 || t_len == 3 || t_len == 2) {
        atom.bonding_idx = (t_len == 3) ? 2 : 3;
      } else if (t_len == 6) {
        atom.bonding_idx = 5;
      } else if (t_len == 1) {
        atom.bonding_idx = 3;
      }
    } else if (atom.el == El::Br) {
      if (t_len == 3)
        atom.bonding_idx = 3;
    } else if (atom.el == El::H || atom.el == El::D) {
      atom.bonding_idx = 0;  // H has SP0 hybridization (no p orbitals)
    }
  }

  for (auto& atom : atoms) {
    int t_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
    }

    if (atom.el == El::O) {
      if (t_len == 2 && atom.par_charge == 0.0f) {
        bool l_sp2 = false;
        for (int nb : neighbors[atom.index]) {
          if (atoms[nb].bonding_idx == 2) {
            l_sp2 = true;
            break;
          }
        }
        if (l_sp2)
          atom.bonding_idx = 2;
      }
    }
  }

  std::map<int, int> pre_bonding;
  for (const auto& atom : atoms)
    pre_bonding[atom.index] = atom.bonding_idx;

  for (auto& atom : atoms) {
    int t_len = 0;
    for (int nb : neighbors[atom.index]) {
      if (!atoms[nb].is_metal)
        t_len++;
    }
    if (atom.el == El::N || atom.el == El::As) {
      if (t_len == 3) {
        if (atom.charge == 0.0f) {
          bool l_sp2 = false;
          for (int nb : neighbors[atom.index]) {
            if (pre_bonding[nb] == 2 && atoms[nb].el != El::O) {
              l_sp2 = true;
              break;
            }
          }
          if (l_sp2) {
            if (num_conn_map[atom.index][1] != 0)
              atom.bonding_idx = 3;
            else
              atom.bonding_idx = 2;
          } else {
            atom.bonding_idx = 3;
          }
        } else if (atom.charge == 1.0f) {
          atom.bonding_idx = 2;
        }
      }
    }
    if (atom.el == El::S) {
      if (t_len == 2 && atom.charge == 0.0f) {
        bool l_sp2 = false;
        for (int nb : neighbors[atom.index]) {
          if (pre_bonding[nb] == 2 && atoms[nb].el != El::O) {
            l_sp2 = true;
            break;
          }
        }
        if (l_sp2)
          atom.bonding_idx = 2;
      }
    }
    if (atom.el == El::C) {
      if (t_len == 3 && atom.charge == -1.0f) {
        std::vector<int> sp2_set;
        for (int nb : neighbors[atom.index]) {
          if (atoms[nb].bonding_idx == 2)
            sp2_set.push_back(nb);
        }
        if (sp2_set.size() == 2)
          atom.bonding_idx = 2;
        else
          atom.bonding_idx = 3;
      }
    }
  }
}

int AcedrgTables::get_num_oxy_connect(const std::vector<CodAtomInfo>& atoms,
                                             const CodAtomInfo& atom,
                                             const std::vector<std::vector<int>>& neighbors) const {
  int nO = 0;
  for (int nb : neighbors[atom.index])
    if (atoms[nb].el == El::O)
      nO++;
  return nO;
}

int AcedrgTables::get_min_ring2_from_cod_class(const std::string& cod_class) const {
  int r_size = 0;
  if (!cod_class.empty()) {
    std::vector<std::string> tmp1 = split_str(cod_class, '(');
    if (!tmp1.empty()) {
      if (tmp1[0].find('[') != std::string::npos) {
        std::vector<std::string> tmp2 = split_str(tmp1[0], '[');
        if (tmp2.size() > 1) {
          if (tmp2[1].find(',') != std::string::npos) {
            std::vector<std::string> tmp3 = split_str(tmp2[1], ',');
            if (!tmp3.empty()) {
              if (tmp3[0].find('x') != std::string::npos) {
                std::vector<std::string> tmp4 = split_str(tmp3[0], 'x');
                if (tmp4.size() > 1)
                  r_size = str_to_int(tmp4[1]);
              } else {
                r_size = str_to_int(tmp3[0]);
              }
            }
          } else {
            std::vector<std::string> tmp3 = split_str(tmp2[1], ']');
            if (!tmp3.empty()) {
              if (tmp3[0].find('x') != std::string::npos) {
                std::vector<std::string> tmp4 = split_str(tmp3[0], 'x');
                if (tmp4.size() > 1)
                  r_size = str_to_int(tmp4[1]);
              } else {
                r_size = str_to_int(tmp3[0]);
              }
            }
          }
        }
      }
    }
  }
  return r_size;
}

bool AcedrgTables::cod_class_is_aromatic(const std::string& cod_class) const {
  if (cod_class.empty())
    return false;
  std::vector<std::string> secs = split_str(cod_class, '(');
  if (secs.empty())
    return false;
  if (secs[0].find('[') != std::string::npos) {
    std::vector<std::string> rs = split_str(secs[0], '[');
    if (rs.size() > 1)
      return rs[1].find('a') != std::string::npos;
  }
  return false;
}

void AcedrgTables::get_small_family(const std::string& in_str, NB1stFam& fam) const {
  fam.name.clear();
  fam.NB2ndList.clear();
  std::vector<std::string> ch_list = {",", "x"};
  std::string name_str;
  bool l_r = false;
  int n_rep = 1;
  for (size_t i = 0; i < in_str.size(); ++i) {
    char c = in_str[i];
    if (std::isalpha(static_cast<unsigned char>(c))) {
      if (std::toupper(c) == c) {
        if (!name_str.empty()) {
          if (fam.name.empty()) {
            fam.name = name_str;
            n_rep = 1;
          } else {
            for (int l = 0; l < n_rep; ++l)
              fam.NB2ndList.push_back(name_str);
            n_rep = 1;
          }
        }
        name_str = c;
      } else {
        name_str += c;
      }
    } else if (c == '[') {
      name_str += c;
      l_r = true;
    } else if (c == ']') {
      name_str += c;
      l_r = false;
    } else if (std::find(ch_list.begin(), ch_list.end(), std::string(1, c)) != ch_list.end()) {
      name_str += c;
    } else if (std::isdigit(static_cast<unsigned char>(c))) {
      if (l_r) {
        name_str += c;
      } else {
        n_rep = std::stoi(in_str.substr(i, 1));
      }
    }
  }

  if (!name_str.empty()) {
    if (fam.name.empty()) {
      fam.name = name_str;
    } else {
      for (int l = 0; l < n_rep; ++l)
        fam.NB2ndList.push_back(name_str);
    }
  }
}


bool AcedrgTables::are_in_same_ring(const CodAtomInfo& a1,
                                           const CodAtomInfo& a2) const {
  for (int ring_idx : a1.in_rings)
    if (std::find(a2.in_rings.begin(), a2.in_rings.end(), ring_idx) != a2.in_rings.end())
      return true;
  return false;
}

int AcedrgTables::angle_ring_size(const CodAtomInfo& center,
                                         const CodAtomInfo& a1,
                                         const CodAtomInfo& a3) const {
  for (const auto& it : center.ring_rep) {
    if (a1.ring_rep.find(it.first) != a1.ring_rep.end() &&
        a3.ring_rep.find(it.first) != a3.ring_rep.end())
      return it.second;
  }
  return 0;
}

void AcedrgTables::order_bond_atoms(const CodAtomInfo& a1, const CodAtomInfo& a2,
                                           const CodAtomInfo*& first,
                                           const CodAtomInfo*& second) const {
  if (a1.hashing_value < a2.hashing_value) {
    first = &a1;
    second = &a2;
    return;
  }
  if (a1.hashing_value > a2.hashing_value) {
    first = &a2;
    second = &a1;
    return;
  }
  // AceDRG: for equal hash, order by codAtmMain (length, then case-insensitive
  // lexicographic). When codAtmMain is also equal, the atom with the greater
  // name (under compare_no_case2) goes first.
  if (compare_no_case2(a1.cod_main, a2.cod_main)) {
    first = &a1;
    second = &a2;
  } else if (compare_no_case2(a2.cod_main, a1.cod_main)) {
    first = &a2;
    second = &a1;
  } else if (!compare_no_case2(a1.id, a2.id)) {
    first = &a1;
    second = &a2;
  } else {
    first = &a2;
    second = &a1;
  }
}

void AcedrgTables::order_angle_flanks(const CodAtomInfo& a1, const CodAtomInfo& a3,
                                             const CodAtomInfo*& flank1,
                                             const CodAtomInfo*& flank3) const {
  if (a1.hashing_value < a3.hashing_value) {
    flank1 = &a1;
    flank3 = &a3;
    return;
  }
  if (a1.hashing_value > a3.hashing_value) {
    flank1 = &a3;
    flank3 = &a1;
    return;
  }
  if (a1.cod_main.size() > a3.cod_main.size()) {
    flank1 = &a1;
    flank3 = &a3;
    return;
  }
  if (a1.cod_main.size() < a3.cod_main.size()) {
    flank1 = &a3;
    flank3 = &a1;
    return;
  }
  if (compare_no_case2(a1.cod_main, a3.cod_main)) {
    flank1 = &a1;
    flank3 = &a3;
  } else if (compare_no_case2(a3.cod_main, a1.cod_main)) {
    flank1 = &a3;
    flank3 = &a1;
  } else if (!compare_no_case2(a1.id, a3.id)) {
    flank1 = &a1;
    flank3 = &a3;
  } else {
    flank1 = &a3;
    flank3 = &a1;
  }
}

Hybridization AcedrgTables::hybrid_from_bonding_idx(int bonding_idx,
                                                           bool is_metal_atom,
                                                           int connectivity) const {
  if (is_metal_atom) {
    if (connectivity <= 4) return Hybridization::SPD5;
    if (connectivity == 5) return Hybridization::SPD5;
    if (connectivity == 6) return Hybridization::SPD6;
    if (connectivity == 7) return Hybridization::SPD7;
    return Hybridization::SPD8;
  }
  if (bonding_idx <= 0)
    return Hybridization::SP_NON;
  if (bonding_idx == 1)
    return Hybridization::SP1;
  if (bonding_idx == 2)
    return Hybridization::SP2;
  if (bonding_idx == 3)
    return Hybridization::SP3;
  if (bonding_idx == 5)
    return Hybridization::SPD5;
  if (bonding_idx == 6)
    return Hybridization::SPD6;
  if (bonding_idx == 7)
    return Hybridization::SPD7;
  if (bonding_idx >= 8)
    return Hybridization::SPD8;
  return Hybridization::SP_NON;
}

// ============================================================================
// CCP4 atom type assignment (AceDRG)
// ============================================================================

int AcedrgTables::ccp4_material_type(Element el) {
  switch (el.elem) {
    case El::H: case El::D:
      return 1;
    case El::C: case El::N: case El::O: case El::P: case El::S: case El::Se:
      return 2;
    case El::Li: case El::Na: case El::K: case El::Rb: case El::Cs: case El::Fr:
      return 3;
    case El::Be: case El::Mg: case El::Ca: case El::Sr: case El::Ba: case El::Ra:
      return 4;
    case El::Sc: case El::Y: case El::Ti: case El::Zr: case El::Hf: case El::Rf:
    case El::V: case El::Nb: case El::Ta: case El::Db:
    case El::Cr: case El::Mo: case El::W: case El::Sg:
    case El::Mn: case El::Tc: case El::Re: case El::Bh:
    case El::Fe: case El::Ru: case El::Os: case El::Hs:
    case El::Co: case El::Rh: case El::Ir: case El::Mt:
    case El::Ni: case El::Pd: case El::Pt: case El::Ds:
    case El::Cu: case El::Ag: case El::Au: case El::Rg:
    case El::Zn: case El::Cd: case El::Hg: case El::Cn:
      return 5;
    case El::Al: case El::Ga: case El::In: case El::Tl: case El::Sn:
    case El::Pb: case El::Bi:
      return 6;
    case El::B: case El::Si: case El::Ge: case El::As: case El::Sb:
    case El::Te: case El::Po:
      return 7;
    case El::F: case El::Cl: case El::Br: case El::I: case El::At:
      return 8;
    case El::La: case El::Ce: case El::Pr: case El::Nd: case El::Pm:
    case El::Sm: case El::Eu: case El::Gd: case El::Tb: case El::Dy:
    case El::Ho: case El::Er: case El::Tm: case El::Yb: case El::Lu:
    case El::Ac: case El::Th: case El::Pa: case El::U: case El::Np:
    case El::Pu: case El::Am: case El::Cm: case El::Bk: case El::Cf:
    case El::Es: case El::Fm: case El::Md: case El::No: case El::Lr:
      return 9;
    case El::He: case El::Ne: case El::Ar: case El::Kr: case El::Xe: case El::Rn:
      return 10;
    default:
      return 0;
  }
}

std::string AcedrgTables::bond_order_key(BondType type) {
  std::string s = bond_type_to_string(type);
  if (s.empty())
    s = "single";
  s = to_upper(s);
  if (s.size() > 4)
    s.resize(4);
  return s;
}

void AcedrgTables::load_ccp4_bonds(const std::string& path) {
  try {
    cif::Document doc = read_cif_gz(path);
    if (doc.blocks.empty())
      return;
    for (const auto& row : doc.blocks[0].find("_lib_bond.",
                    {"atom_type_1", "atom_type_2", "type",
                     "length", "value_esd"})) {
      std::string type1 = row.str(0);
      std::string type2 = row.str(1);
      std::string order = row.str(2);
      double length = cif::as_number(row[3], NAN);
      double sigma = cif::as_number(row[4], NAN);
      if (type1.empty() || order.empty() || std::isnan(length))
        continue;
      std::string order_key = to_upper(order);
      if (order_key.size() > 4)
        order_key.resize(4);
      ccp4_bonds_[type1][type2][order_key] = {length, sigma};
    }
  } catch (std::exception&) {
    return;
  }
}

bool AcedrgTables::search_ccp4_bond(const std::string& type1,
                                           const std::string& type2,
                                           const std::string& order,
                                           ValueStats& out) const {
  auto it1 = ccp4_bonds_.find(type1);
  if (it1 == ccp4_bonds_.end())
    return false;
  auto it2 = it1->second.find(type2);
  if (it2 == it1->second.end())
    return false;
  auto it3 = it2->second.find(order);
  if (it3 == it2->second.end())
    return false;
  out.value = it3->second.length;
  out.sigma = it3->second.sigma;
  out.count = 1;
  return true;
}

std::vector<std::string> AcedrgTables::compute_ccp4_types(
    const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    const std::vector<std::vector<int>>& neighbors) const {
  std::vector<Ccp4AtomInfo> atoms;
  atoms.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    Ccp4AtomInfo info;
    info.el = cc.atoms[i].el;
    info.chem_type = cc.atoms[i].el.name();
    info.ccp4_type = info.chem_type;
    info.bonding_idx = atom_info[i].bonding_idx;
    info.ring_rep = atom_info[i].ring_rep;
    info.conn_atoms = neighbors[i];
    info.conn_atoms_no_metal.clear();
    for (int nb : neighbors[i]) {
      if (cc.atoms[nb].is_hydrogen())
        info.conn_h_atoms.push_back(nb);
      if (!cc.atoms[nb].el.is_metal())
        info.conn_atoms_no_metal.push_back(nb);
    }
    // AceDRG treats metal-coordinated carbons with 4+ total neighbors as sp3 in CCP4 types.
    if (info.el == El::C && info.bonding_idx == 2) {
      size_t total_conn = info.conn_atoms.size();
      if (total_conn >= 4 && total_conn > info.conn_atoms_no_metal.size())
        info.bonding_idx = 3;
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(cc.atoms[i].charge));
    atoms.emplace_back(std::move(info));
  }

  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);

  std::vector<std::string> out;
  out.reserve(atoms.size());
  for (const auto& atom : atoms)
    out.push_back(atom.ccp4_type);
  return out;
}

void AcedrgTables::set_one_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                            size_t idx) {
  int ntype = ccp4_material_type(atoms[idx].el);
  switch (ntype) {
    case 1:
      set_hydro_ccp4_type(atoms, idx);
      break;
    case 2:
      set_org_ccp4_type(atoms, idx);
      break;
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
      atoms[idx].ccp4_type = atoms[idx].chem_type;
      break;
    default:
      atoms[idx].ccp4_type = atoms[idx].chem_type;
      break;
  }
  atoms[idx].ccp4_type = to_upper(atoms[idx].ccp4_type);
}

void AcedrgTables::set_hydro_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                              size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  atom.ccp4_type = "H";
  if (atom.conn_atoms.size() == 1) {
    int nb = atom.conn_atoms[0];
    if (atoms[nb].chem_type == "S")
      atom.ccp4_type = "HSH1";
  }
}

void AcedrgTables::set_org_ccp4_type(std::vector<Ccp4AtomInfo>& atoms,
                                            size_t idx) {
  Ccp4AtomInfo& atom = atoms[idx];
  int r5 = 0;
  int r6 = 0;
  for (const auto& item : atom.ring_rep) {
    if (item.second == 5)
      r5 += 1;
    if (item.second == 6)
      r6 += 1;
  }

  const size_t nconn = atom.conn_atoms_no_metal.size();
  const size_t nh = atom.conn_h_atoms.size();

  if (atom.chem_type == "C") {
    if (atom.bonding_idx == 2) {
      if (r5 && r6) {
        atom.ccp4_type = "CR56";
      } else if (r5 == 2) {
        atom.ccp4_type = "CR55";
      } else if (r6 == 2) {
        atom.ccp4_type = "CR66";
      } else if (r5 == 1) {
        if (nh == 1) {
          atom.ccp4_type = "CR15";
        } else if (nh == 0) {
          atom.ccp4_type = "CR5";
        }
      } else if (r6 == 1) {
        if (nh == 1) {
          atom.ccp4_type = "CR16";
        } else if (nh == 0) {
          atom.ccp4_type = "CR6";
        }
      } else {
        if (nh == 1) {
          atom.ccp4_type = "C1";
        } else if (nh == 2) {
          atom.ccp4_type = "C2";
        } else if (nh == 0) {
          atom.ccp4_type = "C";
        }
      }
    } else if (atom.bonding_idx == 3) {
      if (nh == 0) {
        atom.ccp4_type = "CT";
      } else if (nh == 1) {
        atom.ccp4_type = "CH1";
      } else if (nh == 2) {
        atom.ccp4_type = "CH2";
      } else if (nh == 3) {
        atom.ccp4_type = "CH3";
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "CSP";
    }
  } else if (atom.chem_type == "N") {
    if (atom.bonding_idx == 2) {
      if (nconn == 3) {
        if (nh == 1) {
          atom.ccp4_type = "NH1";
        } else if (nh == 2) {
          atom.ccp4_type = "NH2";
        } else if (nh == 0) {
          atom.ccp4_type = "NH0";
        } else {
          atom.ccp4_type = "N";
        }
      } else if (nconn == 2) {
        if (nh == 1) {
          atom.ccp4_type = "N21";
        } else if (nh == 0) {
          atom.ccp4_type = "N20";
        } else {
          atom.ccp4_type = "N";
        }
      }
    } else if (atom.bonding_idx == 3) {
      if (nconn == 4) {
        if (nh == 1) {
          atom.ccp4_type = "NT1";
        } else if (nh == 2) {
          atom.ccp4_type = "NT2";
        } else if (nh == 3) {
          atom.ccp4_type = "NT3";
        } else if (nh == 4) {
          atom.ccp4_type = "NT4";
        } else if (nh == 0) {
          atom.ccp4_type = "NT";
        } else {
          atom.ccp4_type = "N";
        }
      } else if (nconn == 3) {
        if (nh == 1) {
          atom.ccp4_type = "N31";
        } else if (nh == 2) {
          atom.ccp4_type = "N32";
        } else if (nh == 3) {
          atom.ccp4_type = "N33";
        } else if (nh == 0) {
          atom.ccp4_type = "N30";
        } else {
          atom.ccp4_type = "N3";
        }
      }
    } else if (atom.bonding_idx == 1) {
      atom.ccp4_type = "NSP";
    }
  } else if (atom.chem_type == "P") {
    atom.ccp4_type = (nconn == 4 ? "P" : "P1");
  } else if (atom.chem_type == "O") {
    bool lP = false;
    bool lS = false;
    bool lB = false;
    for (int nb : atom.conn_atoms) {
      if (atoms[nb].chem_type == "P") lP = true;
      if (atoms[nb].chem_type == "S") lS = true;
      if (atoms[nb].chem_type == "B") lB = true;
    }

    bool has_par_charge = std::fabs(atom.par_charge) > 1e-6f;
    bool has_negative_charge = atom.formal_charge < 0;

    if (atom.bonding_idx == 2) {
      // par_charge check - but only negative charges trigger OP/OS/OB/OC types
      if (has_par_charge && atom.par_charge < 0) {
        if (lP) {
          atom.ccp4_type = "OP";
        } else if (lS) {
          atom.ccp4_type = "OS";
        } else if (lB) {
          atom.ccp4_type = "OB";
        } else {
          atom.ccp4_type = "OC";
        }
      } else if (nconn == 2) {
        if (nh == 1) {
          atom.ccp4_type = "OH1";
        } else if (nh == 2) {
          atom.ccp4_type = "OH2";
        } else {
          atom.ccp4_type = "O";
        }
      } else {
        // Only negative charges trigger OP/OS/OB/OC types
        // Positive charges (like C≡O+ ligands) stay as "O"
        if (has_negative_charge) {
          if (lP) {
            atom.ccp4_type = "OP";
          } else if (lS) {
            atom.ccp4_type = "OS";
          } else if (lB) {
            atom.ccp4_type = "OB";
          } else {
            atom.ccp4_type = "OC";
          }
        } else {
          atom.ccp4_type = "O";
        }
      }
    } else if (atom.bonding_idx == 3) {
      bool lC = false;
      for (int nb : atom.conn_atoms)
        if (atoms[nb].chem_type == "C")
          lC = true;
      if (lC && nh == 1 && nconn == 2) {
        atom.ccp4_type = "OH1";
      } else if (nh == 2) {
        atom.ccp4_type = "OH2";
      } else if (nconn == 2) {
        if (has_par_charge) {
          atom.ccp4_type = "OC2";
        } else if (nh == 1) {
          atom.ccp4_type = "OH1";
        } else {
          atom.ccp4_type = "O2";
        }
      } else if (nconn == 1 && has_negative_charge) {
        // Metal-bonded oxygen with negative formal charge from valence calculation
        if (lP) {
          atom.ccp4_type = "OP";
        } else if (lS) {
          atom.ccp4_type = "OS";
        } else if (lB) {
          atom.ccp4_type = "OB";
        } else {
          atom.ccp4_type = "OC";
        }
      }
    } else if (nconn == 1) {
      if (has_negative_charge) {
        if (lP) {
          atom.ccp4_type = "OP";
        } else if (lS) {
          atom.ccp4_type = "OS";
        } else if (lB) {
          atom.ccp4_type = "OB";
        } else {
          atom.ccp4_type = "OC";
        }
      } else {
        atom.ccp4_type = "O";
      }
    } else {
      atom.ccp4_type = "O";
    }
  } else if (atom.chem_type == "S") {
    if (nconn == 3 || nconn == 4) {
      if (nh == 0) {
        atom.ccp4_type = "S3";
      } else if (nh == 1) {
        atom.ccp4_type = "SH1";
      }
    } else if (nconn == 2) {
      if (nh == 0) {
        atom.ccp4_type = "S2";
      } else {
        atom.ccp4_type = "SH1";
      }
    } else if (nconn == 1) {
      atom.ccp4_type = "S1";
    } else if (nh == 0) {
      atom.ccp4_type = "S";
    } else if (nh == 1) {
      atom.ccp4_type = "SH1";
    } else {
      atom.ccp4_type = "S";
    }
  } else if (atom.chem_type == "Se") {
    // Selenium - similar to sulfur
    if (nconn == 3 || nconn == 4) {
      atom.ccp4_type = "SE";
    } else if (nconn == 2) {
      atom.ccp4_type = "SE";
    } else if (nconn == 1) {
      atom.ccp4_type = "SE";
    } else {
      atom.ccp4_type = "SE";
    }
  } else {
    atom.ccp4_type = atom.chem_type;
  }
}

void AcedrgTables::assign_ccp4_types(ChemComp& cc) const {
  if (cc.atoms.empty())
    return;

  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);
  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);

  std::vector<Ccp4AtomInfo> atoms;
  atoms.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    Ccp4AtomInfo info;
    info.el = cc.atoms[i].el;
    info.chem_type = cc.atoms[i].el.name();
    info.ccp4_type = info.chem_type;
    info.bonding_idx = atom_info[i].bonding_idx;
    info.ring_rep = atom_info[i].ring_rep;
    info.conn_atoms = neighbors[i];
    info.conn_atoms_no_metal.clear();
    for (int nb : neighbors[i]) {
      if (cc.atoms[nb].is_hydrogen())
        info.conn_h_atoms.push_back(nb);
      if (!cc.atoms[nb].el.is_metal())
        info.conn_atoms_no_metal.push_back(nb);
    }
    // AceDRG treats metal-coordinated carbons with 4+ total neighbors as sp3 in CCP4 types.
    if (info.el == El::C && info.bonding_idx == 2) {
      size_t total_conn = info.conn_atoms.size();
      if (total_conn >= 4 && total_conn > info.conn_atoms_no_metal.size())
        info.bonding_idx = 3;
    }
    info.par_charge = cc.atoms[i].charge;
    info.formal_charge = static_cast<int>(std::round(cc.atoms[i].charge));
    atoms.emplace_back(std::move(info));
  }

  // Valence-based charge calculation for atoms bonded to metals.
  // Acedrg computes ligand charges via valence bookkeeping: for non-metal atoms,
  // sum only non-metal bond orders, then set charge = -(expectedValence - sumBo).
  // This ensures metal-bonded oxygens (like O bonded to both C and Hg) get
  // formal_charge=-1 and thus type "OC".
  // Note: atoms bonded ONLY to metals (no non-metal neighbors except H) are excluded
  // from this adjustment - their type stays "O" even with formal charge.
  for (size_t i = 0; i < atoms.size(); ++i) {
    const Ccp4AtomInfo& info = atoms[i];
    // Skip metals, hydrogens, and atoms already charged
    if (info.el.is_metal() || info.el == El::H || info.formal_charge != 0)
      continue;
    // Check if this atom has any metal neighbors
    bool has_metal_neighbor = info.conn_atoms.size() > info.conn_atoms_no_metal.size();
    if (!has_metal_neighbor)
      continue;
    // Skip if bonded only to metals (and possibly H) - no non-metal heavy atom neighbors
    bool has_non_metal_heavy_neighbor = false;
    for (int nb : info.conn_atoms_no_metal) {
      if (!cc.atoms[nb].is_hydrogen()) {
        has_non_metal_heavy_neighbor = true;
        break;
      }
    }
    if (!has_non_metal_heavy_neighbor)
      continue;
    // Get expected valence for common elements
    int expected_valence = 0;
    if (info.el == El::O) expected_valence = 2;
    else if (info.el == El::N) expected_valence = 3;
    else if (info.el == El::S) expected_valence = 2;
    else if (info.el == El::C) expected_valence = 4;
    else continue;  // Skip elements we don't have valence info for
    // Sum bond orders to non-metal neighbors only
    float sum_bo = 0.0f;
    for (const BondInfo& bi : adjacency[i]) {
      if (!cc.atoms[bi.neighbor_idx].el.is_metal())
        sum_bo += order_of_bond_type(bi.type);
    }
    // Calculate remaining valence and set formal charge
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    if (rem_v != 0)
      atoms[i].formal_charge = -rem_v;
  }

  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) != 1)
      set_one_ccp4_type(atoms, i);
  for (size_t i = 0; i < atoms.size(); ++i)
    if (ccp4_material_type(atoms[i].el) == 1)
      set_one_ccp4_type(atoms, i);

  for (size_t i = 0; i < atoms.size(); ++i)
    cc.atoms[i].chem_type = atoms[i].ccp4_type;
}

void AcedrgTables::apply_metal_charge_corrections(ChemComp& cc) const {
  if (cc.atoms.empty())
    return;

  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);

  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    const Element& el = cc.atoms[i].el;
    if (el.is_metal() || el == El::H)
      continue;

    int non_metal_conn = 0;
    for (int nb : neighbors[i]) {
      if (!cc.atoms[nb].el.is_metal())
        ++non_metal_conn;
    }
    if (non_metal_conn == static_cast<int>(neighbors[i].size()))
      continue;  // no metal neighbors

    bool has_non_metal_heavy_neighbor = false;
    for (int nb : neighbors[i]) {
      if (!cc.atoms[nb].el.is_metal() && !cc.atoms[nb].is_hydrogen()) {
        has_non_metal_heavy_neighbor = true;
        break;
      }
    }
    if (!has_non_metal_heavy_neighbor)
      continue;

    int expected_valence = 0;
    if (el == El::O) expected_valence = 2;
    else if (el == El::N) expected_valence = 3;
    else if (el == El::S) expected_valence = 2;
    else if (el == El::C) expected_valence = 4;
    else if (el == El::P) expected_valence = 3;
    else continue;

    float sum_bo = 0.0f;
    for (const BondInfo& bi : adjacency[i]) {
      if (!cc.atoms[bi.neighbor_idx].el.is_metal())
        sum_bo += order_of_bond_type(bi.type);
    }
    int rem_v = expected_valence - static_cast<int>(std::round(sum_bo));
    cc.atoms[i].charge = static_cast<float>(-rem_v);
  }
}

// ============================================================================
// Implementation - Bond search
// ============================================================================

void AcedrgTables::fill_restraints(ChemComp& cc) const {
  if (!tables_loaded_)
    return;

  // Classify atoms
  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);
  std::vector<std::vector<BondInfo>> adjacency = build_adjacency(cc);
  std::vector<std::vector<int>> neighbors = build_neighbors(adjacency);
  std::vector<std::string> ccp4_types;
  if (!ccp4_bonds_.empty())
    ccp4_types = compute_ccp4_types(cc, atom_info, neighbors);

  // Print atom classification if verbose (level 1 = bonds only, level 2+ = atoms too)
  if (verbose >= 1) {
    std::fprintf(stderr, "  Atom classification:\n");
    for (size_t i = 0; i < atom_info.size(); ++i) {
      const auto& a = atom_info[i];
      std::fprintf(stderr, "    %s: el=%s conn=%d ring=%d arom=%d hybr=%s "
                   "bonding_idx=%d hash=%d cod_class=%s\n",
                   a.id.c_str(), element_name(a.el), a.connectivity,
                   a.min_ring_size, a.is_aromatic ? 1 : 0,
                   hybridization_to_string(a.hybrid), a.bonding_idx,
                   a.hashing_value, a.cod_class.c_str());
      if (verbose >= 2)
        std::fprintf(stderr, "        cod_main=%s nb1nb2_sp=%s nb=%s nb2=%s\n",
                     a.cod_main.c_str(), a.nb1nb2_sp.c_str(),
                     a.nb_symb.c_str(), a.nb2_symb.c_str());
    }
  }

  // Fill bonds
  for (auto& bond : cc.rt.bonds) {
    if (std::isnan(bond.value)) {
      int match_level = fill_bond(cc, atom_info, bond);
      int idx1 = cc.find_atom_index(bond.id1.atom);
      int idx2 = cc.find_atom_index(bond.id2.atom);
      auto apply_ccp4 = [&](bool override_existing) {
        if (ccp4_types.empty() || idx1 < 0 || idx2 < 0)
          return false;
        if (!override_existing && !std::isnan(bond.value))
          return false;
        std::string order = bond_order_key(bond.type);
        ValueStats vs;
        if (search_ccp4_bond(ccp4_types[idx1], ccp4_types[idx2], order, vs) ||
            search_ccp4_bond(ccp4_types[idx2], ccp4_types[idx1], order, vs)) {
          bond.value = vs.value;
          bond.esd = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
          return true;
        }
        ValueStats v1, v2;
        bool has1 = search_ccp4_bond(ccp4_types[idx1], ".", order, v1);
        bool has2 = search_ccp4_bond(ccp4_types[idx2], ".", order, v2);
        if (has1 && has2) {
          bond.value = 0.5 * (v1.value + v2.value);
          bond.esd = 0.02;
          return true;
        }
        return false;
      };
      // CCP4 energetic library fallback
      if (std::isnan(bond.value))
        apply_ccp4(false);
      // CCP4 DELO override for explicit delocalized bonds only.
      // Applying this to regular single/double C-O bonds causes mismatches with AceDRG.
      if (!std::isnan(bond.value) && !ccp4_types.empty() &&
          bond.type == BondType::Deloc && match_level < 4) {
        if (idx1 >= 0 && idx2 >= 0) {
          const std::string& t1 = ccp4_types[idx1];
          const std::string& t2 = ccp4_types[idx2];
          // Check for carboxylate C-O bond (C bonded to O or OC)
          bool is_c_o_bond = (t1 == "C" && (t2 == "O" || t2 == "OC")) ||
                             (t2 == "C" && (t1 == "O" || t1 == "OC"));
          if (is_c_o_bond) {
            // Check if this C has two terminal O neighbors (carboxylate pattern)
            int c_idx = (t1 == "C") ? idx1 : idx2;
            int terminal_o_count = 0;
            for (const auto& b : cc.rt.bonds) {
              if (b.id1.atom == cc.atoms[c_idx].id || b.id2.atom == cc.atoms[c_idx].id) {
                const std::string& other_name = (b.id1.atom == cc.atoms[c_idx].id) ? b.id2.atom : b.id1.atom;
                auto it_other = cc.find_atom(other_name);
                if (it_other != cc.atoms.end() && it_other->el == El::O) {
                  // Check if this oxygen is terminal (only bonded to C and possibly H)
                  int heavy_neighbors = 0;
                  for (const auto& b2 : cc.rt.bonds) {
                    if (b2.id1.atom == other_name || b2.id2.atom == other_name) {
                      const std::string& nb = (b2.id1.atom == other_name) ? b2.id2.atom : b2.id1.atom;
                      auto it_nb = cc.find_atom(nb);
                      if (it_nb != cc.atoms.end() && it_nb->el != El::H)
                        ++heavy_neighbors;
                    }
                  }
                  if (heavy_neighbors == 1)  // Only bonded to C (the carboxyl carbon)
                    ++terminal_o_count;
                }
              }
            }
            if (terminal_o_count >= 2) {
              // True carboxylate pattern - use DELO value
              ValueStats vs;
              if (search_ccp4_bond("OC", "C", "DELO", vs) ||
                  search_ccp4_bond("C", "OC", "DELO", vs)) {
                bond.value = vs.value;
                bond.esd = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
              }
            }
          }
        }
      }
    }
  }

  // Populate value_nucleus/esd_nucleus for X-H bonds from prot_hydr_dists
  for (auto& bond : cc.rt.bonds) {
    int idx1 = cc.find_atom_index(bond.id1.atom);
    int idx2 = cc.find_atom_index(bond.id2.atom);
    if (idx1 < 0 || idx2 < 0)
      continue;
    const CodAtomInfo& a1 = atom_info[idx1];
    const CodAtomInfo& a2 = atom_info[idx2];
    if (a1.el == El::H || a2.el == El::H) {
      const CodAtomInfo& h_atom = (a1.el == El::H) ? a1 : a2;
      const CodAtomInfo& heavy_atom = (a1.el == El::H) ? a2 : a1;
      ProtHydrDist phd = search_prot_hydr_dist(h_atom, heavy_atom);
      if (!std::isnan(phd.nucleus_val)) {
        bond.value_nucleus = phd.nucleus_val;
        bond.esd_nucleus = phd.nucleus_sigma;
      }
    }
  }

  // Fill angles
  for (auto& angle : cc.rt.angles) {
    if (std::isnan(angle.value)) {
      fill_angle(cc, atom_info, angle);
    }
  }

  // AceDRG adjustment: enforce planar ring angle sum ((n-2)*180/n) for SP2 rings.
  auto atom_index = make_atom_index(cc);
  std::map<int, std::vector<size_t>> rings;
  for (size_t i = 0; i < atom_info.size(); ++i)
    for (int ring_id : atom_info[i].in_rings)
      rings[ring_id].push_back(i);
  for (const auto& ring : rings) {
    const std::vector<size_t>& ring_atoms = ring.second;
    if (ring_atoms.size() < 3)
      continue;
    bool planar = true;
    for (size_t idx : ring_atoms) {
      if (atom_info[idx].bonding_idx == 3) {
        planar = false;
        break;
      }
    }
    if (!planar)
      continue;
    std::set<std::string> ring_ids;
    for (size_t idx : ring_atoms)
      ring_ids.insert(cc.atoms[idx].id);
    std::vector<size_t> ring_angles;
    ring_angles.reserve(ring_atoms.size());
    for (size_t i = 0; i < cc.rt.angles.size(); ++i) {
      const auto& ang = cc.rt.angles[i];
      if (ring_ids.count(ang.id2.atom) &&
          ring_ids.count(ang.id1.atom) &&
          ring_ids.count(ang.id3.atom))
        ring_angles.push_back(i);
    }
    if (ring_angles.size() != ring_atoms.size())
      continue;
    double sum = 0.0;
    for (size_t idx : ring_angles)
      sum += cc.rt.angles[idx].value;
    double mean = sum / static_cast<double>(ring_angles.size());
    double target = (static_cast<double>(ring_atoms.size()) - 2.0) *
                    180.0 / static_cast<double>(ring_atoms.size());
    double shift = target - mean;
    for (size_t idx : ring_angles)
      cc.rt.angles[idx].value += shift;
  }

  // AceDRG adjustment: enforce 360-degree sum for sp2 centers with 3 angles,
  // keeping ring angles fixed (checkRingAngleConstraints behavior).
  std::map<std::string, std::vector<size_t>> center_angles;
  for (size_t i = 0; i < cc.rt.angles.size(); ++i)
    center_angles[cc.rt.angles[i].id2.atom].push_back(i);
  for (const auto& entry : center_angles) {
    if (entry.second.size() != 3)
      continue;
    auto it = atom_index.find(entry.first);
    if (it == atom_index.end())
      continue;
    size_t center_idx = it->second;
    if (atom_info[center_idx].hybrid != Hybridization::SP2)
      continue;
    std::vector<size_t> fixed, free;
    for (size_t idx_ang : entry.second) {
      const auto& ang = cc.rt.angles[idx_ang];
      auto it1 = atom_index.find(ang.id1.atom);
      auto it3 = atom_index.find(ang.id3.atom);
      if (it1 == atom_index.end() || it3 == atom_index.end())
        continue;
      const CodAtomInfo& a1 = atom_info[it1->second];
      const CodAtomInfo& a3 = atom_info[it3->second];
      if (angle_ring_size(atom_info[center_idx], a1, a3) > 0)
        fixed.push_back(idx_ang);
      else
        free.push_back(idx_ang);
    }
    if (free.empty())
      continue;
    double fixed_sum = 0.0;
    for (size_t idx_ang : fixed)
      fixed_sum += cc.rt.angles[idx_ang].value;
    double free_sum = 0.0;
    for (size_t idx_ang : free)
      free_sum += cc.rt.angles[idx_ang].value;
    double diff = (360.0 - fixed_sum - free_sum) / static_cast<double>(free.size());
    if (std::fabs(diff) > 0.01) {
      double new_sum = 0.0;
      for (size_t idx_ang : free) {
        cc.rt.angles[idx_ang].value += diff;
        new_sum += cc.rt.angles[idx_ang].value;
      }
      cc.rt.angles[free[0]].value += (360.0 - fixed_sum - new_sum);
    } else {
      cc.rt.angles[free[0]].value += diff;
    }
  }
}

int AcedrgTables::fill_bond(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Bond& bond) const {

  int idx1 = cc.find_atom_index(bond.id1.atom);
  int idx2 = cc.find_atom_index(bond.id2.atom);
  if (idx1 < 0 || idx2 < 0)
    return 0;

  const CodAtomInfo& a1 = atom_info[idx1];
  const CodAtomInfo& a2 = atom_info[idx2];

  const char* source = "no_match";

  // Check for metal bond
  if (a1.is_metal || a2.is_metal) {
    const CodAtomInfo& metal = a1.is_metal ? a1 : a2;
    const CodAtomInfo& ligand = a1.is_metal ? a2 : a1;

    std::string mkey = metal.el.name();
    std::string lkey = ligand.el.name();
    for (char& c : mkey)
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    for (char& c : lkey)
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    auto it_m = covalent_radii_.find(mkey);
    auto it_l = covalent_radii_.find(lkey);
    if (it_m != covalent_radii_.end() && it_l != covalent_radii_.end()) {
      bond.value = it_m->second + it_l->second;
      bond.esd = 0.04;
      source = "metal_cova";
      if (verbose)
        std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f)\n",
                     bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                     a1.hashing_value, a2.hashing_value,
                     hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                     source, bond.value, bond.esd);
      return 10;
    }

    ValueStats vs = search_metal_bond(metal, ligand, atom_info);
    if (vs.count > 0) {
      bond.value = vs.value;
      double sigma = std::isnan(vs.sigma) ? 0.02 : vs.sigma;
      bond.esd = std::max(0.02, clamp_bond_sigma(sigma));
      source = "metal";
      if (verbose)
        std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                     bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                     a1.hashing_value, a2.hashing_value,
                     hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                     source, bond.value, bond.esd, vs.count);
      return 10;  // metal bond - treat as high-specificity match
    }
  }

  // Try both multilevel and HRS
  bool same_ring = are_in_same_ring(a1, a2);
  ValueStats vs_ml = search_bond_multilevel(a1, a2);
  ValueStats vs_hrs = search_bond_hrs(a1, a2, same_ring);

  // Acedrg's logic: use multilevel when it matches with sufficient threshold.
  // Acedrg iterates from start_level upward and uses the first level that meets threshold.
  // The threshold check is done inside search_bond_multilevel (entry count for levels 1-8,
  // and no threshold for HRS-like levels 9-11). When it returns with level >= 0, the threshold was met.
  // There's no special preference for HRS over type-based (levels 0-2) matches.
  ValueStats vs;
  vs.level = -1;  // no match by default; prevents false acceptance below
  if (vs_ml.level >= 0) {
    // Multilevel matched with threshold met - use it
    vs = vs_ml;
    source = "multilevel";
  } else if (vs_hrs.count > 0) {
    // Multilevel didn't match - try HRS
    vs = vs_hrs;
    source = "HRS";
  }

  bool is_hrs = (source && std::strcmp(source, "HRS") == 0);
  bool accept = false;
  if (is_hrs || vs.level >= 9) {
    accept = vs.count > 0;
  } else {
    // For multilevel levels 0-8, thresholding is handled inside search_bond_multilevel().
    accept = vs.level >= 0;
  }
  if (accept) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    if (verbose)
      std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                   bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                   a1.hashing_value, a2.hashing_value,
                   hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                   source, bond.value, bond.esd, vs.count);
    return vs.level;  // return match level for multilevel, or 0 for HRS
  }

  // Try element+hybridization fallback
  vs = search_bond_en(a1, a2);
  if (vs.count > 0) {
    bond.value = vs.value;
    bond.esd = clamp_bond_sigma(vs.sigma);
    source = "EN";
    if (verbose)
      std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s → %s (%.3f, %.3f, n=%d)\n",
                   bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                   a1.hashing_value, a2.hashing_value,
                   hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                   source, bond.value, bond.esd, vs.count);
    return 0;  // fallback - no good type-specific match
  }

  // No match found - leave bond.value as NaN so CCP4 fallback in fill_restraints
  // can be applied. AceDRG's search order is: multilevel -> HRS -> EN -> CCP4.
  if (verbose)
    std::fprintf(stderr, "  bond %s-%s: hash %d-%d hybr %s-%s -> %s\n",
                 bond.id1.atom.c_str(), bond.id2.atom.c_str(),
                 a1.hashing_value, a2.hashing_value,
                 hybridization_to_string(a1.hybrid), hybridization_to_string(a2.hybrid),
                 source);
  return 0;  // no type-specific match
}

ValueStats AcedrgTables::search_bond_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& a2) const {

  if (verbose >= 2)
    std::fprintf(stderr, "      search_bond_multilevel: input a1=%s(hash=%d) a2=%s(hash=%d)\n",
                 a1.id.c_str(), a1.hashing_value, a2.id.c_str(), a2.hashing_value);

  const CodAtomInfo* left = nullptr;
  const CodAtomInfo* right = nullptr;
  order_bond_atoms(a1, a2, left, right);

  if (verbose >= 2)
    std::fprintf(stderr, "      after order: left=%s(hash=%d) right=%s(hash=%d)\n",
                 left->id.c_str(), left->hashing_value, right->id.c_str(), right->hashing_value);

  // Build lookup keys
  int ha1 = left->hashing_value;
  int ha2 = right->hashing_value;

  std::string h1 = hybridization_to_string(left->hybrid);
  std::string h2 = hybridization_to_string(right->hybrid);
  if (h1 > h2) std::swap(h1, h2);
  std::string hybr_comb = h1 + "_" + h2;

  std::string in_ring = are_in_same_ring(a1, a2) ? "Y" : "N";
  // AceDRG: if requested ring key is missing, fall back to the other (Y/N).
  auto has_ring_key = [&](const std::string& key) -> bool {
    auto it1 = bond_idx_full_.find(ha1);
    if (it1 == bond_idx_full_.end()) return false;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return false;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return false;
    return it3->second.find(key) != it3->second.end();
  };
  if (!has_ring_key(in_ring)) {
    std::string alt = (in_ring == "Y") ? "N" : "Y";
    if (has_ring_key(alt)) {
      if (verbose >= 2)
        std::fprintf(stderr, "      ring key '%s' missing, using '%s'\n",
                     in_ring.c_str(), alt.c_str());
      in_ring = alt;
    }
  }

  // Use neighbor symbols as additional keys.
  std::string a1_nb2 = left->nb2_symb;
  std::string a2_nb2 = right->nb2_symb;
  std::string a1_nb = left->nb1nb2_sp;
  std::string a2_nb = right->nb1nb2_sp;

  // Use COD main type as atom type (for 1D lookup).
  std::string a1_type = left->cod_main;
  std::string a2_type = right->cod_main;

  // Full COD class keys (as generated, AceDRG uses codClass directly).
  std::string a1_class = left->cod_class_no_charge;
  std::string a2_class = right->cod_class_no_charge;

  int num_th = min_observations_bond;
  if (a1.el == El::As || a2.el == El::As || a1.el == El::Ge || a2.el == El::Ge)
    num_th = 1;

  if (verbose >= 2)
    std::fprintf(stderr, "      lookup: hash=%d/%d hybr=%s ring=%s nb1nb2_sp=%s/%s nb2=%s/%s type=%s/%s\n",
                 ha1, ha2, hybr_comb.c_str(), in_ring.c_str(),
                 a1_nb.c_str(), a2_nb.c_str(), a1_nb2.c_str(), a2_nb2.c_str(),
                 a1_type.c_str(), a2_type.c_str());

  // AceDRG first tries the exact full codClass match (approx level 0) before
  // entering inter-level fallback.
  bool has_a1_class_only = false;  // a1 class exists, a2 class missing
  {
    auto it1 = bond_idx_full_.find(ha1);
    if (it1 != bond_idx_full_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(hybr_comb);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(in_ring);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_nb2);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_nb2);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a1_nb);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a2_nb);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a1_type);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a2_type);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(a1_class);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(a2_class);
                          if (it12 != it11->second.end()) {
                            if (it12->second.count >= num_th) {
                              ValueStats exact = it12->second;
                              exact.level = 0;
                              if (verbose >= 2)
                                std::fprintf(stderr,
                                             "      matched exact codClass level=0: value=%.3f sigma=%.3f count=%d\n",
                                             exact.value, exact.sigma, exact.count);
                              return exact;
                            }
                          } else {
                            // Special AceDRG fallback only when a2 class is absent.
                            has_a1_class_only = true;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  const auto* map_1d = [&]() -> const std::map<std::string, std::map<std::string, std::vector<ValueStats>>>* {
    auto it1 = bond_idx_1d_.find(ha1);
    if (it1 == bond_idx_1d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    auto it5 = it4->second.find(a1_nb2);
    if (it5 == it4->second.end()) return nullptr;
    auto it6 = it5->second.find(a2_nb2);
    if (it6 == it5->second.end()) return nullptr;
    auto it7 = it6->second.find(a1_nb);
    if (it7 == it6->second.end()) return nullptr;
    auto it8 = it7->second.find(a2_nb);
    if (it8 == it7->second.end()) return nullptr;
    return &it8->second;
  }();

  const auto* map_2d = [&]() -> const std::map<std::string, std::map<std::string, std::vector<ValueStats>>>* {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 == bond_idx_2d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    auto it5 = it4->second.find(a1_nb2);
    if (it5 == it4->second.end()) return nullptr;
    auto it6 = it5->second.find(a2_nb2);
    if (it6 == it5->second.end()) return nullptr;
    return &it6->second;
  }();

  // Keep the signal only for level-gating compatibility with AceDRG's
  // class-missing branches. Do not short-circuit to direct 2D stats here:
  // AceDRG still proceeds through start-level logic in these cases.
  (void) has_a1_class_only;

  const auto* map_nb2 = [&]() -> const std::map<std::string, std::map<std::string,
                                      std::map<std::string, std::map<std::string, std::vector<ValueStats>>>>>* {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 == bond_idx_2d_.end()) return nullptr;
    auto it2 = it1->second.find(ha2);
    if (it2 == it1->second.end()) return nullptr;
    auto it3 = it2->second.find(hybr_comb);
    if (it3 == it2->second.end()) return nullptr;
    auto it4 = it3->second.find(in_ring);
    if (it4 == it3->second.end()) return nullptr;
    return &it4->second;
  }();

  // Determine key availability for gating (mirror AceDRG's presence checks).
  bool has_hybr = false;
  bool has_in_ring = false;
  bool has_a1_nb2 = false;
  bool has_a2_nb2 = false;
  bool has_a1_nb = false;
  bool has_a2_nb = false;
  bool has_a1_type = false;
  bool has_a2_type = false;
  {
    auto it1 = bond_idx_2d_.find(ha1);
    if (it1 != bond_idx_2d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(hybr_comb);
        if (it3 != it2->second.end()) {
          has_hybr = true;
          auto it4 = it3->second.find(in_ring);
          if (it4 != it3->second.end()) {
            has_in_ring = true;
            auto it5 = it4->second.find(a1_nb2);
            if (it5 != it4->second.end()) {
              has_a1_nb2 = true;
              auto it6 = it5->second.find(a2_nb2);
              if (it6 != it5->second.end()) {
                has_a2_nb2 = true;
                auto it7 = it6->second.find(a1_nb);
                if (it7 != it6->second.end()) {
                  has_a1_nb = true;
                  auto it8 = it7->second.find(a2_nb);
                  if (it8 != it7->second.end())
                    has_a2_nb = true;
                }
              }
            }
          }
        }
      }
    }
  }
  if (map_1d) {
    auto it1 = map_1d->find(a1_type);
    if (it1 != map_1d->end()) {
      has_a1_type = true;
      if (it1->second.find(a2_type) != it1->second.end())
        has_a2_type = true;
    }
  }

  // Determine start level based on available keys (matching acedrg's dynamic logic from codClassify.cpp)
  // acedrg iterates from start_level upward (0→1→2→...) until threshold is met
  int start_level = 0;
  if (!has_a2_type) start_level = std::max(start_level, 1);   // a2M missing
  if (!has_a1_type) start_level = std::max(start_level, 2);   // a1M missing
  if (!has_a2_nb) start_level = std::max(start_level, 4);     // a2NB missing
  if (!has_a1_nb) start_level = std::max(start_level, 5);     // a1NB missing
  if (!has_a2_nb2) start_level = std::max(start_level, 7);    // a2NB2 missing
  if (!has_a1_nb2) start_level = std::max(start_level, 8);    // a1NB2 missing
  if (!has_in_ring) start_level = std::max(start_level, 10);
  if (!has_hybr) start_level = std::max(start_level, 11);

  if (verbose >= 2)
    std::fprintf(stderr, "      has: hybr=%d ring=%d a1_nb2=%d a2_nb2=%d a1_nb=%d a2_nb=%d a1_type=%d a2_type=%d start_level=%d\n",
                 has_hybr, has_in_ring, has_a1_nb2, has_a2_nb2, has_a1_nb, has_a2_nb, has_a1_type, has_a2_type, start_level);

  // AceDRG-like multilevel fallback: iterate from start_level upward until threshold is met.
  // For the special "a2 class missing" branch, AceDRG continues with post-search overrides,
  // so we must not return early here.
  ValueStats matched_vs;
  bool have_matched_vs = false;
  for (int level = start_level; level < 12; level++) {
    // Skip levels that require keys we don't have
    if (level <= 2 && !map_1d) continue;  // type levels need map_1d
    if (level >= 3 && level <= 6 && !map_2d) continue;  // nb levels need map_2d
    if ((level == 7 || level == 8) && !map_nb2) continue;  // nb2 levels need map_nb2
    if (level == 9 && !has_in_ring) continue;
    if (level == 10 && !has_hybr) continue;
    ValueStats vs;
    int values_size = 0;

    if (level == 0) {
      if (map_1d) {
        auto it1 = map_1d->find(a1_type);
        if (it1 != map_1d->end()) {
          auto it2 = it1->second.find(a2_type);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            vs = it2->second.front();
            values_size = 1;
          }
        }
      }
    } else if (level == 1) {
      // Level 1: a1_type matches exactly, aggregate over all a2_types
      if (map_1d) {
        std::vector<ValueStats> values;
        auto it1 = map_1d->find(a1_type);
        if (it1 != map_1d->end()) {
          for (const auto& it2 : it1->second) {
            if (!it2.second.empty())
              values.push_back(it2.second.front());
          }
        }
        // Add entries where a2_type matches but a1_type differs (AceDRG tLev==1)
        for (const auto& it1b : *map_1d) {
          if (it1b.first == a1_type)
            continue;
          auto it2 = it1b.second.find(a2_type);
          if (it2 != it1b.second.end() && !it2->second.empty())
            values.push_back(it2->second.front());
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 2) {
      // Level 2: a2_type matches exactly, aggregate over different a1_types
      if (map_1d) {
        std::vector<ValueStats> values;
        for (const auto& it1 : *map_1d) {
          if (it1.first != a1_type) {
            auto it2 = it1.second.find(a2_type);
            if (it2 != it1.second.end() && !it2->second.empty())
              values.push_back(it2->second.front());
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 3) {
      if (map_2d) {
        auto it1 = map_2d->find(a1_nb);
        if (it1 != map_2d->end()) {
          auto it2 = it1->second.find(a2_nb);
          if (it2 != it1->second.end() && !it2->second.empty()) {
            values_size = static_cast<int>(it2->second.size());
            vs = aggregate_stats(it2->second);
            if (verbose >= 2)
              std::fprintf(stderr, "      level 3: found %d entries, value=%.3f count=%d (need %d)\n",
                           values_size, vs.value, vs.count, num_th);
          }
        }
        if (verbose >= 2 && values_size == 0)
          std::fprintf(stderr, "      level 3 miss: a1_nb='%s' a2_nb='%s'\n", a1_nb.c_str(), a2_nb.c_str());
      }
    } else if (level == 4 || level == 5) {
      if (map_2d) {
        std::vector<ValueStats> values;
        // AceDRG level 4: combine entries with a1NB plus entries with a2NB but not a1NB.
        // AceDRG level 5: only entries with a2NB but not a1NB.
        if (level == 4) {
          auto it1 = map_2d->find(a1_nb);
          if (it1 != map_2d->end()) {
            for (const auto& it2 : it1->second)
              for (const auto& v : it2.second)
                values.push_back(v);
          }
        }
        for (const auto& it1 : *map_2d) {
          if (it1.first == a1_nb)
            continue;
          auto it2 = it1.second.find(a2_nb);
          if (it2 != it1.second.end()) {
            for (const auto& v : it2->second)
              values.push_back(v);
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 6) {
      if (map_2d) {
        std::vector<ValueStats> values;
        for (const auto& it1 : *map_2d)
          for (const auto& it2 : it1.second)
            for (const auto& v : it2.second)
              values.push_back(v);
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 7 || level == 8) {
      if (map_nb2) {
        std::vector<ValueStats> values;
        if (level == 7) {
          auto it1 = map_nb2->find(a1_nb2);
          if (it1 != map_nb2->end()) {
            for (const auto& it2 : it1->second)
              for (const auto& it3 : it2.second)
                for (const auto& it4 : it3.second)
                  for (const auto& v : it4.second)
                    values.push_back(v);
          }
        }
        for (const auto& it2 : *map_nb2) {
          if (it2.first == a1_nb2)
            continue;
          auto it3 = it2.second.find(a2_nb2);
          if (it3 != it2.second.end()) {
            for (const auto& it4 : it3->second)
              for (const auto& it5 : it4.second)
                for (const auto& v : it5.second)
                  values.push_back(v);
          }
        }
        values_size = static_cast<int>(values.size());
        if (!values.empty())
          vs = aggregate_stats(values);
      }
    } else if (level == 9) {
      auto it1 = bond_hasp_2d_.find(ha1);
      if (it1 != bond_hasp_2d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end()) {
          auto it3 = it2->second.find(hybr_comb);
          if (it3 != it2->second.end()) {
            auto it4 = it3->second.find(in_ring);
            if (it4 != it3->second.end() && !it4->second.empty()) {
              values_size = static_cast<int>(it4->second.size());
              vs = aggregate_stats(it4->second);
            }
          }
        }
      }
    } else if (level == 10) {
      auto it1 = bond_hasp_1d_.find(ha1);
      if (it1 != bond_hasp_1d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end()) {
          auto it3 = it2->second.find(hybr_comb);
          if (it3 != it2->second.end() && !it3->second.empty()) {
            values_size = static_cast<int>(it3->second.size());
            vs = aggregate_stats(it3->second);
          }
        }
      }
    } else if (level == 11) {
      auto it1 = bond_hasp_0d_.find(ha1);
      if (it1 != bond_hasp_0d_.end()) {
        auto it2 = it1->second.find(ha2);
        if (it2 != it1->second.end() && !it2->second.empty()) {
          values_size = static_cast<int>(it2->second.size());
          vs = aggregate_stats(it2->second);
        }
      }
    }

    // Check threshold (matching acedrg's logic)
    // Levels 1-8 use entry count threshold, level 0 and others use observation count
    bool threshold_met = false;
    if (level >= 1 && level <= 8) {
      threshold_met = values_size >= num_th;  // entry count
    } else if (level >= 9) {
      threshold_met = values_size > 0;        // HRS entries only require presence
    } else {
      threshold_met = vs.count >= num_th;     // observation count
    }
    if (values_size > 0 && threshold_met) {
      if (verbose >= 2)
        std::fprintf(stderr, "      matched: level=%d\n", level);
      vs.level = level;
      if (!has_a1_class_only)
        return vs;
      matched_vs = vs;
      have_matched_vs = true;
      break;
    }
  }

  // AceDRG branch (searchCodOrgBonds2_2, a2C missing):
  // after inter-level search, prefer exact (a1M,a2M)[0] if enough observations,
  // otherwise use exact (a1NB1NB2,a2NB1NB2) 2D stats regardless threshold.
  if (has_a1_class_only) {
    if (map_1d) {
      auto it1 = map_1d->find(a1_type);
      if (it1 != map_1d->end()) {
        auto it2 = it1->second.find(a2_type);
        if (it2 != it1->second.end() && !it2->second.empty()) {
          ValueStats exact_1d = it2->second.front();
          if (exact_1d.count >= num_th) {
            exact_1d.level = 0;
            return exact_1d;
          }
        }
      }
    }
    if (map_2d) {
      auto it1 = map_2d->find(a1_nb);
      if (it1 != map_2d->end()) {
        auto it2 = it1->second.find(a2_nb);
        if (it2 != it1->second.end() && !it2->second.empty()) {
          ValueStats exact_2d = aggregate_stats(it2->second);
          exact_2d.level = 3;
          return exact_2d;
        }
      }
    }
    if (have_matched_vs)
      return matched_vs;
  }

  if (have_matched_vs)
    return matched_vs;

  // No match found in detailed tables
  if (verbose >= 2)
    std::fprintf(stderr, "      matched: level=none (no multilevel match)\n");
  ValueStats no_match;
  no_match.level = -1;  // sentinel for no match
  return no_match;
}

ValueStats AcedrgTables::search_bond_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& a2, bool in_ring) const {

  const CodAtomInfo* left = nullptr;
  const CodAtomInfo* right = nullptr;
  order_bond_atoms(a1, a2, left, right);

  BondHRSKey key;
  key.hash1 = left->hashing_value;
  key.hash2 = right->hashing_value;

  std::string h1 = hybridization_to_string(left->hybrid);
  std::string h2 = hybridization_to_string(right->hybrid);
  if (h1 > h2) std::swap(h1, h2);
  key.hybrid_pair = h1 + "_" + h2;
  key.in_ring = in_ring ? "Y" : "N";

  auto it = bond_hrs_.find(key);
  if (it != bond_hrs_.end()) {
    return it->second;
  }

  return ValueStats();
}

ValueStats AcedrgTables::search_bond_en(const CodAtomInfo& a1,
    const CodAtomInfo& a2) const {

  std::string elem1 = a1.el.name();
  std::string elem2 = a2.el.name();
  std::string sp1 = hybridization_to_string(a1.hybrid);
  std::string sp2 = hybridization_to_string(a2.hybrid);

  // Canonicalize order
  if (elem1 > elem2 || (elem1 == elem2 && sp1 > sp2)) {
    std::swap(elem1, elem2);
    std::swap(sp1, sp2);
  }

  auto it1 = en_bonds_.find(elem1);
  if (it1 == en_bonds_.end()) return ValueStats();

  auto it2 = it1->second.find(sp1);
  if (it2 == it1->second.end()) return ValueStats();

  auto it3 = it2->second.find(elem2);
  if (it3 == it2->second.end()) return ValueStats();

  auto it4 = it3->second.find(sp2);
  if (it4 == it3->second.end()) return ValueStats();

  if (it4->second.empty()) return ValueStats();

  return aggregate_stats(it4->second);
}

ProtHydrDist AcedrgTables::search_prot_hydr_dist(const CodAtomInfo& /* h_atom */,
    const CodAtomInfo& heavy_atom) const {
  // Build type key: H_sp[1|2|3][_arom]_ELEM
  // Examples: H_sp3_C, H_sp2_arom_C, H_sp2_N
  std::string hybr;
  switch (heavy_atom.hybrid) {
    case Hybridization::SP1: hybr = "sp1"; break;
    case Hybridization::SP2: hybr = heavy_atom.is_aromatic ? "sp2_arom" : "sp2"; break;
    case Hybridization::SP3: hybr = "sp3"; break;
    default: return ProtHydrDist();
  }

  std::string type_key = std::string("H_") + hybr + "_" + heavy_atom.el.name();

  auto it = prot_hydr_dists_.find(type_key);
  if (it != prot_hydr_dists_.end()) {
    if (verbose >= 2)
      std::fprintf(stderr, "      prot_hydr_dists: found %s -> electron=%.4f nucleus=%.4f\n",
                   type_key.c_str(), it->second.electron_val, it->second.nucleus_val);
    return it->second;
  }

  // Try without aromatic qualifier if sp2_arom not found
  if (heavy_atom.is_aromatic && heavy_atom.hybrid == Hybridization::SP2) {
    type_key = std::string("H_sp2_") + heavy_atom.el.name();
    it = prot_hydr_dists_.find(type_key);
    if (it != prot_hydr_dists_.end()) {
      if (verbose >= 2)
        std::fprintf(stderr, "      prot_hydr_dists: found %s (fallback from arom) -> electron=%.4f nucleus=%.4f\n",
                     type_key.c_str(), it->second.electron_val, it->second.nucleus_val);
      return it->second;
    }
  }

  return ProtHydrDist();
}

ValueStats AcedrgTables::search_metal_bond(const CodAtomInfo& metal,
    const CodAtomInfo& ligand, const std::vector<CodAtomInfo>& atoms) const {
  int ha1 = metal.connectivity;
  int ha2 = static_cast<int>(ligand.conn_atoms_no_metal.size());

  std::vector<std::string> nb_ids;
  nb_ids.reserve(ligand.conn_atoms_no_metal.size());
  for (int nb : ligand.conn_atoms_no_metal) {
    if (nb >= 0 && nb < static_cast<int>(atoms.size()))
      nb_ids.push_back(atoms[nb].el.name());
  }
  std::sort(nb_ids.begin(), nb_ids.end(), compare_no_case);

  std::string ligand_class;
  for (const auto& id : nb_ids)
    ligand_class.append("_" + id);

  const MetalBondEntry* class_entry = nullptr;
  const MetalBondEntry* pre_entry = nullptr;

  for (const auto& entry : metal_bonds_) {
    if (entry.metal != metal.el || entry.ligand != ligand.el)
      continue;
    if (entry.metal_coord != ha1 || entry.ligand_coord != ha2)
      continue;
    if (!pre_entry)
      pre_entry = &entry;
    if (!ligand_class.empty() && entry.ligand_class == ligand_class)
      class_entry = &entry;
  }

  if (class_entry) {
    if (class_entry->count > metal_class_min_count)
      return ValueStats(class_entry->value, class_entry->sigma,
                        class_entry->count);
    if (pre_entry)
      return ValueStats(pre_entry->pre_value, pre_entry->pre_sigma,
                        pre_entry->pre_count);
  }

  if (pre_entry)
    return ValueStats(pre_entry->pre_value, pre_entry->pre_sigma,
                      pre_entry->pre_count);

  return ValueStats();
}

bool AcedrgTables::lookup_pep_tors(const std::string& a1,
    const std::string& a2, const std::string& a3, const std::string& a4,
    TorsionEntry& out) const {
  auto it = pep_tors_.find(a1 + "_" + a2 + "_" + a3 + "_" + a4);
  if (it == pep_tors_.end())
    return false;
  out = it->second;
  return true;
}

// ============================================================================
// Implementation - Angle search
// ============================================================================

ValueStats AcedrgTables::search_angle_multilevel(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3) const {

  // Build lookup keys - canonicalize flanking atoms
  // Table format: ha1=left_flank, ha2=center, ha3=right_flank
  const CodAtomInfo *flank1, *flank3;
  order_angle_flanks(a1, a3, flank1, flank3);
  int ha1 = flank1->hashing_value;
  int ha2 = center.hashing_value;
  int ha3 = flank3->hashing_value;

  // Build hybridization tuple - table uses hash order: minFlank_center_maxFlank
  std::string h1 = hybridization_to_string(flank1->hybrid);  // min hash
  std::string h2 = hybridization_to_string(center.hybrid);   // center
  std::string h3 = hybridization_to_string(flank3->hybrid);  // max hash
  std::string hybr_tuple = h1 + "_" + h2 + "_" + h3;

  // Build valueKey (ring:hybr_tuple)
  int ring_val = angle_ring_size(center, *flank1, *flank3);
  std::string value_key = std::to_string(ring_val) + ":" + hybr_tuple;

  // Get neighbor symbols - table format: a1=flank1, a2=center, a3=flank3
  std::string a1_nb2 = flank1->nb2_symb;
  std::string a2_nb2 = center.nb2_symb;
  std::string a3_nb2 = flank3->nb2_symb;
  std::string a1_root = flank1->cod_root;
  std::string a2_root = center.cod_root;
  std::string a3_root = flank3->cod_root;
  std::string a1_nb = flank1->nb_symb;
  std::string a2_nb = center.nb_symb;
  std::string a3_nb = flank3->nb_symb;
  std::string a1_type = flank1->cod_main;
  std::string a2_type = center.cod_main;
  std::string a3_type = flank3->cod_main;

  if (verbose >= 2) {
    std::fprintf(stderr,
                 "      angle key %s-%s-%s: hash=%d/%d/%d value=%s "
                 "nb2=%s/%s/%s nb=%s/%s/%s type=%s/%s/%s\n",
                 a1.id.c_str(), center.id.c_str(), a3.id.c_str(),
                 ha1, ha2, ha3, value_key.c_str(),
                 a1_nb2.c_str(), a2_nb2.c_str(), a3_nb2.c_str(),
                 a1_nb.c_str(), a2_nb.c_str(), a3_nb.c_str(),
                 a1_type.c_str(), a2_type.c_str(), a3_type.c_str());
  }

  // Level 1D: Try exact match with atom types
  {
    auto it1 = angle_idx_1d_.find(ha1);
    if (it1 != angle_idx_1d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(a1_nb);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(a2_nb);
                          if (it12 != it11->second.end()) {
                            auto it13 = it12->second.find(a3_nb);
                            if (it13 != it12->second.end()) {
                              auto it14 = it13->second.find(a1_type);
                              if (it14 != it13->second.end()) {
                                auto it15 = it14->second.find(a2_type);
                                if (it15 != it14->second.end()) {
                                  auto it16 = it15->second.find(a3_type);
                                  if (it16 != it15->second.end() && !it16->second.empty()) {
                                    ValueStats vs = it16->second.front();
                                    if (vs.count >= min_observations_angle) {
                                      if (verbose >= 2)
                                        std::fprintf(stderr,
                                                     "      matched angle %s-%s-%s: level=1D\n",
                                                     a1.id.c_str(), center.id.c_str(),
                                                     a3.id.c_str());
                                      return vs;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 2D: No atom types
  {
    auto it1 = angle_idx_2d_.find(ha1);
    if (it1 != angle_idx_2d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end()) {
                        auto it11 = it10->second.find(a1_nb);
                        if (it11 != it10->second.end()) {
                          auto it12 = it11->second.find(a2_nb);
                          if (it12 != it11->second.end()) {
                            auto it13 = it12->second.find(a3_nb);
                            if (it13 != it12->second.end() && !it13->second.empty()) {
                              ValueStats vs = it13->second.front();
                              if (vs.count >= min_observations_angle_fallback) {
                                if (verbose >= 2)
                                  std::fprintf(stderr,
                                               "      matched angle %s-%s-%s: level=2D\n",
                                               a1.id.c_str(), center.id.c_str(),
                                               a3.id.c_str());
                                return vs;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 3D: Hash + valueKey + NB2 only
  {
    auto it1 = angle_idx_3d_.find(ha1);
    if (it1 != angle_idx_3d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end()) {
                  auto it8 = it7->second.find(a1_nb2);
                  if (it8 != it7->second.end()) {
                    auto it9 = it8->second.find(a2_nb2);
                    if (it9 != it8->second.end()) {
                      auto it10 = it9->second.find(a3_nb2);
                      if (it10 != it9->second.end() && !it10->second.empty()) {
                        ValueStats vs = it10->second.front();
                        if (vs.count >= min_observations_angle_fallback) {
                          if (verbose >= 2)
                            std::fprintf(stderr,
                                         "      matched angle %s-%s-%s: level=3D\n",
                                         a1.id.c_str(), center.id.c_str(),
                                         a3.id.c_str());
                          return vs;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 4D: Hash + valueKey only
  {
    auto it1 = angle_idx_4d_.find(ha1);
    if (it1 != angle_idx_4d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end()) {
            auto it5 = it4->second.find(a1_root);
            if (it5 != it4->second.end()) {
              auto it6 = it5->second.find(a2_root);
              if (it6 != it5->second.end()) {
                auto it7 = it6->second.find(a3_root);
                if (it7 != it6->second.end() && !it7->second.empty()) {
                  ValueStats vs = it7->second.front();
                  if (vs.count >= min_observations_angle_fallback) {
                    if (verbose >= 2)
                      std::fprintf(stderr,
                                   "      matched angle %s-%s-%s: level=4D (value=%.3f sigma=%.3f count=%d)\n",
                                   a1.id.c_str(), center.id.c_str(),
                                   a3.id.c_str(), vs.value, vs.sigma, vs.count);
                    return vs;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Level 5D: Hash + valueKey only
  {
    auto it1 = angle_idx_5d_.find(ha1);
    if (it1 != angle_idx_5d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end()) {
          auto it4 = it3->second.find(value_key);
          if (it4 != it3->second.end() && !it4->second.empty()) {
            ValueStats vs = it4->second.front();
            if (vs.count >= min_observations_angle_fallback) {
              if (verbose >= 2)
                std::fprintf(stderr,
                             "      matched angle %s-%s-%s: level=5D\n",
                             a1.id.c_str(), center.id.c_str(),
                             a3.id.c_str());
              return vs;
            }
          }
        }
      }
    }
  }

  // Level 6D: Hash only
  {
    auto it1 = angle_idx_6d_.find(ha1);
    if (it1 != angle_idx_6d_.end()) {
      auto it2 = it1->second.find(ha2);
      if (it2 != it1->second.end()) {
        auto it3 = it2->second.find(ha3);
        if (it3 != it2->second.end() && !it3->second.empty()) {
          if (verbose >= 2)
            std::fprintf(stderr, "      matched angle %s-%s-%s: level=6D\n",
                         a1.id.c_str(), center.id.c_str(), a3.id.c_str());
          return it3->second.front();
        }
      }
    }
  }

  // No match found in detailed tables
  return ValueStats();
}

void AcedrgTables::fill_angle(const ChemComp& cc,
    const std::vector<CodAtomInfo>& atom_info,
    Restraints::Angle& angle) const {

  int idx1 = cc.find_atom_index(angle.id1.atom);
  int idx2 = cc.find_atom_index(angle.id2.atom);  // center
  int idx3 = cc.find_atom_index(angle.id3.atom);
  if (idx1 < 0 || idx2 < 0 || idx3 < 0)
    return;

  const CodAtomInfo& a1 = atom_info[idx1];
  const CodAtomInfo& center = atom_info[idx2];
  const CodAtomInfo& a3 = atom_info[idx3];

  // Check for metal center
  if (center.is_metal) {
    std::vector<double> ideal_angles = get_metal_angles(center.el,
                                                        center.connectivity);
    if (!ideal_angles.empty()) {
      // Use the most common angle for this geometry
      angle.value = ideal_angles[0];
      angle.esd = clamp_angle_sigma(3.0);
      return;
    }
  }

  int ring_size = angle_ring_size(center, a1, a3);

  // Try detailed multilevel search first
  ValueStats vs = search_angle_multilevel(a1, center, a3);
  if (vs.count >= min_observations_angle) {
    angle.value = vs.value;
    angle.esd = clamp_angle_sigma(vs.sigma);
    return;
  }

  // Try HRS table
  vs = search_angle_hrs(a1, center, a3, ring_size);
  if (vs.count >= min_observations_angle) {
    angle.value = vs.value;
    angle.esd = clamp_angle_sigma(vs.sigma);
    return;
  }

  // Fallback based on hybridization
  switch (center.hybrid) {
    case Hybridization::SP1:
      angle.value = 180.0;
      break;
    case Hybridization::SP2:
      angle.value = 120.0;
      break;
    case Hybridization::SP3:
      angle.value = 109.47;  // tetrahedral
      break;
    case Hybridization::SPD5:
    case Hybridization::SPD6:
      angle.value = 90.0;    // octahedral
      break;
    default:
      angle.value = 109.47;  // default tetrahedral
      break;
  }
  angle.esd = upper_angle_sigma;
}

ValueStats AcedrgTables::search_angle_hrs(const CodAtomInfo& a1,
    const CodAtomInfo& center, const CodAtomInfo& a3, int ring_size) const {

  AngleHRSKey key;
  // Loaded data has: hash1=min(flank), hash2=center, hash3=max(flank)
  key.hash1 = std::min(a1.hashing_value, a3.hashing_value);
  key.hash2 = center.hashing_value;
  key.hash3 = std::max(a1.hashing_value, a3.hashing_value);

  // Build hybrid tuple - table uses hash order: first_center_third
  // where first/third are determined by min/max of flank hash values
  std::string h1, h2, h3;
  h2 = hybridization_to_string(center.hybrid);
  if (a1.hashing_value <= a3.hashing_value) {
    h1 = hybridization_to_string(a1.hybrid);
    h3 = hybridization_to_string(a3.hybrid);
  } else {
    h1 = hybridization_to_string(a3.hybrid);
    h3 = hybridization_to_string(a1.hybrid);
  }
  std::string hybr_tuple = h1 + "_" + h2 + "_" + h3;
  key.value_key = std::to_string(ring_size) + ":" + hybr_tuple;

  auto it = angle_hrs_.find(key);
  if (it != angle_hrs_.end()) {
    return it->second;
  }

  return ValueStats();
}

std::vector<double> AcedrgTables::get_metal_angles(Element metal,
    int coord_number) const {

  std::vector<double> angles;

  // First check if we have specific angles in the table
  for (const auto& entry : metal_angles_) {
    if (entry.metal == metal && entry.coord_number == coord_number) {
      angles.push_back(entry.angle);
    }
  }

  if (!angles.empty()) {
    return angles;
  }

  // Fall back to geometry defaults
  auto geo_it = metal_coord_geo_.find(metal);
  if (geo_it != metal_coord_geo_.end()) {
    auto cn_it = geo_it->second.find(coord_number);
    if (cn_it != geo_it->second.end()) {
      switch (cn_it->second) {
        case CoordGeometry::LINEAR:
          angles.push_back(180.0);
          break;
        case CoordGeometry::TRIGONAL_PLANAR:
          angles.push_back(120.0);
          break;
        case CoordGeometry::T_SHAPED:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::TETRAHEDRAL:
          angles.push_back(109.47);
          break;
        case CoordGeometry::SQUARE_PLANAR:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::TRIGONAL_BIPYRAMIDAL:
          angles.push_back(90.0);
          angles.push_back(120.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::SQUARE_PYRAMIDAL:
          angles.push_back(90.0);
          angles.push_back(180.0);
          break;
        case CoordGeometry::OCTAHEDRAL:
          angles.push_back(90.0);
          break;
        default:
          break;
      }
    }
  }

  // Ultimate fallback based on coordination number
  if (angles.empty()) {
    switch (coord_number) {
      case 2: angles.push_back(180.0); break;
      case 3: angles.push_back(120.0); break;
      case 4: angles.push_back(109.47); break;
      case 5: angles.push_back(90.0); angles.push_back(120.0); break;
      case 6: angles.push_back(90.0); break;
      default: angles.push_back(90.0); break;
    }
  }

  return angles;
}

// ============================================================================
// Implementation - Statistical helpers
// ============================================================================

ValueStats AcedrgTables::aggregate_stats(
    const std::vector<ValueStats>& values) const {

  if (values.empty())
    return ValueStats();

  if (values.size() == 1)
    return values[0];

  // Match AceDRG setValueSet() aggregation (weighted mean + pooled sigma).
  double sum_val = 0.0;
  double sum1 = 0.0;
  int total_count = 0;

  for (const auto& v : values) {
    if (v.count > 0) {
      sum_val += v.value * v.count;
      sum1 += (v.count - 1) * v.sigma * v.sigma + v.count * v.value * v.value;
      total_count += v.count;
    }
  }

  if (total_count == 0)
    return ValueStats();

  double mean = sum_val / total_count;
  double sum2 = mean * sum_val;

  double sigma = 0.0;
  if (total_count > 1)
    sigma = std::sqrt(std::fabs(sum1 - sum2) / (total_count - 1));
  else
    sigma = std::sqrt(std::fabs(sum1 - sum2) / total_count);

  return ValueStats(mean, sigma, total_count);
}

int AcedrgTables::str_to_int(const std::string& s) {
  std::istringstream iss(s);
  int value = 0;
  iss >> value;
  return value;
}

bool AcedrgTables::compare_no_case(const std::string& first,
                                   const std::string& second) {
  size_t i = 0;
  while (i < first.length() && i < second.length()) {
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(first[i])));
    char b = static_cast<char>(std::toupper(static_cast<unsigned char>(second[i])));
    if (a < b)
      return true;
    if (a > b)
      return false;
    ++i;
  }
  return first.length() > second.length();
}

bool AcedrgTables::compare_no_case2(const std::string& first,
                                    const std::string& second) {
  if (first.length() > second.length())
    return true;
  if (first.length() < second.length())
    return false;
  for (size_t i = 0; i < first.length() && i < second.length(); ++i) {
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(first[i])));
    char b = static_cast<char>(std::toupper(static_cast<unsigned char>(second[i])));
    if (a < b)
      return true;
    if (a > b)
      return false;
  }
  return false;
}

bool AcedrgTables::desc_sort_map_key(const SortMap& a, const SortMap& b) {
  return a.key.length() > b.key.length();
}

bool AcedrgTables::desc_sort_map_key2(const SortMap2& a, const SortMap2& b) {
  if (a.key.length() > b.key.length())
    return true;
  if (a.key.length() == b.key.length())
    return a.nNB > b.nNB;
  return false;
}

// Helper to compact element list: ["H","H","H"] -> "H3"
// Only compacts runs of 3+ identical elements (acedrg convention)
static std::string compact_element_list(const std::vector<std::string>& elems) {
  std::string result;
  for (size_t i = 0; i < elems.size(); ) {
    size_t count = 1;
    while (i + count < elems.size() && elems[i + count] == elems[i])
      ++count;
    if (count >= 3) {
      // Compact runs of 3+
      result += elems[i];
      result += std::to_string(count);
    } else {
      // Don't compact runs of 1-2, just repeat the element
      for (size_t j = 0; j < count; ++j)
        result += elems[i];
    }
    i += count;
  }
  return result;
}

std::string AcedrgTables::compute_acedrg_type(
    const CodAtomInfo& atom,
    const std::vector<CodAtomInfo>& atoms,
    const std::vector<std::vector<int>>& neighbors) const {
  std::string result = element_name(atom.el);

  // Build neighbor descriptions: neighbor_element + neighbor's_other_neighbors
  // Note: AceDRG excludes metals from the atom typing system entirely
  std::map<std::string, int> neighbor_groups;
  for (int nb_idx : neighbors[atom.index]) {
    if (atoms[nb_idx].is_metal)  // Skip metals in neighbor description
      continue;
    std::string desc = element_name(atoms[nb_idx].el);
    // Collect second neighbors (excluding the central atom and metals)
    std::vector<std::string> second_nbs;
    for (int nb2_idx : neighbors[nb_idx]) {
      if (nb2_idx != atom.index && !atoms[nb2_idx].is_metal)
        second_nbs.push_back(element_name(atoms[nb2_idx].el));
    }
    // Sort alphabetically
    std::sort(second_nbs.begin(), second_nbs.end());
    // Compact runs: H,H,H -> H3
    desc += compact_element_list(second_nbs);
    neighbor_groups[desc]++;
  }

  // Sort groups: by length desc, then alphabetically
  std::vector<std::pair<std::string, int>> sorted_groups(
      neighbor_groups.begin(), neighbor_groups.end());
  std::sort(sorted_groups.begin(), sorted_groups.end(),
    [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
      if (a.first.size() != b.first.size())
        return a.first.size() > b.first.size();  // longer first
      return a.first < b.first;  // then alphabetically
    });

  for (const auto& group : sorted_groups) {
    result += "(" + group.first + ")";
    if (group.second > 1)
      result += std::to_string(group.second);
  }
  return result;
}

std::vector<std::string> AcedrgTables::compute_acedrg_types(const ChemComp& cc) const {
  std::vector<CodAtomInfo> atom_info = classify_atoms(cc);

  std::vector<std::string> types;
  types.reserve(cc.atoms.size());
  for (size_t i = 0; i < cc.atoms.size(); ++i) {
    types.push_back(atom_info[i].cod_class);
  }
  return types;
}

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

  tables.apply_metal_charge_corrections(cc);

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
