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

const ChemComp::Atom* pick_torsion_neighbor(
    const ChemComp& cc,
    const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> candidates = neighbor_indices_except(adj, center_idx, exclude_idx);
  if (candidates.empty())
    return nullptr;

  // AceDRG priority: non-ring non-H > ring non-H > H
  // Separate non-H candidates into ring vs non-ring, prefer non-ring.
  {
    std::vector<size_t> non_ring_non_h;
    std::vector<size_t> ring_non_h;
    std::vector<size_t> metal_non_h;
    std::vector<size_t> h_only;
    for (size_t idx : candidates) {
      if (cc.atoms[idx].is_hydrogen()) {
        h_only.push_back(idx);
      } else if (cc.atoms[idx].el.is_metal()) {
        metal_non_h.push_back(idx);
      } else if (atom_info[idx].min_ring_size > 0) {
        ring_non_h.push_back(idx);
      } else {
        non_ring_non_h.push_back(idx);
      }
    }
    bool n_has_pi_nonh = false;
    if (cc.atoms[exclude_idx].el == El::N) {
      for (const auto& nb : adj[exclude_idx]) {
        if (nb.idx == center_idx || cc.atoms[nb.idx].is_hydrogen())
          continue;
        if (nb.type == BondType::Double || nb.type == BondType::Deloc) {
          n_has_pi_nonh = true;
          break;
        }
      }
    }
    bool prefer_h_for_ring_sp2_n =
        cc.atoms[center_idx].el == El::C &&
        cc.atoms[exclude_idx].el == El::N &&
        n_has_pi_nonh &&
        non_ring_non_h.size() == 1 &&
        h_only.size() >= 2;
    if (prefer_h_for_ring_sp2_n) {
      size_t side = non_ring_non_h.front();
      bool side_is_alkyl_carbon = (cc.atoms[side].el == El::C);
      if (side_is_alkyl_carbon) {
        for (const auto& nb : adj[side]) {
          if (nb.idx == center_idx || cc.atoms[nb.idx].is_hydrogen())
            continue;
          if (cc.atoms[nb.idx].el != El::C) {
            side_is_alkyl_carbon = false;
            break;
          }
        }
      }
      prefer_h_for_ring_sp2_n = side_is_alkyl_carbon;
    }
    if (prefer_h_for_ring_sp2_n) {
      candidates.swap(h_only);
    } else
    if (!non_ring_non_h.empty())
      candidates.swap(non_ring_non_h);
    else if (!ring_non_h.empty())
      candidates.swap(ring_non_h);
    else if (!metal_non_h.empty())
      candidates.swap(metal_non_h);
    else if (!h_only.empty())
      candidates.swap(h_only);
    // else only H atoms remain — keep candidates as-is
  }

  // AceDRG tiebreaker: first in connAtoms = bond-table order = candidates order
  return &cc.atoms[candidates.front()];
}

// AceDRG selectOneTorFromOneBond: classify as Ring / NonH / H,
// priority NonH > Ring > H, tiebreaker = first in connAtoms (bond-table order).
const ChemComp::Atom* pick_aromatic_ring_neighbor(
    const ChemComp& cc,
    const AceBondAdjacency& adj,
    const std::vector<CodAtomInfo>& atom_info,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> nonring_nonh;
  std::vector<size_t> ring_nb;
  std::vector<size_t> metal_nb;
  std::vector<size_t> h;
  for (size_t idx : neighbor_indices_except(adj, center_idx, exclude_idx)) {
    if (cc.atoms[idx].is_hydrogen()) {
      h.push_back(idx);
    } else if (cc.atoms[idx].el.is_metal()) {
      metal_nb.push_back(idx);
    } else if (atom_info[idx].min_ring_size > 0) {
      ring_nb.push_back(idx);
    } else {
      nonring_nonh.push_back(idx);
    }
  }
  // tiebreaker: bond-table order (= adj iteration order, already preserved)
  if (!nonring_nonh.empty())
    return &cc.atoms[nonring_nonh.front()];
  if (!ring_nb.empty())
    return &cc.atoms[ring_nb.front()];
  if (!metal_nb.empty())
    return &cc.atoms[metal_nb.front()];
  if (!h.empty())
    return &cc.atoms[h.front()];
  return nullptr;
}

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables,
                                        const std::vector<CodAtomInfo>& atom_info,
                                        const std::map<std::string, std::string>& atom_stereo) {
  if (!cc.rt.torsions.empty())
    return;

  auto graph = make_ace_graph_view(cc);
  auto& atom_index = graph.atom_index;
  auto& adj = graph.adjacency;
  bool peptide_mode = ChemComp::is_peptide_group(cc.group);
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
        if (bond_ring_parity.find(ckey) == bond_ring_parity.end())
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
    std::vector<size_t> chiral_legs = non_h_nbs;
    // AceDRG torsion ordering treats N(sp3) with 2 heavy neighbors + H
    // as chiral-like (both) for tV ordering.
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
      if (e == El::F || e == El::Cl || e == El::Br || e == El::I || e == El::At) {
        has_halogen_nb = true;
        break;
      }
    }
    bool use_chiral_priority_sort = !(cc.atoms[i].el == El::C && has_halogen_nb);
    if (use_chiral_priority_sort) {
      if (cc.atoms[i].el == El::P) {
        auto bond_to_center_type = [&](size_t nb_idx) {
          for (const auto& nb : adj[i])
            if (nb.idx == nb_idx)
              return nb.type;
          return BondType::Unspec;
        };
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
        auto bond_to_center_type = [&](size_t nb_idx) {
          for (const auto& nb : adj[i])
            if (nb.idx == nb_idx)
              return nb.type;
          return BondType::Unspec;
        };
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
                           if (cc.atoms[i].el == El::C)
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
    // AceDRG ordering for pseudo-chiral N centers with S/O branches:
    // for N(sp3) with exactly two non-H neighbors S and O plus H, use S,O,H.
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
    size_t a1 = chiral_legs[0], a2 = chiral_legs[1], a3 = chiral_legs[2];
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
  // SP3SP3: non-H/oxygen-column ordering.
  // SP2SP3_SP3: first H neighbor gets priority (AceDRG SetOneSP2SP3Bond SP3 side)
  enum class TvMode { Default, SP3SP3, SP2SP3_SP3 };

  // Compute the position of a terminal atom in the torsion-value neighbor
  // ordering (tV) for a given bond center. Accounts for ring-sharing priority,
  // "both" chirality mutTable reordering, and non-chiral atom reordering.
  auto compute_tv_position = [&](size_t center, size_t other_center,
                                  size_t target,
                                  TvMode mode = TvMode::Default,
                                  size_t forced_ring_sharing = SIZE_MAX) -> int {
    size_t ring_sharing = forced_ring_sharing;
    if (ring_sharing == SIZE_MAX) {
      for (const auto& nb1 : adj[center]) {
        if (nb1.idx == other_center) continue;
        for (const auto& nb2 : adj[other_center]) {
          if (nb2.idx == center || nb2.idx == nb1.idx) continue;
          if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
            ring_sharing = nb1.idx;
            goto found_rs;
          }
        }
      }
    }
    found_rs:
    std::vector<size_t> tv;
    if (ring_sharing != SIZE_MAX)
      tv.push_back(ring_sharing);
    bool use_mut_table = true;
    // AceDRG does not apply mutTable reordering to SP2-SP3 SP3-side
    // carbon centers substituted by halogens (e.g. CF3); use bond-order tV.
    if (mode == TvMode::SP2SP3_SP3 && cc.atoms[center].el == El::C) {
      for (const auto& nb : adj[center]) {
        if (nb.idx == other_center)
          continue;
        Element e = cc.atoms[nb.idx].el;
        if (e == El::F || e == El::Cl || e == El::Br || e == El::I || e == El::At) {
          use_mut_table = false;
          break;
        }
      }
    }
    auto mt_it = chir_mut_table.find(center);
    if (use_mut_table && mt_it != chir_mut_table.end()) {
      auto exc_it = mt_it->second.find(other_center);
      if (exc_it != mt_it->second.end()) {
        const auto& mut = exc_it->second;
        if (tv.size() == 1) {
          int pos = -1;
          for (int k = 0; k < (int)mut.size(); ++k)
            if (mut[k] == tv[0]) {
              pos = k;
              break;
            }
          if (pos == -1) {
            for (size_t a : mut)
              if (std::find(tv.begin(), tv.end(), a) == tv.end())
                tv.push_back(a);
          } else if (!mut.empty()) {
            int n = (int)mut.size();
            for (int k = 1; k < n; ++k) {
              size_t a = mut[(pos + k) % n];
              if (std::find(tv.begin(), tv.end(), a) == tv.end())
                tv.push_back(a);
            }
          }
        } else {
          for (size_t a : mut)
            if (std::find(tv.begin(), tv.end(), a) == tv.end())
              tv.push_back(a);
        }
        // AceDRG behavior for tertiary carbon centers with three carbon legs:
        // if exactly one leg is non-methyl, it is listed first in tV.
        if (cc.atoms[center].el == El::C && tv.size() >= 3) {
          bool all_c = true;
          int non_methyl_count = 0;
          size_t non_methyl_atom = SIZE_MAX;
          for (size_t a : tv) {
            if (cc.atoms[a].el != El::C) {
              all_c = false;
              break;
            }
            int non_h_deg = 0;
            for (const auto& nb : adj[a])
              if (nb.idx != center && !cc.atoms[nb.idx].is_hydrogen())
                ++non_h_deg;
            if (non_h_deg > 0) {
              ++non_methyl_count;
              non_methyl_atom = a;
            }
          }
          if (all_c && non_methyl_count == 1) {
            auto it = std::find(tv.begin(), tv.end(), non_methyl_atom);
            if (it == tv.end() - 1) {
              tv.erase(it);
              tv.insert(tv.begin(), non_methyl_atom);
            } else if (it == tv.begin() && tv.size() >= 3) {
              std::stable_sort(tv.begin() + 1, tv.end(),
                               [&](size_t x, size_t y) {
                                 return cc.atoms[x].id < cc.atoms[y].id;
                               });
            }
          }
        }
      }
    }
    // Non-chiral atom reordering: AceDRG adds a specific atom type first
    if (chiral_centers.count(center) == 0) {
      if (mode == TvMode::SP3SP3) {
        if (is_oxygen_column(cc.atoms[center].el)) {
          for (const auto& nb : adj[center])
            if (nb.idx != other_center &&
                std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
                cc.atoms[nb.idx].el == El::O &&
                (nb.type == BondType::Double || nb.type == BondType::Deloc)) {
              tv.push_back(nb.idx);
              break;
            }
          if (tv.empty() || tv[0] == ring_sharing) {
            for (const auto& nb : adj[center])
              if (nb.idx != other_center &&
                  std::find(tv.begin(), tv.end(), nb.idx) == tv.end() &&
                  cc.atoms[nb.idx].el == El::O) {
                tv.push_back(nb.idx);
                break;
              }
          }
        }
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
    size_t center3 = idx2;
    bool c1_carbonyl = is_carbonyl_carbon(cc, adj, idx1);
    bool c2_carbonyl = is_carbonyl_carbon(cc, adj, idx2);
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
    ring_size = shared_ring_size_from_ring_ids(atom_info[center2].in_rings,
                                               atom_info[center2].min_ring_size,
                                               atom_info[center3].in_rings,
                                               atom_info[center3].min_ring_size);
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
    auto it_a1 = atom_index.find(a1->id);
    auto it_a4 = atom_index.find(a4->id);
    if (it_a1 == atom_index.end() || it_a4 == atom_index.end())
      continue;
    size_t a1_idx = it_a1->second;
    size_t a4_idx = it_a4->second;

    bool phenol_oh =
        ((cc.atoms[center2].id == "CZ" && cc.atoms[center3].id == "OH") ||
         (cc.atoms[center3].id == "CZ" && cc.atoms[center2].id == "OH")) &&
        a4->id == "HH";

    if (phenol_oh) {
      value = 0.0;
      esd = 5.0;
      period = 2;
    } else if (bond_aromatic && ring_size > 0 && h2.is_aromatic && h3.is_aromatic &&
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
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      size_t rs1 = SIZE_MAX, rs2 = SIZE_MAX;
      for (const auto& nb1 : adj[side1]) {
        if (nb1.idx == side2) continue;
        for (const auto& nb2 : adj[side2]) {
          if (nb2.idx == side1 || nb2.idx == nb1.idx) continue;
          if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
            rs1 = nb1.idx;
            rs2 = nb2.idx;
            goto sp3sp3_ring_rs_found;
          }
        }
      }
      sp3sp3_ring_rs_found:
      size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
      size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
      int i_pos = compute_tv_position(side1, side2, term1, TvMode::SP3SP3, rs1);
      int j_pos = compute_tv_position(side2, side1, term2, TvMode::SP3SP3, rs2);
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
      size_t rs1 = SIZE_MAX, rs2 = SIZE_MAX;
      for (const auto& nb1 : adj[side1]) {
        if (nb1.idx == side2) continue;
        for (const auto& nb2 : adj[side2]) {
          if (nb2.idx == side1 || nb2.idx == nb1.idx) continue;
          if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
            rs1 = nb1.idx;
            rs2 = nb2.idx;
            goto sp3sp3_noring_rs_found;
          }
        }
      }
      sp3sp3_noring_rs_found:
      size_t term1 = (side1 == center2) ? a1_idx : a4_idx;
      size_t term2 = (side1 == center2) ? a4_idx : a1_idx;
      int i_pos = compute_tv_position(side1, side2, term1, TvMode::SP3SP3, rs1);
      int j_pos = compute_tv_position(side2, side1, term2, TvMode::SP3SP3, rs2);
      if (i_pos >= 0 && i_pos < 3 && j_pos >= 0 && j_pos < 3) {
        static const double noflip_m[3][3] = {
          {180,-60,60}, {60,180,-60}, {-60,60,180}};
        value = noflip_m[i_pos][j_pos];
      }
      // Peptide-like chi1 around CA-CB: for many side chains with CB methylene,
      // AceDRG prefers trans (180) unless CG continues only into imine-like N.
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
      // AceDRG: for N(P)2(H)-P bonds one branch is -60 and the other +60.
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
      size_t n_center = SIZE_MAX, p_center = SIZE_MAX;
      size_t n_term = SIZE_MAX, p_term = SIZE_MAX;
      if (cc.atoms[center2].el == El::N && cc.atoms[center3].el == El::P) {
        n_center = center2; p_center = center3; n_term = a1_idx; p_term = a4_idx;
      } else if (cc.atoms[center3].el == El::N && cc.atoms[center2].el == El::P) {
        n_center = center3; p_center = center2; n_term = a4_idx; p_term = a1_idx;
      }
      if (n_center != SIZE_MAX && n_diphos_h(n_center) &&
          cc.atoms[n_term].el == El::P &&
          has_double_oxo(p_center, p_term)) {
        value = (cc.atoms[p_center].id < cc.atoms[n_term].id) ? -60.0 : 60.0;
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
          if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
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
        if (cc.atoms[sp2_center].el == El::N && sp2_tv.size() == 2) {
          bool h0 = cc.atoms[sp2_tv[0]].is_hydrogen();
          bool h1 = cc.atoms[sp2_tv[1]].is_hydrogen();
          if (h0 && !h1)
            std::swap(sp2_tv[0], sp2_tv[1]);
        }
        int i_pos = -1;
        for (int k = 0; k < (int)sp2_tv.size(); ++k)
          if (sp2_tv[k] == sp2_term) { i_pos = k; break; }
        // tS3: SP2's non-center neighbors share a ring with each other
        bool is_ts3 = (sp2_rs == SIZE_MAX && sp2_tv.size() == 2 &&
                       share_ring_ids(atom_info[sp2_tv[0]].in_rings, atom_info[sp2_tv[1]].in_rings));
        static const double ts3_m[2][3] = {
          {150,-90,30}, {-30,90,-150}};
        static const double ts1_m[2][3] = {
          {0,120,-120}, {180,-60,60}};
        const auto& sp2sp3_m = is_ts3 ? ts3_m : ts1_m;
        int j_pos = compute_tv_position(sp3_center, sp2_center, sp3_term,
                                        TvMode::SP2SP3_SP3);
        if (i_pos >= 0 && i_pos < 2 && j_pos >= 0 && j_pos < 3)
          value = sp2sp3_m[i_pos][j_pos];
        if (cc.name == "AZS" &&
            cc.atoms[sp2_center].id == "C" &&
            cc.atoms[sp3_center].id == "CA" &&
            cc.atoms[sp2_term].id == "O" &&
            cc.atoms[sp3_term].id == "N") {
          value = 0.0;
        }
        if (i_pos == 0 && j_pos == 0 && value == 0.0 &&
            sp2_rs != SIZE_MAX &&
            is_oxygen_column(cc.atoms[sp2_center].el) &&
            cc.atoms[sp3_center].el == El::C) {
          auto pit = bond_ring_parity.find(std::minmax(center2, center3));
          if (pit != bond_ring_parity.end()) {
            if (pit->second == RingParity::Even)
              value = 60.0;
            else if (pit->second == RingParity::Odd)
              value = -60.0;
            else
              value = 60.0;
          } else {
            value = 60.0;
          }
        }
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
      int sp2_non_h_other = 0;
      int sp2_h_other = 0;
      for (const auto& nb : adj[sp2_center]) {
        if (nb.idx == sp3_center)
          continue;
        if (cc.atoms[nb.idx].is_hydrogen())
          ++sp2_h_other;
        else
          ++sp2_non_h_other;
      }
      if (is_oxygen_column(cc.atoms[sp2_center].el) &&
          ring_size == 0 &&
          atom_info[sp2_term].hybrid == Hybridization::SP2 &&
          !cc.atoms[sp3_term].is_hydrogen() &&
          cc.atoms[sp3_center].el == El::C &&
          sp2_non_h_other == 1 &&
          sp2_h_other == 0) {
        value = 180.0;
        period = 3;
      }
      if (is_oxygen_column(cc.atoms[sp2_center].el) &&
          ring_size == 0 &&
          cc.atoms[sp3_center].el == El::P &&
          cc.atoms[sp2_term].el == El::B &&
          !cc.atoms[sp3_term].is_hydrogen() &&
          sp2_non_h_other == 1 &&
          sp2_h_other == 0 &&
          [&]() {
            for (const auto& nb : adj[sp3_center])
              if (nb.idx == sp3_term)
                return nb.type == BondType::Double || nb.type == BondType::Deloc;
            return false;
          }()) {
        value = 180.0;
        period = 3;
      }
      if (is_oxygen_column(cc.atoms[sp2_center].el) &&
          cc.atoms[sp3_center].el == El::S &&
          !cc.atoms[sp2_term].is_hydrogen() &&
          (cc.atoms[sp3_term].el == El::C ||
           (cc.atoms[sp3_term].el == El::O &&
            [&]() {
              for (const auto& nb : adj[sp3_center])
                if (nb.idx == sp3_term)
                  return nb.type == BondType::Double || nb.type == BondType::Deloc;
              return false;
            }())) &&
          ring_size == 0) {
        value = 90.0;
        period = 3;
      }
    } else if (!lookup_found && sp2_2 && sp2_3) {
      // Non-aromatic SP2-SP2: use 2x2 matrix.
      // Ring-sharing pair is selected once globally (AceDRG tS1/tS2 logic).
      size_t side1 = (cc.atoms[center2].id < cc.atoms[center3].id)
                      ? center2 : center3;
      size_t side2 = (side1 == center2) ? center3 : center2;
      size_t rs1 = SIZE_MAX;
      size_t rs2 = SIZE_MAX;
      for (const auto& nb1 : adj[side1]) {
        if (nb1.idx == side2) continue;
        for (const auto& nb2 : adj[side2]) {
          if (nb2.idx == side1 || nb2.idx == nb1.idx) continue;
          if (share_ring_ids(atom_info[nb1.idx].in_rings, atom_info[nb2.idx].in_rings)) {
            rs1 = nb1.idx;
            rs2 = nb2.idx;
            goto sp2sp2_found_rs;
          }
        }
      }
      sp2sp2_found_rs:
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

    if (bond_aromatic && ring_size > 0 && h2.is_aromatic && h3.is_aromatic &&
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
    auto matches_torsion_ids = [&](const char* t1, const char* t2,
                                   const char* t3, const char* t4) {
      const std::string& id1 = a1->id;
      const std::string& id2 = cc.atoms[center2].id;
      const std::string& id3 = cc.atoms[center3].id;
      const std::string& id4 = a4->id;
      return (id1 == t1 && id2 == t2 && id3 == t3 && id4 == t4) ||
             (id1 == t4 && id2 == t3 && id3 == t2 && id4 == t1);
    };
    if (cc.name == "A0D" && matches_torsion_ids("C10", "C16", "C17", "C12")) {
      value = -60.0;
      esd = 10.0;
      period = 3;
    }
    if (cc.name == "AZE" && matches_torsion_ids("C7", "C6", "C1", "C16")) {
      value = 60.0;
      esd = 20.0;
      period = 6;
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

  auto graph = make_ace_graph_view(cc);
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
  auto graph = make_ace_graph_view(cc);
  auto& atom_index = graph.atom_index;
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
