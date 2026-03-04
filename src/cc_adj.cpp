// Chemical adjustment phase split from ace_cc.cpp

#include "gemmi/cc_adj.hpp"
#include "gemmi/ace_graph.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <unordered_set>

namespace gemmi {
namespace {

bool ace_trace_mode() {
  if (const char* env = std::getenv("GEMMI_ACE_TRACE"))
    return env[0] != '\0' && env[0] != '0';
  return false;
}

struct AceRuleStats {
  int atom_count = 0;
  int bond_count = 0;
  int angle_count = 0;
  int torsion_count = 0;
  int chir_count = 0;
  int plane_count = 0;
  double charge_sum = 0.0;
};

AceRuleStats collect_rule_stats(const ChemComp& cc) {
  AceRuleStats s;
  s.atom_count = static_cast<int>(cc.atoms.size());
  s.bond_count = static_cast<int>(cc.rt.bonds.size());
  s.angle_count = static_cast<int>(cc.rt.angles.size());
  s.torsion_count = static_cast<int>(cc.rt.torsions.size());
  s.chir_count = static_cast<int>(cc.rt.chirs.size());
  s.plane_count = static_cast<int>(cc.rt.planes.size());
  for (const auto& atom : cc.atoms)
    s.charge_sum += atom.charge;
  return s;
}

bool rule_stats_changed(const AceRuleStats& before, const AceRuleStats& after) {
  if (before.atom_count != after.atom_count) return true;
  if (before.bond_count != after.bond_count) return true;
  if (before.angle_count != after.angle_count) return true;
  if (before.torsion_count != after.torsion_count) return true;
  if (before.chir_count != after.chir_count) return true;
  if (before.plane_count != after.plane_count) return true;
  return std::fabs(before.charge_sum - after.charge_sum) > 1e-6;
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
      if (bond.type == BondType::Double || is_aromatic_or_deloc(bond.type))
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
    bool has_ester_o = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom != c_neighbor && bond.id2.atom != c_neighbor)
        continue;
      std::string other = (bond.id1.atom == c_neighbor) ? bond.id2.atom : bond.id1.atom;
      int other_idx = cc.find_atom_index(other);
      if (other_idx >= 0 && cc.atoms[other_idx].el == El::O) {
        ++o_count;
        if (bond.type == BondType::Double || bond.type == BondType::Deloc) {
          has_double_o = true;
        } else {
          auto o_it = neighbors.find(other);
          if (o_it != neighbors.end()) {
            for (const std::string& o_nb : o_it->second) {
              if (o_nb == c_neighbor)
                continue;
              int o_nb_idx = cc.find_atom_index(o_nb);
              if (o_nb_idx >= 0 && cc.atoms[o_nb_idx].el != El::H) {
                has_ester_o = true;
                break;
              }
            }
          }
        }
      }
    }
    if (o_count >= 2 && has_double_o && !has_ester_o)
      return true;
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
  if (o_count < 2 || !has_double_o || has_n_neighbor)
    return;

  oxt_it->charge = -1.0f;
  remove_atom_by_id(cc, "HXT");
}

} // namespace

bool add_n_terminal_h3(ChemComp& cc) {
  if (cc.find_atom("H3") != cc.atoms.end())
    return false;
  if (cc.find_atom("H") == cc.atoms.end() || cc.find_atom("H2") == cc.atoms.end())
    return false;
  auto n_it = cc.find_atom("N");
  if (n_it == cc.atoms.end())
    return false;
  if (std::fabs(n_it->charge - 1.0f) > 0.6f)
    return false;
  int n_neighbors = 0;
  int n_h = 0;
  std::string ca_atom;
  for (const auto& bond : cc.rt.bonds) {
    std::string other;
    if (bond.id1.atom == "N")
      other = bond.id2.atom;
    else if (bond.id2.atom == "N")
      other = bond.id1.atom;
    else
      continue;
    ++n_neighbors;
    if (other == "H" || other == "H2")
      ++n_h;
    else if (other == "CA")
      ca_atom = other;
    else {
      auto it = cc.find_atom(other);
      if (it != cc.atoms.end() && it->el != El::H)
        return false;
    }
  }
  if (n_neighbors != 3 || n_h != 2 || ca_atom.empty())
    return false;
  if (atom_has_unsaturated_bond(cc, ca_atom))
    return false;
  auto neighbors = make_neighbor_names(cc);
  if (!has_carboxylic_acid_neighbor(cc, neighbors, ca_atom))
    return false;

  cc.atoms.push_back(ChemComp::Atom{"H3", "", El::H, 0.0f, "H", "", Position()});
  cc.rt.bonds.push_back({{1, "N"}, {1, "H3"}, BondType::Single, false,
                         NAN, NAN, NAN, NAN});

  auto add_angle_if_present = [&](const std::string& a1,
                                  const std::string& center,
                                  const std::string& a3) {
    if (cc.find_atom(a1) == cc.atoms.end() || cc.find_atom(center) == cc.atoms.end() ||
        cc.find_atom(a3) == cc.atoms.end())
      return;
    cc.rt.angles.push_back({{1, a1}, {1, center}, {1, a3}, NAN, NAN});
  };
  add_angle_if_present("CA", "N", "H3");
  add_angle_if_present("H", "N", "H3");
  add_angle_if_present("H2", "N", "H3");
  return true;
}

namespace {

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

} // namespace

void sync_n_terminal_h3_angles(ChemComp& cc) {
  if (cc.find_atom("H3") == cc.atoms.end() || cc.find_atom("N") == cc.atoms.end())
    return;
  auto copy_angle = [&](const std::string& a1, const std::string& a3,
                        const std::string& src1, const std::string& src3) {
    Restraints::Angle* target = find_angle(cc, "N", a1, a3);
    if (!target)
      return;
    Restraints::Angle* source = find_angle(cc, "N", src1, src3);
    if (!source)
      return;
    target->value = source->value;
    target->esd = source->esd;
  };
  {
    Restraints::Angle* ca_h3 = find_angle(cc, "N", "CA", "H3");
    if (ca_h3 && std::isnan(ca_h3->value)) {
      Restraints::Angle* ca_h = find_angle(cc, "N", "CA", "H");
      Restraints::Angle* ca_h2 = find_angle(cc, "N", "CA", "H2");
      if (ca_h && std::isfinite(ca_h->value)) {
        ca_h3->value = ca_h->value;
        ca_h3->esd = ca_h->esd;
      } else if (ca_h2 && std::isfinite(ca_h2->value)) {
        ca_h3->value = ca_h2->value;
        ca_h3->esd = ca_h2->esd;
      } else {
        ca_h3->value = 109.5;
        ca_h3->esd = 3.0;
      }
    }
  }
  copy_angle("H", "H3", "H", "H2");
  copy_angle("H2", "H3", "H", "H2");
}

namespace {

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

void adjust_nitro_group(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);
  for (auto& n_atom : cc.atoms) {
    if (n_atom.el != El::N)
      continue;
    std::vector<std::string> o_neighbors;
    int c_single = 0;
    bool other = false;
    for (const std::string& nb : neighbors[n_atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0)
        continue;
      Element el = cc.atoms[idx].el;
      BondType bt = get_bond_type(cc, n_atom.id, nb);
      if (el == El::O) {
        if (bt == BondType::Double || bt == BondType::Deloc)
          o_neighbors.push_back(nb);
        else
          other = true;
      } else if (el == El::C) {
        if (bt == BondType::Single)
          ++c_single;
        else
          other = true;
      } else {
        other = true;
      }
    }
    if (other || c_single != 1 || o_neighbors.size() != 2)
      continue;

    std::string o_single = (o_neighbors[0] < o_neighbors[1]) ? o_neighbors[0] : o_neighbors[1];
    for (auto& bond : cc.rt.bonds) {
      if ((bond.id1.atom == n_atom.id && bond.id2.atom == o_single) ||
          (bond.id2.atom == n_atom.id && bond.id1.atom == o_single)) {
        bond.type = BondType::Single;
        break;
      }
    }
    int o_idx = cc.find_atom_index(o_single);
    if (o_idx >= 0)
      cc.atoms[o_idx].charge = -1.0f;
    n_atom.charge = 1.0f;
  }
}

void adjust_single_bond_oxide(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    if (std::fabs(atom.charge) > 0.5f)
      continue;
    int n_h = 0, n_heavy = 0;
    std::string heavy_id;
    for (const std::string& nb : neighbors[atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0)
        continue;
      if (cc.atoms[idx].el == El::H) {
        ++n_h;
      } else {
        ++n_heavy;
        heavy_id = nb;
      }
    }
    if (n_h != 0 || n_heavy != 1)
      continue;
    if (get_bond_type(cc, atom.id, heavy_id) != BondType::Single)
      continue;
    atom.charge = -1.0f;
  }
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
    std::string atom_id = atom.id;
    cc.atoms.push_back({new_h, "", El::H, 0.0f, "H", "", Position()});
    cc.rt.bonds.push_back({{1, atom_id}, {1, new_h}, BondType::Single, false, NAN, NAN, NAN, NAN});
    for (const std::string& nb : neighbors[atom_id])
      cc.rt.angles.push_back({{1, nb}, {1, atom_id}, {1, new_h}, NAN, NAN});
  }
}

void adjust_carboxy_asp(ChemComp& cc) {
  auto neighbors = make_neighbor_names(cc);
  std::set<std::string> matched_atoms;
  for (const auto& atom : cc.atoms) {
    if (atom.el != El::C)
      continue;
    std::vector<std::string> o_neighbors;
    std::string c_neighbor;
    bool ok = true;
    for (const std::string& nid : neighbors[atom.id]) {
      int idx = cc.find_atom_index(nid);
      if (idx < 0)
        continue;
      if (cc.atoms[idx].el == El::O)
        o_neighbors.push_back(nid);
      else if (cc.atoms[idx].el == El::C)
        c_neighbor = nid;
      else
        ok = false;
    }
    if (!ok || o_neighbors.size() != 2 || c_neighbor.empty())
      continue;
    bool c_is_sp3 = true;
    if (atom_has_unsaturated_bond(cc, c_neighbor))
      c_is_sp3 = false;
    if (!c_is_sp3)
      continue;
    matched_atoms.insert(o_neighbors[0]);
    matched_atoms.insert(o_neighbors[1]);
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

std::string acedrg_h_name(const ChemComp& cc, const std::string& n_id,
                          const std::vector<std::string>& h_on_n) {
  std::set<std::string> used_names;
  for (const auto& atom : cc.atoms)
    used_names.insert(atom.id);
  std::string root;
  for (char c : n_id)
    if (std::isalpha(static_cast<unsigned char>(c)))
      root += c;
  std::string h_root = (root.size() == 2) ? "H" + root.substr(1) : "H";
  if (h_root.size() >= 2 && used_names.find(h_root) == used_names.end())
    return h_root;
  if (!(h_root == "H" && !h_on_n.empty()) && used_names.find(h_root) == used_names.end())
    return h_root;
  auto numeric_suffix = [](const std::string& h) -> int {
    size_t d = 0;
    for (; d < h.size(); ++d)
      if (std::isdigit(static_cast<unsigned char>(h[d])))
        break;
    if (d >= h.size())
      return 0;
    const std::string suffix = h.substr(d);
    if (!std::all_of(suffix.begin(), suffix.end(),
                     [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }))
      return 0;
    return std::stoi(suffix);
  };
  int idx_max = 0;
  for (const std::string& h : h_on_n)
    idx_max = std::max(idx_max, numeric_suffix(h));
  std::string cand = cat(h_root, (idx_max > 0) ? idx_max + 1 : 2);
  if (used_names.find(cand) == used_names.end())
    return cand;
  for (int n = 2; n < 10000; ++n) {
    cand = cat(h_root, n);
    if (used_names.find(cand) == used_names.end())
      return cand;
  }
  return {};
}

void add_hydrogen_with_restraints(ChemComp& cc,
                                  std::map<std::string, std::vector<std::string>>& neighbors,
                                  const std::string& center_id,
                                  const std::string& new_h_id,
                                  const std::vector<std::string>& angle_flanks) {
  cc.atoms.push_back(ChemComp::Atom{new_h_id, "", El::H, 0.0f, "H", "", Position()});
  cc.rt.bonds.push_back({{1, center_id}, {1, new_h_id}, BondType::Single, false, NAN, NAN, NAN, NAN});
  for (const std::string& flank_id : angle_flanks)
    cc.rt.angles.push_back({{1, flank_id}, {1, center_id}, {1, new_h_id}, NAN, NAN});
  neighbors[center_id].push_back(new_h_id);
  neighbors[new_h_id].push_back(center_id);
}

void adjust_guanidinium_group(ChemComp& cc,
                              std::map<std::string, std::vector<std::string>>& neighbors) {
  for (auto& atom : cc.atoms) {
    if (atom.el != El::C)
      continue;
    std::string c_id = atom.id;
    const auto& nb = neighbors[c_id];
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
          if (onb_idx < 0 || cc.atoms[onb_idx].el == El::H)
            continue;
          Element onb_el = cc.atoms[onb_idx].el;
          if (onb_el != El::C && onb_el != El::Si && onb_el != El::Ge) {
            bool allow_hydroxyl = false;
            if (onb_el == El::O && get_bond_type(cc, other_n_id, onb) == BondType::Single) {
              for (const std::string& ob_nb : neighbors[onb]) {
                if (ob_nb == other_n_id)
                  continue;
                int ob_idx = cc.find_atom_index(ob_nb);
                if (ob_idx >= 0 && cc.atoms[ob_idx].is_hydrogen()) {
                  allow_hydroxyl = true;
                  break;
                }
              }
            }
            if (!allow_hydroxyl) {
              other_n_has_unsat_sub = true;
              break;
            }
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
      if (has_substituent || h_neighbors.size() != 1)
        continue;
      int n_idx = cc.find_atom_index(n_id);
      if (n_idx < 0)
        continue;
      cc.atoms[n_idx].charge = 1.0f;
      std::string new_h_id = acedrg_h_name(cc, n_id, h_neighbors);
      if (new_h_id.empty())
        continue;
      add_hydrogen_with_restraints(cc, neighbors, n_id, new_h_id, {c_id, h_neighbors[0]});
    }
  }
}

void adjust_amino_ter_amine(ChemComp& cc,
                            std::map<std::string, std::vector<std::string>>& neighbors) {
  for (auto& n1 : cc.atoms) {
    if (n1.el != El::N || std::fabs(n1.charge) > 0.5f)
      continue;
    std::vector<std::string> h_ids;
    std::string c1_id;
    for (const std::string& nb : neighbors[n1.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0) continue;
      if (cc.atoms[idx].el == El::H)
        h_ids.push_back(nb);
      else if (cc.atoms[idx].el == El::C)
        c1_id = nb;
    }
    if (h_ids.size() != 2 || c1_id.empty())
      continue;
    bool c1_has_double = false;
    for (const auto& bond : cc.rt.bonds)
      if ((bond.id1.atom == c1_id || bond.id2.atom == c1_id) &&
          bond.type == BondType::Double) {
        c1_has_double = true;
        break;
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
            if (amide_n_nb == c2_id) continue;
            int an_idx = cc.find_atom_index(amide_n_nb);
            if (an_idx >= 0 && cc.atoms[an_idx].el != El::H) {
              has_amide_n = true;
              break;
            }
          }
        }
      }
      if (!(has_carbonyl && has_amide_n))
        continue;
      n1.charge = 1.0f;
      std::string n1_id = n1.id;
      std::string new_h = acedrg_h_name(cc, n1_id, h_ids);
      if (new_h.empty())
        break;
      std::vector<std::string> angle_flanks = h_ids;
      angle_flanks.push_back(c1_id);
      add_hydrogen_with_restraints(cc, neighbors, n1_id, new_h, angle_flanks);
      break;
    }
  }
}

void adjust_terminal_amine(ChemComp& cc,
                           std::map<std::string, std::vector<std::string>>& neighbors) {
  for (auto& n_atom : cc.atoms) {
    if (n_atom.el != El::N || std::fabs(n_atom.charge) > 0.5f)
      continue;
    int n_carbon = 0, n_hydrogen = 0, n_other = 0;
    std::string alpha_carbon_id;
    std::vector<std::string> h_ids;
    for (const std::string& nb : neighbors[n_atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0) continue;
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
    if (n_carbon != 1 || n_other != 0 || n_hydrogen != 2)
      continue;
    if (atom_has_unsaturated_bond(cc, alpha_carbon_id) ||
        atom_has_unsaturated_bond(cc, n_atom.id) ||
        !has_carboxylic_acid_neighbor(cc, neighbors, alpha_carbon_id))
      continue;
    n_atom.charge = 1.0f;
    std::string n_atom_id = n_atom.id;
    std::string new_h = acedrg_h_name(cc, n_atom_id, h_ids);
    if (new_h.empty())
      continue;
    std::vector<std::string> angle_flanks = h_ids;
    angle_flanks.push_back(alpha_carbon_id);
    add_hydrogen_with_restraints(cc, neighbors, n_atom_id, new_h, angle_flanks);
  }
}

void adjust_protonated_amide_n(ChemComp& cc,
                               std::map<std::string, std::vector<std::string>>& neighbors) {
  auto is_carbonyl = [&](const std::string& c_id) {
    int c_idx = cc.find_atom_index(c_id);
    if (c_idx < 0 || cc.atoms[c_idx].el != El::C)
      return false;
    bool has_double_o = false;
    for (const auto& bond : cc.rt.bonds) {
      if (bond.id1.atom != c_id && bond.id2.atom != c_id)
        continue;
      if (bond.type != BondType::Double && bond.type != BondType::Deloc)
        continue;
      const std::string& other = (bond.id1.atom == c_id) ? bond.id2.atom : bond.id1.atom;
      int other_idx = cc.find_atom_index(other);
      if (other_idx >= 0 && cc.atoms[other_idx].el == El::O)
        has_double_o = true;
      else
        return false;
    }
    return has_double_o;
  };
  auto carbonyl_has_n_neighbor = [&](const std::string& c_id) {
    if (!is_carbonyl(c_id)) return false;
    auto it = neighbors.find(c_id);
    if (it == neighbors.end()) return false;
    for (const std::string& nb : it->second) {
      int idx = cc.find_atom_index(nb);
      if (idx >= 0 && cc.atoms[idx].el == El::N)
        return true;
    }
    return false;
  };
  auto has_adjacent_carbonyl_with_n = [&](const std::string& c_id) {
    auto it = neighbors.find(c_id);
    if (it == neighbors.end()) return false;
    for (const std::string& nb : it->second) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0 || cc.atoms[idx].el != El::C) continue;
      if (carbonyl_has_n_neighbor(nb))
        return true;
    }
    return false;
  };
  for (auto& n_atom : cc.atoms) {
    if (n_atom.el != El::N || std::fabs(n_atom.charge) > 0.5f ||
        atom_has_unsaturated_bond(cc, n_atom.id))
      continue;
    int n_hydrogen = 0, n_heavy = 0;
    std::string c_id;
    std::vector<std::string> h_ids;
    for (const std::string& nb : neighbors[n_atom.id]) {
      int idx = cc.find_atom_index(nb);
      if (idx < 0) continue;
      if (cc.atoms[idx].el == El::H) {
        ++n_hydrogen;
        h_ids.push_back(nb);
      } else {
        ++n_heavy;
        if (cc.atoms[idx].el == El::C)
          c_id = nb;
      }
    }
    if (n_hydrogen != 2 || n_heavy != 1 || c_id.empty() ||
        !is_carbonyl(c_id) || !has_adjacent_carbonyl_with_n(c_id))
      continue;
    n_atom.charge = 1.0f;
    std::string n_atom_id = n_atom.id;
    std::string new_h = acedrg_h_name(cc, n_atom_id, h_ids);
    if (new_h.empty())
      continue;
    std::vector<std::string> angle_flanks = h_ids;
    angle_flanks.push_back(c_id);
    add_hydrogen_with_restraints(cc, neighbors, n_atom_id, new_h, angle_flanks);
  }
}

} // namespace

void apply_chemical_adjustments(ChemComp& cc) {
  auto n_neighbors = make_neighbor_names(cc);
  struct AceChemRule {
    const char* name;
    std::function<void()> apply;
  };
  std::vector<AceChemRule> rules;
  rules.push_back({"oxoacid_phosphate", [&] { adjust_oxoacid_group(cc, El::P, 3, true); }});
  rules.push_back({"oxoacid_sulfate", [&] { adjust_oxoacid_group(cc, El::S, 4, false); }});
  rules.push_back({"nitro_group", [&] { adjust_nitro_group(cc); }});
  rules.push_back({"single_bond_oxide", [&] { adjust_single_bond_oxide(cc); }});
  rules.push_back({"hexafluorophosphate", [&] { adjust_hexafluorophosphate(cc); }});
  rules.push_back({"carboxy_asp", [&] { adjust_carboxy_asp(cc); }});
  rules.push_back({"terminal_carboxylate", [&] { adjust_terminal_carboxylate(cc); }});
  rules.push_back({"guanidinium", [&] { adjust_guanidinium_group(cc, n_neighbors); }});
  rules.push_back({"amino_ter_amine", [&] { adjust_amino_ter_amine(cc, n_neighbors); }});
  rules.push_back({"terminal_amine", [&] { adjust_terminal_amine(cc, n_neighbors); }});
  rules.push_back({"protonated_amide_n", [&] { adjust_protonated_amide_n(cc, n_neighbors); }});

  bool trace = ace_trace_mode();
  for (const AceChemRule& rule : rules) {
    AceRuleStats before = trace ? collect_rule_stats(cc) : AceRuleStats();
    rule.apply();
    if (trace) {
      AceRuleStats after = collect_rule_stats(cc);
      if (rule_stats_changed(before, after)) {
        std::fprintf(stderr,
                     "[ace-rule %s] atoms %+d bonds %+d angles %+d torsions %+d chirs %+d planes %+d charge %+0.3f\n",
                     rule.name,
                     after.atom_count - before.atom_count,
                     after.bond_count - before.bond_count,
                     after.angle_count - before.angle_count,
                     after.torsion_count - before.torsion_count,
                     after.chir_count - before.chir_count,
                     after.plane_count - before.plane_count,
                     after.charge_sum - before.charge_sum);
      }
    }
  }
}

} // namespace gemmi
