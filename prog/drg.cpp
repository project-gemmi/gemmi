// Copyright 2025 Global Phasing Ltd.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <set>
#include "gemmi/read_cif.hpp"     // for read_cif_gz
#include "gemmi/chemcomp.hpp"     // for ChemComp, make_chemcomp_from_block
#include "gemmi/acedrg_tables.hpp" // for AcedrgTables
#include "gemmi/to_cif.hpp"       // for write_cif_to_stream
#include "gemmi/to_chemcomp.hpp"  // for add_chemcomp_to_block
#include "gemmi/fstream.hpp"      // for Ofstream

#define GEMMI_PROG drg
#include "options.h"
#include "timer.h"

using namespace gemmi;

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

void adjust_terminal_carboxylate(ChemComp& cc) {
  // Check for standard amino acid C-terminus: OXT-C(=O)-...
  // Don't apply to other groups (e.g., boronic acid B(OH)3)
  auto oxt_it = cc.find_atom("OXT");
  auto hxt_it = cc.find_atom("HXT");
  if (oxt_it == cc.atoms.end() || hxt_it == cc.atoms.end())
    return;
  if (oxt_it->el != El::O || hxt_it->el != El::H)
    return;

  // Find what OXT is bonded to
  std::string oxt_neighbor;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == "OXT" && bond.id2.atom != "HXT")
      oxt_neighbor = bond.id2.atom;
    else if (bond.id2.atom == "OXT" && bond.id1.atom != "HXT")
      oxt_neighbor = bond.id1.atom;
  }
  if (oxt_neighbor.empty())
    return;

  // Check that OXT is bonded to a carbon (not boron, etc.)
  auto neighbor_it = cc.find_atom(oxt_neighbor);
  if (neighbor_it == cc.atoms.end() || neighbor_it->el != El::C)
    return;

  // Verify it's a carboxyl: the carbon should have another oxygen neighbor
  int o_count = 0;
  for (const auto& bond : cc.rt.bonds) {
    std::string other;
    if (bond.id1.atom == oxt_neighbor)
      other = bond.id2.atom;
    else if (bond.id2.atom == oxt_neighbor)
      other = bond.id1.atom;
    else
      continue;
    auto other_it = cc.find_atom(other);
    if (other_it != cc.atoms.end() && other_it->el == El::O)
      ++o_count;
  }
  if (o_count < 2)
    return;

  // It's a carboxylate - deprotonate OXT and remove HXT
  oxt_it->charge = -1.0f;
  remove_atom_by_id(cc, "HXT");
}

bool add_n_terminal_h3(ChemComp& cc) {
  if (cc.find_atom("N") == cc.atoms.end())
    return false;
  if (cc.find_atom("H3") != cc.atoms.end())
    return false;
  if (cc.find_atom("H") == cc.atoms.end() || cc.find_atom("H2") == cc.atoms.end())
    return false;

  bool n_has_h = false;
  bool n_has_h2 = false;
  bool n_has_h3 = false;
  for (const auto& bond : cc.rt.bonds) {
    if (bond.id1.atom == "N" || bond.id2.atom == "N") {
      const std::string& other = (bond.id1.atom == "N") ? bond.id2.atom : bond.id1.atom;
      if (other == "H")
        n_has_h = true;
      else if (other == "H2")
        n_has_h2 = true;
      else if (other == "H3")
        n_has_h3 = true;
    }
  }
  if (!n_has_h || !n_has_h2 || n_has_h3)
    return false;

  cc.atoms.push_back(ChemComp::Atom{"H3", "", El::H, 0.0f, "H", "", Position()});
  // Bond/angle values set to NAN - will be filled by fill_restraints()
  cc.rt.bonds.push_back({{1, "N"}, {1, "H3"}, BondType::Single, false,
                        NAN, NAN, NAN, NAN});

  auto add_angle_if_present = [&](const std::string& a1, const std::string& a2,
                                  const std::string& a3) {
    if (cc.find_atom(a1) == cc.atoms.end() ||
        cc.find_atom(a2) == cc.atoms.end() ||
        cc.find_atom(a3) == cc.atoms.end())
      return;
    cc.rt.angles.push_back({{1, a1}, {1, a2}, {1, a3}, NAN, NAN});  // filled later
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
  // Identify phosphate oxygens with attached H so we can deprotonate them.
  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }

  std::vector<std::string> phos_h;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    std::string h_id;
    bool has_p = false;
    for (const std::string& nid : nb) {
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
      if (el == El::H)
        h_id = nid;
      else if (el == El::P)
        has_p = true;
    }
    if (!has_p || h_id.empty())
      continue;
    atom.charge = -1.0f;
    phos_h.push_back(h_id);
  }
  for (const std::string& atom_id : phos_h)
    remove_atom_by_id(cc, atom_id);
}

void adjust_carboxylate_group(ChemComp& cc) {
  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }

  std::vector<std::string> hydrogens_to_remove;
  for (auto& atom : cc.atoms) {
    if (atom.el != El::O)
      continue;
    const auto& nb = neighbors[atom.id];
    std::string h_id;
    std::string c_id;
    for (const std::string& nid : nb) {
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      Element el = cc.atoms[it->second].el;
      if (el == El::H)
        h_id = nid;
      else if (el == El::C)
        c_id = nid;
    }
    if (h_id.empty() || c_id.empty())
      continue;
    const auto& c_nb = neighbors[c_id];
    int o_count = 0;
    for (const std::string& nid : c_nb) {
      auto it = atom_index.find(nid);
      if (it == atom_index.end())
        continue;
      if (cc.atoms[it->second].el == El::O)
        o_count += 1;
    }
    if (o_count < 2)
      continue;
    atom.charge = -1.0f;
    hydrogens_to_remove.push_back(h_id);
  }

  for (const std::string& h_id : hydrogens_to_remove)
    remove_atom_by_id(cc, h_id);
}

// Protonate guanidinium groups to neutral pH form (+NH2 instead of =NH)
// Guanidinium has pKa ~12.5, so it's protonated at neutral pH.
void adjust_guanidinium_group(ChemComp& cc) {
  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  std::map<std::string, std::vector<std::string>> neighbors;
  for (const auto& bond : cc.rt.bonds) {
    neighbors[bond.id1.atom].push_back(bond.id2.atom);
    neighbors[bond.id2.atom].push_back(bond.id1.atom);
  }

  // Find guanidinium carbons: C bonded to exactly 3 N atoms
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

    // This is a guanidinium carbon - find N with double bond
    for (const std::string& n_id : n_neighbors) {
      // Check if this N has a double bond to the guanidinium C
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

      // Count hydrogens on this nitrogen
      const auto& n_nb = neighbors[n_id];
      std::vector<std::string> h_neighbors;
      for (const std::string& hid : n_nb) {
        auto it = atom_index.find(hid);
        if (it != atom_index.end() && cc.atoms[it->second].el == El::H)
          h_neighbors.push_back(hid);
      }

      // If only 1 hydrogen, add another (protonate =NH to +NH2)
      if (h_neighbors.size() == 1) {
        auto n_it = atom_index.find(n_id);
        if (n_it == atom_index.end())
          continue;

        // Set charge on the nitrogen
        cc.atoms[n_it->second].charge = 1.0f;

        // Generate new hydrogen name - acedrg uses sequential H3, H4, etc.
        std::string new_h_id;
        for (int n = 3; n < 100; ++n) {
          new_h_id = "H" + std::to_string(n);
          if (cc.find_atom(new_h_id) == cc.atoms.end())
            break;
        }
        if (cc.find_atom(new_h_id) != cc.atoms.end())
          continue;  // Can't find unique name, skip

        // Add the new hydrogen atom
        cc.atoms.push_back(ChemComp::Atom{new_h_id, "", El::H, 0.0f, "H", "", Position()});
        // Add bond (value will be filled by fill_restraints)
        cc.rt.bonds.push_back({{1, n_id}, {1, new_h_id}, BondType::Single, false,
                              NAN, NAN, NAN, NAN});
        // Update atom_index and neighbors for subsequent iterations
        atom_index[new_h_id] = cc.atoms.size() - 1;
        neighbors[n_id].push_back(new_h_id);
        neighbors[new_h_id].push_back(n_id);
      }
    }
  }
}

void add_angles_from_bonds_if_missing(ChemComp& cc) {
  if (!cc.rt.angles.empty())
    return;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

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
      return na > nb; // prefer higher index for HB1/2/3
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

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables) {
  if (!cc.rt.torsions.empty())
    return;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  auto adj = build_bond_adjacency(cc, atom_index);
  std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);
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
    const AcedrgTables& tables) {
  if (!cc.rt.chirs.empty())
    return;

  if (!atom_stereo.empty()) {
    std::map<std::string, size_t> atom_index;
    for (size_t i = 0; i < cc.atoms.size(); ++i)
      atom_index[cc.atoms[i].id] = i;
    auto adj = build_bond_adjacency(cc, atom_index);
    std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);

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

void add_planes_if_missing(ChemComp& cc, const AcedrgTables& tables) {
  if (!cc.rt.planes.empty())
    return;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;
  auto adj = build_bond_adjacency(cc, atom_index);
  std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);

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

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle, OutputDir, TypeOut
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n " EXE_NAME " [options] --output-dir=DIR INPUT1.cif [INPUT2.cif ...]"
    "\n\nFill missing restraint values (bonds, angles) in a monomer CIF file"
    "\nusing COD/CSD statistical data from AceDRG tables."
    "\nIf OUTPUT.cif is -, the output is printed to stdout."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Tables, 0, "t", "tables", Arg::Required,
    "  -t, --tables=DIR  \tDirectory with AceDRG tables (default: $ACEDRG_TABLES"
    "\n\t\tor $CCP4/share/acedrg/tables)." },
  { OutputDir, 0, "o", "output-dir", Arg::Required,
    "  -o, --output-dir=DIR  \tOutput directory for batch processing." },
  { Sigma, 0, "", "sigma", Arg::Float,
    "  --sigma=NUM  \tMaximum sigma for bond restraints (default: 0.02)." },
  { Timing, 0, "", "timing", Arg::None,
    "  --timing  \tPrint timing information." },
  { CifStyle, 0, "", "style", Arg::CifStyle,
    "  --style=STYLE  \tOutput style: default, pdbx, aligned." },
  { TypeOut, 0, "", "typeOut", Arg::None,
    "  --typeOut  \tWrite _chem_comp_atom.atom_type column." },
  { 0, 0, 0, 0, 0, 0 }
};

std::string get_filename(const std::string& path) {
  size_t pos = path.find_last_of("/\\");
  return pos == std::string::npos ? path : path.substr(pos + 1);
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  bool batch_mode = p.options[OutputDir];
  if (batch_mode) {
    if (p.nonOptionsCount() < 1) {
      std::fprintf(stderr, "ERROR: At least one input file required with --output-dir.\n");
      return 1;
    }
  } else {
    p.require_positional_args(2);
  }

  int verbose = p.options[Verbose].count();

  // Get tables directory
  std::string tables_dir;
  if (p.options[Tables]) {
    tables_dir = p.options[Tables].arg;
  } else {
    // Try ACEDRG_TABLES environment variable
    const char* env = std::getenv("ACEDRG_TABLES");
    if (env) {
      tables_dir = env;
    } else {
      // Fallback to $CCP4/share/acedrg/tables
      const char* ccp4 = std::getenv("CCP4");
      if (ccp4)
        tables_dir = std::string(ccp4) + "/share/acedrg/tables";
    }
  }

  if (tables_dir.empty()) {
    std::fprintf(stderr, "ERROR: No tables directory specified.\n"
                         "Use --tables=DIR or set ACEDRG_TABLES or CCP4 environment variable.\n");
    return 1;
  }

  Timer timer(p.options[Timing]);

  try {
    // Load tables
    if (verbose)
      std::fprintf(stderr, "Loading tables from %s ...\n", tables_dir.c_str());
    timer.start();
    AcedrgTables tables;
    tables.load_tables(tables_dir);
    timer.print("Tables loaded in");

    if (p.options[Sigma])
      tables.lower_bond_sigma = std::strtod(p.options[Sigma].arg, nullptr);

    // Set verbose level: -v=atoms, -vv=lookup, -vvv=1D/2D failures
    tables.verbose = verbose;

    // Build list of (input, output) pairs
    std::vector<std::pair<std::string, std::string>> files;
    if (batch_mode) {
      std::string output_dir = p.options[OutputDir].arg;
      // Ensure output_dir ends with separator
      if (!output_dir.empty() && output_dir.back() != '/' && output_dir.back() != '\\')
        output_dir += '/';
      for (int i = 0; i < p.nonOptionsCount(); ++i) {
        std::string input = p.nonOption(i);
        std::string output = output_dir + get_filename(input);
        files.emplace_back(input, output);
      }
    } else {
      files.emplace_back(p.nonOption(0), p.nonOption(1));
    }

    int total_filled = 0;
    for (const auto& file_pair : files) {
      const std::string& input = file_pair.first;
      const std::string& output = file_pair.second;

      // Read input CIF
      if (verbose)
        std::fprintf(stderr, "Reading %s ...\n", input.c_str());
      timer.start();
      cif::Document doc = read_cif_gz(input);
      timer.print("Input CIF read in");

      int filled_count = 0;
      for (cif::Block& block : doc.blocks) {
        // Skip blocks that don't look like monomer definitions
        if (!block.find_values("_chem_comp_atom.atom_id"))
          continue;

        if (verbose)
          std::fprintf(stderr, "Processing block %s ...\n", block.name.c_str());

        std::map<std::string, std::string> atom_stereo;
        cif::Table stereo_tab = block.find("_chem_comp_atom.",
                                           {"atom_id", "pdbx_stereo_config"});
        if (stereo_tab.ok()) {
          for (auto row : stereo_tab) {
            std::string stereo = row.str(1);
            if (stereo.empty() || stereo[0] == '.' || stereo[0] == '?')
              continue;
            atom_stereo.emplace(row.str(0), std::move(stereo));
          }
        }

        ChemComp cc = make_chemcomp_from_block(block);
        add_angles_from_bonds_if_missing(cc);
        adjust_phosphate_group(cc);
        adjust_carboxylate_group(cc);
        adjust_guanidinium_group(cc);

        // Convert to zwitterionic form (add H3, deprotonate carboxyl) BEFORE
        // classification and fill_restraints. Acedrg tables are built using
        // the zwitterionic form where N has 4 neighbors (CA, H, H2, H3).
        bool added_h3 = add_n_terminal_h3(cc);
        adjust_terminal_carboxylate(cc);

        // Count missing values before
        int missing_bonds = count_missing_values(cc.rt.bonds);
        int missing_angles = count_missing_values(cc.rt.angles);

        bool need_fill = (missing_bonds > 0 || missing_angles > 0);
        bool need_typeout = p.options[TypeOut];

        if (!need_fill && !need_typeout) {
          if (verbose)
            std::fprintf(stderr, "  No missing values.\n");
          continue;
        }

        // Compute acedrg_types if requested (on zwitterionic form with H3)
        std::vector<std::string> acedrg_types;
        if (p.options[TypeOut])
          acedrg_types = tables.compute_acedrg_types(cc);

        // Fill restraints if needed
        if (need_fill) {
          timer.start();
          tables.fill_restraints(cc);
          // Sync H3 angles from H/H2 angles (they use the same values)
          if (added_h3)
            sync_n_terminal_h3_angles(cc);
          tables.assign_ccp4_types(cc);
          add_torsions_from_bonds_if_missing(cc, tables);
          add_chirality_if_missing(cc, atom_stereo, tables);
          add_planes_if_missing(cc, tables);
          timer.print("Restraints filled in");

          // Count filled values
          int filled_bonds = missing_bonds - count_missing_values(cc.rt.bonds);
          int filled_angles = missing_angles - count_missing_values(cc.rt.angles);

          if (verbose)
            std::fprintf(stderr, "  Filled %d/%d bonds, %d/%d angles.\n",
                         filled_bonds, missing_bonds, filled_angles, missing_angles);

          filled_count += filled_bonds + filled_angles;
        } else {
          // H3 and carboxylate already adjusted above, just sync angles if needed
          if (added_h3)
            sync_n_terminal_h3_angles(cc);
        }

        // Update the block with new values
        add_chemcomp_to_block(cc, block, acedrg_types);
      }

      if (verbose)
        std::fprintf(stderr, "Writing %s ...\n", output.c_str());

      timer.start();
      Ofstream os(output, &std::cout);
      write_cif_to_stream(os.ref(), doc, cif_write_options(p.options[CifStyle]));
      timer.print("Output written in");

      if (verbose)
        std::fprintf(stderr, "Done with %s. Filled %d restraint values.\n",
                     input.c_str(), filled_count);
      total_filled += filled_count;
    }

    if (batch_mode && verbose)
      std::fprintf(stderr, "All done. Processed %zu files, filled %d total restraint values.\n",
                   files.size(), total_filled);

  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
