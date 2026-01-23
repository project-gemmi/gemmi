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
  bool has_oxt = false;
  bool has_hxt = false;
  for (const auto& atom : cc.atoms) {
    if (atom.id == "OXT" && atom.el == El::O)
      has_oxt = true;
    else if (atom.id == "HXT" && atom.el == El::H)
      has_hxt = true;
  }
  if (!has_oxt || !has_hxt)
    return;
  for (auto& atom : cc.atoms)
    if (atom.id == "OXT")
      atom.charge = -1.0f;
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

  cc.atoms.push_back(ChemComp::Atom{"H3", "", El::H, 0.0f, "H", Position()});
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
  for (auto& atom : cc.atoms) {
    if (atom.el == El::O &&
        (atom.id == "OP2" || atom.id == "OP3")) {
      atom.charge = -1.0f;
    }
  }
  std::vector<std::string> phos_h;
  for (const auto& atom : cc.atoms) {
    if (atom.el == El::H &&
        (atom.id == "HOP2" || atom.id == "HOP3")) {
      phos_h.push_back(atom.id);
    }
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
    adj[idx1].push_back({idx2, bond.type});
    adj[idx2].push_back({idx1, bond.type});
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

const ChemComp::Atom* pick_torsion_neighbor(
    const ChemComp& cc,
    const std::vector<std::vector<NeighborBond>>& adj,
    size_t center_idx,
    size_t exclude_idx) {
  std::vector<size_t> candidates;
  for (const auto& nb : adj[center_idx])
    if (nb.idx != exclude_idx)
      candidates.push_back(nb.idx);
  if (candidates.empty())
    return nullptr;

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
    int p_idx = element_priority(cc.atoms[idx].el);
    int p_best = element_priority(cc.atoms[best].el);
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

int shared_ring_size(const CodAtomInfo& a, const CodAtomInfo& b) {
  if (a.min_ring_size > 0 && b.min_ring_size > 0)
    return std::min(a.min_ring_size, b.min_ring_size);
  return 0;
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
    int p_idx = element_priority(cc.atoms[idx].el);
    int p_best = element_priority(cc.atoms[best].el);
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

void add_torsions_from_bonds_if_missing(ChemComp& cc, const AcedrgTables& tables) {
  if (!cc.rt.torsions.empty())
    return;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;

  auto adj = build_bond_adjacency(cc, atom_index);
  std::vector<CodAtomInfo> atom_info = tables.classify_atoms(cc);

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
    if (c1_carbonyl != c2_carbonyl) {
      center2 = c1_carbonyl ? idx1 : idx2;
      center3 = c1_carbonyl ? idx2 : idx1;
    } else if (cc.atoms[idx1].id == "CA" || cc.atoms[idx2].id == "CA") {
      center2 = (cc.atoms[idx1].id == "CA") ? idx1 : idx2;
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

    const ChemComp::Atom* a1 = pick_torsion_neighbor(cc, adj, center2, center3);
    const ChemComp::Atom* a4 = pick_torsion_neighbor(cc, adj, center3, center2);
    ring_size = shared_ring_size(atom_info[center2], atom_info[center3]);
    if (ring_size > 0) {
      const ChemComp::Atom* ring_a1 =
          pick_torsion_neighbor_in_ring(cc, adj, atom_info, center2, center3, center3);
      const ChemComp::Atom* ring_a4 =
          pick_torsion_neighbor_in_ring(cc, adj, atom_info, center3, center2, center2);
      if (ring_a1)
        a1 = ring_a1;
      if (ring_a4)
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
    if (tables.lookup_pep_tors(a1->id, cc.atoms[center2].id,
                               cc.atoms[center3].id, a4->id, tors_entry)) {
      cc.rt.torsions.push_back({"auto",
                                {1, a1->id},
                                {1, cc.atoms[center2].id},
                                {1, cc.atoms[center3].id},
                                {1, a4->id},
                                tors_entry.value, 10.0, tors_entry.period});
      continue;
    }

    double value = 180.0;
    double esd = 10.0;
    int period = 3;
    if (ring_size == 5 && sp3_2 && sp3_3) {
      bool ca_n_bond =
          (cc.atoms[center2].id == "CA" && cc.atoms[center3].el == El::N) ||
          (cc.atoms[center3].id == "CA" && cc.atoms[center2].el == El::N);
      if (ca_n_bond) {
        value = 60.0;
      } else if (a1->el == El::N || a4->el == El::N) {
        value = -60.0;
      }
    } else if ((sp2_2 && sp3_3) || (sp3_2 && sp2_3)) {
      value = 0.0;
      esd = 20.0;
      period = 6;
    } else if (sp2_2 && sp2_3) {
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

void add_chirality_if_missing(ChemComp& cc) {
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

void add_planes_if_missing(ChemComp& cc) {
  if (!cc.rt.planes.empty())
    return;

  std::map<std::string, size_t> atom_index;
  for (size_t i = 0; i < cc.atoms.size(); ++i)
    atom_index[cc.atoms[i].id] = i;
  auto adj = build_bond_adjacency(cc, atom_index);

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
    Restraints::Plane plane;
    plane.label = "plan-1";
    plane.esd = 0.02;
    plane.ids.push_back({1, cc.atoms[idx].id});
    plane.ids.push_back({1, other.front()});
    for (const std::string& o : oxy)
      plane.ids.push_back({1, o});
    cc.rt.planes.push_back(std::move(plane));
    break;
  }
}

enum OptionIndex {
  Tables=4, Sigma, Timing, CifStyle, OutputDir
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT.cif OUTPUT.cif"
    "\n " EXE_NAME " [options] --output-dir=DIR INPUT1.cif [INPUT2.cif ...]"
    "\n\nFill missing restraint values (bonds, angles) in a monomer CIF file"
    "\nusing COD/CSD statistical data from AceDRG tables."
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

        ChemComp cc = make_chemcomp_from_block(block);
        adjust_terminal_carboxylate(cc);
        add_angles_from_bonds_if_missing(cc);
        bool added_h3 = add_n_terminal_h3(cc);
        adjust_phosphate_group(cc);
        adjust_carboxylate_group(cc);

        // Count missing values before
        int missing_bonds = count_missing_values(cc.rt.bonds);
        int missing_angles = count_missing_values(cc.rt.angles);

        if (missing_bonds == 0 && missing_angles == 0) {
          if (verbose)
            std::fprintf(stderr, "  No missing values.\n");
          continue;
        }

        // Fill restraints
        timer.start();
        tables.fill_restraints(cc);
        if (added_h3)
          sync_n_terminal_h3_angles(cc);
        tables.assign_ccp4_types(cc);
        add_torsions_from_bonds_if_missing(cc, tables);
        add_chirality_if_missing(cc);
        add_planes_if_missing(cc);
        timer.print("Restraints filled in");

        // Count filled values
        int filled_bonds = missing_bonds - count_missing_values(cc.rt.bonds);
        int filled_angles = missing_angles - count_missing_values(cc.rt.angles);

        if (verbose)
          std::fprintf(stderr, "  Filled %d/%d bonds, %d/%d angles.\n",
                       filled_bonds, missing_bonds, filled_angles, missing_angles);

        filled_count += filled_bonds + filled_angles;

        // Update the block with new values
        add_chemcomp_to_block(cc, block);
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
