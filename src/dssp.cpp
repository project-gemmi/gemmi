// Copyright Global Phasing Ltd.

#include <gemmi/topo.hpp>
#include <gemmi/dssp.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/calculate.hpp>
#include <gemmi/unitcell.hpp>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace gemmi {

namespace {

// Constants from DSSP algorithm
constexpr double COUPLING_CONSTANT = 27.888;  // = 332 * 0.42 * 0.2
constexpr double MIN_ATOM_DISTANCE = 0.5;      // Angstrom
constexpr double MIN_ENERGY = -9.9;            // kcal/mol
constexpr double MAX_HBOND_DISTANCE = 0.35;    // nm for geometry-based hbonds
constexpr double MAX_HBOND_ANGLE = 30.0;       // degrees for geometry-based hbonds

double calculate_atomic_distance(const Atom* atom1, const Atom* atom2) {
  if (!atom1 || !atom2) return 0.0;
  return atom1->pos.dist(atom2->pos);
}

double calculate_dihedral_angle(Atom* a1, Atom* a2, Atom* a3, Atom* a4) {
  if (!a1 || !a2 || !a3 || !a4) return 0.0;
  return calculate_dihedral(a1->pos, a2->pos, a3->pos, a4->pos);
}


} // anonymous namespace

// DsspCalculator implementation
std::string DsspCalculator::calculate_secondary_structure(NeighborSearch& ns, Topo& topo) {
  // Clear previous calculations
  ss_info.clear();
  bridges_.clear();
  res_infos.clear();

  // Build working set of residue info pointers
  for (Topo::ChainInfo& chain_info : topo.chain_infos) {
    for (Topo::ResInfo& res_info : chain_info.res_infos) {
      if (res_info.res && res_info.res->get_ca()) {  // Only include residues with CA atoms
        res_infos.push_back(&res_info);
      }
    }
  }

  // Initialize secondary structure info
  ss_info.resize(res_infos.size());

  // Calculate hydrogen bonds
  calculate_hydrogen_bonds(topo, ns);

  // Find bends and breaks first
  find_bends_and_breaks();

  // Find bridges and strands
  find_bridges_and_strands();

  // Find turns and helices
  find_turns_and_helices();

  // Find polyproline helices if enabled
  if (options.search_polyproline) {
    find_polyproline_helices();
  }

  // Generate final secondary structure string
  return generate_ss_string();
}


void DsspCalculator::calculate_hydrogen_bonds(Topo& topo, NeighborSearch& ns) {
  // Process only residues in our working set
  for (size_t i = 0; i < res_infos.size(); ++i) {
    Topo::ResInfo* res_info = res_infos[i];
    const Atom* ca = res_info->res->get_ca();
    if (!ca || !res_info->res || !res_info->res->get_n())
      continue;

    // Uses neighbor search for efficiency
    auto marks = ns.find_neighbors(*ca, 0, 0);
    for (gemmi::NeighborSearch::Mark* mark : marks) {
      CRA cra = mark->to_cra(*ns.model);
      if (!cra.residue || !cra.residue->get_c() || !cra.residue->find_atom("O", '*', El::O))
        continue;

      // Find the corresponding ResInfo in our working set
      Topo::ResInfo* neighbor_resinfo = nullptr;
      size_t j = 0;
      for (j = 0; j < res_infos.size(); ++j) {
        if (res_infos[j]->res == cra.residue) {
          neighbor_resinfo = res_infos[j];
          break;
        }
      }

      // Skip self-interactions and avoid duplicate processing by only processing when i < j
      if (!neighbor_resinfo || neighbor_resinfo == res_info || j <= i)
        continue;

      if (options.hbond_definition == HBondDefinition::Energy) {
        calculate_hbond_energy(res_info, neighbor_resinfo);
        calculate_hbond_energy(neighbor_resinfo, res_info);
      } else {
        calculate_hbond_geometry(res_info, neighbor_resinfo);
        calculate_hbond_geometry(neighbor_resinfo, res_info);
      }
    }
  }
}

static Position calculate_h_pos(Topo::ResInfo* donor, const Atom* donor_n, HydrogenMode h_mode) {

  Position h_pos = donor_n->pos;
  if (h_mode == HydrogenMode::Existing) {
    h_pos = donor->res->find_atom("H", '*', El::H)->pos;
  } else { // HydrogenMode::Calculate
    if (!donor->prev.empty()) {
      if (Residue* prev  = donor->prev[0].res1) {
        // Calculate hydrogen position from previous residue C-O vector
        const Atom* prev_o = prev->find_atom("O", '*', El::O);
        const Atom* prev_c = prev->get_c();
        if (prev_o && prev_c) {
          // official DSSP implementation assumes that N-H vector is a normalized
          // O-C of the previous residue
          Vec3 prev_co = (prev_o->pos - prev_c->pos).normalized();
          h_pos -= Position(prev_co);
        }
      }
    }
  }
  return h_pos;
}

void DsspCalculator::calculate_hbond_energy(Topo::ResInfo* donor, Topo::ResInfo* acceptor) {
  // Check for proline and required atoms
  if (donor->res->name == "PRO" ||
      !acceptor->res->get_c() || !acceptor->res->get_o() ||
      !donor->res->get_n()) {
    return;
  }

  // Check CA distance first for efficiency
  const Atom* donor_ca = donor->res->get_ca();
  const Atom* acceptor_ca = acceptor->res->get_ca();
  if (donor_ca && acceptor_ca) {
    double ca_dist = calculate_atomic_distance(donor_ca, acceptor_ca);
    if (ca_dist > options.min_ca_distance) {
      return;
    }
  }

  // Calculate distances
  const Atom* donor_n = donor->res->get_n();
  const Atom* acceptor_o = acceptor->res->get_o();
  const Atom* acceptor_c = acceptor->res->get_c();
  double dist_NO = donor_n->pos.dist(acceptor_o->pos);
  double dist_NC = donor_n->pos.dist(acceptor_c->pos);
  Position h_pos = calculate_h_pos(donor, donor_n, options.hydrogen_mode);
  double dist_HO = h_pos.dist(acceptor->res->get_o()->pos);
  double dist_HC = h_pos.dist(acceptor->res->get_c()->pos);

  // Calculate hydrogen bond energy
  double energy = 0;
  if (dist_NO < MIN_ATOM_DISTANCE || dist_HC < MIN_ATOM_DISTANCE ||
      dist_HO < MIN_ATOM_DISTANCE || dist_NC < MIN_ATOM_DISTANCE) {
    energy = 0;
  } else {
    energy = COUPLING_CONSTANT * ((1.0/dist_NO) + (1.0/dist_HC) - (1.0/dist_HO) - (1.0/dist_NC));
  }

  if (energy < MIN_ENERGY) {
    energy = MIN_ENERGY;
  }

  // Store the hydrogen bond information
  if (energy < donor->acceptor_energies[0]) {
    donor->acceptors[1] = donor->acceptors[0];
    donor->acceptor_energies[1] = donor->acceptor_energies[0];
    donor->acceptors[0] = acceptor;
    donor->acceptor_energies[0] = energy;
  } else if (energy < donor->acceptor_energies[1]) {
    donor->acceptors[1] = acceptor;
    donor->acceptor_energies[1] = energy;
  }

  if (energy < acceptor->donor_energies[0]) {
    acceptor->donors[1] = acceptor->donors[0];
    acceptor->donor_energies[1] = acceptor->donor_energies[0];
    acceptor->donors[0] = donor;
    acceptor->donor_energies[0] = energy;
  } else if (energy < acceptor->donor_energies[1]) {
    acceptor->donors[1] = donor;
    acceptor->donor_energies[1] = energy;
  }
}

void DsspCalculator::calculate_hbond_geometry(Topo::ResInfo* donor, Topo::ResInfo* acceptor) {
  if (donor->res->name == "PRO" ||
      !acceptor->res->get_c() || !acceptor->res->get_o() || !donor->res->get_n()) {
    return;
  }

  // Check distance criterion
  const Atom* donor_n = donor->res->get_n();
  const Atom* acceptor_o = acceptor->res->get_o();
  if (!donor_n || !acceptor_o)
    return;
  Position o_pos = acceptor_o->pos;
  Position n_pos = donor_n->pos;
  double dist_NO = (n_pos - o_pos).length();

  if (dist_NO > MAX_HBOND_DISTANCE)
    return;
  Position h_pos = calculate_h_pos(donor, donor_n, options.hydrogen_mode);

  // Check angle criterion
  double angle = (o_pos - n_pos).angle(h_pos - n_pos);
  if (angle <= MAX_HBOND_ANGLE) {
    // Store hydrogen bond (geometry-based uses nullptr for energies)
    if (donor->acceptors[0] == nullptr) {
      donor->acceptors[0] = acceptor;
      acceptor->donors[0] = donor;
    } else if (donor->acceptors[1] == nullptr) {
      donor->acceptors[1] = donor->acceptors[0];
      donor->acceptors[0] = acceptor;
      acceptor->donors[1] = acceptor->donors[0];
      acceptor->donors[0] = donor;
    }
  }
}

bool DsspCalculator::has_hbond_between(Topo::ResInfo* donor, Topo::ResInfo* acceptor) const {
    for (int i = 0; i < 2; ++i) {
        if (donor->acceptors[i] == acceptor && donor->acceptor_energies[i] < options.hbond_energy_cutoff) {
            return true;
        }
    }

    for (int i = 0; i < 2; ++i) {
        if (acceptor->donors[i] == donor && acceptor->donor_energies[i] < options.hbond_energy_cutoff) {
            return true;
        }
    }

    return false;
}

bool DsspCalculator::no_chain_breaks_between(size_t res1_idx, size_t res2_idx) const {
  size_t start = std::min(res1_idx, res2_idx);
  size_t end = std::max(res1_idx, res2_idx);

  for (size_t i = start; i < end; ++i)
    if (ss_info[i].has_break)
      return false;
  return true;
}

BridgeType DsspCalculator::calculate_bridge_type(size_t i, size_t j) const {
  if (i >= res_infos.size() || j >= res_infos.size())
    return BridgeType::None;

  // antiparallel
  bool anti1 = has_hbond_between(res_infos[i], res_infos[j]) && has_hbond_between(res_infos[j], res_infos[i]);
  bool anti2 = (i > 0 && j + 1 < res_infos.size() && j > 0 && i + 1 < res_infos.size()) &&
               has_hbond_between(res_infos[i - 1], res_infos[j + 1]) && has_hbond_between(res_infos[j - 1], res_infos[i + 1]);
  if (anti1 && anti2)
    return BridgeType::AntiParallel;
  // parallel
  bool para1 = (j + 1 < res_infos.size() && i + 1 < res_infos.size()) &&
               has_hbond_between(res_infos[i], res_infos[j + 1]) && has_hbond_between(res_infos[j], res_infos[i + 1]);
  bool para2 = (i > 0 && j > 0) &&
               has_hbond_between(res_infos[i - 1], res_infos[j]) && has_hbond_between(res_infos[j - 1], res_infos[i]);
  if (para1 && para2)
    return BridgeType::Parallel;

  return BridgeType::None;
}

void DsspCalculator::find_bridges_and_strands() {
  // Find bridges
  for (size_t i = 1; i + 1 < res_infos.size(); ++i) {
    for (size_t j = i + 1; j < res_infos.size(); ++j) {
      BridgeType bridge_type = calculate_bridge_type(i, j);
      if (bridge_type != BridgeType::None) {
        bridges_.emplace_back(Bridge{i, j, bridge_type});
      }
    }
  }

  // Find strands from bridge patterns
  for (size_t i = 0; i < bridges_.size(); ++i) {
      for (size_t j = i + 1; j < bridges_.size(); ++j) {
          const Bridge& b1 = bridges_[i];
          const Bridge& b2 = bridges_[j];
          if (b1.type == b2.type && no_chain_breaks_between(b1.partner1, b1.partner2) && no_chain_breaks_between(b2.partner1, b2.partner2)) {
              if (b1.partner1 == b2.partner1 + 1 && b1.partner2 == b2.partner2 - 1 && b1.type == BridgeType::AntiParallel) {
                  ss_info[b1.partner1].ss_type = SecondaryStructure::Strand;
                  ss_info[b1.partner2].ss_type = SecondaryStructure::Strand;
                  ss_info[b2.partner1].ss_type = SecondaryStructure::Strand;
                  ss_info[b2.partner2].ss_type = SecondaryStructure::Strand;
              }
              if (b1.partner1 == b2.partner1 + 1 && b1.partner2 == b2.partner2 + 1 && b1.type == BridgeType::Parallel) {
                  ss_info[b1.partner1].ss_type = SecondaryStructure::Strand;
                  ss_info[b1.partner2].ss_type = SecondaryStructure::Strand;
                  ss_info[b2.partner1].ss_type = SecondaryStructure::Strand;
                  ss_info[b2.partner2].ss_type = SecondaryStructure::Strand;
              }
          }
      }
  }

  // Mark isolated bridges
  for (const auto& bridge : bridges_) {
      if (ss_info[bridge.partner1].ss_type == SecondaryStructure::Loop)
        ss_info[bridge.partner1].ss_type = SecondaryStructure::Bridge;
      if (ss_info[bridge.partner2].ss_type == SecondaryStructure::Loop)
        ss_info[bridge.partner2].ss_type = SecondaryStructure::Bridge;
  }
}

void DsspCalculator::find_turns_and_helices() {
  // First pass: mark all helix positions based on individual H-bonds
  for (int turn_type = 0; turn_type < 3; ++turn_type) { // 3, 4, 5 turns
    size_t stride = turn_type + 3;
    TurnType turn = static_cast<TurnType>(stride);

    for (size_t j = 0; j + stride < res_infos.size(); ++j) {
      if (has_hbond_between(res_infos[j + stride], res_infos[j]) && no_chain_breaks_between(j, j + stride)) {
        // Mark end position
        ss_info[j + stride].set_helix_position(turn, HelixPosition::End);

        // Mark middle positions
        for (size_t k = 1; k < stride; ++k) {
          if (ss_info[j + k].get_helix_position(turn) == HelixPosition::None) {
            ss_info[j + k].set_helix_position(turn, HelixPosition::Middle);
          }
        }

        // Mark start position
        if (ss_info[j].get_helix_position(turn) == HelixPosition::End) {
          ss_info[j].set_helix_position(turn, HelixPosition::StartAndEnd);
        } else {
          ss_info[j].set_helix_position(turn, HelixPosition::Start);
        }
      }
    }
  }

  // Second pass: detect helix patterns from overlapping H-bonds
  // Process in order of precedence: alpha-helix (4-turn), 3-10 helix (3-turn), pi-helix (5-turn)

  // Alpha helix detection (4-turn, highest precedence)
  for (size_t i = 0; i + 4 < res_infos.size(); ++i) {
    // Look for overlapping 4-turn patterns: i→i+4 AND i+1→i+5
    if (has_hbond_between(res_infos[i + 4], res_infos[i]) && has_hbond_between(res_infos[i + 5], res_infos[i + 1]) &&
        no_chain_breaks_between(i, i + 5)) {

      // Found overlapping alpha-helix pattern - mark as helix
      for (size_t k = i + 1; k <= i + 4; ++k) {
        ss_info[k].ss_type = SecondaryStructure::Helix_4;
      }

      // Continue the helix while overlapping patterns exist
      for (size_t j = i + 2; j + 4 < res_infos.size(); ++j) {
        if (has_hbond_between(res_infos[j + 4], res_infos[j]) && no_chain_breaks_between(j, j + 4)) {
          for (size_t k = std::max(j + 1, i + 5); k <= j + 4; ++k) {
            ss_info[k].ss_type = SecondaryStructure::Helix_4;
          }
        } else {
          break;
        }
      }
    }
  }

  // 3-10 helix detection (3-turn)
  for (size_t i = 0; i + 3 < res_infos.size(); ++i) {
    // Look for overlapping 3-turn patterns: i→i+3 AND i+1→i+4
    if (has_hbond_between(res_infos[i + 3], res_infos[i]) && has_hbond_between(res_infos[i + 4], res_infos[i + 1]) &&
        no_chain_breaks_between(i, i + 4)) {

      // Check if positions are available (not already alpha-helix)
      bool can_assign = true;
      for (size_t k = i + 1; k <= i + 3 && can_assign; ++k) {
        can_assign = (ss_info[k].ss_type < SecondaryStructure::Helix_4);
      }

      if (can_assign) {
        // Found overlapping 3-10 helix pattern
        for (size_t k = i + 1; k <= i + 3; ++k) {
          ss_info[k].ss_type = SecondaryStructure::Helix_3;
        }

        // Continue the helix while overlapping patterns exist
        for (size_t j = i + 2; j + 3 < res_infos.size(); ++j) {
          if (has_hbond_between(res_infos[j + 3], res_infos[j]) && no_chain_breaks_between(j, j + 3)) {
            bool can_extend = true;
            for (size_t k = std::max(j + 1, i + 4); k <= j + 3 && can_extend; ++k) {
              can_extend = (ss_info[k].ss_type < SecondaryStructure::Helix_4);
            }
            if (can_extend) {
              for (size_t k = std::max(j + 1, i + 4); k <= j + 3; ++k) {
                ss_info[k].ss_type = SecondaryStructure::Helix_3;
              }
            } else {
              break;
            }
          } else {
            break;
          }
        }
      }
    }
  }

  // Pi-helix detection (5-turn, lowest precedence)
  for (size_t i = 0; i + 5 < res_infos.size(); ++i) {
    // Look for overlapping 5-turn patterns: i→i+5 AND i+1→i+6
    if (has_hbond_between(res_infos[i + 5], res_infos[i]) && has_hbond_between(res_infos[i + 6], res_infos[i + 1]) &&
        no_chain_breaks_between(i, i + 6)) {

      // Check if positions are available
      bool can_assign = true;
      for (size_t k = i + 1; k <= i + 5 && can_assign; ++k) {
        SecondaryStructure current = ss_info[k].ss_type;
        can_assign = (current < SecondaryStructure::Helix_3) ||
                    (options.pi_helix_preference && current == SecondaryStructure::Helix_4);
      }

      if (can_assign) {
        // Found overlapping pi-helix pattern
        for (size_t k = i + 1; k <= i + 5; ++k) {
          ss_info[k].ss_type = SecondaryStructure::Helix_5;
        }

        // Continue the helix while overlapping patterns exist
        for (size_t j = i + 2; j + 5 < res_infos.size(); ++j) {
          if (has_hbond_between(res_infos[j + 5], res_infos[j]) && no_chain_breaks_between(j, j + 5)) {
            bool can_extend = true;
            for (size_t k = std::max(j + 1, i + 6); k <= j + 5 && can_extend; ++k) {
              SecondaryStructure current = ss_info[k].ss_type;
              can_extend = (current < SecondaryStructure::Helix_3) ||
                          (options.pi_helix_preference && current == SecondaryStructure::Helix_4);
            }
            if (can_extend) {
              for (size_t k = std::max(j + 1, i + 6); k <= j + 5; ++k) {
                ss_info[k].ss_type = SecondaryStructure::Helix_5;
              }
            } else {
              break;
            }
          } else {
            break;
          }
        }
      }
    }
  }

  // Third pass: mark remaining isolated turns
  for (size_t i = 1; i + 1 < res_infos.size(); ++i) {
    if (ss_info[i].ss_type == SecondaryStructure::Loop) {
      bool is_turn = false;

      // Check if this residue is part of any isolated turn pattern (not already in helix)
      for (int turn_type = 0; turn_type < 3 && !is_turn; ++turn_type) {
        size_t stride = turn_type + 3;
        TurnType turn = static_cast<TurnType>(stride);

        for (size_t k = 1; k < stride && !is_turn; ++k) {
          if (i >= k) {
            HelixPosition pos = ss_info[i - k].get_helix_position(turn);
            if (pos == HelixPosition::Start || pos == HelixPosition::StartAndEnd) {
              is_turn = true;
            }
          }
        }
      }

      if (is_turn) {
        ss_info[i].ss_type = SecondaryStructure::Turn;
      }
    }
  }
}

void DsspCalculator::find_bends_and_breaks() {
  // Find chain breaks
  for (size_t i = 0; i + 1 < res_infos.size(); ++i) {
    bool has_break = false;

    const Atom* c_atom = res_infos[i]->res->get_c();
    const Atom* n_atom = res_infos[i + 1]->res->get_n();

    if (c_atom && n_atom) {
      double dist = calculate_atomic_distance(c_atom, n_atom);

      if (dist > options.max_peptide_bond_distance) {
        has_break = true;
      }
    } else {
      has_break = true;
    }

    if (has_break) {
      ss_info[i].has_break = true;
    }
  }

  // Find bends
  for (size_t i = 1; i + 3 < res_infos.size(); ++i) {
    if (no_chain_breaks_between(i - 1, i + 4)) {
        const Atom* ca_prev_atom = res_infos[i - 1]->res->get_ca();
        const Atom* ca_curr_atom = res_infos[i + 1]->res->get_ca();
        const Atom* ca_next_atom = res_infos[i + 3]->res->get_ca();

        if (ca_prev_atom && ca_curr_atom && ca_next_atom) {
          Position ca_prev = ca_prev_atom->pos;
          Position ca_curr = ca_curr_atom->pos;
          Position ca_next = ca_next_atom->pos;

          Vec3 v1 = ca_curr - ca_prev;
          Vec3 v2 = ca_next - ca_curr;

          double angle = deg(std::acos(v1.dot(v2) / (v1.length() * v2.length())));

          if (angle > options.bend_angle_min && angle < 360.0) {
            ss_info[i+1].ss_type = SecondaryStructure::Bend;
          }
        }
    }
  }
}

void DsspCalculator::find_polyproline_helices() {
  // TODO: Implement polyproline helix detection based on phi/psi angles
  // This requires dihedral angle calculations
}

std::string DsspCalculator::generate_ss_string() const {
  std::string result;
  result.reserve(ss_info.size());

  for (size_t i = 0; i < ss_info.size(); ++i) {
    result += static_cast<char>(ss_info[i].ss_type);
    if (i + 1 < ss_info.size() && ss_info[i].has_break)
        result += static_cast<char>(SecondaryStructure::Break);
  }

  return result;
}

// Convenience function
std::string calculate_dssp(NeighborSearch& ns, Topo::ChainInfo& cinfo, const DsspOptions& opts) {
  DsspCalculator calculator(opts);
  // Create a temporary topology with just this chain
  Topo topo;
  topo.chain_infos.push_back(cinfo);
  return calculator.calculate_secondary_structure(ns, topo);
}

} // namespace gemmi

