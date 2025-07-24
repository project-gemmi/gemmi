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

double calculate_atomic_distance(Atom* atom1, Atom* atom2) {
  if (!atom1 || !atom2) return 0.0;
  return atom1->pos.dist(atom2->pos);
}

double calculate_dihedral_angle(Atom* a1, Atom* a2, Atom* a3, Atom* a4) {
  if (!a1 || !a2 || !a3 || !a4) return 0.0;
  return calculate_dihedral(a1->pos, a2->pos, a3->pos, a4->pos);
}


} // anonymous namespace

// DsspCalculator implementation
std::string DsspCalculator::calculate_secondary_structure(NeighborSearch& ns, Topo::ChainInfo& cinfo) {
  // Clear previous calculations
  residue_info.clear();
  ss_info.clear();
  bridges_.clear();

  // Setup residue information
  setup_residue_info(cinfo);

  if (residue_info.empty()) {
    return "";
  }

  // Initialize secondary structure info
  ss_info.resize(residue_info.size());

  // Calculate hydrogen bonds
  calculate_hydrogen_bonds(ns);

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

void DsspCalculator::setup_residue_info(Topo::ChainInfo& cinfo) {
  residue_info.reserve(cinfo.res_infos.size());

  for (auto& res_info : cinfo.res_infos) {
    if (!res_info.res) continue;

    ResidueInfo rinfo;
    rinfo.res_info = &res_info;

    // Check if this is proline
    if (res_info.res->name == "PRO") {
      rinfo.is_proline = true;
    }

    // Find backbone atoms
    for (auto& atom : res_info.res->atoms) {
      if (atom.name == "CA") {
        rinfo.set_backbone_atom(ResidueInfo::CA, &atom);
      } else if (atom.name == "C") {
        rinfo.set_backbone_atom(ResidueInfo::C, &atom);
      } else if (atom.name == "O") {
        rinfo.set_backbone_atom(ResidueInfo::O, &atom);
      } else if (atom.name == "N") {
        rinfo.set_backbone_atom(ResidueInfo::N, &atom);
        if (options.hydrogen_mode == HydrogenMode::Calculate) {
          // In calculate mode, use N position for H calculation
          rinfo.set_backbone_atom(ResidueInfo::H, &atom);
        }
      } else if (atom.name == "H" && options.hydrogen_mode == HydrogenMode::Existing) {
        rinfo.set_backbone_atom(ResidueInfo::H, &atom);
      }
    }

    residue_info.push_back(rinfo);
  }

  // Remove incomplete residues if requested
  if (options.clear_defective_residues) {
    auto is_incomplete = [](const ResidueInfo& rinfo) {
      return !rinfo.has_atom(ResidueInfo::CA) ||
             !rinfo.has_atom(ResidueInfo::C) ||
             !rinfo.has_atom(ResidueInfo::O) ||
             !rinfo.has_atom(ResidueInfo::N) ||
             !rinfo.has_atom(ResidueInfo::H);
    };

    residue_info.erase(std::remove_if(residue_info.begin(), residue_info.end(), is_incomplete),
                      residue_info.end());
  }

  // Setup prev/next pointers
  for (size_t i = 0; i < residue_info.size(); ++i) {
    if (i > 0) {
      residue_info[i].prev_residue = &residue_info[i-1];
    }
    if (i < residue_info.size() - 1) {
      residue_info[i].next_residue = &residue_info[i+1];
    }
  }
}

void DsspCalculator::calculate_hydrogen_bonds(NeighborSearch& ) {
  if (options.use_neighbor_search) {
    // Use neighbor search for efficiency
    std::vector<Position> ca_positions;
    ca_positions.reserve(residue_info.size());

    for (const auto& rinfo : residue_info) {
      if (rinfo.has_atom(ResidueInfo::CA)) {
        ca_positions.push_back(rinfo.get_backbone_atom(ResidueInfo::CA)->pos);
      }
    }

    // TODO: Implement neighbor search based hydrogen bond calculation
    // For now, fall back to all-vs-all
  }

  // All-vs-all hydrogen bond calculation
  for (size_t i = 0; i < residue_info.size(); ++i) {
    for (size_t j = i + 1; j < residue_info.size(); ++j) {
      if (options.hbond_definition == HBondDefinition::Energy) {
        calculate_hbond_energy(&residue_info[i], &residue_info[j]);
        calculate_hbond_energy(&residue_info[j], &residue_info[i]);
      } else {
        calculate_hbond_geometry(&residue_info[i], &residue_info[j]);
        calculate_hbond_geometry(&residue_info[j], &residue_info[i]);
      }
    }
  }
}

static Position calculate_h_pos(ResidueInfo* donor, const Atom* donor_n, HydrogenMode h_mode) {

  Position h_pos= donor->get_backbone_atom(ResidueInfo::H)->pos;
  if (h_mode == HydrogenMode::Calculate) {
    // Calculate hydrogen position from previous residue C-O vector
    if (donor->prev_residue &&
        donor->prev_residue->has_atom(ResidueInfo::C) &&
        donor->prev_residue->has_atom(ResidueInfo::O)) {
      // official DSSP implementation assumes that N-H vector is a normalized
      // O-C of the previous residue
      Vec3 prev_o = donor->prev_residue->get_backbone_atom(ResidueInfo::O)->pos;
      Vec3 prev_c = donor->prev_residue->get_backbone_atom(ResidueInfo::C)->pos;
      Vec3 prev_co = (prev_o - prev_c).normalized();
      h_pos = donor_n->pos - Position(prev_co);
    } else {
      // no previous residue, leaving H set to N
    }
  }
  return h_pos;
}

void DsspCalculator::calculate_hbond_energy(ResidueInfo* donor, ResidueInfo* acceptor) {
  if (donor->is_proline ||
      !acceptor->has_atom(ResidueInfo::C) || !acceptor->has_atom(ResidueInfo::O) ||
      !donor->has_atom(ResidueInfo::N) || !donor->has_atom(ResidueInfo::H)) {
    return;
  }

  // Check CA distance first for efficiency
  if (donor->has_atom(ResidueInfo::CA) && acceptor->has_atom(ResidueInfo::CA)) {
    double ca_dist = calculate_atomic_distance(
        donor->get_backbone_atom(ResidueInfo::CA),
        acceptor->get_backbone_atom(ResidueInfo::CA));

   if (ca_dist >= options.min_ca_distance) {
      return;
    }
  }

  // Calculate distances
  const Atom* donor_n = donor->get_backbone_atom(ResidueInfo::N);
  const Atom* acceptor_o = acceptor->get_backbone_atom(ResidueInfo::O);
  const Atom* acceptor_c = acceptor->get_backbone_atom(ResidueInfo::C);
  double dist_NO = donor_n->pos.dist(acceptor_o->pos);
  double dist_NC = donor_n->pos.dist(acceptor_c->pos);
  Position h_pos = calculate_h_pos(donor, donor_n, options.hydrogen_mode);
  double dist_HO = h_pos.dist(acceptor->get_backbone_atom(ResidueInfo::O)->pos);
  double dist_HC = h_pos.dist(acceptor->get_backbone_atom(ResidueInfo::C)->pos);

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
    donor->acceptors[0] = acceptor->res_info;
    donor->acceptor_energies[0] = energy;
  } else if (energy < donor->acceptor_energies[1]) {
    donor->acceptors[1] = acceptor->res_info;
    donor->acceptor_energies[1] = energy;
  }

  if (energy < acceptor->donor_energies[0]) {
    acceptor->donors[1] = acceptor->donors[0];
    acceptor->donor_energies[1] = acceptor->donor_energies[0];
    acceptor->donors[0] = donor->res_info;
    acceptor->donor_energies[0] = energy;
  } else if (energy < acceptor->donor_energies[1]) {
    acceptor->donors[1] = donor->res_info;
    acceptor->donor_energies[1] = energy;
  }
}

void DsspCalculator::calculate_hbond_geometry(ResidueInfo* donor, ResidueInfo* acceptor) {
  if (donor->is_proline ||
      !acceptor->has_atom(ResidueInfo::C) || !acceptor->has_atom(ResidueInfo::O) ||
      !donor->has_atom(ResidueInfo::N) || !donor->has_atom(ResidueInfo::H)) {
    return;
  }

  // Check distance criterion
  const Atom* donor_n = donor->get_backbone_atom(ResidueInfo::N);
  Position n_pos = donor_n->pos;
  Position o_pos = acceptor->get_backbone_atom(ResidueInfo::O)->pos;
  double dist_NO = (n_pos - o_pos).length();

  if (dist_NO > MAX_HBOND_DISTANCE)
    return;
  Position h_pos = calculate_h_pos(donor, donor_n, options.hydrogen_mode);
  /*
  Position h_pos = donor->get_backbone_atom(ResidueInfo::H)->pos;
  if (options.hydrogen_mode == HydrogenMode::Calculate) {
    // Calculate hydrogen position from previous residue C-O vector
    if (donor->prev_residue &&
        donor->prev_residue->has_atom(ResidueInfo::C) &&
        donor->prev_residue->has_atom(ResidueInfo::O)) {
      // official DSSP implementation assumes that N-H vector is a normalized
      // O-C of the previous residue
      Vec3 prev_o = donor->prev_residue->get_backbone_atom(ResidueInfo::O)->pos;
      Vec3 prev_c = donor->prev_residue->get_backbone_atom(ResidueInfo::C)->pos;
      Vec3 prev_co = (prev_o - prev_c).normalized();
      h_pos = donor_n->pos - Position(prev_co);
    } else {
      // no previous residue, leaving H set to N
    }
  }
  */

  // Check angle criterion
  double angle = (o_pos - n_pos).angle(h_pos - n_pos);
  if (angle <= MAX_HBOND_ANGLE) {
    // Store hydrogen bond (geometry-based uses nullptr for energies)
    if (donor->acceptors[0] == nullptr) {
      donor->acceptors[0] = acceptor->res_info;
      acceptor->donors[0] = donor->res_info;
    } else if (donor->acceptors[1] == nullptr) {
      donor->acceptors[1] = donor->acceptors[0];
      donor->acceptors[0] = acceptor->res_info;
      acceptor->donors[1] = acceptor->donors[0];
      acceptor->donors[0] = donor->res_info;
    }
  }
}

bool DsspCalculator::has_hbond_between(size_t donor_idx, size_t acceptor_idx) const {
    const ResidueInfo& donor = residue_info[donor_idx];
    const ResidueInfo& acceptor = residue_info[acceptor_idx];

    for (int i = 0; i < 2; ++i) {
        if (donor.acceptors[i] == acceptor.res_info && donor.acceptor_energies[i] < options.hbond_energy_cutoff) {
            return true;
        }
    }

    for (int i = 0; i < 2; ++i) {
        if (acceptor.donors[i] == donor.res_info && acceptor.donor_energies[i] < options.hbond_energy_cutoff) {
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
  // antiparallel
  bool anti1 = has_hbond_between(i, j) && has_hbond_between(j, i);
  bool anti2 = has_hbond_between(i - 1, j + 1) && has_hbond_between(j - 1, i + 1);
  if (anti1 && anti2)
    return BridgeType::AntiParallel;
  // parallel
  bool para1 = has_hbond_between(i, j + 1) && has_hbond_between(j, i + 1);
  bool para2 = has_hbond_between(i - 1, j) && has_hbond_between(j - 1, i);
  if (para1 && para2)
    return BridgeType::Parallel;

  return BridgeType::None;
}

void DsspCalculator::find_bridges_and_strands() {
  // Find bridges
  for (size_t i = 1; i + 1 < residue_info.size(); ++i) {
    for (size_t j = i + 1; j < residue_info.size(); ++j) {
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

    for (size_t j = 0; j + stride < residue_info.size(); ++j) {
      if (has_hbond_between(j + stride, j) && no_chain_breaks_between(j, j + stride)) {
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
  for (size_t i = 0; i + 4 < residue_info.size(); ++i) {
    // Look for overlapping 4-turn patterns: i→i+4 AND i+1→i+5
    if (has_hbond_between(i + 4, i) && has_hbond_between(i + 5, i + 1) &&
        no_chain_breaks_between(i, i + 5)) {

      // Found overlapping alpha-helix pattern - mark as helix
      for (size_t k = i + 1; k <= i + 4; ++k) {
        ss_info[k].ss_type = SecondaryStructure::Helix_4;
      }

      // Continue the helix while overlapping patterns exist
      for (size_t j = i + 2; j + 4 < residue_info.size(); ++j) {
        if (has_hbond_between(j + 4, j) && no_chain_breaks_between(j, j + 4)) {
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
  for (size_t i = 0; i + 3 < residue_info.size(); ++i) {
    // Look for overlapping 3-turn patterns: i→i+3 AND i+1→i+4
    if (has_hbond_between(i + 3, i) && has_hbond_between(i + 4, i + 1) &&
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
        for (size_t j = i + 2; j + 3 < residue_info.size(); ++j) {
          if (has_hbond_between(j + 3, j) && no_chain_breaks_between(j, j + 3)) {
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
  for (size_t i = 0; i + 5 < residue_info.size(); ++i) {
    // Look for overlapping 5-turn patterns: i→i+5 AND i+1→i+6
    if (has_hbond_between(i + 5, i) && has_hbond_between(i + 6, i + 1) &&
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
        for (size_t j = i + 2; j + 5 < residue_info.size(); ++j) {
          if (has_hbond_between(j + 5, j) && no_chain_breaks_between(j, j + 5)) {
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
  for (size_t i = 1; i + 1 < residue_info.size(); ++i) {
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
  for (size_t i = 0; i + 1 < residue_info.size(); ++i) {
    bool has_break = false;

    if (residue_info[i].has_atom(ResidueInfo::C) &&
        residue_info[i + 1].has_atom(ResidueInfo::N)) {
      double dist = calculate_atomic_distance(
          residue_info[i].get_backbone_atom(ResidueInfo::C),
          residue_info[i + 1].get_backbone_atom(ResidueInfo::N));

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
  for (size_t i = 1; i + 3 < residue_info.size(); ++i) {
    if (no_chain_breaks_between(i - 1, i + 4)) {
        if (residue_info[i - 1].has_atom(ResidueInfo::CA) &&
            residue_info[i + 1].has_atom(ResidueInfo::CA) &&
            residue_info[i + 3].has_atom(ResidueInfo::CA)) {

          Position ca_prev = residue_info[i - 1].get_backbone_atom(ResidueInfo::CA)->pos;
          Position ca_curr = residue_info[i + 1].get_backbone_atom(ResidueInfo::CA)->pos;
          Position ca_next = residue_info[i + 3].get_backbone_atom(ResidueInfo::CA)->pos;

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
  return calculator.calculate_secondary_structure(ns, cinfo);
}

// Legacy function for backwards compatibility
std::vector<HBond> dssp_determine_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& cinfo) {
  std::vector<HBond> hbonds;

  DsspCalculator calculator;
  calculator.calculate_secondary_structure(ns, cinfo);

  // Extract hydrogen bonds from the calculation
  // This is a simplified implementation for backward compatibility
  return hbonds;
}

} // namespace gemmi

