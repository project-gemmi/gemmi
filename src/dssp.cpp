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
constexpr double RAD_TO_DEG = 180.0 / M_PI;

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
        if (j != i + 1) {
          calculate_hbond_energy(&residue_info[j], &residue_info[i]);
        }
      } else {
        calculate_hbond_geometry(&residue_info[i], &residue_info[j]);
        if (j != i + 1) {
          calculate_hbond_geometry(&residue_info[j], &residue_info[i]);
        }
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
      Vec3 prev_co = (prev_c - prev_o).normalized();
      h_pos = donor_n->pos + Position(prev_co);
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
  double energy;
  if (dist_NO < MIN_ATOM_DISTANCE || dist_HC < MIN_ATOM_DISTANCE ||

      dist_HO < MIN_ATOM_DISTANCE || dist_NC < MIN_ATOM_DISTANCE) {
    energy = MIN_ENERGY;
  } else {
    energy = COUPLING_CONSTANT * ((1.0/dist_NO) + (1.0/dist_HC) - (1.0/dist_HO) - (1.0/dist_NC));
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
      Vec3 prev_co = (prev_c - prev_o).normalized();
      h_pos = donor_n->pos + Position(prev_co);
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

  bool has_bond = false;

  if (donor.acceptors[0] == acceptor.res_info) {
    has_bond = (options.hbond_definition == HBondDefinition::Geometry) ||
               (donor.acceptor_energies[0] < options.hbond_energy_cutoff);
  } else if (donor.acceptors[1] == acceptor.res_info) {
    has_bond = (options.hbond_definition == HBondDefinition::Geometry) ||
               (donor.acceptor_energies[1] < options.hbond_energy_cutoff);
  }

  return has_bond;
}

bool DsspCalculator::no_chain_breaks_between(size_t res1_idx, size_t res2_idx) const {
  size_t start = std::min(res1_idx, res2_idx);
  size_t end = std::max(res1_idx, res2_idx);

  for (size_t i = start; i < end; ++i)
    if (ss_info[i].has_break)
      return false;
  return true;
}

BridgeType DsspCalculator::calculate_bridge_type(size_t res1_idx, size_t res2_idx) const {
  if (res1_idx == 0 || res2_idx == 0 ||
      res1_idx + 1 >= residue_info.size() || res2_idx + 1 >= residue_info.size()) {
    return BridgeType::None;
  }

  if (!no_chain_breaks_between(res1_idx - 1, res1_idx + 1) ||
      !no_chain_breaks_between(res2_idx - 1, res2_idx + 1)) {
    return BridgeType::None;
  }

  // Check parallel bridge patterns
  if ((has_hbond_between(res1_idx + 1, res2_idx) && has_hbond_between(res2_idx, res1_idx - 1)) ||
      (has_hbond_between(res2_idx + 1, res1_idx) && has_hbond_between(res1_idx, res2_idx - 1))) {
    return BridgeType::Parallel;
  }

  // Check antiparallel bridge patterns
  if ((has_hbond_between(res1_idx + 1, res2_idx - 1) && has_hbond_between(res2_idx + 1, res1_idx - 1)) ||
      (has_hbond_between(res2_idx, res1_idx) && has_hbond_between(res1_idx, res2_idx))) {
    return BridgeType::AntiParallel;
  }

  return BridgeType::None;
}

void DsspCalculator::find_bridges_and_strands() {
  // Find bridges
  for (size_t i = 1; i + 4 < ss_info.size(); ++i) {
    for (size_t j = i + 3; j + 1 < ss_info.size(); ++j) {
      BridgeType bridge_type = calculate_bridge_type(i, j);
      if (bridge_type != BridgeType::None) {
        ss_info[i].add_bridge(j, bridge_type);
        ss_info[j].add_bridge(i, bridge_type);
      }
    }
  }

  // Find strands from bridge patterns
  for (size_t i = 1; i + 1 < ss_info.size(); ++i) {
    for (size_t j = 1; j < 3 && i + j < ss_info.size(); ++j) {
      for (BridgeType bridge_type : {BridgeType::Parallel, BridgeType::AntiParallel}) {
        if (ss_info[i].has_bridges(bridge_type) && ss_info[i + j].has_bridges(bridge_type) &&
            no_chain_breaks_between(i - 1, i + 1) && no_chain_breaks_between(i + j - 1, i + j + 1)) {

          const auto& i_partners = (bridge_type == BridgeType::Parallel) ?
                                   ss_info[i].parallel_bridges : ss_info[i].antiparallel_bridges;
          const auto& j_partners = (bridge_type == BridgeType::Parallel) ?
                                   ss_info[i + j].parallel_bridges : ss_info[i + j].antiparallel_bridges;

          for (size_t i_partner : i_partners) {
            for (size_t j_partner : j_partners) {
              if (std::abs(static_cast<int>(i_partner) - static_cast<int>(j_partner)) < 6) {
                size_t strand_start = std::min(i_partner, j_partner);
                size_t strand_end = std::max(i_partner, j_partner);

                // Mark second strand as extended
                for (size_t k = strand_start; k <= strand_end; ++k) {
                  ss_info[k].ss_type = SecondaryStructure::Strand;
                }

                // Mark first strand as extended
                for (size_t k = 0; k <= j; ++k) {
                  ss_info[i + k].ss_type = SecondaryStructure::Strand;
                }
              }
            }
          }
        }
      }
    }
  }

  // Mark isolated bridges
  for (size_t i = 1; i + 1 < ss_info.size(); ++i) {
    if (ss_info[i].ss_type != SecondaryStructure::Strand &&
        (ss_info[i].has_bridges(BridgeType::Parallel) || ss_info[i].has_bridges(BridgeType::AntiParallel))) {
      ss_info[i].ss_type = SecondaryStructure::Bridge;
    }
  }
}

void DsspCalculator::find_turns_and_helices() {
  // Find helix patterns for different turn types
  for (TurnType turn_type : {TurnType::Turn_3, TurnType::Turn_4, TurnType::Turn_5}) {
    size_t stride = static_cast<size_t>(turn_type);

    // Find individual turns
    for (size_t i = 0; i + stride < ss_info.size(); ++i) {
      if (has_hbond_between(i + stride, i) && no_chain_breaks_between(i, i + stride)) {
        ss_info[i + stride].set_helix_position(turn_type, HelixPosition::End);

        for (size_t k = 1; k < stride; ++k) {
          if (ss_info[i + k].get_helix_position(turn_type) == HelixPosition::None) {
            ss_info[i + k].set_helix_position(turn_type, HelixPosition::Middle);
          }
        }

        if (ss_info[i].get_helix_position(turn_type) == HelixPosition::End) {
          ss_info[i].set_helix_position(turn_type, HelixPosition::StartAndEnd);
        } else {
          ss_info[i].set_helix_position(turn_type, HelixPosition::Start);
        }
      }
    }
  }

  // Assign helix secondary structure
  for (TurnType turn_type : {TurnType::Turn_4, TurnType::Turn_3, TurnType::Turn_5}) {
    size_t stride = static_cast<size_t>(turn_type);

    for (size_t i = 1; i + stride < ss_info.size(); ++i) {
      if ((ss_info[i - 1].get_helix_position(turn_type) == HelixPosition::Start ||
           ss_info[i - 1].get_helix_position(turn_type) == HelixPosition::StartAndEnd) &&
          (ss_info[i].get_helix_position(turn_type) == HelixPosition::Start ||
           ss_info[i].get_helix_position(turn_type) == HelixPosition::StartAndEnd)) {

        SecondaryStructure helix_type;
        bool can_assign = true;

        switch (turn_type) {
          case TurnType::Turn_3:
            helix_type = SecondaryStructure::Helix_3;
            for (size_t k = 0; k < stride && can_assign; ++k) {
              can_assign = (ss_info[i + k].ss_type <= SecondaryStructure::Helix_3);
            }
            break;
          case TurnType::Turn_5:
            helix_type = SecondaryStructure::Helix_5;
            for (size_t k = 0; k < stride && can_assign; ++k) {
              can_assign = (ss_info[i + k].ss_type <= SecondaryStructure::Helix_5) ||
                          (options.pi_helix_preference && ss_info[i + k].ss_type == SecondaryStructure::Helix_4);
            }
            break;
          default:
            helix_type = SecondaryStructure::Helix_4;
            break;
        }

        if (can_assign || helix_type == SecondaryStructure::Helix_4) {
          for (size_t k = 0; k < stride; ++k) {
            ss_info[i + k].ss_type = helix_type;
          }
        }
      }
    }
  }

  // Assign turns
  for (size_t i = 1; i + 1 < ss_info.size(); ++i) {
    if (ss_info[i].ss_type <= SecondaryStructure::Turn) {
      bool is_turn = false;

      for (TurnType turn_type : {TurnType::Turn_3, TurnType::Turn_4, TurnType::Turn_5}) {
        size_t stride = static_cast<size_t>(turn_type);

        for (size_t k = 1; k < stride && !is_turn; ++k) {
          if (i >= k) {
            is_turn = (ss_info[i - k].get_helix_position(turn_type) == HelixPosition::Start ||
                      ss_info[i - k].get_helix_position(turn_type) == HelixPosition::StartAndEnd);
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
      ss_info[i + 1].has_break = true;
      ss_info[i].break_partners.push_back(&ss_info[i + 1]);
      ss_info[i + 1].break_partners.push_back(&ss_info[i]);
    }
  }

  // Find bends
  for (size_t i = 2; i + 2 < residue_info.size(); ++i) {
    if (ss_info[i - 2].has_break || ss_info[i - 1].has_break ||
        ss_info[i].has_break || ss_info[i + 1].has_break) {
      continue;
    }

    if (residue_info[i - 2].has_atom(ResidueInfo::CA) &&
        residue_info[i].has_atom(ResidueInfo::CA) &&
        residue_info[i + 2].has_atom(ResidueInfo::CA)) {

      Position ca_prev = residue_info[i - 2].get_backbone_atom(ResidueInfo::CA)->pos;
      Position ca_curr = residue_info[i].get_backbone_atom(ResidueInfo::CA)->pos;
      Position ca_next = residue_info[i + 2].get_backbone_atom(ResidueInfo::CA)->pos;

      Vec3 v1 = ca_curr - ca_prev;
      Vec3 v2 = ca_next - ca_curr;

      double angle = std::acos(v1.dot(v2) / (v1.length() * v2.length())) * RAD_TO_DEG;

      if (angle > options.bend_angle_min && angle < 360.0) {
        ss_info[i].ss_type = SecondaryStructure::Bend;
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

  for (const auto& info : ss_info) {
    result += static_cast<char>(info.ss_type);
  }

  // Insert break symbols
  for (size_t i = 0; i < ss_info.size(); ++i) {
    if (ss_info[i].has_break) {
      result.insert(result.begin() + i + 1, static_cast<char>(SecondaryStructure::Break));
    }
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

