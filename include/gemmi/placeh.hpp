// Copyright 2018 Global Phasing Ltd.
//
// Place hydrogens according to bond lengths and angles from monomer library.

#ifndef GEMMI_PLACEH_HPP_
#define GEMMI_PLACEH_HPP_

#include <cmath>         // for sqrt, sin, cos
#include <memory>        // for unique_ptr
#include "model.hpp"     // for Atom
#include "topo.hpp"      // for Topo
#include "chemcomp.hpp"  // for ChemComp
#include "calculate.hpp" // for calculate_angle
#include "modify.hpp"    // for remove_hydrogens

namespace gemmi {

// Assumes no hydrogens in the residue.
// Position and serial number are not assigned for new atoms.
inline void add_hydrogens_without_positions(Topo::ResInfo& ri) {
  Residue& res = *ri.res;
  // Add H atom for each conformation (altloc) of the parent atom.
  for (size_t i = 0, size = res.atoms.size(); i != size; ++i) {
    const ChemComp& cc = ri.get_final_chemcomp(res.atoms[i].altloc);
    for (const Restraints::Bond& bond : cc.rt.bonds) {
      // res.atoms may get re-allocated, so we can't set parent earlier
      const Atom& parent = res.atoms[i];
      assert(!parent.is_hydrogen());
      const Restraints::AtomId* atom_id;
      if (bond.id1 == parent.name)
        atom_id = &bond.id2;
      else if (bond.id2 == parent.name)
        atom_id = &bond.id1;
      else
        continue;
      auto it = cc.find_atom(atom_id->atom);
      if (it == cc.atoms.end())
        fail("inconsistent _chem_comp " + cc.name);
      if (it->is_hydrogen()) {
        gemmi::Atom atom;
        atom.name = it->id;
        atom.altloc = parent.altloc;
        atom.element = it->el;
        // calc_flag will be changed to Calculated when the position is set
        atom.calc_flag = CalcFlag::Dummy;
        atom.occ = parent.occ;
        atom.b_iso = parent.b_iso;
        res.atoms.push_back(atom);
      }
    }
  }
}


// Calculate position using one angle (theta) and one dihedral angle (tau).
// Returns position of x4 in x1-x2-x3-x4, where dist=|x3-x4| and
// theta is angle(x2, x3, x4).
// Based on section 3.3 of Paciorek et al, Acta Cryst. A52, 349 (1996).
inline Position position_from_angle_and_torsion(const Position& x1,
                                                const Position& x2,
                                                const Position& x3,
                                                double dist,  // |x3-x4|
                                                double theta, // angle x2-x3-x4
                                                double tau) { // dihedral angle
  using std::sin;
  using std::cos;
  Vec3 u = x2 - x1;
  Vec3 v = x3 - x2;
  Vec3 e1 = v.normalized();
  double delta = u.dot(e1);
  Vec3 e2 = -(u - delta * e1).normalized();
  Vec3 e3 = e1.cross(e2);
  return x3 + Position(dist * (-cos(theta) * e1 +
                               sin(theta) * (cos(tau) * e2 + sin(tau) * e3)));
}

// Rodrigues' rotation formula, rotate vector v given axis of rotation
// (which must be a unit vector) and angle (in radians).
inline
Vec3 rotate_about_axis(const Vec3& v, const Vec3& axis, double theta) {
  double sin_theta = std::sin(theta);
  double cos_theta = std::cos(theta);
  return v * cos_theta + axis.cross(v) * sin_theta +
         axis * (axis.dot(v) * (1 - cos_theta));
}

inline Vec3 get_vector_to_line(const Position& point,
                               const Position& point_on_the_line,
                               const Vec3& unit_vector) {
  // en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
  // the component of a - p perpendicular to the line is: (a-p) - ((a-p).n)n
  Vec3 ap = point_on_the_line - point;
  return ap - ap.dot(unit_vector) * unit_vector;
}

// If no points satisfy the distances returns a pair of NaNs
inline
std::pair<Position, Position> trilaterate(const Position& p1, double r1sq,
                                          const Position& p2, double r2sq,
                                          const Position& p3, double r3sq) {
  // It was based on https://en.wikipedia.org/wiki/Trilateration
  // but apparently that page has changed in the meantime.
  Vec3 ex = (p2 - p1).normalized();
  double i = ex.dot(p3-p1);
  Vec3 ey = (Vec3(p3) - p1 - i*ex).normalized();
  Vec3 ez = ex.cross(ey);
  double d = (p2-p1).length();
  double j = ey.dot(p3-p1);
  double x = (r1sq - r2sq + d*d) / (2*d);
  double y = (r1sq - r3sq + i*i + j*j) / (2*j) - x*i/j;
  double z2 = r1sq - x*x - y*y;
  double z = std::sqrt(z2);  // may result in NaN
  return std::make_pair(p1 + Position(x*ex + y*ey + z*ez),
                        p1 + Position(x*ex + y*ey - z*ez));
}

// Calculate position using two angles.
// Returns p4. Topology: p1 is bonded to p2, p3 and p4.
inline std::pair<Position, Position>
position_from_two_angles(const Position& p1,
                         const Position& p2,
                         const Position& p3,
                         double dist14,     // |p4-p1|
                         double theta214,   // angle p2-p1-p4
                         double theta314) { // angle p3-p1-p4
  double d12sq = p1.dist_sq(p2);
  double d13sq = p1.dist_sq(p3);
  double d14sq = dist14 * dist14;
  // the law of cosines
  double d24sq = d14sq + d12sq - 2 * std::sqrt(d14sq * d12sq) * cos(theta214);
  double d34sq = d14sq + d13sq - 2 * std::sqrt(d14sq * d13sq) * cos(theta314);
  auto t = trilaterate(p1, d14sq, p2, d24sq, p3, d34sq);
  return t;
}


inline void place_hydrogens(const Topo& topo, const Atom& atom) {
  using Angle = Restraints::Angle;
  struct BondedAtom {
    Atom* ptr;
    Position& pos; // == ptr->pos;
    double dist;
  };

  // put atoms bonded to atom into two lists
  std::vector<BondedAtom> known; // heavy atoms with known positions
  std::vector<BondedAtom> hs;    // H atoms (unknown)
  known.reserve(3);
  hs.reserve(4);

  auto range = topo.bond_index.equal_range(&atom);
  for (auto i = range.first; i != range.second; ++i) {
    const Topo::Bond* t = i->second;
    Atom* other = t->atoms[t->atoms[0] == &atom ? 1 : 0];
    if (other->altloc) {
      if (atom.altloc && atom.altloc != other->altloc)
        continue;
      if (atom.altloc == '\0' &&
          in_vector_f([&](const BondedAtom& a) { return a.ptr->name == other->name; },
                      known))
        continue;
    }
    auto& atom_list = other->is_hydrogen() ? hs : known;
    atom_list.push_back({other, other->pos, t->restr->value});
  }

  if (hs.size() == 0)
    return;

  auto giveup = [&](const std::string& message) {
    for (BondedAtom& bonded_h : hs)
      bonded_h.ptr->occ = 0;
    fail(message);
  };

  for (BondedAtom& bonded_h : hs)
    bonded_h.ptr->calc_flag = CalcFlag::Calculated;

  // ==== only hydrogens ====
  if (known.size() == 0) {
    // we can only arbitrarily pick directions of atoms
    for (BondedAtom& bonded_h : hs)
      bonded_h.ptr->occ = 0;
    hs[0].pos = atom.pos + Position(hs[0].dist, 0, 0);
    if (hs.size() > 1) {
      double theta = pi();
      if (const Angle* ang = topo.take_angle(hs[1].ptr, &atom, hs[0].ptr))
        theta = ang->radians();
      hs[1].pos = atom.pos + Position(hs[1].dist * cos(theta),
                                      hs[1].dist * sin(theta), 0);
    }
    if (hs.size() > 2) {
      if (hs.size() == 3) {
        // for now only NH3 (NH2.cif and NH3.cif) has such configuration,
        // so we are cheating here a little.
        double y = 2 * atom.pos.y - hs[1].pos.y;
        hs[2].pos = Position(hs[1].pos.x, y, hs[1].pos.z);
      } else if (hs.size() == 4) {
        // similarly, only CH4 (CH2.cif) and NH4 (NH4.cif) are handled here
        const Angle* ang1 = topo.take_angle(hs[2].ptr, &atom, hs[0].ptr);
        const Angle* ang2 = topo.take_angle(hs[2].ptr, &atom, hs[1].ptr);
        double theta1 = rad(ang1 ? ang1->value : 109.47122);
        double theta2 = rad(ang2 ? ang2->value : 109.47122);
        auto pos = position_from_two_angles(atom.pos, hs[0].pos, hs[1].pos,
                                            hs[2].dist, theta1, theta2);
        hs[2].pos = pos.first;
        hs[3].pos = pos.second;
      }
    }

  // ==== one heavy atom and hydrogens ====
  } else if (known.size() == 1) {
    const BondedAtom& h = hs[0];
    const BondedAtom& heavy = known[0];
    const Restraints::Angle* angle = topo.take_angle(h.ptr, &atom, heavy.ptr);
    if (!angle)
      giveup("No angle restraint for " + h.ptr->name + ".\n");
    if (std::abs(angle->value - 180.0) < 0.5) {
      Vec3 u = atom.pos - h.pos;
      h.pos = atom.pos + Position(u * (h.dist / u.length()));
      if (hs.size() > 1)
        giveup("Unusual: one of two H atoms has angle restraint 180 deg.");
      return;
    }
    double theta = angle->radians();
    double tau = 0.0;
    int period = 0;
    const Atom* tau_end = nullptr;
    auto plane_range = topo.plane_index.equal_range(&atom);
    for (auto i = plane_range.first; i != plane_range.second; ++i) {
      const Topo::Plane& plane = *i->second;
      // only Topo::Plane with atoms.size() >= 4 is put into planes
      if (plane.has(h.ptr) && plane.has(heavy.ptr)) {
        for (const Atom* a : plane.atoms) {
          if (!a->is_hydrogen() && a != &atom && a != heavy.ptr) {
            tau_end = a;
            break;
          }
        }
        break;
      }
    }
    Atom* torsion_h = nullptr;
    if (!tau_end) {
      // Using one dihedral angle.
      // We don't check here for which hydrogen the torsion angle is defined.
      // If an atom has 2 or 3 hydrogens, the torsion angle may not be given
      // for the first one, but only for the 2nd or 3rd (e.g. HD22 in ASN).
      auto tor_range = topo.torsion_index.equal_range(&atom);
      for (auto i = tor_range.first; i != tor_range.second; ++i) {
        const Topo::Torsion& tor = *i->second;
        if (tor.atoms[1] == &atom && tor.atoms[2] == heavy.ptr &&
            tor.atoms[0]->is_hydrogen() && !tor.atoms[3]->is_hydrogen()) {
          tau = rad(tor.restr->value);
          torsion_h = tor.atoms[0];
          tau_end = tor.atoms[3];
          period = tor.restr->period;
          break;
        } else if (tor.atoms[2] == &atom && tor.atoms[1] == heavy.ptr &&
                   tor.atoms[3]->is_hydrogen() && !tor.atoms[0]->is_hydrogen()) {
          tau = rad(tor.restr->value);
          torsion_h = tor.atoms[3];
          tau_end = tor.atoms[0];
          period = tor.restr->period;
          break;
        }
      }
    }
    h.pos = position_from_angle_and_torsion(
        tau_end ? tau_end->pos : Position(0, 0, 0),
        heavy.pos, atom.pos, h.dist, theta, tau);
    if (hs.size() == 2) {
      // I think we can assume the two hydrogens are symmetric.
      Vec3 axis = (heavy.pos - atom.pos).normalized();
      Vec3 perpendicular = get_vector_to_line(hs[0].pos, atom.pos, axis);
      hs[1].pos = hs[0].pos + Position(2 * perpendicular);
    } else if (hs.size() == 3) {
      // Here we assume the three hydrogens are in the same distance from
      // the parent atom and that they make an equilateral triangle.
      // First check for which hydrogen was the torsion restraint (if any).
      int idx = 0;
      for (int i : {1, 2})
        if (torsion_h == hs[i].ptr) {
          idx = i;
          hs[i].pos = h.pos;
        }
      // Now get the other positions by rotation
      Vec3 axis = (heavy.pos - atom.pos).normalized();
      Vec3 v1 = h.pos - atom.pos;
      Vec3 v2 = rotate_about_axis(v1, axis, rad(120));
      Vec3 v3 = rotate_about_axis(v1, axis, rad(-120));
      hs[(idx+1) % 3].pos = atom.pos + Position(v2);
      hs[(idx+2) % 3].pos = atom.pos + Position(v3);
    } else if (hs.size() >= 4) {
      giveup("Unusual: atom bonded to one heavy atoms and 4+ hydrogens.");
    }
    // somewhat arbitrary rule
    if (!tau_end || (hs.size() == 1 && period > 1))
      for (BondedAtom& bonded_h : hs)
        bonded_h.ptr->occ = 0;
  // ==== two heavy atoms and hydrogens ====
  } else {  // known.size() >= 2
    const Angle* ang1 = topo.take_angle(hs[0].ptr, &atom, known[0].ptr);
    const Angle* ang2 = topo.take_angle(hs[0].ptr, &atom, known[1].ptr);
    const Angle* ang3 = topo.take_angle(known[0].ptr, &atom, known[1].ptr);

    if (!ang1 || !ang2)
      giveup(cat("Missing angle restraint ", hs[0].ptr->name, '-', atom.name,
                   '-', known[ang1 ? 1 : 0].ptr->name, ".\n"));
    double theta1 = ang1->radians();
    double theta2 = ang2->radians();
    if (ang3) {
      // If all atoms are in the same plane (sum of angles is 360 degree)
      // the calculations can be simplified.
      double theta3 = ang3->radians();
      // The sum of ideal angles in a plane is not always exactly 360 deg.
      if (theta1 + theta2 + theta3 > rad(360 - 3)) {
        Vec3 v12 = known[0].pos - atom.pos;
        Vec3 v13 = known[1].pos - atom.pos;
        // theta3 is the ideal restraint value, cur_theta3 is the current value
        double cur_theta3 = v12.angle(v13);
        double ratio = (2 * pi() - cur_theta3) / (theta1 + theta2);
        Vec3 axis = v13.cross(v12).normalized();
        Vec3 v14 = rotate_about_axis(v12, axis, theta1 * ratio);
        hs[0].pos = atom.pos + Position(hs[0].dist / v14.length() * v14);
        if (hs.size() > 1) {
          topo.err("Unhandled topology of " + std::to_string(hs.size()) +
                   " hydrogens bonded to " + atom.name);
          for (size_t i = 1; i < hs.size(); ++i)
            hs[i].ptr->occ = 0;
        }
        return;
      }
    }
    auto pos = position_from_two_angles(atom.pos, known[0].pos, known[1].pos,
                                        hs[0].dist, theta1, theta2);
    switch (hs.size()) {
      case 1:
        if (known.size() == 2) {
          const Topo::Chirality* chir = topo.get_chirality(&atom);
          if (chir && chir->restr->sign != ChiralityType::Both) {
            hs[0].pos = chir->check() ? pos.first : pos.second;
          } else {
            hs[0].pos = pos.first;
            hs[0].ptr->occ = 0;
          }
        } else { // known.size() > 2
          const Atom* a3 = known[2].ptr;
          if (const Angle* a = topo.take_angle(a3, &atom, hs[0].ptr)){
            double val1 = calculate_angle(a3->pos, atom.pos, pos.first);
            double val2 = calculate_angle(a3->pos, atom.pos, pos.second);
            double diff1 = angle_abs_diff(val1, a->value);
            double diff2 = angle_abs_diff(val2, a->value);
            hs[0].pos = diff1 < diff2 ? pos.first : pos.second;
          } else {
            hs[0].pos = pos.first;
            hs[0].ptr->occ = 0;
          }
        }
        break;
      case 2:
        hs[0].pos = pos.first;
        hs[1].pos = pos.second;
        break;
      default:
        giveup("Unusual: atom bonded to 2+ heavy atoms and 3+ hydrogens.");
    }
  }
}

inline void adjust_hydrogen_distances(Topo& topo, Restraints::DistanceOf of,
                                      double default_scale=1.) {
  for (const Topo::Bond& t : topo.bonds) {
    assert(t.atoms[0] != nullptr && t.atoms[1] != nullptr);
    if (t.atoms[0]->is_hydrogen() || t.atoms[1]->is_hydrogen()) {
      Position u = t.atoms[1]->pos - t.atoms[0]->pos;
      double scale = t.restr->distance(of) / u.length();
      if (std::isnan(scale))
        scale = default_scale;
      if (t.atoms[1]->is_hydrogen())
        t.atoms[1]->pos = t.atoms[0]->pos + u * scale;
      else
        t.atoms[0]->pos = t.atoms[1]->pos - u * scale;
    }
  }
}

inline void place_hydrogens_on_all_atoms(Topo& topo) {
  for (Topo::ChainInfo& chain_info : topo.chain_infos)
    for (Topo::ResInfo& ri : chain_info.res_infos)
      for (Atom& atom : ri.res->atoms)
        if (!atom.is_hydrogen()) {
          try {
            place_hydrogens(topo, atom);
          } catch (const std::runtime_error& e) {
            topo.err("Placing of hydrogen bonded to "
                     + atom_str(chain_info.chain_ref, *ri.res, atom)
                     + " failed:\n  " + e.what());
          }
        }
}

inline void remove_hydrogens_from_atom(Topo::ResInfo* ri,
                                       const std::string& atom_name, char alt) {
  if (!ri)
    return;
  std::vector<Atom>& atoms = ri->res->atoms;
  const Restraints& rt = ri->get_final_chemcomp(alt).rt;
  for (auto it = atoms.end(); it-- != atoms.begin(); ) {
    if (it->is_hydrogen()) {
      const Restraints::AtomId* heavy = rt.first_bonded_atom(it->name);
      if (heavy && heavy->atom == atom_name && (it->altloc == alt || it->altloc == '\0'))
        atoms.erase(it);
    }
  }
}

enum class HydrogenChange { NoChange, Shift, Remove, ReAdd, ReAddButWater };

inline std::unique_ptr<Topo>
prepare_topology(Structure& st, MonLib& monlib, size_t model_index,
                 HydrogenChange h_change, bool reorder,
                 std::ostream* warnings=nullptr, bool ignore_unknown_links=false) {
  std::unique_ptr<Topo> topo(new Topo);
  topo->warnings = warnings;
  if (model_index >= st.models.size())
    fail("no such model index: " + std::to_string(model_index));
  topo->initialize_refmac_topology(st, st.models[model_index], monlib, ignore_unknown_links);

  bool keep = (h_change == HydrogenChange::NoChange || h_change == HydrogenChange::Shift);
  if (!keep || reorder) {
    // remove/add hydrogens, sort atoms in residues
    for (Topo::ChainInfo& chain_info : topo->chain_infos) {
      for (Topo::ResInfo& ri : chain_info.res_infos) {
        Residue& res = *ri.res;
        if (!keep) {
          remove_hydrogens(res);
          if (h_change == HydrogenChange::ReAdd ||
              (h_change == HydrogenChange::ReAddButWater && !res.is_water())) {
            add_hydrogens_without_positions(ri);
            if (h_change == HydrogenChange::ReAddButWater) {
              // a special handling of HIS for compatibility with Refmac
              if (res.name == "HIS") {
                for (gemmi::Atom& atom : ri.res->atoms)
                  if (atom.name == "HD1" || atom.name == "HE2")
                    atom.occ = 0;
              }
            }
          }
        } else {
          // Special handling of Deuterium - mostly for Refmac.
          // Note: if the model has deuterium, it gets modified.
          if (replace_deuterium_with_fraction(res)) {
            // deuterium names usually differ from the names in dictionary
            for (Atom& atom : res.atoms)
              if (atom.name[0] == 'D' && atom.fraction != 0) {
                const ChemComp& cc = ri.get_final_chemcomp(atom.altloc);
                if (cc.find_atom(atom.name) == cc.atoms.end())
                  atom.name[0] = 'H';
              }
            st.has_d_fraction = true;
          }
        }
        if (reorder && ri.orig_chemcomp) {
          const ChemComp& cc = *ri.orig_chemcomp;
          for (Atom& atom : res.atoms) {
            auto it = cc.find_atom(atom.name);
            if (it == cc.atoms.end())
              topo->err("definition not found for " +
                        atom_str(chain_info.chain_ref, *ri.res, atom));
            atom.serial = int(it - cc.atoms.begin()); // temporary, for sorting only
          }
          std::sort(res.atoms.begin(), res.atoms.end(), [](const Atom& a, const Atom& b) {
                      return a.serial != b.serial ? a.serial < b.serial
                                                  : a.altloc < b.altloc;
          });
        }
      }
    }
  }

  // for atoms with ad-hoc links, for now we don't want hydrogens
  if (!ignore_unknown_links && h_change != HydrogenChange::NoChange)
    for (const Topo::Link& link : topo->extras) {
      const ChemLink* cl = monlib.get_link(link.link_id);
      if (cl && starts_with(cl->name, "auto-")) {
        const Restraints::Bond& bond = cl->rt.bonds.at(0);
        remove_hydrogens_from_atom(topo->find_resinfo(link.res1), bond.id1.atom, link.alt1);
        remove_hydrogens_from_atom(topo->find_resinfo(link.res2), bond.id2.atom, link.alt2);
      }
    }

  assign_serial_numbers(st.models[model_index]);
  topo->finalize_refmac_topology(monlib);

  // the hydrogens added previously have positions not set
  if (h_change != HydrogenChange::NoChange)
    place_hydrogens_on_all_atoms(*topo);

  return topo;
}

} // namespace gemmi
#endif
