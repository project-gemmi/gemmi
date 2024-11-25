// Copyright 2018-2022 Global Phasing Ltd.

#include <gemmi/riding_h.hpp>

namespace gemmi {

namespace {

struct BondedAtom {
  Atom* ptr;
  Position& pos;  // == ptr->pos;
  double dist;
};

// Calculate position using one angle (theta) and one dihedral angle (tau).
// Returns position of x4 in x1-x2-x3-x4, where dist=|x3-x4| and
// theta is angle(x2, x3, x4).
// Based on section 3.3 of Paciorek et al, Acta Cryst. A52, 349 (1996).
Position position_from_angle_and_torsion(const Position& x1,
                                         const Position& x2,
                                         const Position& x3,
                                         double dist,  // |x3-x4|
                                         double theta, // angle x2-x3-x4
                                         double tau) { // dihedral angle
  Vec3 u = x2 - x1;
  Vec3 v = x3 - x2;
  Vec3 e1 = v.normalized();
  Vec3 e2 = -(u - u.dot(e1) * e1).normalized();
  Vec3 e3 = e1.cross(e2);
  Vec3 e23 = std::cos(tau) * e2 + std::sin(tau) * e3;
  return x3 + Position(dist * (-std::cos(theta) * e1 + std::sin(theta) * e23));
}

// Similar to position_from_angle_and_torsion(), but x1 and tau are not given.
Position arbitrary_position_from_angle(const Position& x2,
                                       const Position& x3,
                                       double dist,     // |x3-x4|
                                       double theta) {  // angle x2-x3-x4
  Vec3 u(1, 0, 0);
  Vec3 v = x3 - x2;
  Vec3 e1 = v.normalized();
  Vec3 e2_ = -(u - u.dot(e1) * e1);
  if (e2_.length_sq() < 1e-6) {
    // u || v, let's take any non-parallel u
    u = Vec3(0, 1, 0);
    e2_ = -(u - u.dot(e1) * e1);
  }
  Vec3 e2 = e2_.normalized();
  return x3 + Position(dist * (-std::cos(theta) * e1 + std::sin(theta) * e2));
}


Vec3 get_vector_to_line(const Position& point,
                        const Position& point_on_the_line,
                        const Vec3& unit_vector) {
  // en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
  // the component of a - p perpendicular to the line is: (a-p) - ((a-p).n)n
  Vec3 ap = point_on_the_line - point;
  return ap - ap.dot(unit_vector) * unit_vector;
}

// If no points satisfy the distances returns a pair of NaNs
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
std::pair<Position, Position>
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
  return trilaterate(p1, d14sq, p2, d24sq, p3, d34sq);
}

// Returns angle between hydrogen and the plane of heavy atoms
// in 2H tetrahedral configuration. theta0 is the angle between heavy atoms.
double calculate_tetrahedral_delta(double theta0, double theta1, double theta2) {
  // simplified trilateration:
  //   auto r = trilaterate(Position(0, 0, 0), 1,
  //                        Position(1, 0, 0), 2 - 2 * std::cos(theta1),
  //                        Position(std::cos(theta0), std::sin(theta0), 0),
  //                        2 - 2 * std::cos(theta2));
  //   return std::asin(std::fabs(r.first.z));
  double x = std::cos(theta1);
  double y = (std::cos(theta2) - x * std::cos(theta0)) / std::sin(theta0);
  double z2 = 1 - x*x - y*y;
  double z = std::sqrt(z2);  // may result in NaN
  return std::asin(z);
}

double angle_in_triangle(double a, double b, double c) {
  return std::acos((a * a + c * c - b * b) / (2 * a * c));
}

// Used in rare cases when one angle is missing in 2H-tetrahedral configuration.
// Currently it can happen only with metal sites.
//   b   c(metal)
//    \ /
//   atom -- h (possibly two H atoms, symmetric wrt plane abc)
// Angle c-atom-h is missing, b-atom-h is equal alpha, b-atom-c is theta.
// For two H atoms we also need angle h-atom-h.
double missing_angle_2H_tetrahedral(const Topo& topo, const Atom& atom,
                                    const std::vector<BondedAtom>& hs,
                                    double alpha, double theta) {
  if (hs.size() == 1)
    return 2 * pi() - alpha - theta;  // make it coplanar
  if (hs.size() == 2) {
    if (const Restraints::Angle* hh = topo.take_angle(hs[0].ptr, &atom, hs[1].ptr)) {
      double xh = std::cos(alpha);
      double zh = std::sin(0.5 * hh->radians());
      double yh = std::sqrt(1 - xh*xh - zh*zh);
      double st = std::sin(theta);
      double ct = std::cos(theta);
      Vec3 ah = Vec3(xh, -yh * st, zh * st);
      double angle = std::acos((ct * ah.x + st * ah.y) / ah.length());
      //printf("missing_angle: %s %g\n", atom.name.c_str(), deg(angle));
      return angle;
    }
  }
  return 0;
}

[[noreturn]]
void giveup(const std::string& message, const std::vector<BondedAtom>& hs) {
  for (const BondedAtom& bonded_h : hs) {
    bonded_h.ptr->occ = 0;
    bonded_h.ptr->calc_flag = CalcFlag::Dummy;
  }
  fail(message);
}

// known and hs are lists of heavy atoms and hydrogens bonded to atom.
// hs is const, but nevertheless atoms it points to are modified.
void place_hydrogens(const Topo& topo, const Atom& atom,
                     const std::vector<BondedAtom>& known,
                     const std::vector<BondedAtom>& hs) {
  using Angle = Restraints::Angle;
  assert(!hs.empty());

  for (const BondedAtom& bonded_h : hs)
    bonded_h.ptr->calc_flag = CalcFlag::Calculated;

  // ==== only hydrogens ====
  if (known.size() == 0) {
    // we can only arbitrarily pick directions of atoms
    for (const BondedAtom& bonded_h : hs)
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
      const Angle* ang0 = topo.take_angle(hs[2].ptr, &atom, hs[0].ptr);
      const Angle* ang1 = topo.take_angle(hs[2].ptr, &atom, hs[1].ptr);
      double theta0 = rad(ang0 ? ang0->value : 109.47122);
      double theta1 = rad(ang1 ? ang1->value : 109.47122);
      auto pos = position_from_two_angles(atom.pos, hs[0].pos, hs[1].pos,
                                          hs[2].dist, theta0, theta1);
      hs[2].pos = pos.first;
      // 3 hydrogens: for now happens only NH3 (NH2.cif and NH3.cif)
      // 4 hydrogens: CH4 (CH2.cif) and NH4 (NH4.cif)
      if (hs.size() == 4)
        hs[3].pos = pos.second;
    }

  // ==== one heavy atom and hydrogens ====
  } else if (known.size() == 1) {
    const BondedAtom& h = hs[0];
    const BondedAtom& heavy = known[0];
    const Restraints::Angle* angle = topo.take_angle(h.ptr, &atom, heavy.ptr);
    if (!angle)
      giveup("No angle restraint for " + h.ptr->name + ".\n", hs);
    if (std::abs(angle->value - 180.0) < 0.5) {
      Vec3 u = atom.pos - h.pos;
      h.pos = atom.pos + Position(u * (h.dist / u.length()));
      if (hs.size() > 1)
        giveup("Unusual: one of two H atoms has angle restraint 180 deg.", hs);
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
        }
        if (tor.atoms[2] == &atom && tor.atoms[1] == heavy.ptr &&
            tor.atoms[3]->is_hydrogen() && !tor.atoms[0]->is_hydrogen()) {
          tau = rad(tor.restr->value);
          torsion_h = tor.atoms[3];
          tau_end = tor.atoms[0];
          period = tor.restr->period;
          break;
        }
      }
    }

    if (tau_end)
      h.pos = position_from_angle_and_torsion(tau_end->pos, heavy.pos, atom.pos,
                                              h.dist, theta, tau);
    else
      h.pos = arbitrary_position_from_angle(heavy.pos, atom.pos, h.dist, theta);
    if (std::isnan(h.pos.x)) {
      h.pos = Position(0, 0, 0);
      giveup("bonded atoms are exactly overlapping (case 1).", hs);
    }
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
      giveup("Unusual: atom bonded to one heavy atoms and 4+ hydrogens.", hs);
    }
    if (!tau_end || period > (int)hs.size())
      for (const BondedAtom& bonded_h : hs)
        bonded_h.ptr->occ = 0;
  // ==== two heavy atoms and hydrogens ====
  } else if (known.size() == 2) {
    if (hs.size() >= 3)
      giveup("Unusual: atom bonded to 2+ heavy atoms and 3+ hydrogens.", hs);
    const Angle* ang1 = topo.take_angle(hs[0].ptr, &atom, known[0].ptr);
    const Angle* ang2 = topo.take_angle(hs[0].ptr, &atom, known[1].ptr);
    const Angle* ang3 = topo.take_angle(known[0].ptr, &atom, known[1].ptr);

    double theta1 = ang1 ? ang1->radians() : 0;
    double theta2 = ang2 ? ang2->radians() : 0;
    double theta3 = ang3 ? ang3->radians() : 0;
    if (!ang3) {
      // If ang3 atoms form a bonded triangle (e.g. C11-C12-PT in DVW),
      // we can calculate theta3 from the bond lengths.
      const Restraints::Bond* aptr = topo.take_bond(&atom, known[1].ptr);
      const Restraints::Bond* bptr = topo.take_bond(known[0].ptr, known[1].ptr);
      const Restraints::Bond* cptr = topo.take_bond(known[0].ptr, &atom);
      if (!aptr || !bptr || !cptr)
        giveup(cat("Missing angle restraint ", known[0].ptr->name, '-', atom.name,
                   '-', known[1].ptr->name, ".\n"), hs);
      theta3 = angle_in_triangle(aptr->value, bptr->value, cptr->value);
    }
    // Some configurations with metals have only distance restraints.
    // In such cases, we assume that hydrogens are as far as possible from the metal.
    if (!ang1 && ang2 && known[0].ptr->element.is_metal())
      theta1 = missing_angle_2H_tetrahedral(topo, atom, hs, theta2, theta3);
    if (!ang2 && ang1 && known[1].ptr->element.is_metal())
      theta2 = missing_angle_2H_tetrahedral(topo, atom, hs, theta1, theta3);
    if (theta1 == 0 || theta2 == 0)
      giveup(cat("Missing angle restraint ", hs[0].ptr->name, '-', atom.name,
                 '-', known[theta1 == 0 ? 0 : 1].ptr->name, ".\n"), hs);

    // Co-planar case. The sum of angles should be 360 degrees,
    // but in some cif files it differs slightly.
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
        topo.logger.err("Unhandled topology of ", hs.size(), " hydrogens bonded to ", atom.name);
        for (size_t i = 1; i < hs.size(); ++i) {
          hs[i].ptr->occ = 0;
          hs[i].ptr->calc_flag = CalcFlag::Dummy;
        }
      }
      return;
    }

    // Tetrahedral or similar configuration.
    double hh_half = 0;
    if (hs.size() == 2)
      if (const Angle* hh = topo.take_angle(hs[0].ptr, &atom, hs[1].ptr))
        hh_half = 0.5 * hh->radians();
    if (hh_half == 0)
      hh_half = calculate_tetrahedral_delta(theta3, theta1, theta2);

    // Based on Liebschner et al (2020) doi:10.1016/bs.mie.2020.01.007
    // sec. 2.3. 2H-tetrahedral configuration
    double c0 = std::cos(theta3);
    double c1 = std::cos(theta1);
    double c2 = std::cos(theta2);
    double den = 1 / (1 - c0*c0);
    double a = den * (c1 - c0 * c2);
    double b = den * (c2 - c0 * c1);
    // I think the paper defines u10 and u20 in the opposite direction,
    // but I had to reverse it somewhere to make it work.
    Vec3 u10 = (known[0].pos - atom.pos).normalized();
    Vec3 u20 = (known[1].pos - atom.pos).normalized();
    Vec3 v = u10.cross(u20);
    if (std::isnan(v.x))
      giveup("bonded atoms are exactly overlapping (case 2).", hs);
    Vec3 d = a * u10 + b * u20;
    double dist_sin = hs[0].dist * std::sin(hh_half);
    double dist_cos = hs[0].dist * std::cos(hh_half);
    Vec3 v0s = v.changed_magnitude(dist_sin);
    Vec3 d0c = d.changed_magnitude(dist_cos);
    hs[0].pos = atom.pos + Position(d0c + v0s);
    Position other_pos = atom.pos + Position(d0c - v0s);

    if (hs.size() == 1) {
      const Topo::Chirality* chir = topo.get_chirality(&atom);
      if (chir && chir->restr->sign != ChiralityType::Both) {
        if (!chir->check())
          hs[0].pos = other_pos;
      } else {
        hs[0].ptr->occ = 0;
      }
    } else {  // hs.size() == 2
      hs[1].pos = other_pos;
    }

  } else {  // known.size() >= 3
    if (hs.size() > 1)
      giveup("Unusual: atom bonded to 3+ heavy atoms and 2+ hydrogens.", hs);
    // Based on Liebschner et al (2020) doi:10.1016/bs.mie.2020.01.007
    // sec. 2.4. 1H-tetrahedral configuration
    Vec3 u10 = (known[0].pos - atom.pos).normalized();
    Vec3 u20 = (known[1].pos - atom.pos).normalized();
    Vec3 u30 = (known[2].pos - atom.pos).normalized();
    auto cos_tetrahedral = [&](size_t n) {
      const Angle* angle = topo.take_angle(known[n].ptr, &atom, hs[0].ptr);
      return angle ? std::cos(angle->radians()) : -1./3.;
    };
    SMat33<double> m{1., 1., 1., u10.dot(u20), u10.dot(u30), u20.dot(u30)};
    Vec3 rhs(cos_tetrahedral(0), cos_tetrahedral(1), cos_tetrahedral(2));
    double det = m.determinant();
    if (std::fabs(det) < 1e-12)
      giveup("tetrahedral configuration with four co-planar atoms.", hs);
    Vec3 abc = m.inverse_(det).multiply(rhs);
    Vec3 h_dir = abc.x * u10 + abc.y * u20 + abc.z * u30;
    hs[0].pos = atom.pos + Position(h_dir.changed_magnitude(hs[0].dist));
  }
}

}  // anonymous namespace

void place_hydrogens_on_all_atoms(Topo& topo) {
  std::vector<BondedAtom> known;
  std::vector<BondedAtom> hs;
  auto filter = [](char alt, const std::vector<BondedAtom>& v) {
    std::vector<BondedAtom> out;
    out.reserve(v.size());
    for (const BondedAtom& ba : v)
      if (ba.ptr->altloc_matches(alt))
        out.push_back(ba);
    return out;
  };
  for (Topo::ChainInfo& chain_info : topo.chain_infos)
    for (Topo::ResInfo& ri : chain_info.res_infos) {
      // If we don't have monomer description from a cif file,
      // only ad-hoc restraints, don't try to place hydrogens.
      if (ri.orig_chemcomp == nullptr)
        continue;
      for (Atom& atom : ri.res->atoms) {
        if (atom.is_hydrogen())
          continue;
        try {
          // gather bonded atoms
          known.clear();
          hs.clear();
          auto range = topo.bond_index.equal_range(&atom);
          for (auto i = range.first; i != range.second; ++i) {
            const Topo::Bond* t = i->second;
            Atom* other = t->atoms[t->atoms[0] == &atom ? 1 : 0];
            if (other->altloc && atom.altloc) {
              // We support links between different altlocs in Topo (e.g. link A-B),
              // although these are rare, special cases.
              // But if we had bonds between atom 1 (A/B) and atom 2 (A/B/C),
              // and we had bonds B-B and B-C, we'd want to use only one of them (B-B).
              // Checking atom's name is not robust, but should suffice here.
              if (atom.altloc != other->altloc &&
                  in_vector_f([&](const BondedAtom& a) { return a.ptr->name == other->name; },
                              known))
                continue;
            }
            auto& atom_list = other->is_hydrogen() ? hs : known;
            atom_list.push_back({other, other->pos, t->restr->value});
          }
          if (hs.size() == 0)
            continue;
          std::string altlocs;
          // In a special case: Hs with altlocs on a parent without altloc,
          // we need to process conformations one by one.
          if (atom.altloc == '\0')
            for (const auto& h : hs) {
              char alt = h.ptr->altloc;
              if (alt && altlocs.find(alt) == std::string::npos) {
                altlocs += alt;
                place_hydrogens(topo, atom, filter(alt, known), filter(alt, hs));
              }
            }
          // In all other cases, all bonded atoms are from the same conformation.
          if (altlocs.empty())
            place_hydrogens(topo, atom, known, hs);
        } catch (const std::runtime_error& e) {
          topo.logger.err("Placing of hydrogen bonded to ",
                          atom_str(chain_info.chain_ref, *ri.res, atom),
                          " failed:\n  ", e.what());
        }
      }
    }
}

}  // namespace gemmi
