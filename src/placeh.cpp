// Copyright 2018 Global Phasing Ltd.

#include "placeh.h"
#include <gemmi/chemcomp.hpp> // for ChemComp

using gemmi::Topo;
using gemmi::Restraints;
using gemmi::Position;

struct BondedAtom {
  gemmi::Atom* ptr;
  gemmi::Position& pos; // == ptr->pos;
  double dist;
};

// Calculate position using one angle (theta) and one dihedral angle (tau).
// Returns position of x4 in x1-x2-x3-x4, where dist=|x3-x4| and
// theta is angle(x2, x3, x4).
// Based on section 3.3 of Paciorek et al, Acta Cryst. A52, 349 (1996).
static
Position position_from_angle_and_torsion(const Position& x1,
                                         const Position& x2,
                                         const Position& x3,
                                         double dist,  // |x3-x4|
                                         double theta, // angle x2-x3-x4
                                         double tau) { // dihedral angle
  gemmi::Vec3 u = x2 - x1;
  gemmi::Vec3 v = x3 - x2;
  gemmi::Vec3 e1 = v.normalized();
  double delta = u.dot(e1);
  gemmi::Vec3 e2 = -(u - delta * e1).normalized();
  gemmi::Vec3 e3 = e1.cross(e2);
  return x3 + Position(dist * (-cos(theta) * e1 +
                               sin(theta) * (cos(tau) * e2 + sin(tau) * e3)));
}

// Calculate position using two angles.
// Based on https://en.wikipedia.org/wiki/Trilateration
// Returns p4. Topology: p1 is bonded to p2, p3 and p4.
static
Position position_from_two_angles(const Position& p1,
                                  const Position& p2,
                                  const Position& p3,
                                  double dist14,   // |p4-p1|
                                  double theta1,   // angle p4-p1-p2
                                  double theta2) { // angle p4-p1-p3
  double dist12 = p1.dist(p2);
  double dist13 = p1.dist(p3);
  return p1;
}

static
void place_hydrogen_1_1(BondedAtom& h, const gemmi::Atom& a,
                        const BondedAtom& b, const Topo& topo) {
  const Restraints::Angle* angle = topo.take_angle(h.ptr, &a, b.ptr);
  if (!angle)
    return;
  if (angle->value == 180.0) {
    gemmi::Vec3 u = a.pos - h.pos;
    h.pos = a.pos + Position(u * (h.dist / u.length()));
    return;
  }
  double theta = gemmi::rad(angle->value);
  double tau = 0.0;
  const gemmi::Atom* c = nullptr;
  for (const Topo::Plane& plane : topo.planes) {
    if (plane.atoms.size() > 3 &&
        plane.has(h.ptr) && plane.has(&a) && plane.has(b.ptr)) {
      for (const gemmi::Atom* maybe_c : plane.atoms) {
        if (maybe_c != h.ptr && maybe_c != &a && maybe_c != b.ptr) {
          c = maybe_c;
          break;
        }
      }
      break;
    }
  }
  if (!c) {
    // using one dihedral angle
    for (const Topo::Torsion& tor : topo.torsions) {
      if (tor.atoms[0] == h.ptr && !tor.atoms[3]->is_hydrogen() &&
          tor.atoms[1] == &a && tor.atoms[2] == b.ptr) {
        tau = gemmi::rad(tor.restr->value);
        c = tor.atoms[3];
        break;
      } else if (tor.atoms[3] == h.ptr && !tor.atoms[0]->is_hydrogen() &&
                 tor.atoms[1] == b.ptr && tor.atoms[2] == &a) {
        tau = gemmi::rad(tor.restr->value);
        c = tor.atoms[0];
        break;
      }
    }
  }
  h.pos = position_from_angle_and_torsion(c ? c->pos : Position(0, 0, 0),
                                          b.pos, a.pos, h.dist, theta, tau);
  h.ptr->occ = 0; // the position is not unique
}

void place_hydrogens(const gemmi::Atom& atom, Topo::ResInfo& ri,
                     const Topo& topo) {
  using Angle = Restraints::Angle;

  // put atoms bonded to atom into two lists
  std::vector<BondedAtom> known; // heavy atoms with known positions
  std::vector<BondedAtom> hs;    // H atoms (unknown)
  for (const Topo::Force& force : ri.forces)
    if (force.rkind == Topo::RKind::Bond) {
      const Topo::Bond& t = topo.bonds[force.index];
      int n = Topo::has_atom(&atom, t);
      if (n == 0 || n == 1) {
        gemmi::Atom* other = t.atoms[1-n];
        auto& atom_list = other->is_hydrogen() ? hs : known;
        atom_list.push_back({other, other->pos, t.restr->value});
      }
    }

  if (hs.size() == 0)
    return;

  if (known.size() == 0) {
    hs[0].pos = atom.pos + Position(hs[0].dist, 0, 0);
    if (hs.size() > 1) {
      double theta = M_PI;
      if (const Angle* ang = topo.take_angle(hs[1].ptr, &atom, hs[0].ptr))
        theta = gemmi::rad(ang->value);
      hs[1].pos = atom.pos + Position(hs[1].dist * cos(theta),
                                      hs[1].dist * sin(theta), 0);
    }
    if (hs.size() > 2) {
      double theta1 = gemmi::rad(180);
      double theta2 = gemmi::rad(180);
      if (const Angle* ang = topo.take_angle(hs[2].ptr, &atom, hs[0].ptr))
        theta1 = gemmi::rad(ang->value);
      if (const Angle* ang = topo.take_angle(hs[2].ptr, &atom, hs[1].ptr))
        theta2 = gemmi::rad(ang->value);
      // TODO 3rd atom
    }
    // TODO 4th atom
    for (BondedAtom& bonded_h : hs)
      bonded_h.ptr->occ = 0;

  } else if (known.size() == 1 && hs.size() == 1) {
    place_hydrogen_1_1(hs[0], atom, known[0], topo);

  } else if (known.size() == 2 && hs.size() == 1) {
    double phi = gemmi::rad(120);
    double psi = gemmi::rad(120);
    if (const Angle* ang = topo.take_angle(hs[0].ptr, &atom, known[0].ptr))
      phi = gemmi::rad(ang->value);
    if (const Angle* ang = topo.take_angle(hs[0].ptr, &atom, known[1].ptr))
      psi = gemmi::rad(ang->value);
    hs[0].pos = position_from_two_angles(known[0].pos, atom.pos, known[1].pos,
                                         hs[0].dist, phi, psi);

    printf("[2xangle] %5s: \n", hs[0].ptr->name.c_str());

  // TODO: all other cases

  } else {
    printf("[not done] %5s: %zd %zd\n",
           atom.name.c_str(), hs.size(), known.size());
    for (BondedAtom& bonded_h : hs)
      bonded_h.ptr->occ = 0;
  }
}

#ifdef PLACEH_MAIN

#include <stdio.h>
#include <gemmi/cif.hpp>   // for read_file
#include <gemmi/util.hpp>  // for fail

using namespace gemmi;


void print_restraint_summary(const std::string& id, const ChemComp& cc) {
  printf("%-5s %-5s ", cc.name.c_str(), id.c_str());
  fflush(stdout);
  std::string heavy_atom;

  for (const Restraints::Bond& bond : cc.rt.bonds) {
    const Restraints::AtomId* other_end = nullptr;
    if (bond.id1 == id)
      other_end = &bond.id2;
    else if (bond.id2 == id)
      other_end = &bond.id1;
    if (other_end) {
      if (!heavy_atom.empty())
        fail("H atom with 2+ bonds");
      heavy_atom = other_end->atom;
    }
  }
  if (heavy_atom.empty())
    fail("non-bonded H");
  printf("%-4s ", heavy_atom.c_str());

  int angle_count = 0;
  for (const Restraints::Angle& angle : cc.rt.angles) {
    const Restraints::AtomId* other_end = nullptr;
    if (angle.id1 == id)
      other_end = &angle.id3;
    else if (angle.id3 == id)
      other_end = &angle.id1;
    if (other_end) {
      if (angle.id2 != heavy_atom)
        fail("_chem_comp_angle.atom_id_2 is not H's heavy atom.");
      //printf("%.1f deg to %s, ", angle.value, other_end->atom.c_str());
      if (!cc.get_atom(other_end->atom).is_hydrogen()) {
        ++angle_count;
        if (angle.value == 180.0)
          angle_count += 100;
      }
    }
  }

  int tor_count = 0;
  for (const Restraints::Torsion& tor : cc.rt.torsions) {
    //if (tor.period > 1) continue;
    const Restraints::AtomId* other[3] = {nullptr, nullptr, nullptr};
    if (tor.id1 == id) {
      other[0] = &tor.id2;
      other[1] = &tor.id3;
      other[2] = &tor.id4;
    } else if (tor.id4 == id) {
      other[0] = &tor.id3;
      other[1] = &tor.id2;
      other[2] = &tor.id1;
    }
    if (other[0]) {
      if (other[0]->atom != heavy_atom)
        fail("_chem_comp_tor atom next to H is not H's heavy atom.");
      if (!cc.get_atom(other[2]->atom).is_hydrogen())
        ++tor_count;
    }
  }

  int chir_count = 0;
  for (const Restraints::Chirality& chir : cc.rt.chirs) {
    if (chir.id1 == id || chir.id2 == id || chir.id3 == id) {
      if (chir.id_ctr != heavy_atom)
        fail("_chem_comp_chir atom next to H is not H's heavy atom.");
      ++chir_count;
    }
  }

  int plane_count = 0;
  for (const Restraints::Plane& plane : cc.rt.planes) {
    const std::vector<Restraints::AtomId>& ids = plane.ids;
    if (std::find(ids.begin(), ids.end(), id) != ids.end()) {
      if (std::find(ids.begin(), ids.end(), heavy_atom) == ids.end())
        fail("H in _chem_comp_plane without its heavy atom.");
      ++plane_count;
    }
  }

  printf("%d angles, %d torsions, %d chiralities, %d planes\n",
         angle_count, tor_count, chir_count, plane_count);
}

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    cif::Document doc = cif::read_file(argv[i]);
    for (const cif::Block& block : doc.blocks)
      if (block.name != "comp_list") {
        ChemComp cc = make_chemcomp_from_block(block);
        for (const ChemComp::Atom& atom : cc.atoms)
          if (atom.el == El::H)
            try {
              print_restraint_summary(atom.id, cc);
            } catch (std::runtime_error& e) {
              fprintf(stderr, "%s %s:%s\n",
                      block.name.c_str(), atom.id.c_str(), e.what());
            }
      }
  }
}
#endif

// vim:sw=2:ts=2:et
