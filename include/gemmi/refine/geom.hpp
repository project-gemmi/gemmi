// Copyright 2023 MRC Laboratory of Molecular Biology
//

#ifndef GEMMI_REFINE_GEOM_HPP_
#define GEMMI_REFINE_GEOM_HPP_

#include <set>
#include "../model.hpp"    // for Structure, Atom
#include "../contact.hpp"  // for NeighborSearch, ContactSearch
#include "../topo.hpp"     // for Topo
#include "../select.hpp"   // for count_atom_sites

namespace gemmi {

Mat33 eigen_decomp_inv(const SMat33<double> &m, double e, bool for_precond) {
  const auto eig = m.calculate_eigenvalues();
  // good e = 1.e-9 for plane and ~1e-6 or 1e-4 for precondition
  auto f = [&](double v){
    if (std::abs(v) < e) return 0.;
    return for_precond ? 1. / std::sqrt(v) : 1. / v;
  };
  const Vec3 l{f(eig[0]), f(eig[1]), f(eig[2])};
  Mat33 Q; // formed by eigenvectors
  for (int j = 0; j < 3; ++j) {
    const auto v = m.calculate_eigenvector(eig[j]);
    for (int i = 0; i < 3; ++i) Q[i][j] = v.at(i);
  }

  if (for_precond)
    return Q.multiply_by_diagonal(l);
  else
    return Q.multiply_by_diagonal(l).multiply(Q.transpose());
}

struct PlaneDeriv {
  PlaneDeriv(const std::vector<Atom*> &atoms)
    : dvmdx(atoms.size(), std::vector<Vec3>(3)), dDdx(atoms.size()) {
    // centroid
    for (const Atom* atom : atoms) xs += atom->pos;
    xs /= (double) atoms.size();

    // from gemmi/calculate.hpp find_best_plane()
    // for plane normal
    SMat33<double> A{0, 0, 0, 0, 0, 0};
    for (const Atom* atom : atoms) {
      const Vec3 p = Vec3(atom->pos) - xs;
      A.u11 += p.x * p.x;
      A.u22 += p.y * p.y;
      A.u33 += p.z * p.z;
      A.u12 += p.x * p.y;
      A.u13 += p.x * p.z;
      A.u23 += p.y * p.z;
    }
    const auto eig = A.calculate_eigenvalues();
    double eigmin = eig[0];
    for (double d : {eig[1], eig[2]})
      if (std::fabs(d) < std::fabs(eigmin))
        eigmin = d;

    vm = A.calculate_eigenvector(eigmin);
    D = xs.dot(vm);

    // derivatives
    const Mat33 pinv = eigen_decomp_inv(SMat33<double>{eigmin,eigmin,eigmin,0,0,0}-A, 1e-9, false);
    for (size_t i = 0; i < atoms.size(); ++i) {
      const Vec3 p = Vec3(atoms[i]->pos) - xs;
      const SMat33<double> dAdx[3] = {{2 * p.x, 0.,      0.,      p.y, p.z, 0.},   // dA/dx
                                      {0.,      2 * p.y, 0.,      p.x, 0.,  p.z},  // dA/dy
                                      {0.,      0.,      2 * p.z, 0.,  p.x, p.y}}; // dA/dz
      for (size_t j = 0; j < 3; ++j)
        dvmdx[i][j] = pinv.multiply(dAdx[j].multiply(vm));
    }
    for (size_t i = 0; i < atoms.size(); ++i) {
      dDdx[i] = vm / (double) atoms.size();
      for (size_t j = 0; j < 3; ++j)
        dDdx[i].at(j) += xs.dot(dvmdx[i][j]);
    }
  }

  void flip() {
    vm *= -1;
    D *= -1;
    for (auto &x : dvmdx)
      for (auto &y : x)
        y *= -1;
    for (auto &x : dDdx)
      x *= -1;
  }

  Vec3 xs; // centroid
  Vec3 vm; // normal vector
  double D; // xs dot vm
  std::vector<std::vector<Vec3>> dvmdx; // derivative of vm wrt positions
  std::vector<Vec3> dDdx; // derivative of D wrt positions
};

struct GeomTarget {
  struct MatPos {
    int ipos;
    int imode;
  };

  void setup(Model &model, bool refine_xyz, int adp_mode) { // call it after setting pairs
    assert(0 <= adp_mode && adp_mode <= 2);
    this->refine_xyz = refine_xyz;
    this->adp_mode = adp_mode;
    const size_t n_atoms = count_atom_sites(model);
    atoms.resize(n_atoms);
    for (CRA cra : model.all())
      atoms[cra.atom->serial-1] = cra.atom;

    target = 0.;
    vn.clear();
    am.clear();
    rest_per_atom.clear();
    rest_pos_per_atom.clear();
    const size_t n_pairs = pairs.size();
    size_t qqm = 0, qqv = 0;
    if (refine_xyz) {
      qqv += 3 * n_atoms;
      qqm += 6 * n_atoms + 9 * n_pairs;
    }
    if (adp_mode == 1) {
      qqv += n_atoms;
      qqm += n_atoms + n_pairs;
    }
    else if (adp_mode == 2) {
      qqv += 6 * n_atoms;
      qqm += 21 * n_atoms + 36 * n_pairs;
    }
    vn.assign(qqv, 0.);
    am.assign(qqm, 0.);

    std::vector<size_t> nrest_per_atom(n_atoms, 0);
    for (auto & p : pairs) {
      ++nrest_per_atom[p.first];
      if (p.first != p.second)
        ++nrest_per_atom[p.second];
    }

    rest_pos_per_atom.assign(n_atoms+1, 0);
    for (size_t ia = 0; ia < n_atoms; ++ia)
      rest_pos_per_atom[ia+1] = rest_pos_per_atom[ia] + nrest_per_atom[ia];

    rest_per_atom.assign(std::accumulate(nrest_per_atom.begin(), nrest_per_atom.end(), 0), 0);
    for (size_t i = 0; i < n_atoms; ++i) nrest_per_atom[i] = 0;
    for (size_t i = 0; i < rest_per_atom.size(); ++i) rest_per_atom[i] = 0;
    for (size_t i = 0; i < pairs.size(); ++i) {
      int ia1 = pairs[i].first, ia2 = pairs[i].second;
      nrest_per_atom[ia1] += 1;
      int ip1 = rest_pos_per_atom[ia1] + nrest_per_atom[ia1] - 1;
      rest_per_atom[ip1] = i;
      if (ia1 != ia2) {
        nrest_per_atom[ia2] += 1;
        int ip2 = rest_pos_per_atom[ia2] + nrest_per_atom[ia2] - 1;
        rest_per_atom[ip2] = i;
      }
    }
    nmpos = 6 * n_atoms;
  }

  //int n_atoms;
  bool refine_xyz;
  int adp_mode; // 0: no refine, 1: iso, 2: aniso
  std::vector<Atom*> atoms;
  double target = 0.; // target function value
  std::vector<double> vn; // first derivatives
  std::vector<double> am; // second derivative sparse matrix
  std::vector<size_t> rest_per_atom;
  std::vector<size_t> rest_pos_per_atom;
  std::vector<std::pair<int,int>> pairs;
  int nmpos;
  size_t n_atoms() const { return atoms.size(); }
  MatPos find_restraint(int ia1, int ia2) const {
    MatPos matpos;
    int idist = -1;
    for (size_t irest = rest_pos_per_atom[ia1]; irest < rest_pos_per_atom[ia1+1]; ++irest) {
      int ir1 = rest_per_atom[irest];
      if (pairs[ir1].first == ia2) {
        // atom ia1 is target atom and atom ia2 is object
        idist = ir1;
        matpos.imode = 0;
        break;
      }
      else if (pairs[ir1].second == ia2) {
        // atom ia1 is object atom and atom ia2 is target
        idist = ir1;
        matpos.imode = 1;
        break;
      }
    }
    if (idist < 0) fail("cannot find atom pair");
    matpos.ipos = nmpos + 9 * idist;
    return matpos;
  }
  void incr_vn(size_t ipos, double w, const Vec3 &deriv) {
    assert(ipos+2 < vn.size());
    vn[ipos]   += w * deriv.x;
    vn[ipos+1] += w * deriv.y;
    vn[ipos+2] += w * deriv.z;
  }
  void incr_am_diag(size_t ipos, double w, const Vec3 &deriv) {
    assert(ipos+5 < am.size());
    am[ipos]   += w * deriv.x * deriv.x;
    am[ipos+1] += w * deriv.y * deriv.y;
    am[ipos+2] += w * deriv.z * deriv.z;
    am[ipos+3] += w * deriv.y * deriv.x;
    am[ipos+4] += w * deriv.z * deriv.x;
    am[ipos+5] += w * deriv.z * deriv.y;
  }
  void incr_am_diag12(size_t ipos, double w, const Vec3 &deriv1, const Vec3 &deriv2) {
    // used when atoms are related with each other for example through symmetry
    assert(ipos+5 < am.size());
    am[ipos]   += w *  deriv1.x * deriv2.x * 2.;
    am[ipos+1] += w *  deriv1.y * deriv2.y * 2.;
    am[ipos+2] += w *  deriv1.z * deriv2.z * 2.;
    am[ipos+3] += w * (deriv1.x * deriv2.y + deriv2.x * deriv1.y);
    am[ipos+4] += w * (deriv1.x * deriv2.z + deriv2.x * deriv1.z);
    am[ipos+5] += w * (deriv1.y * deriv2.z + deriv2.y * deriv1.z);
  }
  void incr_am_ndiag(size_t ipos, double w, const Vec3 &deriv1, const Vec3 &deriv2) {
    assert(ipos+8 < am.size());
    am[ipos]   += w * deriv1.x * deriv2.x;
    am[ipos+1] += w * deriv1.y * deriv2.x;
    am[ipos+2] += w * deriv1.z * deriv2.x;
    am[ipos+3] += w * deriv1.x * deriv2.y;
    am[ipos+4] += w * deriv1.y * deriv2.y;
    am[ipos+5] += w * deriv1.z * deriv2.y;
    am[ipos+6] += w * deriv1.x * deriv2.z;
    am[ipos+7] += w * deriv1.y * deriv2.z;
    am[ipos+8] += w * deriv1.z * deriv2.z;
  }
  void get_am_col_row(int *row, int *col) const {
    size_t i = 0, offset = 0;
    if (refine_xyz) {
      for (size_t j = 0; j < n_atoms(); ++j, i+=6) {
        row[i]   = 3*j;   col[i]   = 3*j;
        row[i+1] = 3*j+1; col[i+1] = 3*j+1;
        row[i+2] = 3*j+2; col[i+2] = 3*j+2;
        row[i+3] = 3*j;   col[i+3] = 3*j+1;
        row[i+4] = 3*j;   col[i+4] = 3*j+2;
        row[i+5] = 3*j+1; col[i+5] = 3*j+2;
      }
      for (size_t j = 0; j < pairs.size(); ++j)
        for (size_t k = 0; k < 3; ++k)
          for (size_t l = 0; l < 3; ++l, ++i) {
            row[i] = 3 * pairs[j].second + l;
            col[i] = 3 * pairs[j].first + k;
          }
      offset = 3 * n_atoms();
    }
    if (adp_mode == 1) {
      for (size_t j = 0; j < n_atoms(); ++j, ++i) {
        row[i] = offset + j;
        col[i] = offset + j;
      }
      for (size_t j = 0; j < pairs.size(); ++j, ++i) {
        row[i] = offset + pairs[j].second;
        col[i] = offset + pairs[j].first;
      }
    } else if (adp_mode == 2) {
      for (size_t j = 0; j < n_atoms(); ++j) {
        for (size_t k = 0; k < 6; ++k, ++i)
          row[i] = col[i] = offset + 6 * j + k;
        for (size_t k = 0; k < 6; ++k)
          for (int l = k + 1; l < 6; ++l, ++i) {
            row[i] = offset + 6 * j + k;
            col[i] = offset + 6 * j + l;
          }
      }
      for (size_t j = 0; j < pairs.size(); ++j) {
        for (size_t k = 0; k < 6; ++k)
          for (size_t l = 0; l < 6; ++l, ++i) {
            row[i] = 3 * pairs[j].second + l;
            col[i] = 3 * pairs[j].first + k;
          }
      }
    }
    assert(i == am.size());
  }
};

struct Geometry {
  struct Reporting;
  struct Bond {
    struct Value {
      Value(double v, double s, double vn, double sn) : value(v), sigma(s), value_nucleus(vn), sigma_nucleus(sn) {}
      double value;
      double sigma;
      double value_nucleus;
      double sigma_nucleus;
      // alpha should be here?
    };
    Bond(Atom* atom1, Atom* atom2) : atoms({atom1, atom2}) {}
    void set_image(const NearestImage& im) {
      sym_idx = im.sym_idx;
      std::copy(std::begin(im.pbc_shift), std::end(im.pbc_shift), std::begin(pbc_shift));
    }
    void swap_atoms() {
      std::reverse(atoms.begin(), atoms.end());
      sym_idx = -sym_idx - 1;
    }
    bool same_asu() const {
      return (sym_idx == 0 || sym_idx == -1) && pbc_shift[0]==0 && pbc_shift[1]==0 && pbc_shift[2]==0;
    }
    std::tuple<int,int,int,int,int,int> sort_key() const {
      return std::tie(atoms[0]->serial, atoms[1]->serial, sym_idx, pbc_shift[0], pbc_shift[1], pbc_shift[2]);
    }
    const Value* find_closest_value(double dist, bool use_nucleus) const {
      double db = std::numeric_limits<double>::infinity();
      const Value* ret = nullptr;
      for (const auto &v : values) {
        double tmp = std::abs((use_nucleus ? v.value_nucleus : v.value) - dist);
        if (tmp < db) {
          db = tmp;
          ret = &v;
        }
      }
      return ret;
    }
    double calc(const UnitCell& cell, bool use_nucleus, double wdskal,
                GeomTarget* target, Reporting *reporting) const;

    int type = 1; // 0-2
    double alpha = 1; // only effective for type=2
    int sym_idx = 0; // if negative, atoms need to be swapped.
    std::array<int, 3> pbc_shift = {{0,0,0}};
    std::array<Atom*, 2> atoms;
    std::vector<Value> values;
  };
  struct Angle {
    struct Value {
      Value(double v, double s) : value(v), sigma(s) {}
      double value;
      double sigma;
    };
    Angle(Atom* atom1, Atom* atom2, Atom* atom3) : atoms({atom1, atom2, atom3}) {}
    void swap_atoms() {
      std::reverse(atoms.begin(), atoms.end());
    }
    std::tuple<int,int,int> sort_key() const {
      return std::tie(atoms[0]->serial, atoms[1]->serial, atoms[2]->serial);
    }
    const Value* find_closest_value(double v) const {
      double db = std::numeric_limits<double>::infinity();
      const Value* ret = nullptr;
      for (const auto &value : values) {
        double tmp = angle_abs_diff(value.value, v);
        if (tmp < db) {
          db = tmp;
          ret = &value;
        }
      }
      return ret;
    }
    double calc(double waskal, GeomTarget* target, Reporting *reporting) const;
    int type = 1; // 0 or not
    std::array<Atom*, 3> atoms;
    std::vector<Value> values;
  };
  struct Torsion {
    struct Value {
      Value(double v, double s, int p): value(v), sigma(s), period(p) {}
      double value;
      double sigma;
      int period;
      std::string label;
    };
    Torsion(Atom* atom1, Atom* atom2, Atom* atom3, Atom* atom4) : atoms({atom1, atom2, atom3, atom4}) {}
    void swap_atoms() {
      std::reverse(atoms.begin(), atoms.end());
    }
    std::tuple<int,int,int,int> sort_key() const {
      return std::tie(atoms[0]->serial, atoms[1]->serial, atoms[2]->serial, atoms[3]->serial);
    }
    const Value* find_closest_value(double v) const {
      double db = std::numeric_limits<double>::infinity();
      const Value* ret = nullptr;
      for (const auto &value : values) {
        double tmp = angle_abs_diff(value.value, v, 360. / std::max(1, value.period));
        if (tmp < db) {
          db = tmp;
          ret = &value;
        }
      }
      return ret;
    }
    double calc(double wtskal, GeomTarget* target, Reporting *reporting) const;
    int type = 1; // 0 or not
    std::array<Atom*, 4> atoms;
    std::vector<Value> values;
  };
  struct Chirality {
    Chirality(Atom* atomc, Atom* atom1, Atom* atom2, Atom* atom3) : atoms({atomc, atom1, atom2, atom3}) {}
    double calc(double wchiral, GeomTarget* target, Reporting *reporting) const;
    double value;
    double sigma;
    ChiralityType sign;
    std::array<Atom*, 4> atoms;
  };
  struct Plane {
    Plane(std::vector<Atom*> a) : atoms(a) {}
    double calc(double wplane, GeomTarget* target, Reporting *reporting) const;
    double sigma;
    std::string label;
    std::vector<Atom*> atoms;
  };
  struct Interval {
    Interval(Atom* atom1, Atom* atom2) : atoms({atom1, atom2}) {}
    double dmin;
    double dmax;
    double smin;
    double smax;
    std::array<Atom*, 2> atoms;
  };
  struct Harmonic {
    Harmonic(Atom* a) : atom(a) {}
    double sigma;
    Atom* atom;
  };
  struct Special {
    Special(Atom* a) : atom(a) {}
    double sigma_t;
    double sigma_u;
    bool u_val_incl;
    Transform trans_t;
    Mat33 mat_u;
    Atom* atom;
  };
  struct Stacking {
    Stacking(std::vector<Atom*> plane1, std::vector<Atom*> plane2) : planes({plane1, plane2}) {}
    double calc(double wstack, GeomTarget* target, Reporting *reporting) const;
    double dist;
    double sd_dist;
    double angle;
    double sd_angle;
    std::array<std::vector<Atom*>, 2> planes;
  };
  struct Vdw {
    Vdw(Atom* atom1, Atom* atom2) : atoms({atom1, atom2}) {}
    void set_image(const NearestImage& im) {
      sym_idx = im.sym_idx;
      std::copy(std::begin(im.pbc_shift), std::end(im.pbc_shift), std::begin(pbc_shift));
    }
    void swap_atoms() {
      std::reverse(atoms.begin(), atoms.end());
      sym_idx = -sym_idx - 1;
    }
    bool same_asu() const {
      return sym_idx == 0 && pbc_shift[0]==0 && pbc_shift[1]==0 && pbc_shift[2]==0;
    }
    double calc(const UnitCell& cell, double wvdw, GeomTarget* target, Reporting *reporting) const;
    int type = 0; // 1: vdw, 2: torsion, 3: hbond, 4: metal, 5: dummy-nondummy, 6: dummy-dummy
    double value; // critical distance
    double sigma;
    int sym_idx = 0; // if negative, atoms need to be swapped.
    std::array<int, 3> pbc_shift = {{0,0,0}};
    std::array<Atom*, 2> atoms;
  };
  struct Reporting {
    using bond_reporting_t = std::tuple<const Bond*, const Bond::Value*, double>;
    using angle_reporting_t = std::tuple<const Angle*, const Angle::Value*, double>;
    using torsion_reporting_t = std::tuple<const Torsion*, const Torsion::Value*, double>;
    using chiral_reporting_t = std::tuple<const Chirality*, double, double>; // delta, ideal
    using plane_reporting_t = std::tuple<const Plane*, std::vector<double>>;
    using stacking_reporting_t = std::tuple<const Stacking*, double, double, double>; // delta_angle, delta_dist1, delta_dist2
    using vdw_reporting_t = std::tuple<const Vdw*, double>;
    std::vector<bond_reporting_t> bonds;
    std::vector<angle_reporting_t> angles;
    std::vector<torsion_reporting_t> torsions;
    std::vector<chiral_reporting_t> chirs;
    std::vector<plane_reporting_t> planes;
    std::vector<stacking_reporting_t> stackings;
    std::vector<vdw_reporting_t> vdws;
  };
  Geometry(Structure& s) : st(s), bondindex(s.first_model()) {}
  void load_topo(const Topo& topo);
  void finalize_restraints(); // sort_restraints?
  void setup_vdw(const EnerLib& ener_lib, float max_distsq_for_adp);
  static Position apply_transform(const UnitCell& cell, int sym_idx, const std::array<int, 3>& pbc_shift, const Position &v) {
    FTransform ft = sym_idx == 0 ? FTransform({}) : cell.images[sym_idx-1];
    ft.vec += Vec3(pbc_shift);
    return Position(cell.orth.combine(ft).combine(cell.frac).apply(v));
  }
  static Transform get_transform(const UnitCell& cell, int sym_idx, const std::array<int, 3>& pbc_shift) {
    FTransform ft = sym_idx == 0 ? FTransform({}) : cell.images[sym_idx-1];
    ft.vec += Vec3(pbc_shift);
    return cell.orth.combine(ft).combine(cell.frac);
  }

  void setup_target(bool refine_xyz, int adp_mode);
  void clear_target() {
    target.target = 0.;
    std::fill(target.vn.begin(), target.vn.end(), 0.);
    std::fill(target.am.begin(), target.am.end(), 0.);
  }
  double calc(bool use_nucleus, bool check_only, double wbond, double wangle, double wtors,
              double wchir, double wplane, double wstack, double wvdw);
  double calc_adp_restraint(bool check_only, double sigma);

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;
  std::vector<Interval> intervals;
  std::vector<Harmonic> harmonics;
  std::vector<Special> specials;
  std::vector<Stacking> stackings;
  std::vector<Vdw> vdws;
  Structure& st;
  BondIndex bondindex;
  std::map<int, std::string> chemtypes;
  std::map<int, char> hbtypes; // hydrogen bond types that override ener_lib
  Reporting reporting;
  GeomTarget target;

private:
  void set_vdw_values(Geometry::Vdw &vdw, int d_1_2, const EnerLib& ener_lib) const;
};

inline void Geometry::load_topo(const Topo& topo) {
  auto add = [&](const Topo::Rule& rule, bool same_asu) {
               if (!same_asu && rule.rkind != Topo::RKind::Bond) return; // not supported
               if (rule.rkind == Topo::RKind::Bond) {
                 const Topo::Bond& t = topo.bonds[rule.index];
                 if (t.restr->esd <= 0) return;
                 bonds.emplace_back(t.atoms[0], t.atoms[1]);
                 bonds.back().values.emplace_back(t.restr->value, t.restr->esd,
                                                  t.restr->value_nucleus, t.restr->esd_nucleus);
                 if (!same_asu) {
                   NearestImage im = st.cell.find_nearest_image(t.atoms[0]->pos, t.atoms[1]->pos,
                                                                Asu::Different);
                     bonds.back().set_image(im);
                 }
               } else if (rule.rkind == Topo::RKind::Angle) {
                 const Topo::Angle& t = topo.angles[rule.index];
                 if (t.restr->esd <= 0) return;
                 angles.emplace_back(t.atoms[0], t.atoms[1], t.atoms[2]);
                 angles.back().values.emplace_back(t.restr->value, t.restr->esd);
               } else if (rule.rkind == Topo::RKind::Torsion) {
                 const Topo::Torsion& t = topo.torsions[rule.index];
                 if (t.restr->esd <= 0) return;
                 torsions.emplace_back(t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]);
                 torsions.back().values.emplace_back(t.restr->value, t.restr->esd, t.restr->period);
                 torsions.back().values.back().label = t.restr->label;
               } else if (rule.rkind == Topo::RKind::Chirality) {
                 const Topo::Chirality& t = topo.chirs[rule.index];
                 const auto val_sigma = topo.ideal_chiral_abs_volume_sigma(t);
                 if (val_sigma.second <= 0) return;
                 chirs.emplace_back(t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]);
                 chirs.back().value = val_sigma.first;
                 chirs.back().sigma = val_sigma.second;
                 chirs.back().sign = t.restr->sign;
               } else if (rule.rkind == Topo::RKind::Plane) {
                 const Topo::Plane& t = topo.planes[rule.index];
                 if (t.restr->esd <= 0) return;
                 planes.emplace_back(t.atoms);
                 planes.back().sigma = t.restr->esd;
                 planes.back().label = t.restr->label;
               }
             };

  auto test_link_r = [&topo](const Topo::Rule& rule, const std::string& link_id) {
                       if (rule.rkind != Topo::RKind::Torsion)
                         return true;
                       const Topo::Torsion& t = topo.torsions[rule.index];
                       if (t.restr->label.find("sp2_sp2") == 0)
                         return true;
                       return ((link_id == "TRANS"   || link_id != "CIS"  ||
                                link_id == "PTRANS"  || link_id != "PCIS" ||
                                link_id == "NMTRANS" || link_id != "NMCIS") &&
                               t.restr->label == "omega");
                     };

  for (const Topo::ChainInfo& chain_info : topo.chain_infos)
    for (const Topo::ResInfo& ri : chain_info.res_infos) {
      // 1. link related
      for (const Topo::Link& prev : ri.prev)
        if (!prev.link_rules.empty())
          for (const Topo::Rule& rule : prev.link_rules)
            if (test_link_r(rule, prev.link_id))
              add(rule, true);

      // 2. monomer related
      if (!ri.monomer_rules.empty())
        for (const Topo::Rule& rule : ri.monomer_rules) {
          if (rule.rkind != Topo::RKind::Torsion ||
              topo.torsions[rule.index].restr->label.find("sp2_sp2") == 0 ||
              (ri.orig_chemcomp && ri.orig_chemcomp->group == ChemComp::Group::Peptide &&
               topo.torsions[rule.index].restr->label.find("chi") == 0))
            add(rule, true);
        }

      // collect chem_types
      for (const Atom& atom : ri.res->atoms) {
        const ChemComp& cc = ri.get_final_chemcomp(atom.altloc);
        auto it = cc.find_atom(atom.name);
        if (it != cc.atoms.end())
          chemtypes.emplace(atom.serial, it->chem_type);
      }
    }

  for (const Topo::Link& extra : topo.extras)
    for (const Topo::Rule& rule : extra.link_rules)
      if (test_link_r(rule, extra.link_id))
        add(rule, extra.asu == Asu::Same);
}

inline void Geometry::finalize_restraints() {
  for (const auto& b : bonds)
    if (b.type < 2)
      bondindex.add_link(*b.atoms[0], *b.atoms[1], b.same_asu());

  // sort out type 0 or 1 bonds
  // type = 0: replace it
  // type = 1: add it
  // check type 2. remove if type 0 (or 1?) bonds defined

  for (auto& b : bonds)
    if (b.atoms[0]->serial > b.atoms[1]->serial)
      b.swap_atoms();
  if (bonds.size() > 1)
    std::stable_sort(bonds.begin(), bonds.end(),
                     [](const Bond& l, const Bond& r) { return l.sort_key() < r.sort_key(); });
  // remove duplicated type 0 bonds
  // remove type 2 bonds if bond/angles defined
  std::vector<size_t> to_remove;
  for (size_t i = 0; i < bonds.size(); ++i) {
    if (bonds[i].type == 2 && bondindex.graph_distance(*bonds[i].atoms[0], *bonds[i].atoms[1], bonds[i].same_asu()) < 3) {
      //std::cout << "remove type 2: " << bonds[i].atoms[0]->name << "-" <<  bonds[i].atoms[1]->name << "\n";
      to_remove.emplace_back(i);
    } else if (i < bonds.size() - 1 && bonds[i].sort_key() == bonds[i+1].sort_key()) {
      if (bonds[i+1].type > 0) // merge. overwrite if next is type 0.
        bonds[i+1].values.insert(bonds[i+1].values.end(), bonds[i].values.begin(), bonds[i].values.end());
      //std::cout << "remove/merge: " << bonds[i].atoms[0]->name << "-" <<  bonds[i].atoms[1]->name << " d="
      //          << bonds[i].atoms[0]->pos.dist(bonds[i].atoms[1]->pos)
      //          << " target0=" << bonds[i].values[0].value << "\n";
      to_remove.emplace_back(i);
    }
  }
  for (auto it = to_remove.rbegin(); it != to_remove.rend(); ++it)
    bonds.erase(bonds.begin() + (*it));

  // sort angles
  for (auto& t : angles)
    if (t.atoms[0]->serial > t.atoms[2]->serial)
      t.swap_atoms();
  if (angles.size() > 1) {
    std::stable_sort(angles.begin(), angles.end(),
                     [](const Angle& l, const Angle& r) { return l.sort_key() < r.sort_key(); });
    // remove duplicated angles
    to_remove.clear();
    for (size_t i = 0; i < angles.size() - 1; ++i)
      if (angles[i].sort_key() == angles[i+1].sort_key()) {
        if (angles[i+1].type > 0) // should we always do this?
          angles[i+1].values.insert(angles[i+1].values.end(), angles[i].values.begin(), angles[i].values.end());
        to_remove.emplace_back(i);
      }
    for (auto it = to_remove.rbegin(); it != to_remove.rend(); ++it)
      angles.erase(angles.begin() + (*it));
  }

  // sort torsions
  for (auto& t : torsions)
    if (t.atoms[0]->serial > t.atoms[3]->serial)
      t.swap_atoms();
  if (torsions.size() > 1) {
    std::stable_sort(torsions.begin(), torsions.end(),
                     [](const Torsion& l, const Torsion& r) { return l.sort_key() < r.sort_key(); });
    // remove duplicated torsions
    to_remove.clear();
    for (size_t i = 0; i < torsions.size() - 1; ++i)
      if (torsions[i].sort_key() == torsions[i+1].sort_key()) {
        if (torsions[i+1].type > 0)
          torsions[i+1].values.insert(torsions[i+1].values.end(), torsions[i].values.begin(), torsions[i].values.end());
        to_remove.emplace_back(i);
      }
    for (auto it = to_remove.rbegin(); it != to_remove.rend(); ++it)
      torsions.erase(torsions.begin() + (*it));
  }

  // no care needed for others?
}

inline void Geometry::set_vdw_values(Geometry::Vdw &vdw, int d_1_2, const EnerLib& ener_lib) const {
  //double wvskal          = 1.00;
  double vdw_sdi_vdw     = 0.2; // VDWR SIGM VDW val
  double vdw_sdi_torsion = 0.2; // VDWR SIGM TORS val
  double vdw_sdi_hbond   = 0.2; // VDWR SIGM HBON val
  double vdw_sdi_metal   = 0.2; // VDWR SIGM META val
  double hbond_dinc_ad   = -0.3; // VDWR INCR ADHB val
  double hbond_dinc_ah   = 0.1; // VDWR INCR AHHB val
  //double dinc_torsion    = -0.3; // not used? // // VDWR INCR TORS val
  double dinc_torsion_o  = -0.1;
  double dinc_torsion_n  = -0.1;
  double dinc_torsion_c  = -0.15; // VDWR INCR TORS val (copied)
  double dinc_torsion_all= -0.15; // VDWR INCR TORS val (copied)
  double dinc_dummy      = -0.7; // VDWR INCR DUMM val
  double vdw_sdi_dummy   = 0.3; // VDWR SIGM DUMM val
  //double dvdw_cut_min    = 1.75; // no need? // VDWR VDWC val
  //double dvdw_cut_min_x  = 1.75; // used as twice in fast_hessian_tabulation.f // VDWR VDWC val

  double vdw_rad[2];
  double ion_rad[2];
  char hb_type[2];
  for (int i = 0; i < 2; ++i) {
    const std::string& chem_type = chemtypes.at(vdw.atoms[i]->serial);
    const auto& libatom = ener_lib.atoms.at(chem_type);
    vdw_rad[i] = std::isnan(libatom.vdwh_radius) ? libatom.vdw_radius : libatom.vdwh_radius; // XXX needs switch. check hydrogen is there?
    ion_rad[i] = libatom.ion_radius;
    auto it = hbtypes.find(vdw.atoms[i]->serial);
    hb_type[i] = (it == hbtypes.end()) ? libatom.hb_type : it->second;
  }

  // check torsion related atoms XXX what if within ring? we can remove if torsion.period<3?
  if (d_1_2 == 3) { // for hydrogen also??
    double dinc_curr[2];
    for (int i = 0; i < 2; ++i) {
      switch (vdw.atoms[i]->element) {
      case El::O: dinc_curr[i] = dinc_torsion_o; break;
      case El::N: dinc_curr[i] = dinc_torsion_n; break;
      case El::C: dinc_curr[i] = dinc_torsion_c; break;
      default:    dinc_curr[i] = dinc_torsion_all;
      }
    }
    vdw.type = 2;
    vdw.value = vdw_rad[0] + vdw_rad[1] + dinc_curr[0] + dinc_curr[1];
    vdw.sigma = vdw_sdi_torsion;
    return;
  }

  // check hydrogen bond
  if ((hb_type[0] == 'A' && (hb_type[1] == 'D' || hb_type[1] == 'B')) ||
      (hb_type[0] == 'D' && (hb_type[1] == 'A' || hb_type[1] == 'B')) ||
      (hb_type[0] == 'B' && (hb_type[1] == 'A' || hb_type[1] == 'D' || hb_type[1] == 'B'))) {
    vdw.value = vdw_rad[0] + vdw_rad[1] + hbond_dinc_ad;
    vdw.type = 3;
  }
  else if ((hb_type[0] == 'A' && hb_type[1] == 'H') || // XXX 'H' type must be set beforehand.
           (hb_type[0] == 'B' && hb_type[1] == 'H')) {
    vdw.value = vdw_rad[0] + hbond_dinc_ah;
    vdw.type = 3;
  }
  else if (hb_type[0] == 'H' && (hb_type[1] == 'A' || hb_type[1] == 'B')) {
    vdw.value = vdw_rad[1] + hbond_dinc_ah;
    vdw.type = 3;
  }
  if (vdw.type == 3) {
    vdw.sigma = vdw_sdi_hbond;
    return;
  }

  // check metal bond?
  if (vdw.atoms[0]->element.is_metal() || vdw.atoms[1]->element.is_metal()) { // XXX should be xor?
    vdw.value = ion_rad[0] + ion_rad[1];
    vdw.type = 4;
    vdw.sigma = vdw_sdi_metal;
    if (!std::isnan(vdw.value))
      return;
  }

  // check dummy XXX we should not depend on atom names?
  bool is_dum_1 = vdw.atoms[0]->name.rfind("DUM", 0) == 0;
  bool is_dum_2 = vdw.atoms[1]->name.rfind("DUM", 0) == 0;
  if ((is_dum_1 && !is_dum_2) || (!is_dum_1 && is_dum_2)) {
    vdw.value = std::max(0.7, vdw_rad[0] + vdw_rad[1] + dinc_dummy);
    vdw.type = 5;
    vdw.sigma = vdw_sdi_dummy;
    return;
  }
  if (is_dum_1 && is_dum_2) {
    vdw.value = vdw_rad[0] + vdw_rad[1];
    vdw.type = 6;
    vdw.sigma = vdw_sdi_dummy;
    return;
  }

  // otherwise
  vdw.value = vdw_rad[0] + vdw_rad[1];
  vdw.type = 1;
  vdw.sigma = vdw_sdi_vdw;
}

inline void Geometry::setup_vdw(const EnerLib& ener_lib, float max_distsq_for_adp) {
  // set hbtypes for hydrogen
  if (hbtypes.empty()) {
    for (auto& b : bonds)
      if (b.atoms[0]->is_hydrogen() != b.atoms[1]->is_hydrogen()) {
        int p = b.atoms[0]->is_hydrogen() ? 1 : 0; // parent
        int h = b.atoms[0]->is_hydrogen() ? 0 : 1; // hydrogen
        const std::string& p_chem_type = chemtypes.at(b.atoms[p]->serial);
        const char p_hb_type = ener_lib.atoms.at(p_chem_type).hb_type;
        hbtypes.emplace(b.atoms[h]->serial, p_hb_type == 'D' || p_hb_type == 'B' ? 'H' : 'N');
      }
  }

  vdws.clear();

  // Reference: Refmac vdw_and_contacts.f
  NeighborSearch ns(st.first_model(), st.cell, 4);
  ns.populate();
  float max_vdwr = 4.; // FIXME
  ContactSearch contacts(std::max(std::sqrt(max_distsq_for_adp), max_vdwr * 2));
  contacts.ignore = ContactSearch::Ignore::Nothing;
  contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                    int sym_idx, float dist_sq) {
    // XXX Refmac uses intervals for distances as well? vdw_and_contacts.f remove_bonds_and_angles()
    NearestImage im = st.cell.find_nearest_pbc_image(cra1.atom->pos, cra2.atom->pos, sym_idx);
    int d_1_2 = bondindex.graph_distance(*cra1.atom, *cra2.atom, im.sym_idx == 0 && im.same_asu());
    if (d_1_2 > 2) {
      vdws.emplace_back(cra1.atom, cra2.atom);
      set_vdw_values(vdws.back(), d_1_2, ener_lib);
      assert(!std::isnan(vdws.back().value) && vdws.back().value > 0);
      if (dist_sq > max_distsq_for_adp // quick trick to use ADP restraints. Need better way.
          && vdws.back().value * vdws.back().value < dist_sq * 0.9) {
        vdws.pop_back();
        return;
      }
      vdws.back().set_image(im);
      if (im.sym_idx != 0 || !im.same_asu())
        vdws.back().type += 6;
    }
  });
}

inline void Geometry::setup_target(bool refine_xyz, int adp_mode) {
  std::vector<std::pair<int,int>> tmp;
  for (const auto &t : bonds)
    tmp.emplace_back(t.atoms[0]->serial-1, t.atoms[1]->serial-1);

  for (const auto &t : angles)
    for (int i = 0; i < 2; ++i)
      for (int j = i+1; j < 3; ++j)
        tmp.emplace_back(t.atoms[i]->serial-1, t.atoms[j]->serial-1);

  for (const auto &t : torsions)
    for (int i = 0; i < 3; ++i)
      for (int j = i+1; j < 4; ++j)
        tmp.emplace_back(t.atoms[i]->serial-1, t.atoms[j]->serial-1);

  for (const auto &t : chirs)
    for (int i = 0; i < 3; ++i)
      for (int j = i+1; j < 4; ++j)
        tmp.emplace_back(t.atoms[i]->serial-1, t.atoms[j]->serial-1);

  for (const auto &t : planes)
    for (size_t i = 1; i < t.atoms.size(); ++i)
      for (size_t j = 0; j < i; ++j)
        tmp.emplace_back(t.atoms[i]->serial-1, t.atoms[j]->serial-1);

  for (const auto &t : stackings) {
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 1; j < t.planes[i].size(); ++j)
        for (size_t k = 0; k < j; ++k)
          tmp.emplace_back(t.planes[i][j]->serial-1, t.planes[i][k]->serial-1);

    for (size_t j = 0; j < t.planes[0].size(); ++j)
      for (size_t k = 0; k < t.planes[1].size(); ++k)
        tmp.emplace_back(t.planes[0][j]->serial-1, t.planes[1][k]->serial-1);
  }

  for (const auto &t : vdws)
    tmp.push_back({t.atoms[0]->serial-1, t.atoms[1]->serial-1});

  // sort_and_compress_distances
  for (auto &p : tmp)
    if (p.first > p.second)
      std::swap(p.first, p.second);

  target.pairs.clear();
  if (!tmp.empty()) {
    std::sort(tmp.begin(), tmp.end());
    target.pairs.push_back(tmp[0]);
    for (size_t i = 1; i < tmp.size(); ++i)
      if (tmp[i] != target.pairs.back() && tmp[i].first != tmp[i].second)
        target.pairs.push_back(tmp[i]); // n_target, n_object
  }

  target.setup(st.first_model(), refine_xyz, adp_mode);
}

inline double Geometry::calc(bool use_nucleus, bool check_only,
                             double wbond, double wangle, double wtors,
                             double wchir, double wplane, double wstack,
                             double wvdw) {
  if (check_only)
    reporting = {};
  else
    assert(target.refine_xyz); // otherwise vector and matrix not ready

  auto target_ptr = check_only ? nullptr : &target;
  auto rep_ptr = check_only ? &reporting : nullptr;
  double ret = 0.;

  for (const auto &t : bonds)
    ret += t.calc(st.cell, use_nucleus, wbond, target_ptr, rep_ptr);
  for (const auto &t : angles)
    ret += t.calc(wangle, target_ptr, rep_ptr);
  for (const auto &t : torsions)
    ret += t.calc(wtors, target_ptr, rep_ptr);
  for (const auto &t : chirs)
    ret += t.calc(wchir, target_ptr, rep_ptr);
  for (const auto &t : planes)
    ret += t.calc(wplane, target_ptr, rep_ptr);
  for (const auto &t : stackings)
    ret += t.calc(wstack, target_ptr, rep_ptr);
  for (const auto &t : vdws)
    ret += t.calc(st.cell, wvdw, target_ptr, rep_ptr);

  // TODO intervals, harmonics, specials
  return ret;
}

inline double Geometry::calc_adp_restraint(bool check_only, double sigma) {
  assert(target.adp_mode > 0);
  const int n_pairs = target.pairs.size();
  const int offset_v = target.refine_xyz ? target.n_atoms() * 3 : 0;
  const int offset_a = target.refine_xyz ? target.n_atoms() * 6 + n_pairs * 9 : 0;
  const double weight = 1. / (sigma * sigma);
  double ret = 0.;
  for (int i = 0; i < n_pairs; ++i) {
    const Atom* atom1 = target.atoms[target.pairs[i].first];
    const Atom* atom2 = target.atoms[target.pairs[i].second];
    // calculate minimum distance - expensive?
    const NearestImage im = st.cell.find_nearest_image(atom1->pos, atom2->pos, Asu::Any);
    const double dsq = im.dist_sq;
    const double w_fac = std::exp(-dsq * dsq / (2*4*4*4*std::log(2.))); // arbitrary weighting factor, adjust so that d=4 gives 0.5
    const double w = weight * w_fac;
    if (target.adp_mode == 1) {
      const double f = 0.5 * w * (atom1->b_iso - atom2->b_iso) * (atom1->b_iso - atom2->b_iso);
      ret += f;
      if (!check_only) {
        target.target += f;
        const double df1 = w * (atom1->b_iso - atom2->b_iso);
        target.vn[offset_v + atom1->serial - 1] += df1;
        target.vn[offset_v + atom2->serial - 1] += -df1;
        // diagonal
        target.am[offset_a + atom1->serial - 1] += w;
        target.am[offset_a + atom2->serial - 1] += w;
        // non-diagonal
        target.am[offset_a + target.n_atoms() + i] += -w;
      }
    } else if (target.adp_mode == 2) {
      const auto& a1 = atom1->aniso.scaled(u_to_b()).elements_pdb();
      const auto& a2 = atom2->aniso.scaled(u_to_b()).elements_pdb();
      for (int j = 0; j < 6; ++j) {
        const double f = 0.5 * w * (a1[j] - a2[j]) * (a1[j] - a2[j]);
        ret += f;
        if (!check_only) {
          target.target += f;
          const double df1 = w * (a1[j] - a2[j]);
          target.vn[offset_v + 6 * (atom1->serial-1) + j] += df1;
          target.vn[offset_v + 6 * (atom2->serial-1) + j] += -df1;
          // diagonal block (6 x 6 symmetric)
          target.am[offset_a + 21 * (atom1->serial-1) + j] += w;
          target.am[offset_a + 21 * (atom2->serial-1) + j] += w;
          // non-diagonal block (6 x 6)
          target.am[offset_a + 21 * target.n_atoms() + 36 * i + 6 * j] += -w;
        }
      }
    }
  }
  return ret;
}

inline double Geometry::Bond::calc(const UnitCell& cell, bool use_nucleus, double wdskal,
                                   GeomTarget* target, Reporting *reporting) const {
  assert(!values.empty());
  const bool swapped = sym_idx < 0;
  const Atom* atom1 = atoms[swapped ? 1 : 0];
  const Atom* atom2 = atoms[swapped ? 0 : 1];
  const Transform tr = get_transform(cell, swapped ? -sym_idx - 1 : sym_idx, pbc_shift);
  const Position& x1 = atom1->pos;
  const Position& x2 = same_asu() ? atom2->pos : Position(tr.apply(atom2->pos));
  const double b = x1.dist(x2);
  auto closest = find_closest_value(b, use_nucleus);
  const double ideal = use_nucleus ? closest->value_nucleus : closest->value;
  const double db = b - ideal;
  const double sigma = (use_nucleus ? closest->sigma_nucleus : closest->sigma);
  const double weight = wdskal / sigma;
  const double y = db * weight;
  double ret, dfdy, d2fdy;

  if (type < 2 || std::abs(alpha - 2) < 1e-3) { // least square
    ret = 0.5 * y * y;
    dfdy = y;
    d2fdy = 1.0;
  } else if (std::abs(alpha) < 1e-3) { // cauchy or lorentz
    ret = std::log(0.5 * y * y + 1.0);
    dfdy = y / (0.5 * y * y + 1.0);
    d2fdy = 1.0 / (0.5 * y * y + 1.0);
  } else if (alpha < -1000) { // -inf. welch
    const double expy = std::exp(-0.5 * y * y);
    ret = 1.0 - expy;
    dfdy = y * expy;
    d2fdy = expy;
  } else { // other alpha
    const double alpha2 = std::abs(alpha - 2.0);
    ret = alpha2 / alpha * (std::pow(y * y / alpha2 + 1, 0.5 * alpha) - 1.0);
    dfdy = y * std::pow(y * y / alpha2 + 1, 0.5 * alpha - 1.0);
    d2fdy = std::pow(y * y / alpha2 + 1, 0.5 * alpha - 1.0);
  }

  // note that second derivative is not exact in some alpha
  if (target != nullptr) {
    const Position dydx1 = weight * (x1 - x2) / std::max(b, 0.02);
    const Position dydx2 = same_asu() ? -dydx1 : Position(tr.mat.transpose().multiply(-dydx1));
    const int ia1 = atom1->serial - 1;
    const int ia2 = atom2->serial - 1;
    target->incr_vn(ia1 * 3, dfdy, dydx1);
    target->incr_vn(ia2 * 3, dfdy, dydx2);
    target->incr_am_diag(ia1 * 6, d2fdy, dydx1);
    target->incr_am_diag(ia2 * 6, d2fdy, dydx2);

    if (ia1 != ia2) {
      auto mp = target->find_restraint(ia1, ia2);
      if (mp.imode == 0)
        target->incr_am_ndiag(mp.ipos, d2fdy, dydx1, dydx2);
      else
        target->incr_am_ndiag(mp.ipos, d2fdy, dydx2, dydx1);
    } else
      target->incr_am_diag12(ia1 * 6, d2fdy, dydx1, dydx2);

    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->bonds.emplace_back(this, closest, db);
  return ret;
}

inline double Geometry::Angle::calc(double waskal, GeomTarget* target, Reporting *reporting) const {
  const Position& x1 = atoms[0]->pos;
  const Position& x2 = atoms[1]->pos;
  const Position& x3 = atoms[2]->pos;
  int ia[3];
  for (int i = 0; i < 3; ++i) ia[i] = atoms[i]->serial - 1;
  const Position v1 = x2 - x1;
  const Position v2 = x2 - x3;
  const double v1n = std::max(v1.length(), 0.02);
  const double v2n = std::max(v2.length(), 0.02);
  const double v12 = v1.dot(v2);
  const double cosa = std::min(1., v12 / v1n / v2n);
  const double sina = std::min(1., std::max(std::sqrt(1 - cosa * cosa), 0.1));
  const double a = deg(std::acos(std::max(-1., std::min(1., cosa))));
  auto closest = find_closest_value(a);
  const double da = a - closest->value;
  const double weight = waskal * waskal / (closest->sigma * closest->sigma);
  const double ret = da * da * weight * 0.5;
  if (target != nullptr) {
    Vec3 dadx[3];
    dadx[0] = ((v2 / (v1n * v2n) - v1 * cosa / (v1n * v1n)) / sina) * deg(1);
    dadx[2] = ((v1 / (v1n * v2n) - v2 * cosa / (v2n * v2n)) / sina) * deg(1);
    dadx[1] = -dadx[0] - dadx[2];

    for(int i = 0; i < 3; ++i) {
      target->incr_vn(ia[i] * 3, weight * da, dadx[i]);
      target->incr_am_diag(ia[i] * 6, weight, dadx[i]);
    }

    for (int i = 0; i < 2; ++i)
      for (int j = i+1; j < 3; ++j) {
        auto mp = target->find_restraint(ia[i], ia[j]);
        if (mp.imode == 0) // ia[i] > ia[j]
          target->incr_am_ndiag(mp.ipos, weight, dadx[i], dadx[j]);
        else
          target->incr_am_ndiag(mp.ipos, weight, dadx[j], dadx[i]);
      }
    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->angles.emplace_back(this, closest, da);
  return ret;
}

inline double Geometry::Torsion::calc(double wtskal, GeomTarget* target, Reporting *reporting) const {
  const Position& x1 = atoms[0]->pos;
  const Position& x2 = atoms[1]->pos;
  const Position& x3 = atoms[2]->pos;
  const Position& x4 = atoms[3]->pos;
  int ia[4];
  for (int i = 0; i < 4; ++i) ia[i] = atoms[i]->serial - 1;
  const Vec3 u = x1 - x2;
  const Vec3 v = x4 - x3;
  const Vec3 w = x3 - x2;
  const Vec3 a = u.cross(w);
  const Vec3 b = v.cross(w);
  const double s = a.cross(b).dot(w);
  const double wl = std::max(0.0001, w.length());
  const double t = wl * a.dot(b);
  const double theta = deg(std::atan2(s, t));
  auto closest = find_closest_value(theta);
  const int period = std::max(1, closest->period);
  const double weight = wtskal * wtskal / (closest->sigma * closest->sigma);
  const double dtheta1 = rad(period * (theta - closest->value));
  const double dtheta2 = deg(std::atan2(std::sin(dtheta1), std::cos(dtheta1)));
  const double dtheta = dtheta2 / period;
  const double ret = dtheta * dtheta * weight * 0.5;

  if (target != nullptr) {
    const double denom = rad(std::max(0.0001, s * s + t * t));
    Vec3 dadx[3][3], dbdx[3][3], dwdx[3][2];
    double dwldx[3][2];
    for (int i = 0; i < 3; ++i) {
      Vec3 drdx; drdx.at(i) = 1.;
      const Vec3 d1 = drdx.cross(w);
      const Vec3 d2 = u.cross(drdx);
      const Vec3 d3 = v.cross(drdx);
      dadx[i][0] = d1;     // da/dx1
      dadx[i][1] = -d1-d2; // da/dx2
      dadx[i][2] = d2;     // da/dx3
      dbdx[i][0] = -d3;    // db/dx2
      dbdx[i][1] = d3-d1;  // db/dx3
      dbdx[i][2] = d1;     // db/dx4
      dwdx[i][0] = -drdx;  // dw/dx2
      dwdx[i][1] = drdx;   // dw/dx3
      dwldx[i][1] = w.dot(drdx)/wl; // dwl/dx3
      dwldx[i][0] = -dwldx[i][1];   // dwl/dx2
    }
    Vec3 dthdx[4];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
        double dsdx = 0.;
        double dtdx = 0.;
        if (j != 3) { // only for x1,x2,x3
          dsdx = dadx[i][j].cross(b).dot(w);
          dtdx = dadx[i][j].dot(b) * wl;
          if (j == 0) { // only for x1
            dthdx[j].at(i) = (t * dsdx - s * dtdx) / denom;
            continue;
          }
        }
        // only for x2,x3,x4
        dsdx += a.cross(dbdx[i][j-1]).dot(w);
        dtdx += a.dot(dbdx[i][j-1]) * wl;
        if (j != 3) { // only for x2,x3
          dsdx += a.cross(b).dot(dwdx[i][j-1]);
          dtdx += t / wl * dwldx[i][j-1];
        }
        dthdx[j].at(i) = (t * dsdx - s * dtdx)/denom;
      }
    }

    for(int i = 0; i < 4; ++i) {
      target->incr_vn(ia[i] * 3, dtheta * weight, dthdx[i]);
      target->incr_am_diag(ia[i] * 6, weight, dthdx[i]);
    }

    for (int i = 0; i < 3; ++i)
      for (int j = i+1; j < 4; ++j) {
        auto mp = target->find_restraint(ia[i], ia[j]);
        if (mp.imode == 0)
          target->incr_am_ndiag(mp.ipos, weight, dthdx[i], dthdx[j]);
        else
          target->incr_am_ndiag(mp.ipos, weight, dthdx[j], dthdx[i]);
      }
    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->torsions.emplace_back(this, closest, dtheta);
  return ret;
}

inline double Geometry::Chirality::calc(double wchiral, GeomTarget* target, Reporting *reporting) const {
  const double weight = wchiral * wchiral / (sigma * sigma);
  const Position& xc = atoms[0]->pos;
  const Position& x1 = atoms[1]->pos;
  const Position& x2 = atoms[2]->pos;
  const Position& x3 = atoms[3]->pos;
  int ia[4];
  for (int i = 0; i < 4; ++i) ia[i] = atoms[i]->serial - 1;
  const Vec3 a1 = x1 - xc;
  const Vec3 a2 = x2 - xc;
  const Vec3 a3 = x3 - xc;
  const Vec3 a1xa2 = a1.cross(a2);
  const double v = a1xa2.dot(a3);
  const bool isneg = (sign == ChiralityType::Negative || (sign == ChiralityType::Both && v < 0));
  const double ideal = (isneg ? -1 : 1) * value;
  const double dv = v - ideal;
  const double ret = dv * dv * weight * 0.5;

  if (target != nullptr) {
    Vec3 dcdx[4];
    for (int i = 0; i < 3; ++i) {
      Vec3 drdx; drdx.at(i) = 1.;
      dcdx[1].at(i) = drdx.cross(a2).dot(a3); // atom1
      dcdx[2].at(i) = a1.cross(drdx).dot(a3); // atom2
      dcdx[3].at(i) = a1xa2.dot(drdx);        // atom3
      dcdx[0].at(i) = -dcdx[1].at(i) - dcdx[2].at(i) - dcdx[3].at(i); //atomc
    }

    for(int i = 0; i < 4; ++i) {
      target->incr_vn(ia[i] * 3, dv * weight, dcdx[i]);
      target->incr_am_diag(ia[i] * 6, weight, dcdx[i]);
    }

    for (int i = 0; i < 3; ++i)
      for (int j = i+1; j < 4; ++j) {
        auto mp = target->find_restraint(ia[i], ia[j]);
        if (mp.imode == 0)
          target->incr_am_ndiag(mp.ipos, weight, dcdx[i], dcdx[j]);
        else
          target->incr_am_ndiag(mp.ipos, weight, dcdx[j], dcdx[i]);
      }
    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->chirs.emplace_back(this, dv, ideal);
  return ret;
}

inline double Geometry::Plane::calc(double wplane, GeomTarget* target, Reporting *reporting) const {
  const double weight = wplane * wplane / (sigma * sigma);
  const int natoms = atoms.size();
  const PlaneDeriv pder(atoms);

  double ret = 0.;
  std::vector<double> deltas(natoms);
  for (int j = 0; j < natoms; ++j) {
    deltas[j] = pder.D - pder.vm.dot(atoms[j]->pos);
    ret += deltas[j] * deltas[j] * weight * 0.5;
  }

  if (target != nullptr) {
    for (int j = 0; j < natoms; ++j) {
      const Position &xj = atoms[j]->pos;
      for (int l = 0; l < natoms; ++l) {
        Position dpdx1;
        for (int m = 0; m < 3; ++m)
          dpdx1.at(m) = pder.dDdx[l].at(m) - xj.dot(pder.dvmdx[l][m]) - (j==l ? pder.vm.at(m) : 0);

        target->incr_vn((atoms[l]->serial-1) * 3, deltas[j] * weight, dpdx1);

        for (int k = l; k < natoms; ++k) {
          Position dpdx2;
          for (int m = 0; m < 3; ++m)
            dpdx2.at(m) = pder.dDdx[k].at(m) - xj.dot(pder.dvmdx[k][m]) - (k==j ? pder.vm.at(m) : 0);

          if (k == l)
            target->incr_am_diag((atoms[l]->serial-1) * 6, weight, dpdx1);
          else {
            auto mp = target->find_restraint(atoms[l]->serial-1, atoms[k]->serial-1);
            if (mp.imode == 0)
              target->incr_am_ndiag(mp.ipos, weight, dpdx1, dpdx2);
            else
              target->incr_am_ndiag(mp.ipos, weight, dpdx2, dpdx1);
          }
        }
      }
    }
    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->planes.emplace_back(this, deltas);
  return ret;
}

inline double Geometry::Stacking::calc(double wstack, GeomTarget* target, Reporting *reporting) const {
  double ret = 0;
  PlaneDeriv pder[2] = {planes[0], planes[1]};
  double vm1vm2 = pder[0].vm.dot(pder[1].vm);
  if (vm1vm2 < 0) {
    pder[1].flip();
    vm1vm2 *= -1;
  }

  // angle
  const double wa = wstack * wstack / (sd_angle * sd_angle);
  const double cosa = std::min(1., vm1vm2);
  const double a = deg(std::acos(std::max(-1., std::min(1., cosa))));
  const double deltaa = a - angle;
  const double deltaa2 = deltaa * deltaa;
  ret += 0.5 * wa * deltaa2;
  if (target != nullptr) {
    const double inv_sina = 1. / std::min(1., std::max(std::sqrt(1 - cosa * cosa), 0.1));
    std::vector<std::vector<Vec3>> dpdx;
    for (size_t i = 0; i < 2; ++i) { // plane index
      dpdx.emplace_back(planes[i].size());
      for (size_t j = 0; j < planes[i].size(); ++j) { // atom index of plane i
        for (size_t m = 0; m < 3; ++m)
          dpdx[i][j].at(m) = -deg(1) * pder[i].dvmdx[j][m].dot(pder[1-i].vm) * inv_sina;
        target->incr_vn((planes[i][j]->serial-1) * 3, wa * deltaa, dpdx[i][j]);

        // second derivatives in the same plane
        for (size_t k = 0; k <= j; ++k) {
          if (k == j)
            target->incr_am_diag((planes[i][j]->serial-1) * 6, wa, dpdx[i][j]);
          else {
            auto mp = target->find_restraint(planes[i][j]->serial-1, planes[i][k]->serial-1);
            if (mp.imode == 0)
              target->incr_am_ndiag(mp.ipos, wa, dpdx[i][j], dpdx[i][k]);
            else
              target->incr_am_ndiag(mp.ipos, wa, dpdx[i][k], dpdx[i][j]);
          }
        }
      }
    }
    // second derivatives between two planes
    for (size_t j = 0; j < planes[0].size(); ++j)
      for (size_t k = 0; k < planes[1].size(); ++k) {
        auto mp = target->find_restraint(planes[0][j]->serial-1, planes[1][k]->serial-1);
        if (mp.imode == 0)
          target->incr_am_ndiag(mp.ipos, wa, dpdx[0][j], dpdx[1][k]);
        else
          target->incr_am_ndiag(mp.ipos, wa, dpdx[1][k], dpdx[0][j]);
      }
  }

  // distance; turned off by default in Refmac
  double deltad[2] = {0, 0};
  if (dist > 0) { // skip if ideal dist < 0
    const double wd = wstack * wstack / (sd_dist * sd_dist);
    for (size_t i = 0; i < 2; ++i) {
      double d = pder[i].xs.dot(pder[1-i].vm) - pder[1-i].D; // distance from i to the other
      if (d < 0) {
        d *= -1;
        pder[1-i].flip();
      }
      deltad[i] = d - dist;
      const double deltad2 = deltad[i] * deltad[i];
      ret += 0.5 * wd * deltad2; // distance between planes is not symmetric, so we add both
    }
    if (target != nullptr) {
      std::vector<std::vector<Vec3>> dpdx;
      dpdx.emplace_back(planes[0].size());
      dpdx.emplace_back(planes[1].size());
      for (size_t i = 0; i < 2; ++i) {
        // for the atoms of this plane
        for (size_t j = 0; j < planes[i].size(); ++j) {
          dpdx[i][j] = pder[1-i].vm / planes[i].size();
          target->incr_vn((planes[i][j]->serial-1) * 3, wd * deltad[i], dpdx[i][j]);
          // second derivatives
          for (size_t k = 0; k <= j; ++k) {
            if (k == j)
              target->incr_am_diag((planes[i][j]->serial-1) * 6, wd, dpdx[i][j]);
            else {
              auto mp = target->find_restraint(planes[i][j]->serial-1, planes[i][k]->serial-1);
              if (mp.imode == 0)
                target->incr_am_ndiag(mp.ipos, wd, dpdx[i][j], dpdx[i][k]);
              else
                target->incr_am_ndiag(mp.ipos, wd, dpdx[i][k], dpdx[i][j]);
            }
          }
        }
        // for the atoms of the other plane
        for (size_t j = 0; j < planes[1-i].size(); ++j) {
          for (size_t m = 0; m < 3; ++m)
            dpdx[1-i][j].at(m) = pder[1-i].dvmdx[j][m].dot(pder[i].xs) - pder[1-i].dDdx[j].at(m);
          target->incr_vn((planes[1-i][j]->serial-1) * 3, wd * deltad[i], dpdx[1-i][j]);
          // second derivatives
          for (size_t k = 0; k <= j; ++k) {
            if (k == j)
              target->incr_am_diag((planes[1-i][j]->serial-1) * 6, wd, dpdx[1-i][j]);
            else {
              auto mp = target->find_restraint(planes[1-i][j]->serial-1, planes[1-i][k]->serial-1);
              if (mp.imode == 0)
                target->incr_am_ndiag(mp.ipos, wd, dpdx[1-i][j], dpdx[1-i][k]);
              else
                target->incr_am_ndiag(mp.ipos, wd, dpdx[1-i][k], dpdx[1-i][j]);
            }
          }
        }
        // second derivatives between two planes
        for (size_t j = 0; j < planes[0].size(); ++j)
          for (size_t k = 0; k < planes[1].size(); ++k) {
            auto mp = target->find_restraint(planes[0][j]->serial-1, planes[1][k]->serial-1);
            if (mp.imode == 0)
              target->incr_am_ndiag(mp.ipos, wd, dpdx[0][j], dpdx[1][k]);
            else
              target->incr_am_ndiag(mp.ipos, wd, dpdx[1][k], dpdx[0][j]);
          }
      }
    }
  }
  if (target != nullptr)
    target->target += ret;
  if (reporting != nullptr)
    reporting->stackings.emplace_back(this, deltaa, deltad[0], deltad[1]);
  return ret;
}


inline double
Geometry::Vdw::calc(const UnitCell& cell, double wvdw, GeomTarget* target, Reporting *reporting) const {
  const double weight = wvdw * wvdw / (sigma * sigma);
  const bool swapped = sym_idx < 0;
  const Atom& atom1 = *atoms[swapped ? 1 : 0];
  const Atom& atom2 = *atoms[swapped ? 0 : 1];
  const int ia1 = atom1.serial - 1;
  const int ia2 = atom2.serial - 1;
  const Transform tr = get_transform(cell, swapped ? -sym_idx - 1 : sym_idx, pbc_shift);
  const Position& x1 = atom1.pos;
  const Position& x2 = same_asu() ? atom2.pos : Position(tr.apply(atom2.pos));
  const double b = x1.dist(x2);
  const double db = b - value;
  if (db > 0)
    return 0.;

  const double ret = db * db * weight * 0.5;
  if (target != nullptr) {
    const Position dbdx1 = (x1 - x2) / std::max(b, 0.02);
    const Position dbdx2 = same_asu() ? -dbdx1 : Position(tr.mat.transpose().multiply(-dbdx1));
    target->incr_vn(ia1 * 3, weight * db, dbdx1);
    target->incr_vn(ia2 * 3, weight * db, dbdx2);
    target->incr_am_diag(ia1 * 6, weight, dbdx1);
    target->incr_am_diag(ia2 * 6, weight, dbdx2);

    if (ia1 != ia2) {
      auto mp = target->find_restraint(ia1, ia2);
      if (mp.imode == 0)
        target->incr_am_ndiag(mp.ipos, weight, dbdx1, dbdx2);
      else
        target->incr_am_ndiag(mp.ipos, weight, dbdx2, dbdx1);
    } else
      target->incr_am_diag12(ia1*6, weight, dbdx1, dbdx2);

    target->target += ret;
  }
  if (reporting != nullptr)
    reporting->vdws.emplace_back(this, db);
  return ret;
}

} // namespace gemmi
#endif
