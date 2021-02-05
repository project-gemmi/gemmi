// Copyright 2020 Global Phasing Ltd.
//

#ifndef GEMMI_BULKSOL_HPP_
#define GEMMI_BULKSOL_HPP_

#include "asudata.hpp"
#include "levmar.hpp"

namespace gemmi {

template<typename Real>
struct BulkSolvent {
  struct Point {
    Miller hkl;
    double stol2;
    std::complex<Real> fcmol, fmask;
    Real fobs, sigma;

    Miller get_x() const { return hkl; }
    double get_y() const { return fobs; }
    double get_weight() const { return 1.0 / sigma; }
  };

  UnitCell cell;
  CrystalSystem crystal_system;
  // model parameters
  double k_overall = 1.;
  // b_star = F B_cart F^T, where F - fractionalization matrix
  SMat33<double> b_star{0, 0, 0, 0, 0, 0};
  bool use_solvent = false;
  // initialize with average values (Fokine & Urzhumtsev, 2002)
  double k_sol = 0.35;
  double b_sol = 46.0;
  gemmi::AsuData<std::complex<Real>> mask_data;
  std::vector<Point> points;

  // pre: calc and obs are sorted
  BulkSolvent(const UnitCell& cell_, const SpaceGroup* sg)
    : cell(cell_),
      crystal_system(sg ? sg->crystal_system() : CrystalSystem::Triclinic) {}

  void add_solvent_and_scale(AsuData<std::complex<Real>>& asu_data) const {
    bool use_scaling = (k_overall != 1 || !b_star.all_zero());
    for (size_t i = 0; i != asu_data.v.size(); ++i) {
      HklValue<std::complex<Real>>& hv = asu_data.v[i];
      if (use_solvent) {
        assert(hv.hkl == mask_data.v[i].hkl);
        double stol2 = cell.calculate_stol_sq(hv.hkl);
        hv.value += (Real) solvent_scale(stol2) * mask_data.v[i].value;
      }
      if (use_scaling)
        hv.value *= (Real) get_scale_factor_aniso(hv.hkl);
    }
  }

  std::vector<double> get_parameters() const {
    switch (crystal_system) {
      case CrystalSystem::Triclinic:
      case CrystalSystem::Monoclinic: // ignoring that two (e.g. U13 and U23) are 0
        return {k_overall,
                b_star.u11, b_star.u22, b_star.u33,
                b_star.u12, b_star.u13, b_star.u23};
      case CrystalSystem::Orthorhombic:
        return {k_overall, b_star.u11, b_star.u22, b_star.u33};
      case CrystalSystem::Tetragonal:
        return {k_overall, b_star.u11, b_star.u33};
      case CrystalSystem::Trigonal:
        return {k_overall, b_star.u11, b_star.u12};
      case CrystalSystem::Hexagonal:
        return {k_overall, b_star.u11, b_star.u33, b_star.u12};
      case CrystalSystem::Cubic:
        return {k_overall, b_star.u11};
    }
    unreachable();
  }

  void set_parameters(const std::vector<double>& p) {
    k_overall = p[0];
    switch (crystal_system) {
      case CrystalSystem::Triclinic:
      case CrystalSystem::Monoclinic:
        b_star = {p[1], p[2], p[3], p[4], p[5], p[6]}; break;
      case CrystalSystem::Orthorhombic:
        b_star = {p[1], p[2], p[3], 0., 0., 0.}; break;
      case CrystalSystem::Tetragonal:
        b_star = {p[1], p[1], p[2], 0., 0., 0.}; break;
      case CrystalSystem::Trigonal:
        b_star = {p[1], p[1], p[1], p[2], p[2], p[2]}; break;
      case CrystalSystem::Hexagonal:
        b_star = {p[1], p[1], p[2], p[3], 0., 0.}; break;
      case CrystalSystem::Cubic:
        b_star = {p[1], p[1], p[1], 0., 0., 0.}; break;
    }
  }

  int count_parameters() const {
    switch (crystal_system) {
      case CrystalSystem::Triclinic: return 1 + 6;
      // for consistency - ignoring that two parameters (e.g. U13 and U23) are 0
      case CrystalSystem::Monoclinic: return 1 + 6;
      case CrystalSystem::Orthorhombic: return 1 + 3;
      case CrystalSystem::Tetragonal: return 1 + 2;
      case CrystalSystem::Trigonal: return 1 + 2;
      case CrystalSystem::Hexagonal: return 1 + 3;
      case CrystalSystem::Cubic: return 1 + 1;
    }
    unreachable();
  }

  void prepare_points(const AsuData<std::complex<Real>>& calc,
                      const AsuData<std::array<Real,2>>& obs) {
    if (use_solvent && mask_data.size() != calc.size())
      fail("prepare_points(): mask data not prepared");
    std::complex<Real> fmask;
    points.reserve(std::min(calc.size(), obs.size()));
    auto c = calc.v.begin();
    for (const HklValue<std::array<Real,2>>& o : obs.v) {
      if (c->hkl != o.hkl) {
        while (*c < o.hkl) {
          ++c;
          if (c == calc.v.end())
            return;
        }
        if (c->hkl != o.hkl)
          continue;
      }
      if (use_solvent) {
        const HklValue<std::complex<Real>>& m = mask_data.v[c - calc.v.begin()];
        if (m.hkl != c->hkl)
          fail("prepare_points(): unexpected data");
        fmask = m.value;
      }
      double stol2 = cell.calculate_stol_sq(o.hkl);
      points.push_back({o.hkl, stol2, c->value, fmask, o.value[0], o.value[1]});
      ++c;
      if (c == calc.v.end())
        break;
    }
  }


  double solvent_scale(double stol2) const {
    return k_sol * std::exp(-b_sol * stol2);
  }

  double get_scale_factor_iso(const Miller& hkl) const {
    return k_overall * std::exp(-b_star.u11 * Vec3(hkl).length_sq());
  }

  void set_b_overall(const SMat33<double>& b_overall) {
    b_star = b_overall.transformed_by(cell.frac.mat);
  }
  SMat33<double> get_b_overall() const {
    return b_star.transformed_by(cell.orth.mat);
  }

  double get_scale_factor_aniso(const Miller& hkl) const {
    return k_overall * std::exp(-0.25 * b_star.r_u_r(hkl));
  }

  // quick linear fit (ignoring sigma) to get initial parameters
  void fit_isotropic_b_approximately() {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (const Point& p : points) {
      double x = p.stol2;
      double fcalc = std::abs(p.fcmol + (Real)solvent_scale(x) * p.fmask);
      double y = std::log(static_cast<float>(p.fobs / fcalc));
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
    }
    double n = (double) points.size();
    double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double intercept = (sy - slope * sx) / n;
    double b_iso = -slope;
    k_overall = exp(intercept);
    set_b_overall({b_iso, b_iso, b_iso, 0, 0, 0});
  }

  void fit_parameters() {
    LevMar levmar;
    levmar.fit(*this, points);
  }


  // interface for fitting
  std::vector<double> compute_values() const {
    std::vector<double> values;
    values.reserve(points.size());
    for (const Point& p : points) {
      double fcalc = std::abs(p.fcmol + (Real)solvent_scale(p.stol2) * p.fmask);
      values.push_back(fcalc * (Real) get_scale_factor_aniso(p.hkl));
    }
    return values;
  }

  void compute_values_and_derivatives(std::vector<double>& yy,
                                      std::vector<double>& dy_da) const {
    assert(yy.size() == points.size());
    int npar = count_parameters();
    assert(dy_da.size() == npar * points.size());
    for (size_t i = 0; i != points.size(); ++i) {
      Vec3 h(points[i].hkl);
      double fcalc;
      if (use_solvent) {
        double solv_scale = k_sol * std::exp(-b_sol * points[i].stol2);
        fcalc = std::abs(points[i].fcmol + (Real)solv_scale * points[i].fmask);
      } else {
        fcalc = std::abs(points[i].fcmol);
      }
      double fe = fcalc * std::exp(-0.25 * b_star.r_u_r(h));
      yy[i] = k_overall * fe;
      dy_da[i * npar + 0] = fe; // k_overall
      double du[6] = {
        -0.25 * yy[i] * (h.x * h.x),
        -0.25 * yy[i] * (h.y * h.y),
        -0.25 * yy[i] * (h.z * h.z),
        -0.5 * yy[i] * (h.x * h.y),
        -0.5 * yy[i] * (h.x * h.z),
        -0.5 * yy[i] * (h.y * h.z),
      };
      double* dy_db = &dy_da[i * npar + 1];
      switch (crystal_system) {
        case CrystalSystem::Triclinic:
        case CrystalSystem::Monoclinic:
          for (int j = 0; j < 6; ++j)
            dy_db[j] = du[j];
          break;
        case CrystalSystem::Orthorhombic:
          for (int j = 0; j < 3; ++j)
            dy_db[j] = du[j];
          break;
        case CrystalSystem::Tetragonal:
          dy_db[0] = du[0] + du[1];
          dy_db[1] = du[2];
          break;
        case CrystalSystem::Trigonal:
          dy_db[0] = du[0] + du[1] + du[2];
          dy_db[1] = du[3] + du[4] + du[5];
          break;
        case CrystalSystem::Hexagonal:
          dy_db[0] = du[0] + du[1];
          dy_db[1] = du[2];
          dy_db[2] = du[3];
          break;
        case CrystalSystem::Cubic:
          dy_db[0] = du[0] + du[1] + du[2];
          break;
      }
    }
  }
};

} // namespace gemmi
#endif
