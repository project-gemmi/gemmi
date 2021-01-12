// Copyright 2020 Global Phasing Ltd.
//

#ifndef GEMMI_BULKSOL_HPP_
#define GEMMI_BULKSOL_HPP_

#include "asudata.hpp"
#include "levmar.hpp"

namespace gemmi {

template<typename Real>
struct FcFo {
  Miller hkl;
  Real fcalc, fobs, sigma;

  Miller get_x() const { return hkl; }
  double get_y() const { return fobs; }
  double get_weight() const { return 1.0 / sigma; }
};

template<typename Real>
std::vector<FcFo<Real>> prepare_fc_fo(const AsuData<std::complex<Real>>& calc,
                                      const AsuData<std::array<Real,2>>& obs) {
  std::vector<FcFo<Real>> data;
  data.reserve(std::min(calc.size(), obs.size()));
  auto c = calc.v.begin();
  for (const HklValue<std::array<Real,2>>& o : obs.v) {
    if (c->hkl != o.hkl) {
      while (*c < o.hkl) {
        ++c;
        if (c == calc.v.end())
          return data;
      }
      if (c->hkl != o.hkl)
        continue;
    }
    data.push_back({o.hkl, std::abs(c->value), o.value[0], o.value[1]});
    ++c;
    if (c == calc.v.end())
      break;
  }
  return data;
}

template<typename Real>
struct BulkSolvent {
  using Point = FcFo<Real>;
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

  // pre: calc and obs are sorted
  BulkSolvent(const UnitCell& cell_, const SpaceGroup* sg)
    : cell(cell_),
      crystal_system(sg ? sg->crystal_system() : CrystalSystem::Triclinic) {}

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
  void quick_iso_fit(const std::vector<Point>& data) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (const Point& p : data) {
      double x = cell.calculate_stol_sq(p.hkl);
      double y = std::log(static_cast<float>(p.fobs / p.fcalc));
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
    }
    double n = (double) data.size();
    double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double intercept = (sy - slope * sx) / n;
    double b_iso = -slope;
    k_overall = exp(intercept);
    set_b_overall({b_iso, b_iso, b_iso, 0, 0, 0});
  }

  void aniso_fit(const std::vector<Point>& data) {
    LevMar levmar;
    levmar.fit(*this, data);
  }


  // interface for fitting
  std::vector<double> compute_values(const std::vector<Point>& data) const {
    std::vector<double> values(data.size());
    for (size_t i = 0; i != data.size(); ++i)
      values[i] = data[i].fcalc * get_scale_factor_aniso(data[i].hkl);
    return values;
  }

  void compute_values_and_derivatives(const std::vector<Point>& data,
                                      std::vector<double>& yy,
                                      std::vector<double>& dy_da) const {
    assert(data.size() == yy.size());
    int npar = count_parameters();
    assert(dy_da.size() == npar * data.size());
    for (size_t i = 0; i != data.size(); ++i) {
      Vec3 h(data[i].hkl);
      double fe = data[i].fcalc * std::exp(-0.25 * b_star.r_u_r(h));
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
