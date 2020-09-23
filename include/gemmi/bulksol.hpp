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
};

template<typename Real>
struct BulkSolvent {
  using Point = FcFo<Real>;
  UnitCell cell;
  CrystalSystem crystal_system;
  // model parameters
  double k_overall = 1.;
  SMat33<double> B_aniso{0, 0, 0, 0, 0, 0};

  // pre: calc and obs are sorted
  BulkSolvent(const UnitCell& cell_, const SpaceGroup* sg)
    : cell(cell_),
      crystal_system(sg ? sg->crystal_system() : CrystalSystem::Triclinic) {}

  std::vector<double> get_parameters() const {
    switch (crystal_system) {
      case CrystalSystem::Triclinic:
      case CrystalSystem::Monoclinic: // ignoring that two (e.g. U13 and U23) are 0
        return {k_overall,
                B_aniso.u11, B_aniso.u22, B_aniso.u33,
                B_aniso.u12, B_aniso.u13, B_aniso.u23};
      case CrystalSystem::Orthorhombic:
        return {k_overall, B_aniso.u11, B_aniso.u22, B_aniso.u33};
      case CrystalSystem::Tetragonal:
        return {k_overall, B_aniso.u11, B_aniso.u33};
      case CrystalSystem::Trigonal:
        return {k_overall, B_aniso.u11, B_aniso.u12};
      case CrystalSystem::Hexagonal:
        return {k_overall, B_aniso.u11, B_aniso.u33, B_aniso.u12};
      case CrystalSystem::Cubic:
        return {k_overall, B_aniso.u11};
    }
    unreachable();
  }

  void set_parameters(const std::vector<double>& p) {
    k_overall = p[0];
    switch (crystal_system) {
      case CrystalSystem::Triclinic:
      case CrystalSystem::Monoclinic:
        B_aniso = {p[1], p[2], p[3], p[4], p[5], p[6]}; break;
      case CrystalSystem::Orthorhombic:
        B_aniso = {p[1], p[2], p[3], 0., 0., 0.}; break;
      case CrystalSystem::Tetragonal:
        B_aniso = {p[1], p[1], p[2], 0., 0., 0.}; break;
      case CrystalSystem::Trigonal:
        B_aniso = {p[1], p[1], p[1], p[2], p[2], p[2]}; break;
      case CrystalSystem::Hexagonal:
        B_aniso = {p[1], p[1], p[2], p[3], 0., 0.}; break;
      case CrystalSystem::Cubic:
        B_aniso = {p[1], p[1], p[1], 0., 0., 0.}; break;
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

  double get_scale_factor_iso(const Miller& hkl) const {
    Real stol2 = cell.calculate_stol_sq(hkl);
    return k_overall * std::exp(-B_aniso.u11 * stol2);
  }

  double get_scale_factor_aniso(const Miller& hkl) const {
    double arh = cell.ar * hkl[0];
    double brk = cell.br * hkl[1];
    double crl = cell.cr * hkl[2];
    double sbs = B_aniso.u11 * arh * arh +
                 B_aniso.u22 * brk * brk +
                 B_aniso.u33 * crl * crl +
             2 * (B_aniso.u12 * arh * brk * cell.cos_gammar +
                  B_aniso.u13 * arh * crl * cell.cos_betar +
                  B_aniso.u23 * brk * crl * cell.cos_alphar);
    return k_overall * std::exp(-0.25 * sbs);
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
    double n = data.size();
    double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double intercept = (sy - slope * sx) / n;
    B_aniso.u11 = B_aniso.u22 = B_aniso.u33 = -slope;
    B_aniso.u12 = B_aniso.u13 = B_aniso.u23 = 0;
    k_overall = exp(intercept);
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
      const Miller& hkl = data[i].hkl;
      double arh = cell.ar * hkl[0];
      double brk = cell.br * hkl[1];
      double crl = cell.cr * hkl[2];
      double sbs = arh * arh * B_aniso.u11 +
                   brk * brk * B_aniso.u22 +
                   crl * crl * B_aniso.u33 +
               2 * (arh * brk * cell.cos_gammar * B_aniso.u12 +
                    arh * crl * cell.cos_betar * B_aniso.u13 +
                    brk * crl * cell.cos_alphar * B_aniso.u23);
      double fe = data[i].fcalc * std::exp(-0.25 * sbs);
      yy[i] = k_overall * fe;
      dy_da[i * npar + 0] = fe; // k_overall
      double du[6] = {
        -0.25 * yy[i] * (arh * arh),
        -0.25 * yy[i] * (brk * brk),
        -0.25 * yy[i] * (crl * crl),
        -0.5 * yy[i] * (arh * brk * cell.cos_gammar),
        -0.5 * yy[i] * (arh * crl * cell.cos_betar),
        -0.5 * yy[i] * (brk * crl * cell.cos_alphar),
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
