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
  // model parameters
  double k_overall = 1.;
  SMat33<double> B_aniso{0, 0, 0, 0, 0, 0};

  // pre: calc and obs are sorted
  BulkSolvent(const UnitCell& cell_) : cell(cell_) {}

  std::vector<double> get_parameters() const {
    return {k_overall,
            B_aniso.u11, B_aniso.u22, B_aniso.u33,
            B_aniso.u12, B_aniso.u13, B_aniso.u23};
  }

  void set_parameters(const std::vector<double>& p) {
    k_overall = p[0];
    B_aniso = {p[1], p[2], p[3], p[4], p[5], p[6]};
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
    assert(dy_da.size() == 7 * data.size());
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
      dy_da[i*7 + 0] = fe; // k_overall
      dy_da[i*7 + 1] = -0.25 * yy[i] * (arh * arh);  // u11
      dy_da[i*7 + 2] = -0.25 * yy[i] * (brk * brk);  // u22
      dy_da[i*7 + 3] = -0.25 * yy[i] * (crl * crl);  // u33
      dy_da[i*7 + 4] = -0.5 * yy[i] * (arh * brk * cell.cos_gammar);  // u12
      dy_da[i*7 + 5] = -0.5 * yy[i] * (arh * crl * cell.cos_betar);   // u13
      dy_da[i*7 + 6] = -0.5 * yy[i] * (brk * crl * cell.cos_alphar);  // u23
    }
  }
};

} // namespace gemmi
#endif
