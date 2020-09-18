// Copyright 2020 Global Phasing Ltd.
//

#ifndef GEMMI_BULKSOL_HPP_
#define GEMMI_BULKSOL_HPP_

#include "asudata.hpp"

namespace gemmi {

template<typename Real>
struct BulkSolvent {
  struct Point {
    Real stol2, fcalc, fobs, sigma;
  };
  std::vector<Point> data;
  UnitCell cell;
  // model parameters
  double k_overall = 1.;
  double B_overall = 0.;
  SMat33<double> B_aniso{0, 0, 0, 0, 0, 0};

  // pre: calc and obs are sorted
  BulkSolvent(AsuData<std::complex<Real>>& calc,
              const AsuData<std::array<Real,2>>& obs) {
    data.reserve(std::min(calc.size(), obs.size()));
    cell = calc.unit_cell();
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
      Real stol2 = cell.calculate_stol_sq(o.hkl);
      data.push_back({stol2, std::abs(c->value), o.value[0], o.value[1]});
      ++c;
      if (c == calc.v.end())
        return;
    }
  }

  double get_scale_factor_iso(const Miller& hkl) const {
    Real stol2 = cell.calculate_stol_sq(hkl);
    return k_overall * std::exp(-B_overall * stol2);
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
  void quick_iso_fit() {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (const Point& p : data) {
      double x = p.stol2;
      double y = std::log(static_cast<float>(p.fobs / p.fcalc));
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
    }
    double n = data.size();
    double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double intercept = (sy - slope * sx) / n;
    B_overall = -slope;
    k_overall = exp(intercept);
  }

  void aniso_fit() {
    B_aniso = {B_overall, B_overall, B_overall, 0., 0., 0.};
    // TODO fitting
    B_overall = (1./3.) * B_aniso.trace();
  }
};


} // namespace gemmi
#endif
