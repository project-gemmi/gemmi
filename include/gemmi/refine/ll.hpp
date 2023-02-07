// Copyright 2023 MRC Laboratory of Molecular Biology
//

#ifndef GEMMI_REFINE_LL_HPP_
#define GEMMI_REFINE_LL_HPP_

#include <vector>
#include <gemmi/grid.hpp>
#include <gemmi/it92.hpp>
#include <gemmi/dencalc.hpp>

namespace gemmi {

// from Refmac subroutine SMOOTH_GAUSS_D in extra_eigen.f
std::pair<double, std::vector<double>> smooth_gauss_d(double kernel_width,
                                                      const std::vector<double> &x_points,
                                                      const std::vector<double> &y_points,
                                                      double x_current) {
  assert(x_points.size() == y_points.size());
  const int n_points = x_points.size();
  if (n_points == 1)
    return std::make_pair(y_points[0], std::vector<double>(1, 1.));

  const double kernel_width2 = kernel_width * kernel_width * 2.0;

  // return values
  double y_current = 0;
  std::vector<double> y_derivs(n_points);

  if (x_current <= x_points.back() && x_current >= x_points.front()) {
    double an = 0.0;
    double fn = 0.0;
    double dx = 0, dx0 = 1.0, dx1 = 0;
    for (int i = 0; i < n_points; ++i) {
      dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2;
      dx0 = std::min(std::abs(dx), dx0);
    }
    for (int i = 0; i < n_points; ++i) {
      dx = (x_current - x_points[i]) * (x_current - x_points[i]) / kernel_width2;
      dx1 = dx - dx0;
      if (dx1 <= 120.0) {
        const double expdx = std::exp(-dx1);
        an += expdx;
        fn += y_points[i] * expdx;
        y_derivs[i] = expdx;
      }
    }
    if (an <= 0.0)
      fail("===> Error in smooth gauss. Width might be too small: 1",
                  std::to_string(n_points), " ", std::to_string(kernel_width), " ",
                  std::to_string(dx), " ", std::to_string(dx0), " ", std::to_string(dx1));

    y_current = fn / an;

    // calculate derivatives
    for (int i = 0; i < n_points; ++i) {
      const double dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2;
      const double dx1 = dx - dx0;
      if (dx1 <= 120.0)
        y_derivs[i] /= an;
      else
        y_derivs[i] = 0.0;
    }
  } else if (x_current > x_points.back()) {
    const double dx1 = (x_current - x_points.back())*(x_current - x_points.back()) / kernel_width2;
    double an = 1.0;
    double fn = y_points.back();
    double dx;
    y_derivs.back() = 1.0;
    for (int i = 0; i < n_points-1; ++i) {
      dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2 - dx1;
      if (dx <= 120.0) {
        const double expdx = std::exp(-dx);
        an += expdx;
        fn += y_points[i] * expdx;
        y_derivs[i] = expdx;
      }
    }
    if (an <= 0.0)
      fail("===> Error in smooth gauss. Width might be too small: 2 ",
                  std::to_string(n_points), " ", std::to_string(kernel_width), " ",
                  std::to_string(dx));

    y_current = fn / an;

    // calculate derivatives
    for (int i = 0; i < n_points; ++i) {
      const double dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2 - dx1;
      if (dx <= 120.0)
        y_derivs[i] /= an;
      else
        y_derivs[i] = 0.0;
    }
  } else if (x_current < x_points.front()) {
    const double dx1 = (x_current-x_points.front())*(x_current-x_points.front()) / kernel_width2;
    double an = 1.0;
    double fn = y_points.front();
    double dx;
    y_derivs[0] = 1.0;
    for (int i = 1; i < n_points; ++i) {
      dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2 - dx1;
      if (dx <= 120.0) {
        const double expdx = std::exp(-dx);
        an += expdx;
        fn = y_points[i]*expdx;
        y_derivs[i] = expdx;
      }
    }
    if (an <= 0.0)
      fail("===> Error in smooth gauss. Width might be too small: 3 ",
                  std::to_string(n_points), " ", std::to_string(kernel_width), " ",
                  std::to_string(dx));

    y_current = fn/an;

    // calculate derivatives
    for (int i = 0; i < n_points; ++i) {
      const double dx = (x_current - x_points[i])*(x_current - x_points[i]) / kernel_width2 - dx1;
      if (dx <= 120.0)
        y_derivs[i] /= an;
      else
        y_derivs[i] = 0.0;
    }
  }
  return std::make_pair(y_current, y_derivs);
}

struct TableS3 {
  double s_min, s_max;
  double delta_s3;
  int n_points;
  std::vector<double> s3_values;
  std::vector<double> y_values;

  TableS3(double d_min, double d_max) {
    s_min = 1. / d_max;
    s_max = 1. / d_min;
    const double smax3 = s_max*s_max*s_max, smin3 = s_min*s_min*s_min;
    delta_s3 = 0.0005;
    n_points = (int) ((smax3 - smin3) / delta_s3);
    n_points = std::max(20, std::min(2000, n_points));
    delta_s3 = (smax3 - smin3) / n_points;
    s3_values.reserve(n_points+1);
    y_values.reserve(n_points+1);
    for (int i = 0; i <= n_points; ++i)
      s3_values.push_back(smin3 + i * delta_s3);
  }

  // From Refmac SUBROUTINE D2DA_RADIAL_BIN in hkon_secder_tch.f
  void make_table(const std::vector<double> &svals, const std::vector<double> &yvals) {
    assert(svals.size() == yvals.size());

    // define bin
    const double smax_ml = *std::max_element(svals.begin(), svals.end()) * 1.0001;
    const double smin_ml = *std::min_element(svals.begin(), svals.end()) * 0.9999;
    const double smin_ml3 = smin_ml * smin_ml * smin_ml;
    const double smax_ml3 = smax_ml * smax_ml * smax_ml;

    // from Refmac SUBROUTINE DEFINE_BINS_FOR_ML in oppro_allocate.f
    double binsize_ml = 0.0005;
    int nbin_ml =  std::max(1, std::min(500,
                                        (int)((smax_ml*smax_ml - smin_ml*smin_ml)/binsize_ml)+1
                                        ));
    binsize_ml = (smax_ml*smax_ml - smin_ml*smin_ml) / nbin_ml;
    int nbin_rad = nbin_ml;
    const double ds3 = (smax_ml3 - smin_ml3) / nbin_rad;

    std::vector<double> smeanb_rad(nbin_rad+1);
    for (int i = 0; i <= nbin_rad; ++i)
      smeanb_rad[i] = smin_ml3 + i * ds3;

    std::vector<double> sec_der_bin(nbin_rad+1);
    std::vector<int> nref_sec_bin(nbin_rad+1);
    for (size_t i = 0; i < svals.size(); ++i) {
      const double s = svals[i];
      const int ibin = std::upper_bound(smeanb_rad.begin(), smeanb_rad.end(), s*s*s) - smeanb_rad.begin() - 1;
      // int ibin;
      // for (int j = 0; j < nbin_rad; ++j)
      //   if (s3 > smeanb_rad[j] && s3 <= smeanb_rad[j+1]) {
      //     ibin = j;
      //     int ibin2 = std::upper_bound(smeanb_rad.begin(), smeanb_rad.end(), s*s*s) - smeanb_rad.begin();
      //     std::cout << "ibin= " << j << " diff= " << ibin2-j << "\n";
      //     break;
      //   }

      sec_der_bin[ibin] += yvals[i];
      nref_sec_bin[ibin] += 1;
    }

    for (int i = 0; i < nbin_rad; ++i)
      if (sec_der_bin[i] > 0 && nref_sec_bin[i] > 0)
        sec_der_bin[i] = std::log(sec_der_bin[i] / nref_sec_bin[i]);

    // sort out bins with no (usable) reflections
    for (int i = 1; i < nbin_rad; ++i)
      if (nref_sec_bin[i] <= 0 && nref_sec_bin[i-1] > 0)
        sec_der_bin[i] = sec_der_bin[i-1];
    for (int i = nbin_rad-2; i >= 0; --i)
      if (nref_sec_bin[i] <= 0 && nref_sec_bin[i+1] > 0)
        sec_der_bin[i] = sec_der_bin[i+1];

    // Start refining parameters of smoothining function
    double tmp1 = sec_der_bin[0];
    sec_der_bin[nbin_rad] = sec_der_bin[nbin_rad-1];
    for (int i = 1; i < nbin_rad; ++i) {
      double tmp2 = sec_der_bin[i];
      sec_der_bin[i] = (tmp1 + tmp2)/2.0;
      tmp1 = tmp2;
    }

    // TODO smooth curve here?

    const double kernel_g_rad = (smeanb_rad[1] - smeanb_rad[0]) / 2.0;
    for (int i = 0; i <= n_points; ++i) {
      const double yval = smooth_gauss_d(kernel_g_rad, smeanb_rad, sec_der_bin, s3_values[i]).first;
      y_values.push_back(std::exp(yval));
    }
  }

  double get_value(double s) const {
    const double s3 = s*s*s;
    const int i =  std::max(0, std::min(n_points,
                                        (int)std::round((s3 - s3_values.front()) / delta_s3)));
    return y_values[i];
  }
};

// crystal is not supported yet
// only works with isotropic ADPs
template <typename Table>
struct LL{
  std::vector<Atom*> atoms;
  UnitCell cell;
  SpaceGroup *sg;
  std::vector<Transform> ncs;
  bool mott_bethe;

  // table (distances x b values)
  std::vector<double> table_bs;
  std::vector<std::vector<double>> pp1; // for x-x diagonal
  std::vector<std::vector<double>> bb;  // for B-B diagonal

  LL(UnitCell cell, SpaceGroup *sg, const std::vector<Atom*> &atoms, bool mott_bethe)
    : atoms(atoms), cell(cell), sg(sg), mott_bethe(mott_bethe) {
    set_ncs({});
  }
  void set_ncs(const std::vector<Transform> &trs) {
    ncs.clear();
    ncs.push_back(Transform()); // make sure first is the identity op
    for (const auto &tr : trs)
      if (!tr.is_identity())
        ncs.push_back(tr);
  }

  // FFT-based gradient calculation: Murshudov et al. (1997) 10.1107/S0907444996012255
  // only assumes cryo-EM SPA
  // den is the Fourier transform of (dLL/dAc-i dLL/dBc)*mott_bethe_factor/s^2
  std::vector<double> calc_grad(Grid<float> &den, bool refine_xyz, bool refine_adp) { // needs <double>?
    const size_t n_atoms = atoms.size();
    const size_t n_v = n_atoms * ((refine_xyz ? 3 : 0) + (refine_adp ? 1 : 0)); // only isotropic ADP for now
    std::vector<double> vn(n_v, 0.);
    for (size_t i = 0; i < n_atoms; ++i) {
      for (const Transform &tr : ncs) { //TODO to use cell images?
        const Atom &atom = *atoms[i];
        const Fractional fpos = cell.fractionalize(Position(tr.apply(atom.pos)));
        const Element &el = atom.element;
        const auto coef = Table::get(el);
        const auto precal = coef.precalculate_density_iso(atom.b_iso,
                                                          mott_bethe ? -el.atomic_number() : 0.);
        // is it ok to use a radius based on ADP? need to apply blur to den?
        const double radius = determine_cutoff_radius(it92_radius_approx(atom.b_iso),
                                                      precal, 1e-7); // TODO cutoff?
        const int N = sizeof(precal.a) / sizeof(precal.a[0]);
        const int du = (int) std::ceil(radius / den.spacing[0]);
        const int dv = (int) std::ceil(radius / den.spacing[1]);
        const int dw = (int) std::ceil(radius / den.spacing[2]);
        Position gx;
        double gb = 0.;
        den.template use_points_in_box<true>(fpos, du, dv, dw,
                                             [&](float& point, const Position& delta, int, int, int) {
                                               const double r2 = delta.length_sq();
                                               if (r2 > radius * radius) return;
                                               double for_x = 0., for_b = 0.;
                                               for (int j = 0; j < N; ++j) {
                                                 const double tmp = precal.a[j] * std::exp(precal.b[j] * r2) * precal.b[j];
                                                 for_x += tmp;
                                                 if (refine_adp) for_b += tmp * (1.5 + r2 * precal.b[j]);
                                               }
                                               gx += for_x * 2 * delta * point; // -1 for flipping delta?
                                               if (refine_adp) gb += for_b * point;
                                             });
        gx *= atom.occ;
        gb *= atom.occ * 0.25 / pi() / pi();
        if (mott_bethe) {
          gx *= -1;
          gb *= -1;
        }
        const auto gx2 = tr.mat.transpose().multiply(gx);
        if (refine_xyz) {
          vn[3*i  ] += gx2.x;
          vn[3*i+1] += gx2.y;
          vn[3*i+2] += gx2.z;
        }
        if (refine_adp)
          vn[(refine_xyz ? n_atoms * 3 : 0) + i] += gb;
      }
    }
    for (auto &v : vn) // to match scale of hessian
      v /= ncs.size();
    return vn;
  }

  /*
  void add_fisher_diagonal_naive(const Vec3 &svec, double d2ll) {
    // TODO symmetry
    // TODO aniso
    const double s2 = svec.length_sq();
    for (size_t i = 0; i < atoms.size(); ++i) {
      const Atom &atom = *atoms[i];
      const auto coef = IT92<double>::get(atom.element); // as Table::?
      const double f = (atom.element.atomic_number() - coef.calculate_sf(s2/4))/s2; // only for mott_bethe
      const double w = atom.occ*atom.occ*f*f*std::exp(-atom.b_iso*s2/2)*d2ll * 4 * pi()*pi();
      const int ipos = i*6;
      for (int sign = 0; sign < 2; ++sign) { // Friedel pair
        const Vec3 s = sign ? -svec : svec;
        am[ipos]   += w * s.x * s.x;
        am[ipos+1] += w * s.y * s.y;
        am[ipos+2] += w * s.z * s.z;
        am[ipos+3] += w * s.y * s.x;
        am[ipos+4] += w * s.z * s.x;
        am[ipos+5] += w * s.z * s.y;
      }
    }
  }
  */

  // preparation for fisher_diag_from_table()
  // Steiner et al. (2003) doi: 10.1107/S0907444903018675
  void make_fisher_table_diag_fast(double b_min, double b_max,
                                   const TableS3 &d2dfw_table) {
    pp1.resize(1);
    bb.resize(1);
    const double b_step = 5;
    const double s_min = d2dfw_table.s_min, s_max = d2dfw_table.s_max;
    const double s_dim = 120; // actually +1 is allocated
    int b_dim = static_cast<int>((b_max - b_min) / b_step) + 1;
    if (b_dim % 2 == 0) ++b_dim; // TODO: need to set maximum b_dim?
    pp1[0].resize(b_dim);
    bb[0].resize(b_dim);

    const double s_step = (s_max - s_min) / s_dim;

    table_bs.clear();
    table_bs.reserve(b_dim);

    // only for D = 0 (same atoms) for now
    for (int ib = 0; ib < b_dim; ++ib) {
      const double b = b_min + b_step * ib;
      table_bs.push_back(b);

      std::vector<double> tpp(s_dim+1), tbb(s_dim+1);
      for (int i = 0; i <= s_dim; ++i) {
        const double s = s_min + s_step * i;
        const double w_c = d2dfw_table.get_value(s); // average of weight
        const double w_c_ft_c = w_c * std::exp(-b*s*s/4.);
        tpp[i] = 16. * pi() * pi() * pi() * w_c_ft_c / 3.; // (2pi)^2 * 4pi/3
        tbb[i] = pi() / 4 * w_c_ft_c * s * s; // 1/16 * 4pi
        if (!mott_bethe) {
          tpp[i] *= s*s*s*s;
          tbb[i] *= s*s*s*s;
        }
      }

      // Numerical integration by Simpson's rule
      double sum_tpp1 = 0, sum_tpp2 = 0, sum_tbb1 = 0, sum_tbb2 = 0;
      for (int i = 1; i < s_dim; i+=2) {
        sum_tpp1 += tpp[i];
        sum_tbb1 += tbb[i];
      }
      for (int i = 2; i < s_dim; i+=2) {
        sum_tpp2 += tpp[i];
        sum_tbb2 += tbb[i];
      }

      pp1[0][ib] = (tpp[0] + tpp.back() + 4 * sum_tpp1 + 2 * sum_tpp2) * s_step / 3.;
      bb[0][ib] = (tbb[0] + tbb.back() + 4 * sum_tbb1 + 2 * sum_tbb2) * s_step / 3.;
    }
  }

  // from Refmac SUBROUTINE LINTER_VALUE2
  // no need to be a member of this class
  double interp_1d(const std::vector<double> &x_points,
                   const std::vector<double> &y_points,
                   double x) const {
    assert(x_points.size() == y_points.size());
    assert(!x_points.empty());

    if (x_points.size() == 1)
      return y_points.front();

    const int k1 = std::min((size_t) (std::lower_bound(x_points.begin(), x_points.end(), x) - x_points.begin()),
                            x_points.size()-2);

    // calculate value of function at the given point
    double b = y_points[k1];
    double a = (y_points[k1+1] - y_points[k1]) / (x_points[k1+1] - x_points[k1]);
    double dx = x - x_points[k1];
    double y = a * dx + b;
    if (x < x_points.front())
      return std::max(0.1 * y_points.front(), std::min(10.0 * y_points.front(), y));
    else if (x > x_points.back())
      return std::max(0.1 * y_points.back(), std::min(10.0 * y_points.back(), y));
    else
      return y;
  }

  std::vector<double> fisher_diag_from_table (bool refine_xyz, bool refine_adp) {
    const size_t n_atoms = atoms.size();
    const size_t n_a = n_atoms * ((refine_xyz ? 3 : 0) + (refine_adp ? 1 : 0)); // only isotropic ADP for now
    const int N = Table::Coef::ncoeffs;
    std::vector<double> am(n_a, 0.);
    for (size_t i = 0; i < n_atoms; ++i) {
      const Atom &atom = *atoms[i];
      const auto coef = Table::get(atom.element);
      const double w = atom.occ * atom.occ;
      const double c = mott_bethe ? coef.c() - atom.element.atomic_number(): coef.c();
      double fac_x = 0., fac_b = 0.;

      // TODO can be reduced for the same elements
      for (int j = 0; j < N + 1; ++j)
        for (int k = 0; k < N + 1; ++k) {
          // * -1 is needed for mott_bethe case, but we only need aj * ak so they cancel.
          const double aj = j < N ? coef.a(j) : c;
          const double ak = k < N ? coef.a(k) : c;
          const double b = 2 * atom.b_iso + (j < N ? coef.b(j) : 0) + (k < N ? coef.b(k) : 0);
          fac_x += aj * ak * interp_1d(table_bs, pp1[0], b);
          fac_b += aj * ak * interp_1d(table_bs, bb[0], b);
        }

      const int ipos = i*3;
      if (refine_xyz) am[ipos] = am[ipos+1] = am[ipos+2] = w * fac_x;
      if (refine_adp) am[(refine_xyz ? n_atoms * 3 : 0) + i] = w * fac_b;
    }
    return am;
  }
};



} // namespace gemmi
#endif
