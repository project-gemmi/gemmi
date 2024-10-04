// Copyright 2020 Global Phasing Ltd.
//
// Anisotropic scaling of data (includes scaling of bulk solvent parameters).

#ifndef GEMMI_SCALING_HPP_
#define GEMMI_SCALING_HPP_

#include "asudata.hpp"
#include "levmar.hpp"
#if WITH_NLOPT
# include <nlopt.h>
#endif

namespace gemmi {

using Vec6 = std::array<double, 6>;

inline double vec6_dot(const Vec6& a, const SMat33<double>& s) {
  return a[0] * s.u11 + a[1] * s.u22 + a[2] * s.u33
       + a[3] * s.u12 + a[4] * s.u13 + a[5] * s.u23;
}

/// Symmetry constraints of ADP.
/// The number of rows is the number of independent coefficients in U.
/// For example, for tetragonal crystal returns two normalized Vec6 vectors
/// in directions [1 1 0 0 0 0] and [0 0 1 0 0 0].
inline std::vector<Vec6> adp_symmetry_constraints(const SpaceGroup* sg) {
  auto constraints = [](const std::initializer_list<int>& l) {
    constexpr double K2 = 0.70710678118654752440;  // sqrt(1/2)
    constexpr double K3 = 0.57735026918962576451;  // sqrt(1/3)
    static const Vec6 vv[10] = {
      {1, 0, 0, 0, 0, 0},  // 0
      {0, 1, 0, 0, 0, 0},  // 1
      {0, 0, 1, 0, 0, 0},  // 2
      {0, 0, 0, 1, 0, 0},  // 3
      {0, 0, 0, 0, 1, 0},  // 4
      {0, 0, 0, 0, 0, 1},  // 5
      {K2, K2, 0, 0, 0, 0},  // 6
      {2/3., 2/3., 0, 1/3., 0, 0},  // 7
      {K3, K3, K3, 0, 0, 0},  // 8
      {0, 0, 0, K3, K3, K3}   // 9
    };
    std::vector<Vec6> c;
    c.reserve(l.size());
    for (int i : l)
      c.push_back(vv[i]);
    return c;
  };
  CrystalSystem cr_system = sg ? sg->crystal_system() : CrystalSystem::Triclinic;
  switch (cr_system) {
    case CrystalSystem::Triclinic:
      return constraints({0, 1, 2, 3, 4, 5});
    case CrystalSystem::Monoclinic:
      // the last index is: a->5, b->4, c->3
      return constraints({0, 1, 2, 3 + 'c' - sg->monoclinic_unique_axis()});
    case CrystalSystem::Orthorhombic:
      return constraints({0, 1, 2});
    case CrystalSystem::Tetragonal:
      return constraints({6, 2});
    case CrystalSystem::Hexagonal:
    case CrystalSystem::Trigonal:
      if (sg->ext == 'R')
        return constraints({8, 9});
      return constraints({7, 2});
    case CrystalSystem::Cubic:
      return constraints({8});
  }
  unreachable();
}

template<typename Real>
struct Scaling {
  struct Point {
    Miller hkl;
    double stol2;
    std::complex<Real> fcmol, fmask;
    Real fobs, sigma;

    Miller get_x() const { return hkl; }
    double get_y() const { return fobs; }
    double get_weight() const { return 1.0 /* / sigma*/; }
  };

  UnitCell cell;
  // model parameters
  double k_overall = 1.;
  // b_star = F B_cart F^T, where F - fractionalization matrix
  SMat33<double> b_star{0, 0, 0, 0, 0, 0};
  std::vector<Vec6> constraint_matrix;
  bool use_solvent = false;
  bool fix_k_sol = false;
  bool fix_b_sol = false;
  // initialize with average values (Fokine & Urzhumtsev, 2002)
  double k_sol = 0.35;
  double b_sol = 46.0;
  std::vector<Point> points;

  Scaling(const UnitCell& cell_, const SpaceGroup* sg)
      : cell(cell_), constraint_matrix(adp_symmetry_constraints(sg)) {}

  // B_{overall} is stored as B* not B_{cartesian}.
  // Use getter and setter to convert from/to B_{cartesian}.
  void set_b_overall(const SMat33<double>& b_overall) {
    b_star = b_overall.transformed_by(cell.frac.mat);
  }
  SMat33<double> get_b_overall() const {
    return b_star.transformed_by(cell.orth.mat);
  }

  // Scale data, optionally adding bulk solvent correction.
  void scale_data(AsuData<std::complex<Real>>& asu_data,
                  const AsuData<std::complex<Real>>* mask_data) const {
    if (use_solvent && !(mask_data && mask_data->size() == asu_data.size()))
      fail("scale_data(): mask data not prepared");
    bool use_scaling = (k_overall != 1 || !b_star.all_zero());
    for (size_t i = 0; i != asu_data.v.size(); ++i) {
      HklValue<std::complex<Real>>& hv = asu_data.v[i];
      if (use_solvent) {
        if (hv.hkl != mask_data->v[i].hkl)
          fail("scale_data(): data arrays don't match");
        double stol2 = cell.calculate_stol_sq(hv.hkl);
        hv.value += (Real)get_solvent_scale(stol2) * mask_data->v[i].value;
      }
      if (use_scaling)
        hv.value *= (Real) get_overall_scale_factor(hv.hkl);
    }
  }

  std::complex<Real> scale_value(const Miller& hkl, std::complex<Real> f_value,
                                 std::complex<Real> mask_value) {
    if (use_solvent) {
      double stol2 = cell.calculate_stol_sq(hkl);
      f_value += (Real)get_solvent_scale(stol2) * mask_value;
    }
    return f_value * (Real) get_overall_scale_factor(hkl);
  }

  std::vector<double> get_parameters() const {
    std::vector<double> ret;
    ret.push_back(k_overall);
    if (use_solvent) {
      if (!fix_k_sol)
        ret.push_back(k_sol);
      if (!fix_b_sol)
        ret.push_back(b_sol);
    }
    for (const Vec6& v : constraint_matrix)
      ret.push_back(vec6_dot(v, b_star));
    return ret;
  }

  /// set k_overall, k_sol, b_sol, b_star
  void set_parameters(const double* p) {
    k_overall = p[0];
    int n = 0;
    if (use_solvent) {
      if (!fix_k_sol)
        k_sol = p[++n];
      if (!fix_b_sol)
        b_sol = p[++n];
    }
    b_star = {0, 0, 0, 0, 0, 0};
    for (const Vec6& row : constraint_matrix) {
      double d = p[++n];
      b_star.u11 += row[0] * d;
      b_star.u22 += row[1] * d;
      b_star.u33 += row[2] * d;
      b_star.u12 += row[3] * d;
      b_star.u13 += row[4] * d;
      b_star.u23 += row[5] * d;
    }
  }

  void set_parameters(const std::vector<double>& p) {
    set_parameters(p.data());
  }

  // pre: all AsuData args are sorted
  void prepare_points(const AsuData<std::complex<Real>>& calc,
                      const AsuData<ValueSigma<Real>>& obs,
                      const AsuData<std::complex<Real>>* mask_data) {
    if (use_solvent && !(mask_data && mask_data->size() == calc.size()))
      fail("prepare_points(): mask data not prepared");
    std::complex<Real> fmask;
    points.reserve(std::min(calc.size(), obs.size()));
    auto c = calc.v.begin();
    for (const HklValue<ValueSigma<Real>>& o : obs.v) {
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
        const HklValue<std::complex<Real>>& m = mask_data->v[c - calc.v.begin()];
        if (m.hkl != c->hkl)
          fail("prepare_points(): unexpected data");
        fmask = m.value;
      }
      double stol2 = cell.calculate_stol_sq(o.hkl);
      if (!std::isnan(o.value.value) && !std::isnan(o.value.sigma))
        points.push_back({o.hkl, stol2, c->value, fmask, o.value.value, o.value.sigma});
      ++c;
      if (c == calc.v.end())
        break;
    }
  }


  double get_solvent_scale(double stol2) const {
    return k_sol * std::exp(-b_sol * stol2);
  }

  double get_overall_scale_factor(const Miller& hkl) const {
    return k_overall * std::exp(-0.25 * b_star.r_u_r(hkl));
  }

  std::complex<Real> get_fcalc(const Point& p) const {
    if (!use_solvent)
      return p.fcmol;
    return p.fcmol + (Real)get_solvent_scale(p.stol2) * p.fmask;
  }

  // quick linear fit (ignoring sigma) to get initial k_overall and isotropic B
  void fit_isotropic_b_approximately() {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (const Point& p : points) {
      if (p.fobs < 1 || p.fobs < p.sigma)  // skip weak reflections
        continue;
      double fcalc = std::abs(get_fcalc(p));
      double x = p.stol2;
      double y = std::log(static_cast<float>(p.fobs / fcalc));
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
      n += 1;
    }
    if (n <= 5)  // this is not expected to happen
      return;
    double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double intercept = (sy - slope * sx) / n;
    double b_iso = -slope;
    k_overall = std::exp(intercept);
    set_b_overall({b_iso, b_iso, b_iso, 0, 0, 0});
  }

  // least-squares fitting of k_overall only
  double lsq_k_overall() const {
    double sxx = 0, sxy = 0;
    for (const Point& p : points) {
      if (p.fobs < 1 || p.fobs < p.sigma)  // skip weak reflections
        continue;
      double x = std::abs(get_fcalc(p));
      double y = p.fobs;
      sxx += x * x;
      sxy += x * y;
    }
    return sxx != 0. ? sxy / sxx : 1.;
  }

  // For testing only, don't use it.
  // Estimates anisotropic b_star using other parameters (incl. isotropic B),
  // following P. Afonine et al, doi:10.1107/S0907444913000462 sec. 2.1.
  // The symmetry constraints are not implemented - don't use it!
  void fit_b_star_approximately() {
    double b_iso = 1/3. * get_b_overall().trace();
    //size_t nc = constraint_matrix.size();
    double M[36] = {};
    double b[6] = {};
    //std::vector<double> Vc(nc);
    for (const Point& p : points) {
      double fcalc = std::abs(get_fcalc(p));
      // the factor 1 / 2 pi^2 will be taken into account later
      double Z = std::log(p.fobs / (k_overall * fcalc)) + b_iso * p.stol2;
      Vec3 h(p.hkl);
      double V[6] = {h.x * h.x, h.y * h.y, h.z * h.z,
                     2 * h.x * h.y, 2 * h.x * h.z, 2 * h.y * h.z};
      //for (size_t i = 0; i < nc; ++i)
      //  Vc[i] = vec6_dot(constraint_matrix[i], V);
      for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j)
          M[6*i+j] += V[i] * V[j];
        b[i] -= Z * V[i];
      }
    }
    jordan_solve(M, b, 6);
    double b_star_iso = 1/3. * b_star.trace();
    SMat33<double> u_star{b[0], b[1], b[2], b[3], b[4], b[5]};
    // u_to_b() / (2 pi^2) = 8 pi^2 / 2 pi^2 = 4
    b_star = u_star.scaled(4.0).added_kI(b_star_iso);
    //auto e = get_b_overall().elements_pdb();
    //printf("fitted B = {%g %g %g  %g %g %g}\n", e[0], e[1], e[2], e[3], e[4], e[5]);
  }

  double fit_parameters() {
    LevMar levmar;
    return levmar.fit(*this);
  }

  double calculate_r_factor() const {
    double abs_diff_sum = 0;
    double denom = 0;
    for (const Point& p : points) {
      abs_diff_sum += std::fabs(p.fobs - compute_value(p));
      denom += p.fobs;
    }
    return abs_diff_sum / denom;
  }

  // interface for fitting
  double compute_value(const Point& p) const {
    return std::abs(get_fcalc(p)) * (Real) get_overall_scale_factor(p.hkl);
  }

  double compute_value_and_derivatives(const Point& p, std::vector<double>& dy_da) const {
    Vec3 h(p.hkl);
    double kaniso = std::exp(-0.25 * b_star.r_u_r(h));
    double fcalc_abs;
    int n = 1;
    if (use_solvent) {
      double solv_b = std::exp(-b_sol * p.stol2);
      double solv_scale = k_sol * solv_b;
      auto fcalc = p.fcmol + (Real)solv_scale * p.fmask;
      fcalc_abs = std::abs(fcalc);
      double dy_dsol = (fcalc.real() * p.fmask.real() +
                        fcalc.imag() * p.fmask.imag()) / fcalc_abs * k_overall * kaniso;
      if (!fix_k_sol)
        dy_da[n++] = solv_b * dy_dsol;
      if (!fix_b_sol)
        dy_da[n++] = -p.stol2 * solv_scale * dy_dsol;
    } else {
      fcalc_abs = std::abs(p.fcmol);
    }
    double fe = fcalc_abs * kaniso;
    double y = k_overall * fe;
    dy_da[0] = fe; // dy/d k_overall
    SMat33<double> du = {
      -0.25 * y * (h.x * h.x),
      -0.25 * y * (h.y * h.y),
      -0.25 * y * (h.z * h.z),
      -0.5 * y * (h.x * h.y),
      -0.5 * y * (h.x * h.z),
      -0.5 * y * (h.y * h.z),
    };
    for (size_t j = 0; j < constraint_matrix.size(); ++j)
      dy_da[n+j] = vec6_dot(constraint_matrix[j], du);
    return y;
  }
};

// only for testing and evaluation - scaling with NLOpt
#if WITH_NLOPT
namespace impl {

template<typename Real>
double calculate_for_nlopt(unsigned n, const double* x, double* grad, void* data) {
  auto scaling = static_cast<Scaling<Real>*>(data);
  scaling->set_parameters(x);
  if (grad)
    return compute_gradients(*scaling, n, grad);
  else
    return compute_wssr(*scaling);
}

inline const char* nlresult_to_string(nlopt_result r) {
  switch (r) {
    case NLOPT_FAILURE: return "failure";
    case NLOPT_INVALID_ARGS: return "invalid arguments";
    case NLOPT_OUT_OF_MEMORY: return "out of memory";
    case NLOPT_ROUNDOFF_LIMITED: return "roundoff errors limit progress";
    case NLOPT_FORCED_STOP: return "interrupted";
    case NLOPT_SUCCESS: return "success";
    case NLOPT_STOPVAL_REACHED: return "stop-value reached";
    case NLOPT_FTOL_REACHED: return "ftol-value reached";
    case NLOPT_XTOL_REACHED: return "xtol-value reached";
    case NLOPT_MAXEVAL_REACHED: return "max. evaluation number reached";
    case NLOPT_MAXTIME_REACHED: return "max. time reached";
    default: return "<unknown result>";
  }
  unreachable();
}

} // namespace impl

template<typename Real>
double fit_parameters_with_nlopt(Scaling<Real>& scaling, const char* optimizer) {
  std::vector<double> params = scaling.get_parameters();
  nlopt_opt opt = nlopt_create(nlopt_algorithm_from_string(optimizer), params.size());
  {  // prepare bounds
    std::vector<double> lb(params.size());
    std::vector<double> ub(params.size());
    lb[0] = 0.7 * params[0];  // k_ov
    ub[0] = 1.2 * params[0];
    size_t n = 1;
    if (scaling.use_solvent) {
      if (!scaling.fix_k_sol) {
        lb[n] = 0.15;
        ub[n] = 0.5;
        ++n;
      }
      if (!scaling.fix_b_sol) {
        lb[n] = 10;
        ub[n] = 80;
        ++n;
      }
    }
    for (; n < params.size(); ++n) {
      lb[n] = -0.01;
      ub[n] = 0.01;
    }
    nlopt_set_lower_bounds(opt, lb.data());
    nlopt_set_upper_bounds(opt, ub.data());
  }
  nlopt_set_min_objective(opt, impl::calculate_for_nlopt<Real>, &scaling);
  nlopt_set_maxeval(opt, 100);
  double minf = NAN;
  nlopt_result r = nlopt_optimize(opt, &params[0], &minf);
  (void) r;
  //if (r < 0)
  //  printf("NLopt result: %s\n", impl::nlresult_to_string(r));
  //else
  //  printf("NLopt minimum value: %g\n", minf);
  nlopt_destroy(opt);
  scaling.set_parameters(params);
  return minf;
}
#endif  // WITH_NLOPT

} // namespace gemmi
#endif
