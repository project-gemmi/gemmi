// Copyright 2024 Global Phasing Ltd.
//
// testing Scaling with NLOpt

#include "scaling.hpp"
#include <nlopt.h>

namespace gemmi {

namespace impl {

template<typename Real>
double calculate_lsq(unsigned n, const double* x, double* grad, void* data) {
  auto scaling = static_cast<Scaling<Real>*>(data);
  scaling->set_parameters(x);
  double y = 0;
  if (grad)
    for (unsigned i = 0; i < n; ++i)
      grad[i] = 0;
  for (const typename Scaling<Real>::Point& pt : scaling->points) {
    Vec3 h(pt.hkl);
    double kaniso = std::exp(-0.25 * scaling->b_star.r_u_r(h));
    double fcalc_abs;
    double grad_k_sol;
    double grad_b_sol;
    if (scaling->use_solvent) {
      double solv_b = std::exp(-scaling->b_sol * pt.stol2);
      double solv_scale = scaling->k_sol * solv_b;
      auto fcalc = pt.fcmol + (Real) solv_scale * pt.fmask;
      fcalc_abs = std::abs(fcalc);
      if (grad) {
        double dy_dsol = (fcalc.real() * pt.fmask.real() +
                          fcalc.imag() * pt.fmask.imag())
                         / fcalc_abs * scaling->k_overall * kaniso;
        grad_k_sol = solv_b * dy_dsol;
        grad_b_sol = -pt.stol2 * solv_scale * dy_dsol;
      }
    } else {
      fcalc_abs = std::abs(pt.fcmol);
    }
    double fe = fcalc_abs * kaniso;
    double dy = pt.fobs - scaling->k_overall * fe;
    y += dy * dy;
    if (grad) {
      double m2dy_fe = -2 * dy * fe;
      double factor = -0.25 * (scaling->k_overall * m2dy_fe);
      grad[0] += m2dy_fe; // dy/d scaling->k_overall
      unsigned offset = 1;
      if (scaling->use_solvent) {
        if (!scaling->fix_k_sol)
          grad[offset++] += -2 * dy * grad_k_sol;
        if (!scaling->fix_b_sol)
          grad[offset++] += -2 * dy * grad_b_sol;
      }
      SMat33<double> du = {
        factor * (h.x * h.x),
        factor * (h.y * h.y),
        factor * (h.z * h.z),
        factor * 2 * (h.x * h.y),
        factor * 2 * (h.x * h.z),
        factor * 2 * (h.y * h.z),
      };
      for (size_t j = 0; j < scaling->constraint_matrix.size(); ++j)
        grad[offset++] += vec6_dot(scaling->constraint_matrix[j], du);
      assert(offset == n);
    }
  }
#if 0
  if (grad) {
    fprintf(stderr, "nlopt>  y=%g\n", y);
    std::vector<double> x_copy(x, x+n);
    for (unsigned i = 0; i < n; ++i) {
      double h = std::max(std::fabs(x[i]), 1e-6) * 1e-3;
      x_copy[i] = x[i] - h;
      double y_left = calculate_lsq<Real>(n, x_copy.data(), nullptr, data);
      x_copy[i] = x[i] + h;
      double y_right = calculate_lsq<Real>(n, x_copy.data(), nullptr, data);
      double numeric = (y_right - y_left) / (2 * h);
      x_copy[i] = x[i];
      //double m = std::max(std::fabs(grad[i]), std::fabs(numeric));
      //if (m > 1e-3 && std::fabs(grad[i] - numeric) > 0.02 * m)
      fprintf(stderr, "!! grad[%u]: %g vs %g  (value: %g)\n", i, grad[i], numeric, x[i]);
    }
  }
#endif
  return y;
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
  nlopt_set_min_objective(opt, impl::calculate_lsq<Real>, &scaling);
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

} // namespace gemmi
