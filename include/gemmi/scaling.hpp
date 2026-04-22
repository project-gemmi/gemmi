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

/// @brief Type alias for a symmetric 3×3 tensor represented as 6 coefficients.
/// Stores (u11, u22, u33, u12, u13, u23).
using Vec6 = std::array<double, 6>;

/// @brief Dot product of Vec6 with a symmetric 3×3 matrix.
/// Used for tensor contractions in ADP refinement.
/// @param a Vec6 vector (6-element array).
/// @param s Symmetric 3×3 matrix.
/// @return Dot product result.
inline double vec6_dot(const Vec6& a, const SMat33<double>& s) {
  return a[0] * s.u11 + a[1] * s.u22 + a[2] * s.u33
       + a[3] * s.u12 + a[4] * s.u13 + a[5] * s.u23;
}

/// @brief Return the symmetry-adapted constraint vectors for the ADP tensor in a space group.
/// The number and content of returned constraint vectors depend on the crystal system.
/// For example, cubic returns one vector [1 1 1 0 0 0] (isotropic); tetragonal returns
/// two vectors in directions [1 1 0 0 0 0] and [0 0 1 0 0 0] (diagonal u11=u22, u33 free).
/// @param sg Pointer to the space group; if null, triclinic symmetry is assumed (6 free components).
/// @return Vector of normalized Vec6 constraint vectors for independent ADP components.
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

/// @brief Anisotropic scaling of calculated structure factors to observed data.
/// Optionally includes bulk solvent correction: Fc + k_sol·exp(-b_sol·stol²)·Fmask.
/// Parameter refinement uses Levenberg-Marquardt (or NLopt if WITH_NLOPT is defined).
/// @tparam Real Floating-point type (float or double).
template<typename Real>
struct Scaling {
  /// @brief One reflection used in the least-squares fit.
  struct Point {
    Miller hkl;                            ///< Miller indices.
    double stol2;                          ///< (sin θ/λ)² for this reflection.
    std::complex<Real> fcmol, fmask;       ///< Calculated molecular and mask structure factors.
    Real fobs, sigma;                      ///< Observed amplitude and its standard uncertainty.

    Miller get_x() const { return hkl; }
    double get_y() const { return fobs; }
    double get_weight() const { return 1.0 /* / sigma*/; }
  };

  UnitCell cell;                           ///< Unit cell parameters.
  double k_overall = 1.;                   ///< Overall scale factor.
  SMat33<double> b_star{0, 0, 0, 0, 0, 0}; ///< Anisotropic B* tensor in reciprocal space.
  std::vector<Vec6> constraint_matrix;     ///< Symmetry constraints on b_star.
  bool use_solvent = false;                ///< If true, include bulk solvent correction.
  bool fix_k_sol = false;                  ///< If true, do not refine k_sol.
  bool fix_b_sol = false;                  ///< If true, do not refine b_sol.
  double k_sol = 0.35;                     ///< Bulk solvent scale factor.
  double b_sol = 46.0;                     ///< Bulk solvent B-factor.
  std::vector<Point> points;               ///< Reflection data for fitting.

  /// @brief Initialize with unit cell and space group.
  /// Sets up constraint_matrix from crystal system symmetry.
  /// @param cell_ Unit cell parameters.
  /// @param sg Pointer to space group (for symmetry constraints).
  Scaling(const UnitCell& cell_, const SpaceGroup* sg)
      : cell(cell_), constraint_matrix(adp_symmetry_constraints(sg)) {}

  /// @brief Set b_star from a real-space B matrix.
  /// Converts from Cartesian coordinates to reciprocal-space B*.
  /// @param b_overall B matrix in real (Cartesian) space.
  void set_b_overall(const SMat33<double>& b_overall) {
    b_star = b_overall.transformed_by(cell.frac.mat);
  }

  /// @brief Get b_star converted back to real (Cartesian) space.
  /// @return B matrix in Cartesian coordinates.
  SMat33<double> get_b_overall() const {
    return b_star.transformed_by(cell.orth.mat);
  }

  /// @brief Apply scaling to all reflections in-place, optionally adding bulk solvent correction.
  /// @param asu_data ASU data to scale (modified in-place).
  /// @param mask_data Mask (solvent) data; required if use_solvent is true, otherwise may be null.
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

  /// @brief Compute scaled value for one reflection.
  /// Applies overall scale factor and optionally bulk solvent correction.
  /// @param hkl Miller indices.
  /// @param f_value Calculated molecular structure factor.
  /// @param mask_value Mask (solvent) structure factor.
  /// @return Scaled complex structure factor.
  std::complex<Real> scale_value(const Miller& hkl, std::complex<Real> f_value,
                                 std::complex<Real> mask_value) {
    if (use_solvent) {
      double stol2 = cell.calculate_stol_sq(hkl);
      f_value += (Real)get_solvent_scale(stol2) * mask_value;
    }
    return f_value * (Real) get_overall_scale_factor(hkl);
  }

  /// @brief Return current parameters as a flat vector.
  /// Includes k_overall, b_star components (via constraint_matrix), and optionally k_sol and b_sol.
  /// @return Vector of parameters in order: k_overall, [k_sol], [b_sol], b_star_components.
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

  /// @brief Set parameters from pointer or vector.
  /// Updates k_overall, b_star components (via constraint_matrix), and optionally k_sol and b_sol.
  /// @param p Pointer to parameter array (order: k_overall, [k_sol], [b_sol], b_star_components).
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

  /// @brief Set parameters from vector.
  /// @param p Vector of parameters.
  void set_parameters(const std::vector<double>& p) {
    set_parameters(p.data());
  }

  /// @brief Populate points from matching reflections in calc and obs datasets.
  /// Precondition: all AsuData arguments must be sorted by Miller indices.
  /// @param calc Calculated structure factors.
  /// @param obs Observed amplitudes and sigmas.
  /// @param mask_data Mask structure factors; required if use_solvent is true.
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


  /// @brief Compute k_sol * exp(-b_sol * stol²) for bulk solvent correction.
  /// @param stol2 (sin θ/λ)² value.
  /// @return Solvent scale factor.
  double get_solvent_scale(double stol2) const {
    return k_sol * std::exp(-b_sol * stol2);
  }

  /// @brief Compute k_overall * exp(-b_star:hkl) for overall scale factor.
  /// @param hkl Miller indices.
  /// @return Overall scale factor including anisotropic B*.
  double get_overall_scale_factor(const Miller& hkl) const {
    return k_overall * std::exp(-0.25 * b_star.r_u_r(hkl));
  }

  /// @brief Compute total calculated structure factor for point p.
  /// Includes molecular part and optionally bulk solvent correction.
  /// @param p Reflection point.
  /// @return Total Fcalc = Fcmol + k_sol·exp(-b_sol·stol²)·Fmask (or just Fcmol if no solvent).
  std::complex<Real> get_fcalc(const Point& p) const {
    if (!use_solvent)
      return p.fcmol;
    return p.fcmol + (Real)get_solvent_scale(p.stol2) * p.fmask;
  }

  /// @brief Quick linear least-squares fit of isotropic B only (ignoring sigma).
  /// Used for initialization before full refinement.
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

  /// @brief Compute optimal k_overall by linear least squares.
  /// @return Optimal scale factor.
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

  /// @brief Fit anisotropic B* by linear approximation (for testing/initialization only).
  /// DO NOT USE for production: symmetry constraints are not implemented.
  /// Based on P. Afonine et al., doi:10.1107/S0907444913000462, section 2.1.
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

  /// @brief Full Levenberg-Marquardt optimization of all parameters.
  /// @return Final R factor.
  double fit_parameters() {
    LevMar levmar;
    return levmar.fit(*this);
  }

  /// @brief Compute R-factor: Σ|Fobs - Fcalc| / Σ|Fobs|.
  /// @return R-factor value.
  double calculate_r_factor() const {
    double abs_diff_sum = 0;
    double denom = 0;
    for (const Point& p : points) {
      abs_diff_sum += std::fabs(p.fobs - compute_value(p));
      denom += p.fobs;
    }
    return abs_diff_sum / denom;
  }

  /// @brief Compute Fcalc for point p (interface for fitting).
  /// @param p Reflection point.
  /// @return Scaled calculated structure factor amplitude.
  double compute_value(const Point& p) const {
    return std::abs(get_fcalc(p)) * (Real) get_overall_scale_factor(p.hkl);
  }

  /// @brief Compute Fcalc and its derivatives with respect to all parameters.
  /// @param p Reflection point.
  /// @param dy_da Output vector to fill with derivatives [dy/dk_overall, dy/dk_sol, dy/db_sol, dy/db_star_components].
  /// @return Calculated structure factor amplitude.
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

/// @brief Fit scaling parameters using NLopt with the named optimizer.
/// Only available when compiled with WITH_NLOPT.
/// @tparam Real Floating-point type (float or double).
/// @param scaling The Scaling object to optimize.
/// @param optimizer Name of the optimizer algorithm (NLopt algorithm string).
/// @return Final residual value (WSSR).
/// @note This function is for testing and evaluation purposes.
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
