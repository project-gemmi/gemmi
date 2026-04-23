// Copyright 2020 Global Phasing Ltd.
//
// Least-squares fitting - Levenberg-Marquardt method.
//
// Based on the code from fityk (but here it's under MPL 2.0).

#ifndef GEMMI_LEVMAR_HPP_
#define GEMMI_LEVMAR_HPP_

#include <cassert>
#include <cmath>      // for fabs
#include <algorithm>  // for min
#include <vector>
#include "fail.hpp"   // for fail
#include "math.hpp"   // for sq

//#define GEMMI_DEBUG_LEVMAR

namespace gemmi {

/// This function solves a set of linear algebraic equations using
/// Gauss-Jordan elimination with partial pivoting.
///
/// A * x = b
///
/// a is n x n matrix (in vector)
/// b is vector of length n,
/// This function returns vector x[] in b[], and 1-matrix in a[].
/// @brief Solve a linear system using Gauss-Jordan elimination with partial pivoting.
/// @details
/// Solves A * x = b by reducing A to the identity matrix via row operations.
/// Handles singular matrices by skipping zero rows/columns.
/// @param a n x n matrix stored in row-major order (modified to identity matrix).
/// @param b Right-hand side vector of length n (modified to solution x).
/// @param n System size.
inline void jordan_solve(double* a, double* b, int n) {
  for (int i = 0; i < n; i++) {
    // looking for a pivot element
    int maxnr = -1;
    double amax = 0;
    for (int j = i; j < n; j++) {
      double aji = std::fabs(a[n * j + i]);
      if (aji > amax) {
        maxnr = j;
        amax = aji;
      }
    }
    // handle singular matrix
    if (maxnr == -1) {
      // i-th column has only zeros.
      // If it's the same about i-th row, and b[i]==0, let x[i]==0.
      for (int j = i; j < n; j++)
        if (a[n * i + j] != 0. || b[i] != 0.)
          fail("Trying to reverse singular matrix. Column ", std::to_string(i), " is zeroed.");
      continue; // x[i]=b[i], b[i]==0
    }
    // interchanging rows
    if (maxnr != i) {
      for (int j = i; j < n; j++)
        std::swap(a[n * maxnr + j], a[n * i + j]);
      std::swap(b[i], b[maxnr]);
    }
    // divide by a_ii -- to get a_ii=1
    double c = 1.0 / a[i * n + i];
    for (int j = i; j < n; j++)
      a[i * n + j] *= c;
    b[i] *= c;
    // subtract -- to zero all remaining elements of this row
    for (int k = 0; k < n; k++)
      if (k != i) {
        double d = a[k * n + i];
        for (int j = i; j < n; j++)
          a[k * n + j] -= a[i * n + j] * d;
        b[k] -= b[i] * d;
      }
  }
}

/// @brief Solve a linear system stored in vector containers.
/// @param a n x n matrix as flat vector (n^2 elements).
/// @param b Right-hand side vector of length n.
inline void jordan_solve(std::vector<double>& a, std::vector<double>& b) {
  assert(a.size() == b.size() * b.size());
  jordan_solve(a.data(), b.data(), (int)b.size());
}

/// @brief Print parameters to stderr (debug only).
/// @param name Label for the parameter set.
/// @param a Vector of parameters to print.
inline void print_parameters(const std::string& name, std::vector<double> &a) {
  fprintf(stderr, " %s:", name.c_str());
  for (double& x : a)
    fprintf(stderr, " %g", x);
  fprintf(stderr, "\n");
}

/// @brief Compute weighted sum of squared residuals for a target.
/// @details
/// Uses long double for accumulation to improve numerical accuracy.
/// Assumes Target provides: points container, get_weight(), get_y(), compute_value().
/// @tparam Target Fitting target type.
/// @param target The fitting target.
/// @return Sum of weighted squared residuals.
template<typename Target>
double compute_wssr(const Target& target) {
  long double wssr = 0; // long double here notably increases the accuracy
  for (const auto& p : target.points)
    wssr += sq(p.get_weight() * (p.get_y() - target.compute_value(p)));
  return (double) wssr;
}

/// @brief Compute function value, residuals, and gradients with respect to parameters.
/// @details
/// Assumes Target provides: points container, get_weight(), get_y(),
/// compute_value_and_derivatives(point, dy_da_vector).
/// @tparam Target Fitting target type.
/// @param target The fitting target.
/// @param n Number of parameters.
/// @param grad Output array of size n for partial derivatives d(wssr)/da.
/// @return Weighted sum of squared residuals.
template<typename Target>
double compute_gradients(const Target& target, unsigned n, double* grad) {
  double wssr = 0;
  for (unsigned i = 0; i < n; ++i)
    grad[i] = 0;
  std::vector<double> dy_da(n);
  for (const auto& p : target.points) {
    double y = target.compute_value_and_derivatives(p, dy_da);
    double dy = p.get_weight() * (p.get_y() - y);
    wssr += sq(dy);
    for (unsigned i = 0; i < n; ++i)
      grad[i] += -2 * dy * dy_da[i];
  }
#if 0
  // calculate numerical derivatives to check analytical formulas
  fprintf(stderr, ">>  y=%g\n", wssr);
  std::vector<double> x = target.get_parameters();
  assert(x.size() == n);
  for (unsigned i = 0; i < n; ++i) {
    double x_orig = x[i];
    double h = std::max(std::fabs(x[i]), 1e-6) * 1e-3;
    x[i] = x_orig - h;
    const_cast<Target&>(target).set_parameters(x);
    double y_left = compute_wssr(target);
    x[i] = x_orig + h;
    const_cast<Target&>(target).set_parameters(x);
    double y_right = compute_wssr(target);
    double numeric = (y_right - y_left) / (2 * h);
    x[i] = x_orig;
    double m = std::max(std::fabs(grad[i]), std::fabs(numeric));
    if (m > 1e-3 && std::fabs(grad[i] - numeric) > 0.02 * m)
      fprintf(stderr, "!! grad[%u]: %g vs %g  (value: %g)\n", i, grad[i], numeric, x[i]);
  }
  const_cast<Target&>(target).set_parameters(x);
#endif
  return wssr;
}

/// @brief Compute Jacobian-based matrices for Levenberg-Marquardt algorithm.
/// @details
/// Computes the normal equations: alpha = J^T*J (approximates Hessian),
/// beta = J^T*residual. These are the building blocks for iterative refinement.
/// Alpha is initially undamped; the LM algorithm applies the damping factor
/// (1 + lambda) to diagonal elements. Both matrices use only the lower triangle
/// and are symmetrized after computation.
/// Assumes Target provides: points container, get_weight(), get_y(),
/// compute_value_and_derivatives(point, dy_da_vector).
/// @tparam Target Fitting target type.
/// @param target The fitting target.
/// @param alpha Output: n x n matrix (stored as flat vector) = J^T*J.
/// @param beta Output: n-element vector = J^T*residual.
/// @return Weighted sum of squared residuals.
template<typename Target>
double compute_lm_matrices(const Target& target,
                           std::vector<double>& alpha,
                           std::vector<double>& beta) {
  assert(!beta.empty());
  assert(alpha.size() == beta.size() * beta.size());
  long double wssr = 0; // long double here notably increases the accuracy
  size_t na = beta.size();
  std::fill(alpha.begin(), alpha.end(), 0.0);
  std::fill(beta.begin(), beta.end(), 0.0);
  std::vector<double> dy_da(na);
  for (const auto& p : target.points) {
    double y = target.compute_value_and_derivatives(p, dy_da);
    double weight = p.get_weight();
    double dy_sig = weight * (p.get_y() - y);
    for (size_t j = 0; j != na; ++j) {
      if (dy_da[j] != 0) {
        dy_da[j] *= weight;
        for (size_t k = j+1; k-- != 0;)
          alpha[na * j + k] += dy_da[j] * dy_da[k];
        beta[j] += dy_sig * dy_da[j];
      }
    }
    wssr += sq(dy_sig);
  }

  // Only half of the alpha matrix was filled above. Fill the rest.
  for (size_t j = 1; j < na; j++)
    for (size_t k = 0; k < j; k++)
      alpha[na * k + j] = alpha[na * j + k];
  return (double) wssr;
}

/// @brief Levenberg-Marquardt non-linear least-squares optimization.
/// @details
/// Implements the Levenberg-Marquardt algorithm for fitting model parameters
/// to minimize the sum of weighted squared residuals.
/// The algorithm adjusts a damping factor (lambda) to interpolate between
/// gradient descent (large lambda) and Newton's method (small lambda),
/// automatically selecting the step size that gives best improvement.
struct LevMar {
  /// @brief Maximum number of function evaluations before terminating.
  int eval_limit = 100;
  /// @brief Stop optimization if damping factor lambda exceeds this value.
  double lambda_limit = 1e+15;
  /// @brief Stop if relative change in WSSR falls below this for two consecutive iterations.
  double stop_rel_change = 1e-5;

  /// @brief Factor by which lambda is multiplied if fit worsens.
  double lambda_up_factor = 10;
  /// @brief Factor by which lambda is multiplied if fit improves.
  double lambda_down_factor = 0.1;
  /// @brief Initial damping factor (typically 0.001).
  double lambda_start = 0.001;

  /// @brief Initial weighted sum of squared residuals (set by fit()).
  double initial_wssr = NAN;
  /// @brief Number of function evaluations performed (set by fit()).
  int eval_count = 0;

  /// @brief Jacobian-based normal equations matrix (J^T*J), size n*n.
  std::vector<double> alpha;
  /// @brief Jacobian-based normal equations vector (J^T*residual), size n.
  std::vector<double> beta;
  /// @brief Working copies of alpha and beta during refinement.
  std::vector<double> temp_alpha, temp_beta;

  /// @brief Run Levenberg-Marquardt optimization.
  /// @details
  /// Iteratively refines parameters to minimize WSSR. At each iteration,
  /// solves (J^T*J + lambda*diag(J^T*J)) * da = J^T*residual for step da,
  /// then updates parameters and checks for improvement.
  /// Terminates when eval_limit is reached, lambda exceeds lambda_limit,
  /// or relative improvement drops below stop_rel_change for two iterations.
  /// @tparam Target Fitting target type (must implement: get_parameters(),
  ///         set_parameters(), points container, and
  ///         compute_value_and_derivatives()).
  /// @param target The fitting target (modified in place).
  /// @return Final weighted sum of squared residuals. initial_wssr and eval_count
  ///         are also set and can be inspected after fit() returns.
  template<typename Target>
  double fit(Target& target) {
    std::vector<double> initial_a = target.get_parameters();
#ifdef GEMMI_DEBUG_LEVMAR
    print_parameters("ini", initial_a);
#endif
    std::vector<double> best_a = initial_a;
    size_t na = initial_a.size();

    double lambda = lambda_start;
    alpha.resize(na * na);
    beta.resize(na);

    initial_wssr = compute_lm_matrices(target, alpha, beta);
    double wssr = initial_wssr;

    int small_change_counter = 0;
    eval_count = 1;  // number of function evaluations so far
    for (int iter = 0; ; iter++) {
      if (eval_limit > 0 && eval_count >= eval_limit)
        break;

      // prepare next parameters -> temp_beta
      temp_alpha = alpha;
      // Using '*=' not '+=' below applies the dampling factor as:
      // J^T J + lambda * diag(J^T J); not ... + lambda * I.
      for (size_t j = 0; j < na; j++)
        temp_alpha[na * j + j] *= (1.0 + lambda);
      temp_beta = beta;

      // Matrix solution (Ax=b)  temp_alpha * da == temp_beta
      jordan_solve(temp_alpha, temp_beta);

      for (size_t i = 0; i < na; i++)
        // put new a[] into temp_beta[]
        temp_beta[i] += best_a[i];

      target.set_parameters(temp_beta);
      double new_wssr = compute_wssr(target);
      ++eval_count;
#ifdef GEMMI_DEBUG_LEVMAR
      fprintf(stderr, " #%d WSSR=%.8g %+g%% (%+.4g%%) lambda=%g\n",
              iter, new_wssr, 100. * (new_wssr / initial_wssr - 1.),
              100. * (new_wssr / wssr - 1.), lambda);
      if (new_wssr < wssr)
        print_parameters("", temp_beta);
#else
      (void) iter;
#endif
      if (new_wssr < wssr) {
        double rel_change = (wssr - new_wssr) / wssr;
        wssr = new_wssr;
        best_a = temp_beta;

        if (wssr == 0)
          break;
        // termination criterion: negligible change of wssr
        if (rel_change < stop_rel_change) {
          if (++small_change_counter >= 2)
            break;
        } else {
          small_change_counter = 0;
        }
        compute_lm_matrices(target, alpha, beta);
        ++eval_count;
        lambda *= lambda_down_factor;
      } else { // worse fitting
        if (lambda > lambda_limit) // termination criterion: large lambda
          break;
        lambda *= lambda_up_factor;
      }
    }

    target.set_parameters(wssr < initial_wssr ? best_a : initial_a);
    return wssr;
  }
};

} // namespace gemmi
#endif
