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

inline void jordan_solve(std::vector<double>& a, std::vector<double>& b) {
  assert(a.size() == b.size() * b.size());
  jordan_solve(a.data(), b.data(), (int)b.size());
}


#ifdef GEMMI_DEBUG_LEVMAR
inline void print_parameters(const std::string& name, std::vector<double> &a) {
  fprintf(stderr, " %s:", name.c_str());
  for (double& x : a)
    fprintf(stderr, " %g", x);
  fprintf(stderr, "\n");
}
#endif

struct LevMar {
  // termination criteria
  int eval_limit = 100;
  double lambda_limit = 1e+15;
  double stop_rel_change = 1e-5;

  // adjustable parameters (normally the default values work fine)
  double lambda_up_factor = 10;
  double lambda_down_factor = 0.1;
  double lambda_start = 0.001;

  // values set in fit() that can be inspected later
  double initial_wssr = NAN;
  int eval_count = 0;  // number of function evaluations

  // arrays used during refinement
  std::vector<double> alpha; // matrix
  std::vector<double> beta;  // vector
  std::vector<double> temp_alpha, temp_beta; // working arrays


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

    initial_wssr = this->compute_derivatives(target);
    double wssr = initial_wssr;

    int small_change_counter = 0;
    eval_count = 1;  // number of function evaluations so far
    for (int iter = 0; ; iter++) {
      if (eval_limit > 0 && eval_count >= eval_limit)
        break;

      // prepare next parameters -> temp_beta
      temp_alpha = alpha;
      for (size_t j = 0; j < na; j++)
        temp_alpha[na * j + j] *= (1.0 + lambda);
      temp_beta = beta;

      // Matrix solution (Ax=b)  temp_alpha * da == temp_beta
      jordan_solve(temp_alpha, temp_beta);

      for (size_t i = 0; i < na; i++)
        // put new a[] into temp_beta[]
        temp_beta[i] += best_a[i];

      target.set_parameters(temp_beta);
      double new_wssr = this->compute_wssr(target);
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
        this->compute_derivatives(target);
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

private:
  template<typename Target>
  double compute_derivatives(const Target& target) {
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

  template<typename Target>
  double compute_wssr(const Target& target) {
    long double wssr = 0; // long double here notably increases the accuracy
    for (const auto& p : target.points)
      wssr += sq(p.get_weight() * (p.get_y() - target.compute_value(p)));
    return (double) wssr;
  }
};

} // namespace gemmi
#endif
