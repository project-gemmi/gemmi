// Copyright 2020 Global Phasing Ltd.
//
// Least-squares fitting - Levenberg-Marquardt method.
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
inline void jordan_solve(std::vector<double>& a, std::vector<double>& b) {
  assert(a.size() == b.size() * b.size());
  int n = (int) b.size();
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
          fail("Trying to reverse singular matrix. Column ", i, " is zeroed.");
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
  double initial_wssr;
  int eval_count;  // number of function evaluations

  // arrays used during refinement
  std::vector<double> alpha; // matrix
  std::vector<double> beta;  // vector
  std::vector<double> temp_alpha, temp_beta; // working arrays


  template<typename Model, typename Point>
  double fit(Model& model, const std::vector<Point>& data) {
    eval_count = 0;
    std::vector<double> initial_a = model.get_parameters();
    initial_wssr = this->compute_wssr(model.compute_values(data), data);
    std::vector<double> best_a = initial_a;
    size_t na = initial_a.size();

    double lambda = lambda_start;
    alpha.resize(na * na);
    beta.resize(na);

    double wssr = initial_wssr;
    this->compute_derivatives(model, data);

    int small_change_counter = 0;
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

      for (size_t i = 0; i < na; i++) {
#ifdef GEMMI_DEBUG_LEVMAR
        fprintf(stderr, "%+g ", temp_beta[i]);
#endif
        // put new a[] into temp_beta[]
        temp_beta[i] = best_a[i] + temp_beta[i];
      }

      model.set_parameters(temp_beta);
      double new_wssr = this->compute_wssr(model.compute_values(data), data);
#ifdef GEMMI_DEBUG_LEVMAR
      fprintf(stderr, "\nWSSR=%g (%g%%) lambda=%g\n",
              new_wssr, 100. * new_wssr / initial_wssr, lambda);
#endif
      if (new_wssr < wssr) {
        double rel_change = (wssr - new_wssr) / wssr;
        wssr = new_wssr;
        best_a = temp_beta;

        if (wssr == 0)
          break;
        // termination criterium: negligible change of wssr
        if (rel_change < stop_rel_change) {
          if (++small_change_counter >= 2)
            break;
        } else {
          small_change_counter = 0;
        }
        this->compute_derivatives(model, data);
        lambda *= lambda_down_factor;
      } else { // worse fitting
        if (lambda > lambda_limit) // termination criterium: large lambda
          break;
        lambda *= lambda_up_factor;
      }
    }

    model.set_parameters(wssr < initial_wssr ? best_a : initial_a);
    return wssr;
  }

private:
  template<typename Model, typename Point>
  void compute_derivatives(const Model& model, const std::vector<Point>& data) {
    assert(alpha.size() == beta.size() * beta.size());
    int na = (int)beta.size();
    assert(na != 0);
    fill(alpha.begin(), alpha.end(), 0.0);
    fill(beta.begin(), beta.end(), 0.0);
    // Iterating over points is tiled to limit memory usage. It's also a little
    // faster than a single loop over all points for large number of points.
    const int kMaxTileSize = 1024;
    std::vector<double> dy_da;
    int n = (int)data.size();
    for (int tstart = 0; tstart < n; tstart += kMaxTileSize) {
      const int der_size = na; // + 1;
      int tsize = std::min(n - tstart, kMaxTileSize);
      std::vector<typename Model::Point> xx(data.begin() + tstart,
                                            data.begin() + tstart + tsize);
      std::vector<double> yy(tsize, 0.);
      dy_da.resize(tsize * der_size);
      fill(dy_da.begin(), dy_da.end(), 0.);
      model.compute_values_and_derivatives(xx, yy, dy_da);
      for (int i = 0; i != tsize; ++i) {
        double weight = data[tstart + i].get_weight();
        double dy_sig = weight * (data[tstart + i].get_y() - yy[i]);
        double* t = &dy_da[i * der_size];
        for (int j = 0; j != na; ++j) {
          if (t[j] != 0) {
            t[j] *= weight;
            for (int k = j; k != -1; --k)
              alpha[na * j + k] += t[j] * t[k];
            beta[j] += dy_sig * t[j];
          }
        }
      }
    }

    // Only half of the alpha matrix was filled above. Fill the rest.
    for (int j = 1; j < na; j++)
      for (int k = 0; k < j; k++)
        alpha[na * k + j] = alpha[na * j + k];
  }

  template<typename Point>
  double compute_wssr(const std::vector<double>& yy, const std::vector<Point>& data) {
    // long double here notably increases the accuracy of calculations
    long double wssr = 0;
    for (int j = 0; j < (int)yy.size(); j++) {
      double dy = data[j].get_weight() * (data[j].get_y() - yy[j]);
      wssr += dy * dy;
    }
    ++eval_count;
    return (double) wssr;
  }
};

} // namespace gemmi
#endif
