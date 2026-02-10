//! @file
//! @brief Statistics utilities: classes Variance, Covariance, Correlation, DataStats.
//!
//! Provides single-pass algorithms for calculating statistical measures from
//! streaming data. Useful for accumulating statistics without storing all data points.

// Copyright 2018 Global Phasing Ltd.

#ifndef GEMMI_STATS_HPP_
#define GEMMI_STATS_HPP_

#include <cstddef>  // for size_t
#include <cmath>    // for sqrt, NAN, INFINITY
#include <vector>

namespace gemmi {

//! @brief Single-pass variance calculation using Welford's algorithm.
//!
//! Popular single-pass algorithm for calculating variance and mean without storing all values.
//! Uses numerically stable update formulas. Can calculate both sample and population variance.
struct Variance {
  int n = 0;         //!< Number of data points
  double sum_sq = 0.;//!< Sum of squared deviations
  double mean_x = 0.;//!< Running mean

  Variance() = default;

  //! @brief Construct and initialize from iterator range.
  //! @tparam T Iterator type
  //! @param begin Start iterator
  //! @param end End iterator
  template <typename T> Variance(T begin, T end) : Variance() {
    for (auto i = begin; i != end; ++i)
      add_point(*i);
  }

  //! @brief Add data point to accumulator.
  //! @param x Data value
  void add_point(double x) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    sum_sq += dx * (x - mean_x);
  }

  //! @brief Calculate sample variance (divides by n-1).
  //! @return Sample variance
  double for_sample() const { return sum_sq / (n - 1); }

  //! @brief Calculate population variance (divides by n).
  //! @return Population variance
  double for_population() const { return sum_sq / n; }
};

//! @brief Single-pass covariance calculation.
//!
//! Extends Variance to calculate covariance between two variables.
//! Accumulates data for both x and y, computing covariance in one pass.
struct Covariance : Variance {
  double mean_y = 0.;  //!< Running mean of y values

  //! @brief Add data point pair to accumulator.
  //! @param x First variable value
  //! @param y Second variable value
  void add_point(double x, double y) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    mean_y += (y - mean_y) / n;
    sum_sq += dx * (y - mean_y);
  }
};

//! @brief Single-pass correlation and regression calculation.
//!
//! Calculates Pearson correlation coefficient, variances, covariance, and
//! linear regression line (slope/intercept) in a single pass through the data.
struct Correlation {
  int n = 0;          //!< Number of data points
  double sum_xx = 0.; //!< Sum of squared x deviations
  double sum_yy = 0.; //!< Sum of squared y deviations
  double sum_xy = 0.; //!< Sum of cross products
  double mean_x = 0.; //!< Running mean of x
  double mean_y = 0.; //!< Running mean of y

  //! @brief Add data point pair to accumulator.
  //! @param x First variable value
  //! @param y Second variable value
  void add_point(double x, double y) {
    ++n;
    double weight = (double)(n - 1) / n;
    double dx = x - mean_x;
    double dy = y - mean_y;
    sum_xx += weight * dx * dx;
    sum_yy += weight * dy * dy;
    sum_xy += weight * dx * dy;
    mean_x += dx / n;
    mean_y += dy / n;
  }

  //! @brief Calculate Pearson correlation coefficient.
  //! @return Correlation coefficient (-1 to 1)
  double coefficient() const { return sum_xy / std::sqrt(sum_xx * sum_yy); }

  //! @brief Calculate variance of x.
  //! @return Variance of x
  double x_variance() const { return sum_xx / n; }

  //! @brief Calculate variance of y.
  //! @return Variance of y
  double y_variance() const { return sum_yy / n; }

  //! @brief Calculate covariance.
  //! @return Covariance of x and y
  double covariance() const { return sum_xy / n; }

  //! @brief Calculate ratio of means.
  //! @return mean_y / mean_x
  double mean_ratio() const { return mean_y / mean_x; }

  //! @brief Calculate slope of regression line y = slope * x + intercept.
  //! @return Slope
  double slope() const { return sum_xy / sum_xx; }

  //! @brief Calculate intercept of regression line y = slope * x + intercept.
  //! @return Intercept
  double intercept() const { return mean_y - slope() * mean_x; }
};

//! @brief Combine two correlation accumulators.
//! @param a First correlation
//! @param b Second correlation
//! @return Combined correlation
//!
//! Merges statistics from two independent Correlation objects into a single one.
//! Useful for parallel computation or combining statistics from different data subsets.
inline Correlation combine_two_correlations(const Correlation& a, const Correlation& b) {
  auto sq = [](double x) { return x * x; };
  Correlation r;
  r.n = a.n + b.n;
  r.mean_x = (a.n * a.mean_x + b.n * b.mean_x) / r.n;
  r.mean_y = (a.n * a.mean_y + b.n * b.mean_y) / r.n;
  r.sum_xx = a.sum_xx + a.n * sq(a.mean_x - r.mean_x)
           + b.sum_xx + b.n * sq(b.mean_x - r.mean_x);
  r.sum_yy = a.sum_yy + a.n * sq(a.mean_y - r.mean_y)
           + b.sum_yy + b.n * sq(b.mean_y - r.mean_y);
  r.sum_xy = a.sum_xy + a.n * (a.mean_x - r.mean_x) * (a.mean_y - r.mean_y)
           + b.sum_xy + b.n * (b.mean_x - r.mean_x) * (b.mean_y - r.mean_y);
  return r;
}

//! @brief Combine multiple correlation accumulators.
//! @param cors Vector of correlations to combine
//! @return Combined correlation
//!
//! Successively combines all correlations in the vector.
inline Correlation combine_correlations(const std::vector<Correlation>& cors) {
  Correlation result;
  for (const Correlation& cor : cors)
    result = combine_two_correlations(result, cor);
  return result;
}

//! @brief Summary statistics for a dataset.
//!
//! Contains min, max, mean, RMS, and count of NaN values.
struct DataStats {
  double dmin = NAN;       //!< Minimum value
  double dmax = NAN;       //!< Maximum value
  double dmean = NAN;      //!< Mean value
  double rms = NAN;        //!< Root mean square
  size_t nan_count = 0;    //!< Number of NaN values
};

//! @brief Calculate summary statistics for a vector of data.
//! @tparam T Data type (must be convertible to double)
//! @param data Vector of data values
//! @return DataStats with summary statistics
//!
//! Computes min, max, mean, and RMS in a single pass. NaN values are excluded
//! from calculations and counted separately.
template<typename T>
DataStats calculate_data_statistics(const std::vector<T>& data) {
  DataStats stats;
  double sum = 0;
  double sq_sum = 0;
  stats.dmin = INFINITY;
  stats.dmax = -INFINITY;
  for (double d : data) {
    if (std::isnan(d)) {
      stats.nan_count++;
      continue;
    }
    sum += d;
    sq_sum += d * d;
    if (d < stats.dmin)
      stats.dmin = d;
    if (d > stats.dmax)
      stats.dmax = d;
  }
  if (stats.nan_count != data.size()) {
    stats.dmean = sum / (data.size() - stats.nan_count);
    stats.rms = std::sqrt(sq_sum / (data.size() - stats.nan_count) - stats.dmean * stats.dmean);
  } else {
    stats.dmin = NAN;
    stats.dmax = NAN;
  }
  return stats;
}

} // namespace gemmi
#endif
