// Copyright 2018 Global Phasing Ltd.
//
// Statistics utilities: classes Covariance, Correlation, DataStats

#ifndef GEMMI_STATS_HPP_
#define GEMMI_STATS_HPP_

#include <cstddef>  // for size_t
#include <cmath>    // for sqrt, NAN, INFINITY
#include <vector>

namespace gemmi {

// popular single-pass algorithm for calculating variance and mean
struct Variance {
  int n = 0;
  double sum_sq = 0.;
  double mean_x = 0.;

  Variance() = default;
  template <typename T> Variance(T begin, T end) : Variance() {
    for (auto i = begin; i != end; ++i)
      add_point(*i);
  }
  void add_point(double x) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    sum_sq += dx * (x - mean_x);
  }
  double for_sample() const { return sum_sq / (n - 1); }
  double for_population() const { return sum_sq / n; }
};

struct Covariance : Variance {
  double mean_y = 0.;
  void add_point(double x, double y) {
    ++n;
    double dx = x - mean_x;
    mean_x += dx / n;
    mean_y += (y - mean_y) / n;
    sum_sq += dx * (y - mean_y);
  }
};

struct Correlation {
  int n = 0;
  double sum_xx = 0.;
  double sum_yy = 0.;
  double sum_xy = 0.;
  double mean_x = 0.;
  double mean_y = 0.;
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
  double coefficient() const { return sum_xy / std::sqrt(sum_xx * sum_yy); }
  double x_variance() const { return sum_xx / n; }
  double y_variance() const { return sum_yy / n; }
  double covariance() const { return sum_xy / n; }
  double mean_ratio() const { return mean_y / mean_x; }
  // the regression line
  double slope() const { return sum_xy / sum_xx; }
  double intercept() const { return mean_y - slope() * mean_x; }
};

struct DataStats {
  double dmin = NAN;
  double dmax = NAN;
  double dmean = NAN;
  double rms = NAN;
  size_t nan_count = 0;
};

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
