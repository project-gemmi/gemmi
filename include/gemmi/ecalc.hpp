// Copyright 2023 Global Phasing Ltd.
//
// Normalization of amplitudes F->E ("Karle" approach, similar to CCP4 ECALC).

#ifndef GEMMI_ECALC_HPP_
#define GEMMI_ECALC_HPP_

#include <cassert>
#include "binner.hpp"

namespace gemmi {

template<typename DataProxy>
std::vector<double> calculate_amplitude_normalizers(const DataProxy& data, int fcol_idx,
                                                    const Binner& binner) {
  struct CountAndSum {
    int n = 0;
    double sum = 0.;
  };
  int nreflections = data.size() / data.stride();
  std::vector<double> multipliers(nreflections, NAN);
  if (data.spacegroup() == nullptr)
    gemmi::fail("unknown space group in the data file");
  GroupOps gops = data.spacegroup()->operations();
  std::vector<double> inv_d2(multipliers.size());
  for (size_t i = 0, n = 0; n < data.size(); n += data.stride(), ++i)
    inv_d2[i] = data.unit_cell().calculate_1_d2(data.get_hkl(n));
  std::vector<int> bin_index = binner.get_bins_from_1_d2(inv_d2);
  std::vector<CountAndSum> stats(binner.size());
  for (size_t i = 0, n = 0; n < data.size(); n += data.stride(), i++) {
    Miller hkl = data.get_hkl(n);
    double f = data.get_num(n + fcol_idx);
    if (!std::isnan(f)) {
      int epsilon = gops.epsilon_factor(hkl);
      double inv_epsilon = 1.0 / epsilon;
      double f2 = f * f * inv_epsilon;
      multipliers[i] = std::sqrt(inv_epsilon);
      CountAndSum& cs = stats[bin_index[i]];
      cs.n++;
      cs.sum += f2;
    }
  }

  // simple smoothing with kernel [0.75 1 0.75]
  std::vector<double> smoothed(stats.size());
  {
    const double k = 0.75;
    smoothed[0] = (stats[0].sum + k * stats[1].sum) / (stats[0].n + k * stats[1].n);
    size_t n = stats.size() - 1;
    for (size_t i = 1; i < n; ++i)
      smoothed[i] = (stats[i].sum + k * (stats[i-1].sum + stats[i+1].sum))
                  / (stats[i].n + k * (stats[i-1].n + stats[i+1].n));
    smoothed[n] = (stats[n].sum + k * stats[n-1].sum) / (stats[n].n + k * stats[n-1].n);
  }

#if 0
  {
    // print shell statistics
    std::vector<int> refl_counts(binner.size());
    printf(" shell\t    #F\t    d\t <F^2>\tsmoothd\t  #refl\t mid d\n");
    for (int idx : bin_index)
      ++refl_counts[idx];
    for (size_t i = 0; i < binner.size(); ++i) {
      double d = 1 / std::sqrt(binner.limits[i]);
      double mid_d = 1 / std::sqrt(binner.mids[i]);
      double avg_f2 = stats[i].sum / stats[i].n;
      printf("%6zu\t%6d\t%7.3f\t%7.0f\t%7.0f\t%6d\t%7.3f\n",
             i+1, stats[i].n, d, avg_f2, smoothed[i], refl_counts[i], mid_d);
    }
    printf("\n");
  }
#endif

  for (double& x : smoothed)
    x = std::sqrt(x);
  for (size_t i = 0; i < multipliers.size(); ++i) {
    double x = inv_d2[i];
    int bin = bin_index[i];
    double rms = smoothed[bin];
    if (x > binner.mids.front() && x < binner.mids.back()) {
      // linear interpolation in 1/d^2
      if (x > binner.mids[bin])
        ++bin;
      double x0 = binner.mids[bin - 1];
      double x1 = binner.mids[bin];
      double y0 = smoothed[bin - 1];
      double y1 = smoothed[bin];
      assert(x0 <= x && x <= x1);
      rms = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }
    multipliers[i] /= rms;
  }
  return multipliers;
}

} // namespace gemmi
#endif
