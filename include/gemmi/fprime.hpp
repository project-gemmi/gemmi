// Cromer-Liberman calculation of anomalous scattering factors,
// with corrections from Kissel & Pratt.

// Read the header comment in src/fprime.cpp for the details.

#ifndef GEMMI_FPRIME_HPP_
#define GEMMI_FPRIME_HPP_

#include "fail.hpp"

namespace gemmi {

/// Cromer-Liberman calculation of anomalous scattering factors.
/// input:
///   z      - atomic number
///   npts   - array length
///   energy - energies in eV
/// output:
///   fp     - f' (real part of anomalous scattering)
///   fpp    - f" (imaginary part of anomalous scattering)
GEMMI_DLL void cromer_liberman_for_array(int z, int npts, const double* energy,
                                         double* fp, double* fpp);

/// returns fp, fpp is returned through the last argument.
inline double cromer_liberman(int z, double energy, double* fpp) {
  double fp = 0., fpp_ = 0.;
  cromer_liberman_for_array(z, 1, &energy, &fp, &fpp_);
  if (fpp)
    *fpp = fpp_;
  return fp;
}

} // namespace gemmi
#endif
