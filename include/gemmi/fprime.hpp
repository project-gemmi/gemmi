// Cromer-Liberman calculation of anomalous scattering factors,
// with corrections from Kissel & Pratt.

// Read the header comment in src/fprime.cpp for the details.

#ifndef GEMMI_FPRIME_HPP_
#define GEMMI_FPRIME_HPP_

#include "fail.hpp"

namespace gemmi {

/// @brief Cromer-Liberman calculation of anomalous scattering factors for an array
///   of energies, with Kissel-Pratt corrections.
/// @param z atomic number
/// @param npts array length
/// @param energy energies in eV
/// @param fp output: f' (real part of anomalous scattering)
/// @param fpp output: f" (imaginary part of anomalous scattering)
/// @par References
/// Cromer, D.T. & Liberman, D.A. (1994). Anomalous dispersion calculations
/// near to and on the long-wavelength side of an absorption edge.
/// Acta Cryst. A51, 416. https://doi.org/10.1107/S0108767394013292
///
/// Kissel, L. & Pratt, R.H. (1990). Rayleigh scattering — elastic photon
/// scattering by bound electrons: status and perspectives.
/// Acta Cryst. A46, 170–175. https://doi.org/10.1107/S0108767389010718
///
/// Brennan, S. & Cowan, P.L. (1992). A suite of programs for calculating
/// x-ray absorption, reflection, and diffraction performance.
/// Rev. Sci. Instrum. 63, 850. https://doi.org/10.1063/1.1142625
GEMMI_DLL void cromer_liberman_for_array(int z, int npts, const double* energy,
                                         double* fp, double* fpp);

/// @brief Single-energy wrapper for cromer_liberman_for_array.
/// @param z atomic number
/// @param energy energy in eV
/// @param fpp output: may be nullptr (optional f" value)
/// @return f' value
inline double cromer_liberman(int z, double energy, double* fpp) {
  double fp = 0., fpp_ = 0.;
  cromer_liberman_for_array(z, 1, &energy, &fp, &fpp_);
  if (fpp)
    *fpp = fpp_;
  return fp;
}

} // namespace gemmi
#endif
