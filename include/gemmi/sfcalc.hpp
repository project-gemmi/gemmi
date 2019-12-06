// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors
// (see rhogrid.hpp + fourier.hpp for FFT-based calculation)

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "model.hpp"   // for Structure, ...

namespace gemmi {

// calculate part of the structure factor: exp(2 pi i r * s)
inline std::complex<double> structure_factor_part(Fractional fpos,
                                                  const Miller& hkl) {
  double arg = 2 * pi() * (hkl[0] * fpos.x + hkl[1] * fpos.y + hkl[2] * fpos.z);
  return std::complex<double>{std::cos(arg), std::sin(arg)};
}

template <typename Table>
std::complex<double> calculate_structure_factor(const Model& model,
                                                const UnitCell& cell,
                                                const Miller& hkl) {
  std::complex<double> sf = 0.;
  std::vector<double> scattering_factors((int) El::END, 0.);
  double stol2 = 0.25 * cell.calculate_1_d2(hkl[0], hkl[1], hkl[2]);
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        El el = atom.element.elem;
        double& sfactor = scattering_factors[(int)el];
        if (sfactor == 0.) {
          if (!Table::has(el))
            fail("Missing scattering factor for ", element_name(el));
          sfactor = Table::get(el).calculate_sf(stol2);
        }
        Fractional fpos = cell.fractionalize(atom.pos);
        std::complex<double> part = structure_factor_part(fpos, hkl);
        for (const FTransform& image : cell.images)
          part += structure_factor_part(image.apply(fpos), hkl);
        sf += atom.occ * sfactor * std::exp(-atom.b_iso * stol2) * part;
      }
  return sf;
}

} // namespace gemmi
#endif
