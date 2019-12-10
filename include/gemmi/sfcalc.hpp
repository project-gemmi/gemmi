// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors
// (see rhogrid.hpp + fourier.hpp for FFT-based calculation)

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "model.hpp"   // for Structure, ...
#include "smodel.hpp"  // for AtomicStructure

namespace gemmi {

template <typename Table>
class StructureFactorCalculator {
public:
  StructureFactorCalculator(const UnitCell& cell) : cell_(cell) {}

  void set_hkl(const Miller& hkl) {
    hkl_ = hkl;
    stol2_ = 0.25 * cell_.calculate_1_d2(hkl);
    scattering_factors.clear();
    scattering_factors.resize((int) El::END, 0.);
  }

  double get_scattering_factor(Element element, double fprim = 0.) {
    double& sfactor = scattering_factors[(int)element.elem];
    if (sfactor == 0.) {
      if (!Table::has(element.elem))
        fail("Missing scattering factor for ", element.name());
      sfactor = Table::get(element.elem).calculate_sf(stol2_) + fprim;
    }
    return sfactor;
  }

  std::complex<double> get_contribution(Element el, const Fractional& fpos,
                                        double b_iso, double fprim=0.) {
    double scat_factor = get_scattering_factor(el, fprim);
    std::complex<double> part = calculate_sf_part(fpos);
    for (const FTransform& image : cell_.images)
      part += calculate_sf_part(image.apply(fpos));
    return scat_factor * std::exp(-b_iso * stol2_) * part;
  }

  std::complex<double> calculate_sf_from_model(const Model& model,
                                               const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_hkl(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& a : res.atoms)
          sf += (double)a.occ * get_contribution(a.element,
                                                 cell_.fractionalize(a.pos),
                                                 a.b_iso);
    return sf;
  }

  // The occupancy is assumed to take into account symmetry,
  // i.e. to be fractional if the atom is on special position.
  std::complex<double> calculate_sf_from_atomic_structure(
                            const AtomicStructure& atomic, const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_hkl(hkl);
    for (const AtomicStructure::Site& site : atomic.sites) {
      double b_iso = 8 * pi() * pi() * site.u_iso;
      sf += site.occ * get_contribution(site.element, site.fract, b_iso);
    }
    return sf;
  }

private:
  // calculate part of the structure factor: exp(2 pi i r * s)
  std::complex<double> calculate_sf_part(const Fractional& fpos) {
    double arg = 2 * pi() * (hkl_[0]*fpos.x + hkl_[1]*fpos.y + hkl_[2]*fpos.z);
    return std::complex<double>{std::cos(arg), std::sin(arg)};
  }

  const UnitCell& cell_;
  Miller hkl_;
  double stol2_;
  std::vector<double> scattering_factors;
};


} // namespace gemmi
#endif
