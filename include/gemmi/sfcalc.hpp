// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors.
// It does not use optimizations described in the literature,
// cf. Bourhis et al (2014) https://doi.org/10.1107/S2053273314022207,
// because direct calculations are not used in MX if performance is important.
// For FFT-based calculations see rhogrid.hpp + fourier.hpp.

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "model.hpp"   // for Structure, ...
#include "small.hpp"   // for SmallStructure
#include "fprime.hpp"  // for cromer_libermann

namespace gemmi {

// calculate exp(-2 pi^2 s.U.s)
template<typename Site>
double calculate_aniso_part(const UnitCell& cell, const Site& site,
                            const Vec3& hkl) {
  double arh = cell.ar * hkl.x;
  double brk = cell.br * hkl.y;
  double crl = cell.cr * hkl.z;
  double sus = arh * arh * site.aniso.u11 +
               brk * brk * site.aniso.u22 +
               crl * crl * site.aniso.u33 +
               2 * (arh * brk * site.aniso.u12 +
                    arh * crl * site.aniso.u13 +
                    brk * crl * site.aniso.u23);
  return std::exp(-2 * pi() * pi() * sus);
}

// calculate part of the structure factor: exp(2 pi i r * s)
std::complex<double> calculate_sf_part(const Fractional& fpos,
                                       const Miller& hkl) {
  double arg = 2 * pi() * (hkl[0]*fpos.x + hkl[1]*fpos.y + hkl[2]*fpos.z);
  return std::complex<double>{std::cos(arg), std::sin(arg)};
}

template <typename Table>
class StructureFactorCalculator {
public:
  StructureFactorCalculator(const UnitCell& cell) : cell_(cell) {}

  void set_stol2_and_scattering_factors(const Miller& hkl) {
    stol2_ = 0.25 * cell_.calculate_1_d2(hkl);
    scattering_factors_.clear();
    scattering_factors_.resize((int) El::END, 0.);
    for (auto const& fprime : fprimes_) {
      El el = fprime.first;
      if (Table::has(el)) {
        double sf = Table::get(el).calculate_sf(stol2_) + fprime.second;
        scattering_factors_[(int)el] = sf;
      }
    }
  }

  double get_scattering_factor(Element element) {
    double& sfactor = scattering_factors_[element.ordinal()];
    if (sfactor == 0.) {
      if (!Table::has(element.elem))
        fail("Missing scattering factor for ", element.name());
      sfactor = Table::get(element.elem).calculate_sf(stol2_);
    }
    return sfactor;
  }

  std::complex<double> calculate_sf_from_model(const Model& model,
                                               const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms) {
          Fractional fract = cell_.fractionalize(site.pos);
          double oc_sf = site.occ * get_scattering_factor(site.element);
          std::complex<double> factor = calculate_sf_part(fract, hkl);
          if (!site.aniso.nonzero()) {
            for (const FTransform& image : cell_.images)
              factor += calculate_sf_part(image.apply(fract), hkl);
            sf += oc_sf * std::exp(-site.b_iso * stol2_) * factor;
          } else {
            Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
            factor *= calculate_aniso_part(cell_, site, vhkl);
            for (const FTransform& image : cell_.images) {
              factor += calculate_sf_part(image.apply(fract), hkl) *
                        calculate_aniso_part(cell_, site,
                                             image.mat.left_multiply(vhkl));
            }
            sf += oc_sf * factor;
          }
        }
    return sf;
  }

  // The occupancy is assumed to take into account symmetry,
  // i.e. to be fractional if the atom is on special position.
  std::complex<double>
  calculate_sf_from_small_structure(const SmallStructure& small,
                                    const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const SmallStructure::Site& site : small.sites) {
      double oc_sf = site.occ * get_scattering_factor(site.element);
      std::complex<double> factor = calculate_sf_part(site.fract, hkl);
      if (!site.aniso.nonzero()) {
        for (const FTransform& image : cell_.images)
          factor += calculate_sf_part(image.apply(site.fract), hkl);
        sf += oc_sf * std::exp(-8 * pi() * pi() * stol2_ * site.u_iso) * factor;
      } else {
        Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
        factor *= calculate_aniso_part(cell_, site, vhkl);
        for (const FTransform& image : cell_.images) {
          factor += calculate_sf_part(image.apply(site.fract), hkl) *
                    calculate_aniso_part(cell_, site,
                                         image.mat.left_multiply(vhkl));
        }
        sf += oc_sf * factor;
      }
    }
    return sf;
  }

  void set_fprime(El el, double val) { fprimes_[el] = val; }
  void set_fprime_if_not_set(El el, double val) { fprimes_.emplace(el, val); }
  const std::map<El, double>& fprimes() const { return fprimes_; }

private:
  const UnitCell& cell_;
  double stol2_;
  std::map<El, double> fprimes_;
  std::vector<double> scattering_factors_;
};

} // namespace gemmi
#endif
