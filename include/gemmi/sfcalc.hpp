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

  // Calculation of Debye-Waller factor with isotropic ADPs
  double dwf_iso(const SmallStructure::Site& site) const {
    return std::exp(-8 * pi() * pi() * stol2_ * site.u_iso);
  }
  double dwf_iso(const Atom& atom) const {
    return std::exp(-stol2_ * atom.b_iso);
  }

  // Calculation of Debye-Waller factor exp(-2 pi^2 s.U.s)
  // cf. B. Rupp's book, p. 641 or RWGK & Adams 2002, J. Appl. Cryst. 35, 477
  // Small molecule and macromolecular anisotropic U's are defined differently,
  // so we have two functions.
  double dwf_aniso(const SmallStructure::Site& site, const Vec3& hkl) const {
    Vec3 arh(cell_.ar * hkl.x, cell_.br * hkl.y, cell_.cr * hkl.z);
    return std::exp(-2 * pi() * pi() * site.aniso.r_u_r(arh));
  }
  double dwf_aniso(const Atom& atom, const Vec3& hkl) const {
    return std::exp(-2 * pi() * pi() *
                    atom.aniso.transformed_by(cell_.frac.mat).r_u_r(hkl));
  }

  template<typename Site>
  std::complex<double> calculate_sf_from_atom(const Fractional& fract,
                                              const Site& site,
                                              const Miller& hkl) {
    double oc_sf = site.occ * get_scattering_factor(site.element);
    std::complex<double> sum = calculate_sf_part(fract, hkl);
    if (!site.aniso.nonzero()) {
      for (const FTransform& image : cell_.images)
        sum += calculate_sf_part(image.apply(fract), hkl);
      return oc_sf * dwf_iso(site) * sum;
    } else {
      Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
      sum *= dwf_aniso(site, vhkl);
      for (const FTransform& image : cell_.images)
        sum += calculate_sf_part(image.apply(fract), hkl) *
               dwf_aniso(site, image.mat.left_multiply(vhkl));
      return oc_sf * sum;
    }
  }

  std::complex<double> calculate_sf_from_model(const Model& model,
                                               const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms)
          sf += calculate_sf_from_atom(cell_.fractionalize(site.pos),
                                       site, hkl);
    return sf;
  }

  // The occupancy is assumed to take into account symmetry,
  // i.e. to be fractional if the atom is on special position.
  std::complex<double>
  calculate_sf_from_small_structure(const SmallStructure& small,
                                    const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const SmallStructure::Site& site : small.sites)
      sf += calculate_sf_from_atom(site.fract, site, hkl);
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
