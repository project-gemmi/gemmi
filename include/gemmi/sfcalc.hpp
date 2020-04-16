// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors
// (see rhogrid.hpp + fourier.hpp for FFT-based calculation)

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "model.hpp"   // for Structure, ...
#include "small.hpp"   // for SmallStructure
#include "fprime.hpp"  // for cromer_libermann

namespace gemmi {

// calculate s.U.s
template<typename Site>
double calculate_s_u_s(const UnitCell& cell, const Site& site, const Vec3& hkl) {
  double arh = cell.ar * hkl.x;
  double brk = cell.br * hkl.y;
  double crl = cell.cr * hkl.z;
  return arh * arh * site.u11 + brk * brk * site.u22 + crl * crl * site.u33 +
    2 * (arh * brk * site.u12 + arh * crl * site.u13 + brk * crl * site.u23);
}

template <typename Table>
class StructureFactorCalculator {
public:
  StructureFactorCalculator(const UnitCell& cell) : cell_(cell) {}

  void set_hkl(const Miller& hkl) {
    hkl_ = hkl;
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
    set_hkl(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& a : res.atoms) {
          Fractional fpos = cell_.fractionalize(a.pos);
          std::complex<double> part = calculate_sf_part(fpos);
          for (const FTransform& image : cell_.images)
            part += calculate_sf_part(image.apply(fpos));
          sf += (double)a.occ * get_scattering_factor(a.element) *
                std::exp(-a.b_iso * stol2_) * part;
        }
    return sf;
  }

  // The occupancy is assumed to take into account symmetry,
  // i.e. to be fractional if the atom is on special position.
  std::complex<double>
  calculate_sf_from_small_structure(const SmallStructure& small,
                                    const Miller& hkl) {
    constexpr double mtwo_pi2 = -2 * pi() * pi();
    std::complex<double> sf = 0.;
    set_hkl(hkl);
    for (const SmallStructure::Site& site : small.sites) {
      double oc_sf = site.occ * get_scattering_factor(site.element);
      std::complex<double> factor = calculate_sf_part(site.fract);
      if (site.u11 == 0.) {
        for (const FTransform& image : cell_.images)
          factor += calculate_sf_part(image.apply(site.fract));
        sf += oc_sf * std::exp(4 * mtwo_pi2 * stol2_ * site.u_iso) * factor;
      } else {
        Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
        factor *= std::exp(mtwo_pi2 * calculate_s_u_s(cell_, site, vhkl));
        for (const FTransform& image : cell_.images) {
          factor += calculate_sf_part(image.apply(site.fract)) *
                    std::exp(mtwo_pi2 * calculate_s_u_s(cell_, site,
                                               image.mat.left_multiply(vhkl)));
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
  // calculate part of the structure factor: exp(2 pi i r * s)
  std::complex<double> calculate_sf_part(const Fractional& fpos) const {
    double arg = 2 * pi() * (hkl_[0]*fpos.x + hkl_[1]*fpos.y + hkl_[2]*fpos.z);
    return std::complex<double>{std::cos(arg), std::sin(arg)};
  }

  const UnitCell& cell_;
  Miller hkl_;
  double stol2_;
  std::map<El, double> fprimes_;
  std::vector<double> scattering_factors_;
};


} // namespace gemmi
#endif
