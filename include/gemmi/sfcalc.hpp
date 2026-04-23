// Copyright 2019 Global Phasing Ltd.
//
// Direct calculation of structure factors.
//
// It does not use optimizations described in the literature,
// cf. Bourhis et al (2014) https://doi.org/10.1107/S2053273314022207,
// because direct calculations are not used in MX if performance is important.
// For FFT-based calculations see dencalc.hpp + fourier.hpp.

#ifndef GEMMI_SFCALC_HPP_
#define GEMMI_SFCALC_HPP_

#include <complex>
#include "addends.hpp" // for Addends
#include "model.hpp"   // for Structure, ...
#include "small.hpp"   // for SmallStructure

namespace gemmi {

/// @brief Compute exp(2πi r·s) for one atom at fractional position for a reflection.
/// @param fpos Fractional coordinates of the atom.
/// @param hkl Miller indices of the reflection.
/// @return Complex phase factor.
inline std::complex<double> calculate_sf_part(const Fractional& fpos,
                                              const Miller& hkl) {
  double arg = 2 * pi() * (hkl[0]*fpos.x + hkl[1]*fpos.y + hkl[2]*fpos.z);
  return std::complex<double>{std::cos(arg), std::sin(arg)};
}

/// @brief Calculates structure factors by direct summation over all atoms.
/// @details Simple direct summation; for optimised FFT-based calculations see
/// dencalc.hpp + fourier.hpp.
/// @tparam Table Scattering factor table type (e.g., IT92, WK95, ElectronTable).
///         Must have a nested Coef type with coef_type member.
/// @par References
/// Bourhis, L.J., Dolomanov, O.V., Gildea, R.J., Howard, J.A.K. &
/// Puschmann, H. (2015). The anatomy of a comprehensive constrained,
/// restrained refinement program for the modern computing environment —
/// Olex2 dissected. Acta Cryst. A70, 300–311.
/// https://doi.org/10.1107/S2053273314022207
template <typename Table>
class StructureFactorCalculator {
public:
  using coef_type = typename Table::Coef::coef_type;

  /// @brief Initialize with the unit cell.
  /// @param cell Unit cell; a reference is stored for use in subsequent calculations.
  StructureFactorCalculator(const UnitCell& cell) : cell_(cell) {}

  /// @brief Cache (sin θ/λ)² for hkl and clear the scattering factor cache.
  /// Must be called before computing structure factors for a new reflection.
  /// @param hkl Miller indices of the reflection.
  void set_stol2_and_scattering_factors(const Miller& hkl) {
    stol2_ = (coef_type) cell_.calculate_stol_sq(hkl);
    scattering_factors_.clear();
    scattering_factors_.resize(addends.size(), 0.);
  }

  /// @brief Return scattering factor for element at current (sin θ/λ)².
  /// Lazily computes and caches the result.
  /// @param element The chemical element.
  /// @param charge Formal charge (signed char; 0 = neutral atom).
  /// @return Scattering factor value.
  double get_scattering_factor(Element element, signed char charge) {
    double& sfactor = scattering_factors_[element.ordinal()];
    if (sfactor == 0.) {
      if (!Table::has(element.elem))
        fail("Missing scattering factor for ", element.name());
      sfactor = Table::get(element.elem, charge).calculate_sf(stol2_) + addends.get(element);
    }
    return sfactor;
  }

  /// @brief Isotropic Debye-Waller factor exp(-B·stol²) for a small-molecule site.
  /// @param site Small-molecule site with u_iso field.
  /// @return Debye-Waller factor value.
  double dwf_iso(const SmallStructure::Site& site) const {
    return std::exp(-u_to_b() * stol2_ * site.u_iso);
  }

  /// @brief Isotropic Debye-Waller factor exp(-B·stol²) for a macromolecular atom.
  /// @param atom Macromolecular atom with b_iso field.
  /// @return Debye-Waller factor value.
  double dwf_iso(const Atom& atom) const {
    return std::exp(-stol2_ * atom.b_iso);
  }

  /// @brief Anisotropic Debye-Waller factor exp(-2π²·s·U·s) for a small-molecule site.
  /// Uses site.aniso, where the anisotropic U tensor is in crystallographic coordinates.
  /// @param site Small-molecule site with aniso field.
  /// @param hkl Miller indices of the reflection.
  /// @return Debye-Waller factor value.
  double dwf_aniso(const SmallStructure::Site& site, const Vec3& hkl) const {
    Vec3 arh(cell_.ar * hkl.x, cell_.br * hkl.y, cell_.cr * hkl.z);
    return std::exp(-2 * pi() * pi() * site.aniso.r_u_r(arh));
  }

  /// @brief Anisotropic Debye-Waller factor exp(-2π²·s·U·s) for a macromolecular atom.
  /// Uses atom.aniso in orthogonal coordinates.
  /// @param atom Macromolecular atom with aniso field.
  /// @param hkl Miller indices of the reflection.
  /// @return Debye-Waller factor value.
  double dwf_aniso(const Atom& atom, const Vec3& hkl) const {
    return std::exp(-2 * pi() * pi() *
                    atom.aniso.transformed_by<>(cell_.frac.mat).r_u_r(hkl));
  }

  /// @brief Contribution of one atom to structure factor, given its scattering factor.
  /// Accounts for Debye-Waller factor, occupancy, and phase.
  /// @tparam Site SmallStructure::Site or Atom.
  /// @param fract Fractional coordinates of the atom.
  /// @param site Site/atom record containing occupancy and aniso data.
  /// @param hkl Miller indices of the reflection.
  /// @param sf Precomputed scattering factor for this atom.
  /// @return Complex structure factor contribution.
  template<typename Site>
  std::complex<double> calculate_sf_from_atom_sf(const Fractional& fract,
                                                 const Site& site,
                                                 const Miller& hkl,
                                                 double sf) {
    double oc_sf = site.occ * sf;
    std::complex<double> sum = calculate_sf_part(fract, hkl);
    if (!site.aniso.nonzero()) {
      for (const FTransform& image : cell_.images)
        sum += calculate_sf_part(image.apply(fract), hkl);
      return oc_sf * dwf_iso(site) * sum;
    }
    Vec3 vhkl(hkl[0], hkl[1], hkl[2]);
    sum *= dwf_aniso(site, vhkl);
    for (const FTransform& image : cell_.images)
      sum += calculate_sf_part(image.apply(fract), hkl) *
             dwf_aniso(site, image.mat.left_multiply(vhkl));
    return oc_sf * sum;
  }

  /// @brief Contribution of one atom to structure factor, with automatic scattering factor lookup.
  /// Accounts for Debye-Waller factor, occupancy, and phase.
  /// @tparam Site SmallStructure::Site or Atom.
  /// @param fract Fractional coordinates of the atom.
  /// @param site Site/atom record containing element and charge information.
  /// @param hkl Miller indices of the reflection.
  /// @return Complex structure factor contribution.
  template<typename Site>
  std::complex<double> calculate_sf_from_atom(const Fractional& fract,
                                              const Site& site,
                                              const Miller& hkl) {
    double atom_sf = get_scattering_factor(site.element, site.charge);
    return calculate_sf_from_atom_sf(fract, site, hkl, atom_sf);
  }

  /// @brief Sum contributions from all atoms in a Model for the given reflection.
  /// @param model The macromolecular model.
  /// @param hkl Miller indices of the reflection.
  /// @return Total structure factor.
  std::complex<double> calculate_sf_from_model(const Model& model, const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms)
          sf += calculate_sf_from_atom(cell_.fractionalize(site.pos), site, hkl);
    return sf;
  }

  /// @brief Compute Z (atomic number sum) component for Mott-Bethe conversion.
  /// Used when a different model is needed for the Z calculation.
  /// @param model The macromolecular model.
  /// @param hkl Miller indices of the reflection.
  /// @param only_h If true, restrict calculation to hydrogen atoms.
  /// @return Z component of structure factor.
  std::complex<double> calculate_mb_z(const Model& model, const Miller& hkl, bool only_h) {
    std::complex<double> sf = 0.;
    stol2_ = (coef_type) cell_.calculate_stol_sq(hkl);
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& site : res.atoms)
          if (!only_h || site.element.is_hydrogen())
            sf += calculate_sf_from_atom_sf(cell_.fractionalize(site.pos), site, hkl,
                                            -site.element.atomic_number());
    return sf;
  }

  /// @brief Return the Mott-Bethe prefactor at current (sin θ/λ)².
  /// Used to convert X-ray structure factors to electron structure factors.
  /// @return Mott-Bethe factor: -mott_bethe_const() / (4·stol²).
  double mott_bethe_factor() const {
    return -mott_bethe_const() / 4 / stol2_;
  }

  /// @brief Sum contributions from all sites in a SmallStructure for the given reflection.
  /// The occupancy is assumed to account for symmetry (i.e., fractional if the atom
  /// is on a special position).
  /// @param small_st The small-molecule structure.
  /// @param hkl Miller indices of the reflection.
  /// @return Total structure factor.
  std::complex<double> calculate_sf_from_small_structure(const SmallStructure& small_st,
                                                         const Miller& hkl) {
    std::complex<double> sf = 0.;
    set_stol2_and_scattering_factors(hkl);
    for (const SmallStructure::Site& site : small_st.sites)
      sf += calculate_sf_from_atom(site.fract, site, hkl);
    return sf;
  }

private:
  const UnitCell& cell_;
  coef_type stol2_;
  std::vector<double> scattering_factors_;
public:
  /// Addends object for anomalous scattering corrections (usually f' for X-rays).
  /// Public and mutable by the caller.
  Addends addends;
};

} // namespace gemmi
#endif
