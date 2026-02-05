//! @file
//! @brief Calculate electron density from molecular models.
//!
//! Tools to prepare a grid with values of electron density of a model.

// Copyright 2019 Global Phasing Ltd.
//
// Tools to prepare a grid with values of electron density of a model.

#ifndef GEMMI_DENCALC_HPP_
#define GEMMI_DENCALC_HPP_

#include "addends.hpp"  // for Addends
#include "formfact.hpp" // for ExpSum
#include "grid.hpp"     // for Grid
#include "model.hpp"    // for Structure, ...
#include "calculate.hpp" // for calculate_b_aniso_range

namespace gemmi {

//! @brief Determine radius where density falls below cutoff.
//! @tparam N Number of Gaussian terms
//! @tparam Real Floating-point type
//! @param x1 Initial radius estimate
//! @param precal Precalculated exponential sum
//! @param cutoff_level Density threshold
//! @return Radius where density equals cutoff_level
template<int N, typename Real>
Real determine_cutoff_radius(Real x1, const ExpSum<N, Real>& precal, Real cutoff_level) {
  Real y1, dy;
  std::tie(y1, dy) = precal.calculate_with_derivative(x1);
  // Generally, density is supposed to decrease with radius.
  // But if we have addends (in particular -Z for Mott-Bothe),
  // it can first rise, then decrease. We want to be after the maximum.
  while (std::copysign(dy, y1 * dy) > 0) { // unlikely
    x1 += 1.0f;
    std::tie(y1, dy) = precal.calculate_with_derivative(x1);
  }
  Real x2 = x1;
  Real y2 = y1;
  if (std::fabs(y1) < cutoff_level) {
    while (std::fabs(y1) < cutoff_level) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5f;
      std::tie(y1, dy) = precal.calculate_with_derivative(x1);
      // with addends it's possible to land on the left side of the maximum
      if (std::copysign(dy, y1 * dy) > 0) { // unlikely
        while (std::copysign(dy, y1 * dy) > 0 && x1 + 0.1f < x2) {
          x1 += 0.1f;
          std::tie(y1, dy) = precal.calculate_with_derivative(x1);
        }
        if (std::fabs(y1) < cutoff_level)
          return x1;
        break;
      }
      if (x1 < 0) { // unlikely
        x1 = 0;
        y1 = precal.calculate(x1 * x1);
        break;
      }
    }
  } else {
    while (std::fabs(y2) > cutoff_level) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5f;
      y2 = precal.calculate(x2 * x2);
    }
  }

  return x1 + (x1 - x2) / (y1 - y2) * (cutoff_level - y1);
}

//! @brief Approximate density radius for IT92 scattering factors.
//! @tparam Real Floating-point type
//! @param b B-factor value
//! @return Approximate radius above cutoff=1e-5 (calibrated for carbon)
//!
//! approximated radius of electron density (IT92) above cutoff=1e-5 for C
template <typename Real>
Real it92_radius_approx(Real b) {
  return (8.5f + 0.075f * b) / (2.4f + 0.0045f * b);
}

//! @brief Calculate electron density on a grid.
//! @tparam Table Scattering factor table type (e.g., IT92)
//! @tparam GReal Grid data type (float or double)
//!
//! Usual usage:
//! - set d_min and optionally also other parameters,
//! - set addends to f' values for your wavelength (see fprime.hpp)
//! - use grid.setup_from() to set grid's unit cell and space group
//! - check that Table has SF coefficients for all elements that are to be used
//! - call put_model_density_on_grid()
//! - do FFT using transform_map_to_f_phi()
//! - if blur is used, multiply the SF by reciprocal_space_multiplier()
template <typename Table, typename GReal>
struct DensityCalculator {
  // GReal = type of grid; CReal = type of coefficients in Table
  using CReal = typename Table::Coef::coef_type;
  Grid<GReal> grid;  //!< Output density grid
  double d_min = 0.;  //!< Minimum resolution (Angstroms)
  double rate = 1.5;  //!< Oversampling rate
  double blur = 0.;  //!< Additional B-factor blur
  float cutoff = 1e-5f;  //!< Density cutoff threshold
#if GEMMI_COUNT_DC
  size_t atoms_added = 0;
  size_t density_computations = 0;
#endif
  Addends addends;  //!< Anomalous scattering corrections (f')

  //! @brief Get requested grid spacing.
  //! @return Grid spacing (d_min / (2 * rate))
  double requested_grid_spacing() const { return d_min / (2 * rate); }

  //! @brief Set blur to match Refmac FFT map calculations.
  //! @param model Model to analyze
  //! @param allow_negative Allow negative blur values
  void set_refmac_compatible_blur(const Model& model, bool allow_negative=false) {
    double spacing = requested_grid_spacing();
    if (spacing <= 0)
      spacing = std::min(std::min(grid.spacing[0], grid.spacing[1]), grid.spacing[2]);
    double b_min = calculate_b_aniso_range(model).first;
    blur = u_to_b() / 1.1 * sq(spacing) - b_min;
    if (!allow_negative && blur < 0)
      blur = 0.;
  }

  //! @brief Add single atom's density to grid.
  //! @param atom Atom to add
  //!
  //! pre: check if Table::has(atom.element)
  void add_atom_density_to_grid(const Atom& atom) {
    Element el = atom.element;
    const auto& coef = Table::get(el, atom.charge, atom.serial);
    do_add_atom_density_to_grid(atom, coef, addends.get(el));
  }

  //! @brief Add constant density contribution for atom.
  //! @param atom Atom position
  //! @param c Constant factor value
  //!
  //! Parameter c is a constant factor and has the same meaning as either addend
  //! or c in scattering factor coefficients (a1, b1, ..., c).
  void add_c_contribution_to_grid(const Atom& atom, float c) {
    do_add_atom_density_to_grid(atom, GaussianCoef<0, 1, CReal>{0}, c);
  }

  template<int N>
  CReal estimate_radius(const ExpSum<N, CReal>& precal, CReal b) const {
    if (N == 1)
      return std::sqrt(std::log(cutoff / std::abs(precal.a[0])) / precal.b[0]);
    CReal x1 = it92_radius_approx(b);
    return determine_cutoff_radius(x1, precal, (CReal)cutoff);
  }

  template<typename Coef>
  void do_add_atom_density_to_grid(const Atom& atom, const Coef& coef, float addend) {
#if GEMMI_COUNT_DC
    ++atoms_added;
#endif
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    if (!atom.aniso.nonzero()) {
      // isotropic
      CReal b = static_cast<CReal>(atom.b_iso + blur);
      auto precal = coef.precalculate_density_iso(b, addend);
      CReal radius = estimate_radius(precal, b);
      grid.template use_points_around<true>(fpos, radius, [&](GReal& point, double r2) {
          point += GReal(atom.occ * precal.calculate((CReal)r2));
#if GEMMI_COUNT_DC
          ++density_computations;
#endif
      }, /*fail_on_too_large_radius=*/false);
    } else {
      // anisotropic
      auto aniso_b = atom.aniso.scaled(CReal(u_to_b())).added_kI(CReal(blur));
      // rough estimate, so we don't calculate eigenvalues
      CReal b_max = std::max(std::max(aniso_b.u11, aniso_b.u22), aniso_b.u33);
      auto precal_iso = coef.precalculate_density_iso(b_max, addend);
      double radius = estimate_radius(precal_iso, b_max);
      auto precal = coef.precalculate_density_aniso_b(aniso_b, addend);
      int du = (int) std::ceil(radius / grid.spacing[0]);
      int dv = (int) std::ceil(radius / grid.spacing[1]);
      int dw = (int) std::ceil(radius / grid.spacing[2]);
      grid.template use_points_in_box<true>(
          fpos, du, dv, dw,
          [&](GReal& point, double, const Position& delta, int, int, int) {
            point += GReal(atom.occ * precal.calculate(delta));
#if GEMMI_COUNT_DC
            ++density_computations;
#endif
          },
          false, radius);
    }
  }

  //! @brief Initialize grid with appropriate size and zero values.
  //! @throws std::runtime_error if d_min not set and grid not configured
  void initialize_grid() {
    grid.data.clear();
    double spacing = requested_grid_spacing();
    if (spacing > 0)
      grid.set_size_from_spacing(spacing, GridSizeRounding::Up);
    else if (grid.point_count() > 0)
      // d_min not set, but a custom grid has been setup by the user
      grid.fill(0.);
    else
      fail("initialize_grid(): d_min is not set");
  }

  //! @brief Add all atoms from model to existing grid.
  //! @param model Model with atoms to add
  void add_model_density_to_grid(const Model& model) {
    grid.check_not_empty();
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          add_atom_density_to_grid(atom);
  }

  //! @brief Calculate complete density map for model.
  //! @param model Model to calculate density for
  //!
  //! Initializes grid, adds all atoms, and symmetrizes.
  void put_model_density_on_grid(const Model& model) {
    initialize_grid();
    add_model_density_to_grid(model);
    grid.symmetrize_sum();
  }

  // deprecated, use directly grid.setup_from(st)
  void set_grid_cell_and_spacegroup(const Structure& st) {
    grid.setup_from(st);
  }

  //! @brief Get blur correction in reciprocal space.
  //! @param inv_d2 Inverse d-spacing squared (1/d^2)
  //! @return Blur multiplier exp(blur * 0.25 * 1/d^2)
  //!
  //! The argument is 1/d^2 - as outputted by unit_cell.calculate_1_d2(hkl).
  double reciprocal_space_multiplier(double inv_d2) const {
    return std::exp(blur * 0.25 * inv_d2);
  }

  //! @brief Get Mott-Bethe correction factor for reflection.
  //! @param hkl Miller indices
  //! @return Mott-Bethe factor with optional blur correction
  double mott_bethe_factor(const Miller& hkl) const {
    double inv_d2 = grid.unit_cell.calculate_1_d2(hkl);
    double factor = -mott_bethe_const() / inv_d2;
    return blur == 0 ? factor : factor * reciprocal_space_multiplier(inv_d2);
  }
};

} // namespace gemmi
#endif
