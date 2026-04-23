/// @file dencalc.hpp
/// @brief Electron density calculation from atomic coordinates using Gaussian density distributions.
///
/// Provides DensityCalculator struct for placing atomic Gaussian electron density onto a 3D grid,
/// with support for isotropic and anisotropic B-factors, occupancy, X-ray scattering factors,
/// and electron scattering (Coulomb potential). Used to generate synthetic electron density maps
/// from an atomic model for comparison with experimental maps or map correlation calculations.

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

/// @brief Find the radius at which a radial density function falls below a cutoff level.
/// @tparam N Number of Gaussian components in the exponential sum
/// @tparam Real Floating-point type (float or double)
/// @param x1 Initial search radius (in Angstroms or grid units)
/// @param precal Precalculated exponential sum object providing calculate() and derivative
/// @param cutoff_level Density threshold; radius is where |density| equals this value
/// @return Radius (distance) at which the density function equals the cutoff level
///
/// Uses binary search with special handling for addends (like Mott-Bethe factor) that
/// may cause the function to rise then fall, rather than monotonically decreasing.
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

/// @brief Approximate radial extent of electron density for a given B-factor.
/// @tparam Real Floating-point type
/// @param b Isotropic B-factor (in Angstrom^2)
/// @return Approximate radius (in Angstroms) at which density falls to ~1e-5 (International Tables Vol. C formula)
///
/// Used as initial guess for cutoff radius search; applicable to X-ray scattering factors (IT92).
template <typename Real>
Real it92_radius_approx(Real b) {
  return (8.5f + 0.075f * b) / (2.4f + 0.0045f * b);
}

/// @brief Calculate electron density from an atomic model on a regular grid.
/// @tparam Table Scattering factor coefficients table (e.g., IT92, Cromer-Mann, electron scattering)
/// @tparam GReal Type for grid values (float or double)
///
/// Places Gaussian electron density contributions from atoms onto a 3D grid.
/// Supports isotropic and anisotropic B-factors, occupancy, X-ray and electron scattering.
///
/// Typical workflow:
/// - Construct DensityCalculator; set grid via grid.setup_from(structure)
/// - Set d_min (resolution) and other parameters (blur, cutoff)
/// - Optionally set addends (wavelength-dependent f' anomalous corrections)
/// - Call put_model_density_on_grid(model) to populate the grid
/// - Use transform_map_to_f_phi() for FFT to reciprocal space
/// - If blur > 0, multiply reciprocal-space data by reciprocal_space_multiplier()
template <typename Table, typename GReal>
struct DensityCalculator {
  /// @brief Type of grid (provided as template parameter)
  using CReal = typename Table::Coef::coef_type;

  /// @brief Output electron density grid (unit cell and space group set via setup_from())
  Grid<GReal> grid;

  /// @brief Target d_min resolution (Angstroms); grid sampling = d_min / (2 * rate)
  double d_min = 0.;

  /// @brief Oversampling rate relative to d_min (default 1.5 = 50% oversampling)
  double rate = 1.5;

  /// @brief Additional blur (B-factor) to apply to all atoms (Angstrom^2); default 0
  double blur = 0.;

  /// @brief Density cutoff for determining atom-to-grid interaction radius (default 1e-5)
  float cutoff = 1e-5f;

#if GEMMI_COUNT_DC
  /// @brief Count of atoms added (debugging; only if GEMMI_COUNT_DC defined)
  size_t atoms_added = 0;
  /// @brief Count of density calculations (debugging; only if GEMMI_COUNT_DC defined)
  size_t density_computations = 0;
#endif

  /// @brief Wavelength-dependent f' corrections (additive anomalous factors)
  Addends addends;

  /// @brief Compute grid spacing (Angstroms per voxel) based on d_min and rate.
  /// @return Grid spacing; 0 if d_min not set
  double requested_grid_spacing() const { return d_min / (2 * rate); }

  /// @brief Set blur to match Refmac5 conventions: depends on existing grid and model B-factors.
  /// @param model Atomic model providing B-factor statistics
  /// @param allow_negative If true, blur can be negative (default false = clamp to 0)
  ///
  /// Calculates blur such that the effective B-factor on the model matches Refmac defaults.
  void set_refmac_compatible_blur(const Model& model, bool allow_negative=false) {
    double spacing = requested_grid_spacing();
    if (spacing <= 0)
      spacing = std::min(std::min(grid.spacing[0], grid.spacing[1]), grid.spacing[2]);
    double b_min = calculate_b_aniso_range(model).first;
    blur = u_to_b() / 1.1 * sq(spacing) - b_min;
    if (!allow_negative && blur < 0)
      blur = 0.;
  }

  /// @brief Add electron density contribution of a single atom to the grid.
  /// @param atom Atom with element, position, B-factor, occupancy, and anisotropic info
  /// @pre Table must have scattering factor coefficients for atom.element
  ///
  /// Places Gaussian density with appropriate radius based on B-factor and scattering factors.
  /// Handles both isotropic and anisotropic B-factors. Occupancy and anomalous addends applied.
  void add_atom_density_to_grid(const Atom& atom) {
    Element el = atom.element;
    const auto& coef = Table::get(el, atom.charge, atom.serial);
    do_add_atom_density_to_grid(atom, coef, addends.get(el));
  }

  /// @brief Add a constant radial density contribution for an atom (for special cases).
  /// @param atom Atom providing position, occupancy
  /// @param c Constant density factor (as in scattering factor coefficients or addends)
  ///
  /// Useful for adding constant density contributions (e.g., Mott-Bethe factor for electron scattering).
  void add_c_contribution_to_grid(const Atom& atom, float c) {
    do_add_atom_density_to_grid(atom, GaussianCoef<0, 1, CReal>{0}, c);
  }

  /// @brief Estimate the interaction radius for density based on precalculated B-factor.
  /// @tparam N Number of exponential components
  /// @param precal Precalculated exponential sum object
  /// @param b Isotropic B-factor (Angstrom^2)
  /// @return Interaction radius (Angstroms) where density falls below cutoff
  ///
  /// For single-Gaussian scattering factors (N=1), computes analytically.
  /// For multi-Gaussian (N>1), uses determine_cutoff_radius() binary search.
  template<int N>
  CReal estimate_radius(const ExpSum<N, CReal>& precal, CReal b) const {
    if (N == 1)
      return std::sqrt(std::log(cutoff / std::abs(precal.a[0])) / precal.b[0]);
    CReal x1 = it92_radius_approx(b);
    return determine_cutoff_radius(x1, precal, (CReal)cutoff);
  }

  /// @brief Internal: place electron density on grid for an atom with given scattering factors.
  /// @tparam Coef Scattering factor coefficients type
  /// @param atom Atom with position, occupancy, B-factor, anisotropic U-tensor
  /// @param coef Precalculated scattering factor coefficients (from Table::get())
  /// @param addend Wavelength-dependent anomalous correction (f') added to constant term
  ///
  /// Handles isotropic B-factor case with radial density sampling and anisotropic case
  /// with box-based sampling respecting the anisotropic U-tensor.
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

  /// @brief Clear grid data and allocate size based on d_min and rate, or use existing grid.
  /// @throws Fails if d_min not set and grid has no existing dimensions
  ///
  /// Sets grid size for FFT-friendly dimensions if d_min > 0.
  /// Otherwise assumes grid already configured by user (via setup_from).
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

  /// @brief Add electron density contributions from all atoms in a model.
  /// @param model Atomic model with chains, residues, atoms
  ///
  /// Iterates through all atoms and calls add_atom_density_to_grid() for each.
  /// Grid must already be initialized.
  void add_model_density_to_grid(const Model& model) {
    grid.check_not_empty();
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          add_atom_density_to_grid(atom);
  }

  /// @brief Initialize grid and add all atom densities from a model.
  /// @param model Atomic model
  ///
  /// Calls initialize_grid(), add_model_density_to_grid(), and symmetrize_sum().
  void put_model_density_on_grid(const Model& model) {
    initialize_grid();
    add_model_density_to_grid(model);
    grid.symmetrize_sum();
  }

  /// @brief Set grid unit cell and space group from a structure.
  /// @param st Structure providing unit cell and space group
  /// @deprecated Use grid.setup_from(st) directly
  void set_grid_cell_and_spacegroup(const Structure& st) {
    grid.setup_from(st);
  }

  /// @brief Compute reciprocal-space multiplier for blur correction factor.
  /// @param inv_d2 Inverse d-spacing squared (1/d^2) from unit_cell.calculate_1_d2(hkl)
  /// @return Exponential factor: exp(blur * 0.25 * inv_d2)
  ///
  /// If blur was applied, multiply structure factors by this to correct for the applied blur.
  double reciprocal_space_multiplier(double inv_d2) const {
    return std::exp(blur * 0.25 * inv_d2);
  }

  /// @brief Compute the Mott-Bethe factor for electron scattering.
  /// @param hkl Miller indices
  /// @return Mott-Bethe correction factor (optionally scaled by reciprocal_space_multiplier if blur > 0)
  ///
  /// For electron scattering, the Mott-Bethe factor approximates the Coulomb potential contribution.
  double mott_bethe_factor(const Miller& hkl) const {
    double inv_d2 = grid.unit_cell.calculate_1_d2(hkl);
    double factor = -mott_bethe_const() / inv_d2;
    return blur == 0 ? factor : factor * reciprocal_space_multiplier(inv_d2);
  }
};

} // namespace gemmi
#endif
