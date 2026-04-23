// Copyright 2021 Global Phasing Ltd.
//
// Flat bulk solvent mask. With helper tools that modify data on grid.

/// @file
/// @brief Solvent masking utilities for crystallographic refinement.

#ifndef GEMMI_SOLMASK_HPP_
#define GEMMI_SOLMASK_HPP_

#include "grid.hpp"      // for Grid
#include "floodfill.hpp" // for FloodFill
#include "model.hpp"     // for Model, Atom, ...

namespace gemmi {

/// Enumeration of atomic radii sets for solvent masking.
/// Determines which radii library is used when computing the protein region.
enum class AtomicRadiiSet {
  VanDerWaals,  ///< Standard van der Waals radii from crystallographic tables
  Cctbx,        ///< CCTBX/CCP4 van der Waals radii
  Refmac,       ///< Refmac radii for bulk solvent correction
  Constant      ///< Constant radius applied to all atoms
};

/// Returns the van der Waals radius from CCTBX/CCP4 library for a given element.
/// Data derived from cctbx/eltbx/van_der_waals_radii.py for compatibility.
/// @param el Chemical element
/// @return Radius in Angstroms
inline float cctbx_vdw_radius(El el) {
  static constexpr float radii[] = {
    /*X*/  1.00f,
    /*H*/  1.20f, /*He*/ 1.40f,
    /*Li*/ 1.82f, /*Be*/ 0.63f, /*B*/  1.75f, /*C*/ 1.775f, /*N*/ 1.50f,
    /*O*/  1.45f, /*F*/  1.47f, /*Ne*/ 1.54f,
    /*Na*/ 2.27f, /*Mg*/ 1.73f, /*Al*/ 1.50f, /*Si*/ 2.10f, /*P*/ 1.90f,
    /*S*/  1.80f, /*Cl*/ 1.75f, /*Ar*/ 1.88f,
    /*K*/  2.75f, /*Ca*/ 1.95f, /*Sc*/ 1.32f, /*Ti*/ 1.95f, /*V*/  1.06f,
    /*Cr*/ 1.13f, /*Mn*/ 1.19f, /*Fe*/ 1.26f, /*Co*/ 1.13f, /*Ni*/ 1.63f,
    /*Cu*/ 1.40f, /*Zn*/ 1.39f, /*Ga*/ 1.87f, /*Ge*/ 1.48f, /*As*/ 0.83f,
    /*Se*/ 1.90f, /*Br*/ 1.85f, /*Kr*/ 2.02f,
    /*Rb*/ 2.65f, /*Sr*/ 2.02f, /*Y*/  1.61f, /*Zr*/ 1.42f, /*Nb*/ 1.33f,
    /*Mo*/ 1.75f, /*Tc*/ 2.00f, /*Ru*/ 1.20f, /*Rh*/ 1.22f, /*Pd*/ 1.63f,
    /*Ag*/ 1.72f, /*Cd*/ 1.58f, /*In*/ 1.93f, /*Sn*/ 2.17f, /*Sb*/ 1.12f,
    /*Te*/ 1.26f, /*I*/  1.98f, /*Xe*/ 2.16f,
    /*Cs*/ 3.01f, /*Ba*/ 2.41f, /*La*/ 1.83f, /*Ce*/ 1.86f, /*Pr*/ 1.62f,
    /*Nd*/ 1.79f, /*Pm*/ 1.76f, /*Sm*/ 1.74f, /*Eu*/ 1.96f, /*Gd*/ 1.69f,
    /*Tb*/ 1.66f, /*Dy*/ 1.63f, /*Ho*/ 1.61f, /*Er*/ 1.59f, /*Tm*/ 1.57f,
    /*Yb*/ 1.54f, /*Lu*/ 1.53f, /*Hf*/ 1.40f, /*Ta*/ 1.22f, /*W*/  1.26f,
    /*Re*/ 1.30f, /*Os*/ 1.58f, /*Ir*/ 1.22f, /*Pt*/ 1.72f, /*Au*/ 1.66f,
    /*Hg*/ 1.55f, /*Tl*/ 1.96f, /*Pb*/ 2.02f, /*Bi*/ 1.73f, /*Po*/ 1.21f,
    /*At*/ 1.12f, /*Rn*/ 2.30f,
    /*Fr*/ 3.24f, /*Ra*/ 2.57f, /*Ac*/ 2.12f, /*Th*/ 1.84f, /*Pa*/ 1.60f,
    /*U*/  1.75f, /*Np*/ 1.71f, /*Pu*/ 1.67f, /*Am*/ 1.66f, /*Cm*/ 1.65f,
    /*Bk*/ 1.64f, /*Cf*/ 1.63f, /*Es*/ 1.62f, /*Fm*/ 1.61f, /*Md*/ 1.60f,
    /*No*/ 1.59f, /*Lr*/ 1.58f, /*Rf*/ 1.00f, /*Db*/ 1.00f, /*Sg*/ 1.00f,
    /*Bh*/ 1.00f, /*Hs*/ 1.00f, /*Mt*/ 1.00f, /*Ds*/ 1.00f, /*Rg*/ 1.00f,
    /*Cn*/ 1.00f, /*Nh*/ 1.00f, /*Fl*/ 1.00f, /*Mc*/ 1.00f, /*Lv*/ 1.00f,
    /*Ts*/ 1.00f, /*Og*/ 1.00f,
    /*D*/  1.20f, /*END*/0.f
  };
  static_assert(ce_almost_eq(radii[static_cast<int>(El::D)], 1.2f), "Hmm");
  static_assert(sizeof(radii) / sizeof(radii[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return radii[static_cast<int>(el)];
}

/// Returns the effective radius for bulk solvent correction from Refmac's ener_lib.cif.
/// Represents ionic radius minus 0.2 Angstroms or vdW radius plus 0.2 Angstroms.
/// For full Refmac compatibility, use @c r_probe=1.0 and @c r_shrink=0.8.
/// @param el Chemical element
/// @return Radius in Angstroms
inline float refmac_radius_for_bulk_solvent(El el) {
#if 0
  static constexpr float radii[] = {
    /*X*/  1.00f,
    /*H*/  1.40f, /*He*/ 1.60f,
    /*Li*/ 0.53f, /*Be*/ 0.21f, /*B*/  0.05f, /*C*/  1.90f, /*N*/  1.12f,
    /*O*/  1.08f, /*F*/  0.99f, /*Ne*/ 0.92f,
    /*Na*/ 0.93f, /*Mg*/ 0.51f, /*Al*/ 0.33f, /*Si*/ 0.20f, /*P*/  0.39f,
    /*S*/  0.20f, /*Cl*/ 1.47f, /*Ar*/ 1.34f,
    /*K*/  1.31f, /*Ca*/ 0.94f, /*Sc*/ 0.69f, /*Ti*/ 0.36f, /*V*/  0.48f,
    /*Cr*/ 0.33f, /*Mn*/ 0.26f, /*Fe*/ 0.48f, /*Co*/ 0.34f, /*Ni*/ 0.43f,
    /*Cu*/ 0.51f, /*Zn*/ 0.54f, /*Ga*/ 0.41f, /*Ge*/ 0.20f, /*As*/ 0.28f,
    /*Se*/ 0.22f, /*Br*/ 0.53f, /*Kr*/ 1.49f,
    /*Rb*/ 1.28f, /*Sr*/ 1.12f, /*Y*/  0.84f, /*Zr*/ 0.53f, /*Nb*/ 0.42f,
    /*Mo*/ 0.35f, /*Tc*/ 0.31f, /*Ru*/ 0.32f, /*Rh*/ 0.49f, /*Pd*/ 0.58f,
    /*Ag*/ 0.61f, /*Cd*/ 0.72f, /*In*/ 0.56f, /*Sn*/ 0.49f, /*Sb*/ 0.70f,
    /*Te*/ 0.37f, /*I*/  0.36f, /*Xe*/ 1.70f,
    /*Cs*/ 1.61f, /*Ba*/ 1.29f, /*La*/ 0.97f, /*Ce*/ 0.81f, /*Pr*/ 0.79f,
    /*Nd*/ 0.92f, /*Pm*/ 0.91f, /*Sm*/ 0.90f, /*Eu*/ 0.89f, /*Gd*/ 0.88f,
    /*Tb*/ 0.70f, /*Dy*/ 0.85f, /*Ho*/ 0.84f, /*Er*/ 0.83f, /*Tm*/ 0.82f,
    /*Yb*/ 0.81f, /*Lu*/ 0.80f, /*Hf*/ 0.52f, /*Ta*/ 0.58f, /*W*/  0.36f,
    /*Re*/ 0.32f, /*Os*/ 0.33f, /*Ir*/ 0.51f, /*Pt*/ 0.51f, /*Au*/ 0.51f,
    /*Hg*/ 0.90f, /*Tl*/ 0.69f, /*Pb*/ 0.59f, /*Bi*/ 0.70f, /*Po*/ 0.61f,
    /*At*/ 0.56f, /*Rn*/ 1.80f,
    /*Fr*/ 1.74f, /*Ra*/ 1.42f, /*Ac*/ 1.06f, /*Th*/ 0.88f, /*Pa*/ 0.72f,
    /*U*/  0.46f, /*Np*/ 0.65f, /*Pu*/ 0.65f, /*Am*/ 0.79f, /*Cm*/ 0.79f,
    /*Bk*/ 0.77f, /*Cf*/ 0.76f, /*Es*/ 1.00f, /*Fm*/ 1.00f, /*Md*/ 1.00f,
    /*No*/ 1.00f, /*Lr*/ 1.00f, /*Rf*/ 1.00f, /*Db*/ 1.00f, /*Sg*/ 1.00f,
    /*Bh*/ 1.00f, /*Hs*/ 1.00f, /*Mt*/ 1.00f, /*Ds*/ 1.00f, /*Rg*/ 1.00f,
    /*Cn*/ 1.00f, /*Nh*/ 1.00f, /*Fl*/ 1.00f, /*Mc*/ 1.00f, /*Lv*/ 1.00f,
    /*Ts*/ 1.00f, /*Og*/ 1.00f,
    /*D*/  1.40f, /*END*/0.f
  };
  static_assert(ce_almost_eq(radii[static_cast<int>(El::D)], 1.40f), "Hmm");
  return radii[static_cast<int>(el)];
#else
  // temporary solution used in Refmac
  switch (el) {
    case El::H: return 1.4f;
    case El::D: return 1.4f;
    case El::O: return 1.08f;
    case El::C: return 2.0f;
    case El::N: return 1.12f;
    default: return 1.6f;
  };
#endif
}

/// Marks grid points within a radius around atoms as masked.
/// Sets all grid points within (atomic radius + probe radius) of each atom to the given value.
/// @tparam T Grid value type (typically int8_t for binary masks or float for weighted masks)
/// @param mask Grid to modify in-place
/// @param model Molecular model containing atoms to mask around
/// @param atomic_radii_set Which atomic radii library to use
/// @param r_probe Probe radius (in Angstroms) added to each atomic radius
/// @param value Value to set for masked grid points
/// @param ignore_hydrogen If true, skip hydrogen atoms
/// @param ignore_zero_occupancy_atoms If true, skip atoms with zero occupancy
template<typename T>
void mask_points_in_radius(Grid<T>& mask, const Model& model,
                           AtomicRadiiSet atomic_radii_set,
                           double r_probe, T value,
                           bool ignore_hydrogen,
                           bool ignore_zero_occupancy_atoms) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        if ((ignore_hydrogen && atom.is_hydrogen()) ||
            (ignore_zero_occupancy_atoms && atom.occ <= 0))
          continue;
        El elem = atom.element.elem;
        double r = r_probe;
        switch (atomic_radii_set) {
          case AtomicRadiiSet::VanDerWaals: r += vdw_radius(elem); break;
          case AtomicRadiiSet::Cctbx: r += cctbx_vdw_radius(elem); break;
          case AtomicRadiiSet::Refmac: r += refmac_radius_for_bulk_solvent(elem); break;
          case AtomicRadiiSet::Constant: /* r is included in r_probe */ break;
        }
        mask.set_points_around(atom.pos, r, value);
      }
}

/// @deprecated Use mask_points_in_radius with AtomicRadiiSet::Constant instead.
/// Marks grid points within a constant radius around atoms.
/// @tparam T Grid value type
/// @param mask Grid to modify
/// @param model Molecular model with atoms to mask around
/// @param radius Constant radius in Angstroms
/// @param value Value to set for masked grid points
/// @param ignore_hydrogen If true, skip hydrogen atoms
/// @param ignore_zero_occupancy_atoms If true, skip atoms with zero occupancy
template<typename T>
void mask_points_in_constant_radius(Grid<T>& mask, const Model& model,
                                    double radius, T value,
                                    bool ignore_hydrogen,
                                    bool ignore_zero_occupancy_atoms) {
  mask_points_in_radius(mask, model, AtomicRadiiSet::Constant, radius, value,
                        ignore_hydrogen, ignore_zero_occupancy_atoms);
}


/// Collects all distinct alternate location indicators from a model.
/// @param model Molecular model to scan
/// @return String containing unique altloc characters found in the model
inline std::string distinct_altlocs(const Model& model) {
  std::string altlocs;
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      add_distinct_altlocs(res, altlocs);
  return altlocs;
}

/// Masks grid points using atom occupancy values to create a weighted mask.
/// Creates a non-binary mask by considering occupancy factors for atoms with alternate locations.
/// Each altloc is processed separately and contributions are accumulated.
/// @param mask Existing float mask to modify in-place (contributions are subtracted)
/// @param model Model to mask; atoms are categorized by alternate location code
/// @param atomic_radii_set Which atomic radii library to use
/// @param r_probe Probe radius (in Angstroms) added to each atomic radius
/// @param ignore_hydrogen If true, skip hydrogen atoms
/// @param ignore_zero_occupancy_atoms If true, skip atoms with zero occupancy
inline
void mask_points_using_occupancy(Grid<float>& mask, const Model& model,
                                  AtomicRadiiSet atomic_radii_set, double r_probe,
                                  bool ignore_hydrogen,
                                  bool ignore_zero_occupancy_atoms) {

  std::string altlocs = distinct_altlocs(model);
  altlocs += '\0';  // no altloc

  gemmi::Grid<float> m;
  m.copy_metadata_from(mask);
  for (char& altloc : altlocs) {
    m.fill(0.0);
    for (const Chain& chain : model.chains) {
      for (const Residue& res : chain.residues) {
        for (const Atom& atom : res.atoms) {
          if ((ignore_hydrogen && atom.is_hydrogen()) ||
              (ignore_zero_occupancy_atoms && atom.occ <= 0))
            continue;
          if (atom.altloc == altloc) {
            El elem = atom.element.elem;
            double r = r_probe;
            switch (atomic_radii_set) {
              case AtomicRadiiSet::VanDerWaals: r += vdw_radius(elem); break;
              case AtomicRadiiSet::Cctbx: r += cctbx_vdw_radius(elem); break;
              case AtomicRadiiSet::Refmac: r += refmac_radius_for_bulk_solvent(elem); break;
              case AtomicRadiiSet::Constant: /* r is included in r_probe */ break;
            }
            Fractional fpos = m.unit_cell.fractionalize(atom.pos);
            m.use_points_around<true>(fpos, r, [&](float& ref, double) {
                ref = std::min(ref, -atom.occ);
            }, false);
          }
        }
      }
    }

    // reduce starting mask
    for (size_t i = 0; i < m.data.size(); ++i) {
      if (m.data[i] < 0.0)
        mask.data[i] += m.data[i];
    }
  }

  // When we reduced the solvent mask by the atom occupancy each time
  // we hit a grid point, we might end up with a value below zero (due
  // to overlapping atoms) and need to limit it here.
  for (float& d : mask.data)
    d = std::max(d, 0.f);
}

/// Creates a margin of points around the boundary of a masked region.
/// Finds grid points at distance <= @p r from points with value equal to @p value,
/// and sets them to @p margin_value. Uses efficient stencil-based neighbor search.
/// @tparam T Grid value type
/// @param mask Grid to modify in-place
/// @param r Distance threshold in Angstroms
/// @param value Boundary value to search for (typically solvent/protein boundary)
/// @param margin_value Value to set for margin points
/// @throws gemmi::Failure if radius exceeds half the unit cell dimensions
template<typename T>
void set_margin_around(Grid<T>& mask, double r, T value, T margin_value) {
  int du = (int) std::floor(r / mask.spacing[0]);
  int dv = (int) std::floor(r / mask.spacing[1]);
  int dw = (int) std::floor(r / mask.spacing[2]);
  double max_spacing2 = sq(std::max(std::max(mask.unit_cell.a / mask.nu,
                                             mask.unit_cell.b / mask.nv),
                                             mask.unit_cell.c / mask.nw)) + 1e-6;
  if (2 * du >= mask.nu || 2 * dv >= mask.nv || 2 * dw >= mask.nw)
    fail("grid operation failed: radius bigger than half the unit cell?");
  std::vector<std::array<int,3>> stencil1;
  std::vector<std::array<int,3>> stencil2;
  for (int w = -dw; w <= dw; ++w)
    for (int v = -dv; v <= dv; ++v)
      for (int u = -du; u <= du; ++u) {
        Fractional fdelta = mask.get_fractional(u, v, w);
        double r2 = mask.unit_cell.orthogonalize_difference(fdelta).length_sq();
        if (r2 <= r * r && r2 != 0.) {
          std::array<int,3> wvu{{w <= 0 ? w : w - mask.nw,
                                 v <= 0 ? v : v - mask.nv,
                                 u <= 0 ? u : u - mask.nu}};
          if (r2 < max_spacing2)
            stencil1.push_back(wvu);
          else
            stencil2.push_back(wvu);
        }
      }
  if (stencil2.empty()) {
    // r is small; it should be faster to go through masked points and
    // for each one check if any point in given radius are unmasked.
    for (typename Grid<T>::Point p : mask)
      if (*p.value < value) {
        for (const auto& wvu : stencil1) {
          size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
          if (mask.data[idx] >= value) {
            *p.value = margin_value;
            break;
          }
        }
      }
  } else {
    // r is large; it should be faster to go through unmasked points
    // and for each one check if it has masked near neighbors (stencil1).
    // These neighbors get marked as margin. If all near neighbors are
    // unmasked, we can skip further neighbors (they will be checked
    // from other points).
    for (typename Grid<T>::Point p : mask) {
      if (*p.value >= value) {
        bool found = false;
        for (const auto& wvu : stencil1) {
          size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
          if (mask.data[idx] < value) {
            mask.data[idx] = margin_value;
            found = true;
          }
        }
        if (found)
          for (const auto& wvu : stencil2) {
            size_t idx = mask.index_near_zero(p.u + wvu[2], p.v + wvu[1], p.w + wvu[0]);
            if (mask.data[idx] < value)
              mask.data[idx] = margin_value;
          }
      }
    }
  }
  //printf("stencil sizes: %zu\n", stencil1.size(), stencil2.size());
  //printf("margin: %zu\n", std::count(mask.data.begin(), mask.data.end(), margin_value));
}

/// Helper class for computing and applying solvent masks to crystallographic grids.
/// Encapsulates parameters and operations for bulk solvent masking, including
/// mask generation, shrinking, inversion, and symmetry handling.
struct SolventMasker {
  AtomicRadiiSet atomic_radii_set;  ///< Which atomic radii library is used
  bool ignore_hydrogen;              ///< If true, hydrogen atoms are skipped
  bool ignore_zero_occupancy_atoms;  ///< If true, atoms with zero occupancy are skipped
  bool use_atom_occupancy = false;   ///< If true, use atom occupancy for weighted masking
  double rprobe;                     ///< Probe radius added to atomic radii (in A)
  double rshrink;                    ///< Shrinking radius applied after masking (in A)
  double island_min_volume;          ///< Minimum volume (as fraction) of protein islands to retain
  double constant_r;                 ///< Constant radius (for AtomicRadiiSet::Constant)
  double requested_spacing = 0.;     ///< Requested grid spacing (0 = auto)

  /// Initialize SolventMasker with a radii set and optional constant radius.
  /// Automatically sets default parameters (rprobe, rshrink) for the chosen set.
  /// @param choice Atomic radii set to use
  /// @param constant_r_ Constant radius (only used if choice is AtomicRadiiSet::Constant)
  SolventMasker(AtomicRadiiSet choice, double constant_r_=0.) {
    set_radii(choice, constant_r_);
  }

  /// Sets the atomic radii set and related parameters.
  /// Updates rprobe, rshrink, and island_min_volume based on the chosen library.
  /// @param choice Atomic radii set to use
  /// @param constant_r_ Constant radius override
  void set_radii(AtomicRadiiSet choice, double constant_r_=0.) {
    atomic_radii_set = choice;
    constant_r = constant_r_;
    ignore_hydrogen = true;
    ignore_zero_occupancy_atoms = true;
    switch (choice) {
      case AtomicRadiiSet::VanDerWaals:
        rprobe = 1.0;
        rshrink = 1.1;
        island_min_volume = 0.;
        break;
      case AtomicRadiiSet::Cctbx:
        rprobe = 1.1;
        rshrink = 0.9;
        island_min_volume = 0.;
        break;
      case AtomicRadiiSet::Refmac:
        rprobe = 1.0;
        rshrink = 0.8;
        island_min_volume = 50;  // the exact value used in Refmac is yet to be found
        break;
      case AtomicRadiiSet::Constant:
        rprobe = 0;
        rshrink = 0;
        island_min_volume = 0.;
        break;
    }
  }

  /// Fills the entire grid with 1 (solvent region).
  /// @tparam T Grid value type
  /// @param grid Grid to fill
  template<typename T> void clear(Grid<T>& grid) const { grid.fill((T)1); }

  /// Sets grid points around atoms to 0 (protein region).
  /// @tparam T Grid value type
  /// @param grid Grid to modify in-place
  /// @param model Molecular model
  template<typename T> void mask_points(Grid<T>& grid, const Model& model) const {
    mask_points_in_radius(grid, model, atomic_radii_set, constant_r + rprobe, (T)0,
                          ignore_hydrogen, ignore_zero_occupancy_atoms);
  }

  /// Sets grid points around atoms to 0, with optional occupancy weighting.
  /// Uses atom occupancy if @c use_atom_occupancy is true; otherwise calls mask_points<float>.
  /// @param grid Grid to modify in-place
  /// @param model Molecular model
  void mask_points(Grid<float>& grid, const Model& model) const {
    if (use_atom_occupancy)
      mask_points_using_occupancy(grid, model, atomic_radii_set, constant_r + rprobe,
                                  ignore_hydrogen, ignore_zero_occupancy_atoms);
    else
      mask_points<float>(grid, model);
  }

  /// Applies space group symmetry to fill the entire grid.
  /// For integer/binary masks, distributes 0-value points via symmetry operators.
  /// For float masks, sets each point to the minimum value across symmetry mates.
  /// @tparam T Grid value type
  /// @param grid Grid to symmetrize in-place
  template<typename T> void symmetrize(Grid<T>& grid) const {
    if (std::is_same<T, std::int8_t>::value) {
      grid.symmetrize([&](T a, T b) { return a == (T)0 || b == (T)0 ? (T)0 : (T)1; });
    }
    else {
      grid.symmetrize([&](T a, T b) { return a < b ? a : b; });
    }
  }

  /// Shrinks the masked (protein) region by marking a margin as solvent.
  /// Marks all points within @c rshrink of the protein-solvent boundary.
  /// @tparam T Grid value type
  /// @param grid Grid to shrink in-place
  template<typename T> void shrink(Grid<T>& grid) const {
    if (rshrink > 0) {
      set_margin_around(grid, rshrink, (T)1, (T)-1);
      grid.change_values((T)-1, (T)1);
    }
  }

  /// Inverts the mask (1 becomes 0, 0 becomes 1).
  /// @tparam T Grid value type
  /// @param grid Grid to invert in-place
  template<typename T> void invert(Grid<T>& grid) const {
    for (auto& v : grid.data)
      v = (T)1 - v;
  }

  /// Removes small disconnected regions (islands) of value 1 using flood fill.
  /// Islands smaller than @c island_min_volume (as a fraction of unit cell volume) are removed.
  /// @tparam T Grid value type (typically int8_t)
  /// @param grid Grid to modify in-place
  /// @return Number of islands removed
  /// @note Does not work with masks generated by mask_points_using_occupancy()
  template<typename T> int remove_islands(Grid<T>& grid) const {
    if (island_min_volume <= 0)
      return 0;
    size_t limit = static_cast<size_t>(island_min_volume * grid.point_count()
                                       / grid.unit_cell.volume);
    int counter = 0;
    FloodFill<T,1> flood_fill{grid};
    flood_fill.for_each_islands([&](typename FloodFill<T,1>::Result& r) {
        //printf("island %d: %zu in %zu (limit: %zu)\n",
        //       counter, r.point_count(), lines.size(), limit);
        if (r.point_count() <= limit) {
          ++counter;
          flood_fill.set_volume_values(r, (T)0);
        }
    });
    return counter;
  }

  /// Generates a complete solvent mask on a grid using the standard pipeline.
  /// Steps: clear grid (1) -> mask atoms (0) -> apply symmetry -> remove islands -> shrink.
  /// @tparam T Grid value type
  /// @param grid Grid to populate with mask
  /// @param model Molecular model to mask
  template<typename T> void put_mask_on_grid(Grid<T>& grid, const Model& model) const {
    clear(grid);
    assert(!grid.data.empty());
    mask_points(grid, model);
    symmetrize(grid);
    remove_islands(grid);
    shrink(grid);
  }

  /// Sets grid points around atoms to 0, applying symmetry without shrinking.
  /// Used to zero out a float density map in the protein region.
  /// @param grid Float grid to modify in-place
  /// @param model Molecular model
  void set_to_zero(Grid<float>& grid, const Model& model) const {
    mask_points(grid, model);
    grid.symmetrize([&](float a, float b) { return b == 0.f ? 0.f : a; });
  }

#if 0
  template<typename T> void put_mask_on_grid(Grid<T>& grid, const Model& model) {
    // use twice finer grid for solvent mask
    Grid<std::int8_t> mask;
    mask.copy_metadata_from(grid);
    mask.set_size(2*grid.nu, 2*grid.nv, 2*grid.nw);
    mask.data.resize(8 * grid.data.size(), 1);
    put_mask_on_grid(mask, model);
    for (int w = 0, idx = 0; w < grid.nw; ++w)
      for (int v = 0; v < grid.nv; ++v)
        for (int u = 0; u < grid.nu; ++u, ++idx) {
          grid.data[idx] = 0;
          for (int wa = 0; wa < 2; ++wa)
            for (int va = 0; va < 2; ++va)
              for (int ua = 0; ua < 2; ++ua)
                grid.data[idx] += mask.get_value_q(2*u + ua, 2*v + va, 2*w + wa);
          grid.data[idx] *= 1. / 8;
        }
#endif
};

/// Information about a grid point's relationship to nearby model atoms.
/// Used internally for interpolation of density maps around model atoms.
struct NodeInfo {
  double dist_sq;  ///< Square of distance from the nearest atom
  bool found = false;  ///< True if a nearby atom was found within radius
  int u = 0, v = 0, w = 0;  ///< Non-normalized near-model grid coordinates of nearest atom
};

/// Populates a grid with NodeInfo for all grid points near model atoms.
/// Finds the nearest atom and its distance for each grid point within @p radius.
/// @param mask NodeInfo grid to populate
/// @param model Molecular model to search
/// @param radius Search radius in Angstroms
inline void mask_with_node_info(Grid<NodeInfo>& mask, const Model& model, double radius) {
  NodeInfo default_ni;
  default_ni.dist_sq = radius * radius;
  mask.fill(default_ni);
  // cf. use_points_around()
  int du = (int) std::ceil(radius / mask.spacing[0]);
  int dv = (int) std::ceil(radius / mask.spacing[1]);
  int dw = (int) std::ceil(radius / mask.spacing[2]);
  mask.template check_size_for_points_in_box<true>(du, dv, dw, false);
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        Fractional frac0 = mask.unit_cell.fractionalize(atom.pos);
        mask.template do_use_points_in_box<true>(
            frac0, du, dv, dw,
            [&](NodeInfo& ni, double d2, const Position&, int u, int v, int w) {
              if (d2 < ni.dist_sq) {
                ni.dist_sq = d2;
                ni.found = true;
                //ni.elem = atom.element;
                ni.u = u;
                ni.v = v;
                ni.w = w;
              }
            },
            radius);
      }
}

/// Removes grid points that are closer to a symmetry mate than to the original model.
/// A grid point is unmasked (marked as @c found=false) if any symmetry image of the
/// nearest atom is closer than that atom. This avoids double-counting in density calculations.
/// Non-crystallographic symmetry (NCS) is ignored.
/// @param mask NodeInfo grid to update in-place
inline void unmask_symmetry_mates(Grid<NodeInfo>& mask) {
  std::vector<GridOp> symmetry_ops = mask.get_scaled_ops_except_id();
  size_t idx = 0;
  for (int w = 0; w != mask.nw; ++w)
    for (int v = 0; v != mask.nv; ++v)
      for (int u = 0; u != mask.nu; ++u, ++idx) {
        NodeInfo& ni = mask.data[idx];
        if (ni.found)
          for (const GridOp& grid_op : symmetry_ops) {
            std::array<int,3> t = grid_op.apply(u, v, w);
            NodeInfo& im_ni = mask.data[mask.index_n(t[0], t[1], t[2])];
            if (im_ni.found && im_ni.dist_sq > ni.dist_sq)
              im_ni.found = false;
          }
      }
}

/// Interpolates grid values from a source grid around atoms in a destination model.
/// Identifies grid points in the destination grid that are near atoms in the destination model,
/// transforms them to source grid coordinates, and interpolates values from the source.
/// Grid points closer to symmetry mates are skipped to avoid double-counting.
/// @tparam T Grid value type
/// @param dest Destination grid to interpolate into (modified in-place)
/// @param src Source grid to interpolate from
/// @param tr Transformation from destination to source
/// @param dest_model Model in destination grid frame (determines which points to interpolate)
/// @param radius Search radius in Angstroms for atoms
/// @param order Interpolation order (1=linear, 3=cubic, default 1)
template<typename T>
void interpolate_grid_around_model(Grid<T>& dest, const Grid<T>& src,
                                   const Transform& tr,
                                   const Model& dest_model, double radius,
                                   int order=1) {
  Grid<NodeInfo> mask;
  mask.copy_metadata_from(dest);
  mask_with_node_info(mask, dest_model, radius);
  unmask_symmetry_mates(mask);
  // Interpolate values for selected nodes.
  FTransform frac_tr = src.unit_cell.frac.combine(tr.combine(dest.unit_cell.orth));
  for (size_t idx = 0; idx != mask.data.size(); ++idx) {
    const NodeInfo& ni = mask.data[idx];
    if (ni.found) {
      Fractional dest_fr = dest.get_fractional(ni.u, ni.v, ni.w);
      Fractional src_fr = frac_tr.apply(dest_fr);
      dest.data[idx] = src.interpolate_value(src_fr, order);
    }
  }
}


/// Adds a smooth transition zone to the boundary of a binary mask.
/// Converts sharp 0/1 boundaries to smooth transitions using a raised cosine function.
/// Grid points at distance < @p width from the boundary are set to cosine-interpolated values.
/// @tparam T Grid value type
/// @param grid Binary mask to smooth in-place (0s become 1s beyond boundary)
/// @param width Width of transition zone in Angstroms
template<typename T>
void add_soft_edge_to_mask(Grid<T>& grid, double width) {
  const double width2 = width * width;
  const int du = (int) std::ceil(width / grid.spacing[0]);
  const int dv = (int) std::ceil(width / grid.spacing[1]);
  const int dw = (int) std::ceil(width / grid.spacing[2]);

  for (int w = 0; w < grid.nw; ++w)
    for (int v = 0; v < grid.nv; ++v)
      for (int u = 0; u < grid.nu; ++u) {
        size_t idx = grid.index_q(u, v, w);
        if (grid.data[idx] >= 1e-3) continue;
        double min_d2 = width2 + 1;
        Fractional fctr = grid.get_fractional(u, v, w);
        grid.template use_points_in_box<true>(
            fctr, du, dv, dw,
            [&](T& point, double d2, const Position&, int, int, int) {
              if (point > 0.999) {
                if (d2 < min_d2)
                min_d2 = d2;
              }
            });
        if (min_d2 < width2)
          grid.data[idx] = T(0.5 + 0.5 * std::cos(pi() * std::sqrt(min_d2) / width));
      }
}

} // namespace gemmi
#endif
