// Copyright 2019 Global Phasing Ltd.
//
// Tools to make a grid with:
// - values of electron density of a model,
// - bulk solvent mask.

#ifndef GEMMI_RHOGRID_HPP_
#define GEMMI_RHOGRID_HPP_

#include <cassert>
#include "grid.hpp"    // for Grid
#include "model.hpp"   // for Structure, ...

namespace gemmi {

// data from cctbx/eltbx/van_der_waals_radii.py used to generate identical mask
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
  static_assert(radii[static_cast<int>(El::D)] == 1.2f, "Hmm");
  static_assert(sizeof(radii) / sizeof(radii[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return radii[static_cast<int>(el)];
}

// mask utilities
template<typename Real>
void mask_points_in_constant_radius(Grid<Real>& mask, const Model& model,
                                    double radius, Real value) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms)
        mask.set_points_around(atom.pos, radius, value);
}

template<typename Real>
void mask_points_in_vdw_radius(Grid<Real>& mask, const Model& model,
                               double r_probe, Real value, bool cctbx_compat) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        double r = cctbx_compat ? cctbx_vdw_radius(atom.element.elem)
                                : vdw_radius(atom.element.elem);
        mask.set_points_around(atom.pos, r + r_probe, value);
      }
}

template<typename Real>
void set_margin_around(Grid<Real>& mask, double r, Real value, Real margin_value) {
  int du = (int) std::floor(r / mask.spacing[0]);
  int dv = (int) std::floor(r / mask.spacing[1]);
  int dw = (int) std::floor(r / mask.spacing[2]);
  if (2 * du >= mask.nu || 2 * dv >= mask.nv || 2 * dw >= mask.nw)
    fail("grid operation failed: radius bigger than half the unit cell?");
  std::vector<std::int8_t> stencil;
  stencil.reserve((2*dw+1) * (2*dv+1) * (2*du+1));
  for (int w = -dw; w <= dw; ++w)
    for (int v = -dv; v <= dv; ++v)
      for (int u = -du; u <= du; ++u) {
        Fractional fdelta{u * (1.0 / mask.nu), v * (1.0 / mask.nv), w * (1.0 / mask.nw)};
        double r2 = mask.unit_cell.orthogonalize_difference(fdelta).length_sq();
        stencil.push_back(r2 <= r * r && r2 != 0.);
      }
  for (int w = 0; w < mask.nw; ++w)
    for (int v = 0; v < mask.nv; ++v)
      for (int u = 0; u < mask.nu; ++u) {
        Real& point = mask.data[mask.index_q(u, v, w)];
        if (point != value) {
          for (int w2 = w-dw, idx = 0; w2 <= w+dw; ++w2)
            for (int v2 = v-dv; v2 <= v+dv; ++v2)
              for (int u2 = u-du; u2 <= u+du; ++u2, ++idx)
                if (stencil[idx] && mask.data[mask.index_n(u2, v2, w2)] == value) {
                  point = margin_value;
                  goto nextpoint;
                }
        }
nextpoint: ;
      }
}


template <typename F>
double determine_cutoff_radius(const F& func, float cutoff_level) {
  float x1 = 3.5f;
  float y1 = func(x1);
  float x2 = x1;
  float y2 = y1;
  if (y1 < cutoff_level)
    while (y1 < cutoff_level) {
      x2 = x1;
      y2 = y1;
      x1 -= 0.5f;
      y1 = func(x1);
    }
  else
    while (y2 > cutoff_level) {
      x1 = x2;
      y1 = y2;
      x2 += 0.5f;
      y2 = func(x2);
    }
  while (x2 - x1 > 0.02f) {
    float new_x = 0.5f * (x2 + x1);
    float new_y = func(new_x);
    if (new_y < cutoff_level) {
      x2 = new_x;
      y2 = new_y;
    } else {
      x1 = new_x;
      y1 = new_y;
    }
  }
  return x2;
}

// Usual usage:
// - set d_min and optionally also other parameters,
// - set fprimes to f' values for your wavelength (see fprime.hpp)
// - use set_grid_cell_and_spacegroup() to set grid's unit cell and space group
// - check that Table has SF coefficients for all elements that are to be used
// - call put_model_density_on_grid()
// - do FFT using transform_map_to_f_phi()
// - if blur is used, multiply the SF by reciprocal_space_multiplier()
template <typename Table, typename Real>
struct DensityCalculator {
  Grid<Real> grid;
  double d_min = 0.;
  double rate = 1.5;
  double blur = 0.;
  float r_cut = 5e-5f;
  std::vector<float> fprimes = std::vector<float>((int)El::END, 0.f);
  // parameters for used only in put_solvent_mask_on_grid()
  double rprobe = 1.0;
  double rshrink = 1.1;
  bool cctbx_vdw = false;

  // pre: check if Table::has(atom.element)
  void add_atom_density_to_grid(const Atom& atom) {
    constexpr double UtoB = 8 * sq(pi());
    auto& scat = Table::get(atom.element);
    float fprime = fprimes[atom.element.ordinal()];
    Fractional fpos = grid.unit_cell.fractionalize(atom.pos);
    if (!atom.aniso.nonzero()) {
      // isotropic
      double b = atom.b_iso + blur;
      auto precal = scat.precalculate_density_iso(b, fprime);
      double radius = determine_cutoff_radius(
          [&](float r) { return (float)precal.calculate(r*r); },
          r_cut);
      grid.use_points_around(fpos, radius, [&](Real& point, double r2) {
          point += Real(atom.occ * precal.calculate((Real)r2));
      }, /*fail_on_too_large_radius=*/false);
    } else {
      // anisotropic
      SMat33<double> aniso_b = atom.aniso.scaled(UtoB).added_kI(blur);
      // rough estimate, so we don't calculate eigenvalues
      double b_max = std::max(std::max(aniso_b.u11, aniso_b.u22), aniso_b.u33);
      auto precal_iso = scat.precalculate_density_iso(b_max, fprime);
      double radius = determine_cutoff_radius(
          [&](float r) { return (float)precal_iso.calculate(r*r); },
          r_cut);
      auto precal = scat.precalculate_density_aniso_b(aniso_b, fprime);
      int du = (int) std::ceil(radius / grid.spacing[0]);
      int dv = (int) std::ceil(radius / grid.spacing[1]);
      int dw = (int) std::ceil(radius / grid.spacing[2]);
      grid.use_points_in_box(fpos, du, dv, dw,
                             [&](Real& point, const Position& delta) {
        if (delta.length_sq() < radius * radius)
          point += Real(atom.occ * precal.calculate(delta));
      }, false);
    }
  }

  void add_model_density_to_grid(const Model& model) {
    for (const Chain& chain : model.chains)
      for (const Residue& res : chain.residues)
        for (const Atom& atom : res.atoms)
          add_atom_density_to_grid(atom);
  }

  void put_model_density_on_grid(const Model& model) {
    grid.data.clear();
    grid.set_size_from_spacing(d_min / (2 * rate), true);
    add_model_density_to_grid(model);
    grid.symmetrize([](Real a, Real b) { return a + b; });
  }

  void set_grid_cell_and_spacegroup(const Structure& st) {
    grid.unit_cell = st.cell;
    grid.spacegroup = st.find_spacegroup();
  }

  // The argument is 1/d^2 - as outputted by unit_cell.calculate_1_d2(hkl).
  double reciprocal_space_multiplier(double inv_d2) {
    return std::exp(blur * 0.25 * inv_d2);
  }

  void put_solvent_mask_on_grid(const Model& model) {
    assert(!grid.data.empty());
    std::fill(grid.data.begin(), grid.data.end(), 1);
    mask_points_in_vdw_radius<Real>(grid, model, rprobe, 0, cctbx_vdw);
    set_margin_around<Real>(grid, rshrink, 1, -1);
    grid.change_values(-1, 1);
  }
};

} // namespace gemmi
#endif
