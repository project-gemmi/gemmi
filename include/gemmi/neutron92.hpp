// Copyright 2022 Global Phasing Ltd.

// Neutron coherent scattering lengths of the elements,
// from Neutron News, Vol. 3, No. 3, 1992.
//
// We use the same data as cctbx/eltbx/neutron.h, which is based on
// https://www.ncnr.nist.gov/resources/n-lengths/list.html
// which in turn is based on Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.

#ifndef GEMMI_NEUTRON92_HPP_
#define GEMMI_NEUTRON92_HPP_

#include "formfact.hpp"  // for GaussianCoef
#include "elem.hpp"      // for El

namespace gemmi {

#if defined(__GNUC__) && __GNUC__-0 > 4
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

template<class Real>
struct Neutron92 {
  using Coef = GaussianCoef<0, 1, Real>;
  static Real data[121];

  static Real& get_(El el) { return data[static_cast<int>(el)]; }
  static bool has(El el) { return get_(el) != 0; }
  static Coef get(El el) { return Coef{{get_(el)}}; }
  //static Real* get_ptr(El el) { return &get_(el); }
};

// real part of the bound coherent neutron scattering length (fm)
template<class Real>
Real Neutron92<Real>::data[121] = {
  /*X*/ 0.0,
  /*H*/ -3.7390, /*He*/ 3.26,
  /*Li*/ -1.90, /*Be*/ 7.79, /*B*/ 5.30, /*C*/ 6.646,
  /*N*/ 9.36, /*O*/ 5.803, /*F*/ 5.654, /*Ne*/ 4.566,
  /*Na*/ 3.63, /*Mg*/ 5.375, /*Al*/ 3.449, /*Si*/ 4.1491,
  /*P*/ 5.13, /*S*/ 2.847, /*Cl*/ 9.577, /*Ar*/ 1.909,
  /*K*/ 3.67, /*Ca*/ 4.70, /*Sc*/ 12.29, /*Ti*/ -3.438,
  /*V*/ -0.3824, /*Cr*/ 3.635, /*Mn*/ -3.73, /*Fe*/ 9.45,
  /*Co*/ 2.49, /*Ni*/ 10.3, /*Cu*/ 7.718, /*Zn*/ 5.68,
  /*Ga*/ 7.288, /*Ge*/ 8.185, /*As*/ 6.58, /*Se*/ 7.97,
  /*Br*/ 6.795, /*Kr*/ 7.81, /*Rb*/ 7.09, /*Sr*/ 7.02,
  /*Y*/ 7.75, /*Zr*/ 7.16, /*Nb*/ 7.054,
  /*Mo*/ 6.715, /*Tc*/ 6.8, /*Ru*/ 7.03, /*Rh*/ 5.88, /*Pd*/ 5.91,
  /*Ag*/ 5.922, /*Cd*/ 4.87, /*In*/ 4.065, /*Sn*/ 6.225,
  /*Sb*/ 5.57, /*Te*/ 5.8, /*I*/ 5.28, /*Xe*/ 4.92,
  /*Cs*/ 5.42, /*Ba*/ 5.07, /*La*/ 8.24, /*Ce*/ 4.84,
  /*Pr*/ 4.58, /*Nd*/ 7.69, /*Pm*/ 12.6, /*Sm*/ 0.8,
  /*Eu*/ 7.22, /*Gd*/ 6.5, /*Tb*/ 7.38, /*Dy*/ 16.9,
  /*Ho*/ 8.01, /*Er*/ 7.79, /*Tm*/ 7.07, /*Yb*/ 12.43,
  /*Lu*/ 7.21, /*Hf*/ 7.7, /*Ta*/ 6.91, /*W*/ 4.86,
  /*Re*/ 9.2, /*Os*/ 10.7, /*Ir*/ 10.6, /*Pt*/ 9.6,
  /*Au*/ 7.63, /*Hg*/ 12.692, /*Tl*/ 8.776,
  /*Pb*/ 9.405, /*Bi*/ 8.532, /*Po*/ 0., /*At*/ 0., /*Rn*/ 0.,
  /*Fr*/ 0., /*Ra*/ 10.0, /*Ac*/ 0., /*Th*/ 10.31, /*Pa*/ 9.1,
  /*U*/ 8.417, /*Np*/ 10.55, /*Pu*/ 0., /*Am*/ 8.3, /*Cm*/ 9.5,
  /*Bk*/ 0., /*Cf*/ 0., /*Es*/ 0., /*Fm*/ 0., /*Md*/ 0.,
  /*No*/ 0., /*Lr*/ 0., /*Rf*/ 0., /*Db*/ 0., /*Sg*/ 0.,
  /*Bh*/ 0., /*Hs*/ 0., /*Mt*/ 0., /*Ds*/ 0., /*Rg*/ 0., /*Cn*/ 0.,
  /*Nh*/ 0., /*Fl*/ 0., /*Mc*/ 0., /*Lv*/ 0., /*Ts*/ 0., /*Og*/ 0.,
  /*D*/ 6.671, /*END*/ 0.
};

#if defined(__GNUC__) && __GNUC__-0 > 4
#pragma GCC diagnostic pop
#endif

} // namespace gemmi
#endif
