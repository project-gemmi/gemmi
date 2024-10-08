// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include "array.h"
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/vector.h>  // for Scaling::get_parameters
#include "gemmi/scaling.hpp"

void add_scaling(nb::module_& m) {
  m.def("adp_symmetry_constraints", &gemmi::adp_symmetry_constraints);
  using Scaling = gemmi::Scaling<float>;
  using FPhiData = gemmi::AsuData<std::complex<float>>;
  nb::class_<Scaling>(m, "Scaling")
    .def(nb::init<const gemmi::UnitCell&, const gemmi::SpaceGroup*>())
    .def_rw("cell", &Scaling::cell)
    .def_rw("k_overall", &Scaling::k_overall)
    .def_prop_rw("b_overall", &Scaling::get_b_overall, &Scaling::set_b_overall)
    .def_rw("use_solvent", &Scaling::use_solvent)
    .def_rw("k_sol", &Scaling::k_sol)
    .def_rw("b_sol", &Scaling::b_sol)
    .def_prop_rw("parameters", &Scaling::get_parameters,
                  (void (Scaling::*)(const std::vector<double>&)) &Scaling::set_parameters)
    .def("prepare_points", &Scaling::prepare_points,
         nb::arg("calc"), nb::arg("obs"), nb::arg("mask")=static_cast<FPhiData*>(nullptr))
    .def("fit_isotropic_b_approximately", &Scaling::fit_isotropic_b_approximately)
    .def("fit_b_star_approximately", &Scaling::fit_b_star_approximately)
    .def("fit_parameters", &Scaling::fit_parameters)
#if WITH_NLOPT
    .def("fit_parameters_with_nlopt", &gemmi::fit_parameters_with_nlopt<float>)
#endif
    .def("calculate_r_factor", &Scaling::calculate_r_factor)
    .def("get_overall_scale_factor", &Scaling::get_overall_scale_factor, nb::arg("hkl"))
    .def("get_overall_scale_factor", [](const Scaling& self, const cpu_miller_array& hkl) {
        return miller_function<double>(self, &Scaling::get_overall_scale_factor, hkl);
    })
    .def("get_solvent_scale", &Scaling::get_solvent_scale, nb::arg("stol2"))
    .def("scale_data", &Scaling::scale_data,
         nb::arg("asu_data"), nb::arg("mask_data")=static_cast<FPhiData*>(nullptr))
    .def("scale_value", &Scaling::scale_value,
         nb::arg("hkl"), nb::arg("f_value"), nb::arg("mask_value"))
    ;
}
