// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>    // for vectorize
#include "gemmi/scaling.hpp"

namespace py = pybind11;

void add_scaling(py::module& m) {
  using Scaling = gemmi::Scaling<float>;
  using FPhiData = gemmi::AsuData<std::complex<float>>;
  py::class_<Scaling>(m, "Scaling")
    .def(py::init<const gemmi::UnitCell&, const gemmi::SpaceGroup*>())
    .def_readwrite("cell", &Scaling::cell)
    .def_readwrite("crystal_system", &Scaling::crystal_system)
    .def_readwrite("k_overall", &Scaling::k_overall)
    .def_property("b_overall", &Scaling::get_b_overall, &Scaling::set_b_overall)
    .def_readwrite("use_solvent", &Scaling::use_solvent)
    .def_readwrite("k_sol", &Scaling::k_sol)
    .def_readwrite("b_sol", &Scaling::b_sol)
    .def("prepare_points", &Scaling::prepare_points,
         py::arg("calc"), py::arg("obs"), py::arg("mask")=FPhiData())
    .def("fit_isotropic_b_approximately", &Scaling::fit_isotropic_b_approximately)
    .def("fit_parameters", &Scaling::fit_parameters)
    .def("get_overall_scale_factor", &Scaling::get_overall_scale_factor, py::arg("hkl"))
    .def("get_overall_scale_factor", [](const Scaling& self, py::array_t<int> hkl) {
        auto h = hkl.unchecked<2>();
        if (h.shape(1) != 3)
          throw std::domain_error("the hkl array must have size N x 3");
        auto len = h.shape(0);
        py::array_t<double> arr(len);
        double* ptr = (double*) arr.request().ptr;
        for (int i = 0; i < len; ++i)
          ptr[i] = self.get_overall_scale_factor({{h(i, 0), h(i, 1), h(i, 2)}});
        return arr;
    })
    .def("get_solvent_scale", py::vectorize(&Scaling::get_solvent_scale), py::arg("stol2"))
    .def("scale_data", &Scaling::scale_data,
         py::arg("asu_data"), py::arg("mask_data")=FPhiData())
    .def("scale_value", &Scaling::scale_value,
         py::arg("hkl"), py::arg("f_value"), py::arg("mask_value"))
    ;
}
