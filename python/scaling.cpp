// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include "gemmi/scaling.hpp"

namespace py = pybind11;

void add_scaling(py::module& m) {
  using Scaling = gemmi::Scaling<double>;
  py::class_<Scaling>(m, "Scaling")
    .def(py::init<const gemmi::UnitCell&, const gemmi::SpaceGroup*>())
    .def_readwrite("cell", &Scaling::cell)
    .def_readwrite("crystal_system", &Scaling::crystal_system)
    .def_readwrite("k_overall", &Scaling::k_overall)
    .def_property("b_overall", &Scaling::get_b_overall, &Scaling::set_b_overall)
    .def_readwrite("use_solvent", &Scaling::use_solvent)
    .def_readwrite("k_sol", &Scaling::k_sol)
    .def_readwrite("b_sol", &Scaling::b_sol)
    .def("prepare_points", &Scaling::prepare_points)
    .def("fit_isotropic_b_approximately", &Scaling::fit_isotropic_b_approximately)
    .def("fit_parameters", &Scaling::fit_parameters)
    .def("scale_data", &Scaling::scale_data)
    ;
}
