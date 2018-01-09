// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/unitcell.hpp"

#include <cstdio>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

void add_mm(py::module& m) {
  py::class_<Position>(m, "Position")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Position::x)
    .def_readwrite("y", &Position::y)
    .def_readwrite("z", &Position::z)
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + std::to_string(self.x) + "," +
               std::to_string(self.y) + "," + std::to_string(self.z) + ")>";
    });
  py::class_<UnitCell>(m, "UnitCell")
    .def(py::init<>())
    // TODO: init(a,b,c,alpha,beta,gamma)
    .def_readonly("a", &UnitCell::a)
    .def_readonly("b", &UnitCell::b)
    .def_readonly("c", &UnitCell::c)
    .def_readonly("alpha", &UnitCell::alpha)
    .def_readonly("beta", &UnitCell::beta)
    .def_readonly("gamma", &UnitCell::gamma)
    .def_readonly("volume", &UnitCell::volume)
    .def("set", &UnitCell::set)
    .def("__repr__", [](const UnitCell& self) {
        char buf[256];
        std::snprintf(buf, 256, "<gemmi.UnitCell(%.9g, %.9g, %.9g, %.9g, %.9g, %.9g)>",
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma);
        return std::string(buf);
    });

  m.def("calculate_dihedral", &calculate_dihedral,
        "Input: four points. Output: dihedral angle in radians.");
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });
}
