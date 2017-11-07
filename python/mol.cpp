// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/unitcell.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

void add_mmread(py::module& mol) {
  py::class_<Position>(mol, "Position")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Position::x)
    .def_readwrite("y", &Position::y)
    .def_readwrite("z", &Position::z)
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + std::to_string(self.x) + "," +
               std::to_string(self.y) + "," + std::to_string(self.z) + ")>";
    });
  mol.def("calculate_dihedral", &calculate_dihedral,
          "Input: four points. Output: dihedral angle in radians.");
  py::class_<Element>(mol, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });
}
