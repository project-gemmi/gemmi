// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_cif(py::module& cif); // cif.cpp

PYBIND11_MODULE(gemmi, mg) {
  using namespace gemmi::mol;

  mg.doc() = "General MacroMolecular I/O";
  py::module cif = mg.def_submodule("cif", "CIF file format");
  init_cif(cif);

  py::module mol = mg.def_submodule("mol", "MacroMolecular models");
  py::class_<Element>(mol, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.mol.Element: " + std::string(self.name()) + ">";
    });
}

// vim:sw=2:ts=2:et
