// Copyright 2018 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/smcif.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

void add_smcif(py::module& m) {
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });

  py::class_<AtomicStructure> atomic_structure(m, "AtomicStructure");
  py::class_<AtomicStructure::Site>(atomic_structure, "Site")
    .def_readonly("label", &AtomicStructure::Site::label)
    .def_readonly("type_symbol", &AtomicStructure::Site::type_symbol)
    .def_readonly("fract", &AtomicStructure::Site::fract)
    .def("__repr__", [](const AtomicStructure::Site& self) {
        return "<gemmi.AtomicStructure.Site " + self.label + ">";
    });

  atomic_structure
    .def(py::init<>())
    .def_readwrite("name", &AtomicStructure::name)
    .def_readwrite("cell", &AtomicStructure::cell)
    .def_readonly("spacegroup_hm", &AtomicStructure::spacegroup_hm)
    .def_readonly("sites", &AtomicStructure::sites)
    .def("get_all_unit_cell_sites", &AtomicStructure::get_all_unit_cell_sites)
    .def("__repr__", [](const AtomicStructure& self) {
        return "<gemmi.AtomicStructure: " + std::string(self.name) + ">";
    });
}
