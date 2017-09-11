// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;
using namespace gemmi::sym;

void init_sym(py::module& sym) {
  py::class_<SymOp>(sym, "SymOp")
    .def(py::init<>())
    .def_readwrite("rot", &SymOp::rot)
    .def_readwrite("tr", &SymOp::tr)
    .def("compose", &SymOp::compose)
    /*
    .def("__repr__", [](const SymOp &self) {
        return "<gemmi.sym.SymOp " + self.rot + self.tr + ">";
    })
    */
    ;

  sym.def("parse_triplet", &parse_triplet, py::arg("s"),
          "Parse coordinate triplet into SymOp.");
  sym.def("parse_triplet_part", &parse_triplet_part, py::arg("s"),
          "Parse one of the three parts of a triplet.");
  sym.def("make_triplet_part", &make_triplet_part,
          "Make one of the three parts of a triplet.");
}
