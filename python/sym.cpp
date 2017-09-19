// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;
using namespace gemmi::sym;

void init_sym(py::module& sym) {
  py::class_<Op>(sym, "Op")
    .def(py::init<>())
    .def(py::init(&parse_triplet))
    .def_readwrite("rot", &Op::rot)
    .def_readwrite("tran", &Op::tran)
    .def("triplet", &Op::triplet)
    .def("det_rot", &Op::det_rot)
    .def("inverted", &Op::inverted)
    .def("negated", &Op::negated)
    .def("__mul__", [](const Op &a, const Op &b) { return combine(a, b); },
         py::is_operator())
    .def("__mul__", [](const Op &a, const std::string &b) {
            return combine(a, parse_triplet(b));
         }, py::is_operator())
    .def("__rmul__", [](const Op &a, const std::string &b) {
            return combine(parse_triplet(b), a);
         }, py::is_operator())
    .def("__eq__", [](const Op &a, const Op &b) { return a == b; },
         py::is_operator())
    .def("__eq__", [](const Op &a, const std::string& b) {
            return a == parse_triplet(b);
         }, py::is_operator())
    .def("__repr__", [](const Op &self) {
        return "<sym.Op(\"" + self.triplet() + "\")>";
    });

  sym.def("parse_triplet", &parse_triplet, py::arg("triplet"),
          "Parse coordinate triplet into sym.Op.");
  sym.def("parse_triplet_part", &parse_triplet_part, py::arg("s"),
          "Parse one of the three parts of a triplet.");
  sym.def("make_triplet_part", &make_triplet_part,
          py::arg("x"), py::arg("y"), py::arg("z"), py::arg("w"),
          "Make one of the three parts of a triplet.");
  sym.def("combine", &combine, py::arg("a"), py::arg("b"),
          "Combine two symmetry operations.");

  py::class_<SymOps>(sym, "SymOps")
    .def(py::init<>())
    .def_readwrite("sym_ops", &SymOps::sym_ops)
    .def_readwrite("cen_ops", &SymOps::cen_ops);

  sym.def("symops_from_hall", &symops_from_hall, py::arg("hall"),
          "Parse Hall notation.");
}
