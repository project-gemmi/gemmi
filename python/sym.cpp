// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi::sym;

void init_sym(py::module& sym) {
  py::class_<Op>(sym, "Op")
    .def(py::init<>(&Op::identity))
    .def(py::init(&parse_triplet))
    .def_property_readonly_static("TDEN",
                                  [](py::object) -> int { return Op::TDEN; })
    .def_readwrite("rot", &Op::rot)
    .def_readwrite("tran", &Op::tran)
    .def("triplet", &Op::triplet)
    .def("det_rot", &Op::det_rot)
    .def("inverted", &Op::inverted)
    .def("negated", &Op::negated)
    .def("translated", &Op::translated)
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
    .def("__hash__", [](const Op &self) { return std::hash<Op>()(self); })
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

  py::class_<GroupOps>(sym, "GroupOps")
    .def(py::init<>())
    .def("__iter__", [](const GroupOps& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def_readwrite("sym_ops", &GroupOps::sym_ops)
    .def_readwrite("cen_ops", &GroupOps::cen_ops)
    .def("change_basis", &GroupOps::change_basis);

  sym.def("generators_from_hall", &generators_from_hall, py::arg("hall"),
          "Parse Hall notation.");
  sym.def("symops_from_hall", &symops_from_hall, py::arg("hall"),
          "Parse Hall notation.");
}
