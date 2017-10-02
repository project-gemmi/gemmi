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
            [](py::object) -> int { return Op::TDEN; },
            "Denominator (integer) for the translation vector.")
    .def_readwrite("rot", &Op::rot, "3x3 integer matrix.")
    .def_readwrite("tran", &Op::tran,
       "Numerators (integers) of the translation vector. Denominator TDEN=24.")
    .def("triplet", &Op::triplet, "Returns coordinate triplet x,y,z.")
    .def("det_rot", &Op::det_rot, "Determinant of the 3x3 matrix.")
    .def("inverse", [](const Op& self) { return self.inverse(); },
         "Returns inverted operator.")
    .def("negated", &Op::negated, "Returns Op with all elements nagated")
    .def("translated", &Op::translated, py::arg("a"), "Adds a to tran")
    .def("wrap", &Op::wrap, "Wrap the translation part to [0,1)")
    .def("seitz", [](const Op& self) {
            auto arr = self.seitz();
            auto mat = py::list();
            py::object fr = py::module::import("fractions").attr("Fraction");
            for (int i = 0; i < 4; ++i) {
              auto row = py::list();
              for (int j = 0; j < 4; ++j) {
                auto v = arr[i][j];
                if (j == 3 && i != 3 && v != 0)
                  row.append(fr(v, (int)Op::TDEN));
                else
                  row.append(v);
              }
              mat.append(row);
            }
            return mat;
         }, "Returns Seitz matrix (integers and fractions)")
    .def("float_seitz", &Op::float_seitz, "Returns Seitz matrix (floats)")
    .def("__mul__", [](const Op &a, const Op &b) {
            return combine(a, b).wrap();
         }, py::is_operator())
    .def("__mul__", [](const Op &a, const std::string &b) {
            return combine(a, parse_triplet(b)).wrap();
         }, py::is_operator())
    .def("__rmul__", [](const Op &a, const std::string &b) {
            return combine(parse_triplet(b), a).wrap();
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
    .def_readwrite("sym_ops", &GroupOps::sym_ops,
               "Symmetry operations (to be combined with centering vectors).")
    .def_readwrite("cen_ops", &GroupOps::cen_ops, "Centering vectors.")
    .def("change_basis", &GroupOps::change_basis, py::arg("cob"),
         "Applies the change-of-basis operator (in place).");

  py::class_<SpaceGroup>(sym, "SpaceGroup")
    .def(py::init<>())
    .def(py::init([](int n) { return get_spacegroup_by_number(n); }),
         py::arg("ccp4"), py::return_value_policy::reference)
    .def(py::init([](const std::string& s) {return get_spacegroup_by_name(s);}),
         py::arg("hm"), py::return_value_policy::reference)
    //.def(py::init(&find_spacegroup_by_name))
    .def("__repr__", [](const SpaceGroup &self) {
        return "<sym.SpaceGroup(\"" + self.xhm() + "\")>";
    })
    .def_readonly("number", &SpaceGroup::number, "number 1-230.")
    .def_readonly("ccp4", &SpaceGroup::ccp4, "ccp4 number")
    .def_readonly("hm", &SpaceGroup::hm, "Hermann–Mauguin name")
    .def_readonly("ext", &SpaceGroup::ext, "Extension (1, 2, H, R or none)")
    .def_readonly("qualifier", &SpaceGroup::qualifier, "e.g. 'cab'")
    .def_readonly("hall", &SpaceGroup::hall, "Hall symbol")
    .def("xhm", &SpaceGroup::xhm, "extended Hermann-Mauguin name")
    .def("operations", &SpaceGroup::operations, "Group of operations");

  sym.def("generators_from_hall", &generators_from_hall, py::arg("hall"),
          "Parse Hall notation.");
  sym.def("symops_from_hall", &symops_from_hall, py::arg("hall"),
          "Parse Hall notation.");
  sym.def("find_spacegroup_by_number", &find_spacegroup_by_number,
          py::arg("ccp4"), py::return_value_policy::reference,
          "Returns space-group of given number.");
  sym.def("find_spacegroup_by_name", &find_spacegroup_by_name, py::arg("hm"),
          py::return_value_policy::reference,
          "Returns space-group with given name.");
  sym.def("find_spacegroup_by_ops", &find_spacegroup_by_ops,
           py::arg("group_ops"), py::return_value_policy::reference,
          "Returns space-group with identical operations.");
}
