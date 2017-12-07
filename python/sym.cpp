// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

void add_symmetry(py::module& m) {
  py::class_<Op>(m, "Op")
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
    .def("combine", &Op::combine, py::arg("b"),
         "Combine two symmetry operations.")
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
    .def("__mul__", [](const Op &a, const Op &b) { return a * b; },
         py::is_operator())
    .def("__mul__", [](const Op &a, const std::string &b) {
            return a * parse_triplet(b);
         }, py::is_operator())
    .def("__rmul__", [](const Op &a, const std::string &b) {
            return parse_triplet(b) * a;
         }, py::is_operator())
    .def("__eq__", [](const Op &a, const Op &b) { return a == b; },
         py::is_operator())
    .def("__eq__", [](const Op &a, const std::string& b) {
            return a == parse_triplet(b);
         }, py::is_operator())
#if PY_MAJOR_VERSION < 3  // in Py3 != is inferred from ==
    .def("__ne__", [](const Op &a, const Op &b) { return a != b; },
         py::is_operator())
    .def("__ne__", [](const Op &a, const std::string& b) {
            return a != parse_triplet(b);
         }, py::is_operator())
#endif
    .def("__hash__", [](const Op &self) { return std::hash<Op>()(self); })
    .def("__repr__", [](const Op &self) {
        return "<gemmi.Op(\"" + self.triplet() + "\")>";
    });

  m.def("parse_triplet", &parse_triplet, py::arg("triplet"),
        "Parse coordinate triplet into gemmi.Op.");
  m.def("parse_triplet_part", &parse_triplet_part, py::arg("s"),
        "Parse one of the three parts of a triplet.");
  m.def("make_triplet_part", &make_triplet_part,
        py::arg("x"), py::arg("y"), py::arg("z"), py::arg("w"),
        "Make one of the three parts of a triplet.");

  py::class_<GroupOps>(m, "GroupOps")
    .def(py::init<>())
    .def("__iter__", [](const GroupOps& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__eq__", [](const GroupOps &a, const GroupOps &b) {
            return a.is_same_as(b);
    }, py::is_operator())
#if PY_MAJOR_VERSION < 3  // in Py3 != is inferred from ==
    .def("__ne__", [](const GroupOps &a, const GroupOps &b) {
            return !a.is_same_as(b);
    }, py::is_operator())
#endif
    .def("__len__", [](const GroupOps& g) {
            return g.sym_ops.size() * g.cen_ops.size();
    })
    .def_readwrite("sym_ops", &GroupOps::sym_ops,
               "Symmetry operations (to be combined with centering vectors).")
    .def_readwrite("cen_ops", &GroupOps::cen_ops, "Centering vectors.")
    .def("find_grid_factors", &GroupOps::find_grid_factors,
         "Minimal multiplicity for real-space grid (e.g. 1,1,6 for P61).")
    .def("change_basis", &GroupOps::change_basis, py::arg("cob"),
         "Applies the change-of-basis operator (in place).");

  py::class_<SpaceGroup>(m, "SpaceGroup")
    .def(py::init<>())
    .def(py::init([](int n) { return get_spacegroup_by_number(n); }),
         py::arg("ccp4"), py::return_value_policy::reference)
    .def(py::init([](const std::string& s) {return get_spacegroup_by_name(s);}),
         py::arg("hm"), py::return_value_policy::reference)
    //.def(py::init(&find_spacegroup_by_name))
    .def("__repr__", [](const SpaceGroup &self) {
        return "<gemmi.SpaceGroup(\"" + self.xhm() + "\")>";
    })
    .def_readonly("number", &SpaceGroup::number, "number 1-230.")
    .def_readonly("ccp4", &SpaceGroup::ccp4, "ccp4 number")
    .def_readonly("hm", &SpaceGroup::hm, "Hermann-Mauguin name")
    .def_readonly("ext", &SpaceGroup::ext, "Extension (1, 2, H, R or none)")
    .def_readonly("qualifier", &SpaceGroup::qualifier, "e.g. 'cab'")
    .def_readonly("hall", &SpaceGroup::hall, "Hall symbol")
    .def("xhm", &SpaceGroup::xhm, "extended Hermann-Mauguin name")
    .def("short_name", &SpaceGroup::short_name,
         "H-M name w/o spaces and with 1's removed in '1 ... 1'.")
    .def("operations", &SpaceGroup::operations, "Group of operations");

  m.def("spacegroup_table", []() {
            return py::make_iterator(spacegroup_tables::main);
        }, py::return_value_policy::reference);
  m.def("generators_from_hall", &generators_from_hall, py::arg("hall"),
        "Parse Hall notation.");
  m.def("symops_from_hall", &symops_from_hall, py::arg("hall"),
        "Parse Hall notation.");
  m.def("find_spacegroup_by_number", &find_spacegroup_by_number,
        py::arg("ccp4"), py::return_value_policy::reference,
        "Returns space-group of given number.");
  m.def("find_spacegroup_by_name", &find_spacegroup_by_name, py::arg("hm"),
        py::return_value_policy::reference,
        "Returns space-group with given name.");
  m.def("find_spacegroup_by_ops", &find_spacegroup_by_ops,
        py::arg("group_ops"), py::return_value_policy::reference,
        "Returns space-group with identical operations.");
}
