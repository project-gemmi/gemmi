// Copyright 2017 Global Phasing Ltd.

#include "cif.hh"
#include "cifgz.hh"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(gemmi) {
  using namespace gemmi::cif;

  py::module m("gemmi", "General MacroMolecular I/O");

  py::class_<Document>(m, "Document")
    .def(py::init<>())
    .def(py::init<const std::string &>())
    .def("__len__", [](const Document& d) { return d.blocks.size(); })
    .def("__iter__", [](const Document& d) {
        return py::make_iterator(d.blocks.begin(), d.blocks.end());
    }, py::keep_alive<0, 1>())
    .def("parse_file", &Document::parse_file, "Parse file")
    .def("sole_block", &Document::sole_block,
         py::return_value_policy::reference,
         "Returns the only block if there is exactly one");

  py::class_<Block>(m, "Block")
    .def(py::init<>())
    .def_readonly("name", &Block::name)
    .def("find_value", &Block::find_value, py::return_value_policy::reference)
    .def("find_loop", &Block::find_loop, py::return_value_policy::reference)
    ;

  py::class_<Loop>(m, "Loop")
    .def(py::init<>())
    .def("width", &Loop::width)
    .def("length", &Loop::length)
    .def("__iter__", [](const Loop& self) {
        return py::make_iterator(self.tags.begin(), self.tags.end());
    }, py::keep_alive<0, 1>())
    .def("val", &Loop::val)
    ;

  py::class_<LoopColumn>(m, "LoopColumn")
    .def(py::init<>())
    .def_readonly("loop", &LoopColumn::loop)
    .def_readwrite("col", &LoopColumn::col)
    .def("__iter__", [](const LoopColumn& self) {
        return py::make_iterator(self.begin(), self.end());
    }, py::keep_alive<0, 1>())
    .def("__bool__", [](const LoopColumn &self) -> bool { return self.loop; });

  m.def("read_any", &read_any, "Reads normal or gzipped cif file.");

  //PYBIND11_MAKE_OPAQUE(std::vector<int>);

  return m.ptr();
}

// vim:sw=2:ts=2:et
