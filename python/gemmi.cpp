// Copyright 2017 Global Phasing Ltd.

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void add_symmetry(py::module& m); // sym.cpp
void add_mm(py::module& m); // mm.cpp
void init_cif(py::module& cif); // cif.cpp

PYBIND11_MODULE(gemmi, mg) {
  mg.doc() = "General MacroMolecular I/O";
  add_symmetry(mg);
  add_mm(mg);
  py::module cif = mg.def_submodule("cif", "CIF file format");
  init_cif(cif);
}

// vim:sw=2:ts=2:et
