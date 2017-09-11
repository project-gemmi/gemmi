// Copyright 2017 Global Phasing Ltd.

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_cif(py::module& cif); // cif.cpp
void init_sym(py::module& sym); // sym.cpp
void init_mol(py::module& mol); // mol.cpp

PYBIND11_MODULE(gemmi, mg) {
  mg.doc() = "General MacroMolecular I/O";
  py::module cif = mg.def_submodule("cif", "CIF file format");
  init_cif(cif);
  py::module sym = mg.def_submodule("sym", "Symmetry, space groups");
  init_sym(sym);
  py::module mol = mg.def_submodule("mol", "MacroMolecular models");
  init_mol(mol);
}

// vim:sw=2:ts=2:et
