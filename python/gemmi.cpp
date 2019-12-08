// Copyright 2017 Global Phasing Ltd.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gemmi/version.hpp"
#include "gemmi/dirwalk.hpp"
#include "gemmi/fileutil.hpp"  // for expand_if_pdb_code
#include "gemmi/elem.hpp"
#include "gemmi/it92.hpp"
#include "gemmi/smodel.hpp"    // for AtomicStructure

namespace py = pybind11;

void add_symmetry(py::module& m); // sym.cpp
void add_grid(py::module& m); // grid.cpp
void add_unitcell(py::module& m); // unitcell.cpp
void add_hkl(py::module& m); // hkl.cpp
void add_mol(py::module& m); // mol.cpp
void add_cif(py::module& cif); // cif.cpp
void add_read_structure(py::module& m); // read.cpp
void add_cif_read(py::module& cif); // read.cpp
void add_monlib(py::module& cif); // monlib.cpp

void add_misc(py::module& m) {
  using gemmi::Element;
  using gemmi::AtomicStructure;
  using IT92 = gemmi::IT92<double>;
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def("__eq__",
         [](const Element &a, const Element &b) { return a.elem == b.elem; },
         py::is_operator())
#if PY_MAJOR_VERSION < 3  // in Py3 != is inferred from ==
    .def("__ne__",
         [](const Element &a, const Element &b) { return a.elem != b.elem; },
         py::is_operator())
#endif
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("covalent_r", &Element::covalent_r)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def_property_readonly("it92", [](const Element& self) {
        return IT92::get_ptr(self.elem);
    }, py::return_value_policy::reference_internal)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });

  py::class_<AtomicStructure> atomic_structure(m, "AtomicStructure");
  py::class_<AtomicStructure::Site>(atomic_structure, "Site")
    .def_readonly("label", &AtomicStructure::Site::label)
    .def_readonly("type_symbol", &AtomicStructure::Site::type_symbol)
    .def_readonly("fract", &AtomicStructure::Site::fract)
    .def_readonly("occ", &AtomicStructure::Site::occ)
    .def_readonly("u_iso", &AtomicStructure::Site::u_iso)
    .def_readonly("element", &AtomicStructure::Site::element)
    .def_readonly("charge", &AtomicStructure::Site::charge)
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

  py::class_<IT92::Coef>(m, "IT92Coef")
    .def("calculate_sf", &IT92::Coef::calculate_sf, py::arg("stol2"))
    .def("calculate_density", &IT92::Coef::calculate_density,
         py::arg("r2"), py::arg("B"))
    ;

  py::class_<gemmi::CifWalk>(m, "CifWalk")
    .def(py::init<const char*>())
    .def("__iter__", [](gemmi::CifWalk& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
  py::class_<gemmi::CoorFileWalk>(m, "CoorFileWalk")
    .def(py::init<const char*>())
    .def("__iter__", [](gemmi::CoorFileWalk& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
  m.def("is_pdb_code", &gemmi::is_pdb_code);
  m.def("expand_pdb_code_to_path", &gemmi::expand_pdb_code_to_path);
  m.def("expand_if_pdb_code", &gemmi::expand_if_pdb_code,
        py::arg("code"), py::arg("filetype")='M');
}

PYBIND11_MODULE(gemmi, mg) {
  mg.doc() = "General MacroMolecular I/O";
  mg.attr("__version__") = GEMMI_VERSION;
  add_symmetry(mg);
  add_grid(mg);
  add_unitcell(mg);
  add_hkl(mg);
  add_misc(mg);
  add_mol(mg);
  add_monlib(mg);
  add_read_structure(mg);
  py::module cif = mg.def_submodule("cif", "CIF file format");
  add_cif(cif);
  add_cif_read(cif);
}
