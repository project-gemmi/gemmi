// Copyright 2017 Global Phasing Ltd.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gemmi/version.hpp"
#include "gemmi/dirwalk.hpp"
#include "gemmi/fileutil.hpp"  // for expand_if_pdb_code
#include "gemmi/elem.hpp"
#include "gemmi/it92.hpp"
#include "gemmi/small.hpp"    // for SmallStructure

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
  using gemmi::SmallStructure;
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

  py::class_<SmallStructure> small_structure(m, "SmallStructure");
  py::class_<SmallStructure::Site>(small_structure, "Site")
    .def_readonly("label", &SmallStructure::Site::label)
    .def_readonly("type_symbol", &SmallStructure::Site::type_symbol)
    .def_readonly("fract", &SmallStructure::Site::fract)
    .def_readonly("occ", &SmallStructure::Site::occ)
    .def_readonly("u_iso", &SmallStructure::Site::u_iso)
    .def_readonly("element", &SmallStructure::Site::element)
    .def_readonly("charge", &SmallStructure::Site::charge)
    .def("__repr__", [](const SmallStructure::Site& self) {
        return "<gemmi.SmallStructure.Site " + self.label + ">";
    });

  using AtomType = SmallStructure::AtomType;
  py::class_<AtomType>(small_structure, "AtomType")
    .def_readonly("symbol", &AtomType::symbol)
    .def_readonly("element", &AtomType::element)
    .def_readwrite("dispersion_real", &AtomType::dispersion_real)
    .def_readwrite("dispersion_imag", &AtomType::dispersion_imag)
    .def("__repr__", [](const AtomType& self) {
        return "<gemmi.SmallStructure.AtomType " + self.symbol + ">";
    });

  small_structure
    .def(py::init<>())
    .def_readwrite("name", &SmallStructure::name)
    .def_readwrite("cell", &SmallStructure::cell)
    .def_readonly("spacegroup_hm", &SmallStructure::spacegroup_hm)
    .def_readonly("sites", &SmallStructure::sites)
    .def_readonly("atom_types", &SmallStructure::atom_types)
    .def("get_all_unit_cell_sites", &SmallStructure::get_all_unit_cell_sites)
    .def("__repr__", [](const SmallStructure& self) {
        return "<gemmi.SmallStructure: " + std::string(self.name) + ">";
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
