// Copyright 2017 Global Phasing Ltd.

#include "common.h"
#include <pybind11/stl.h>
#include "gemmi/version.hpp"
#include "gemmi/dirwalk.hpp"
#include "gemmi/fileutil.hpp"  // for expand_if_pdb_code
#include "gemmi/small.hpp"    // for SmallStructure
#include "gemmi/interop.hpp"  // for atom_to_site, mx_to_sx_structure

namespace py = pybind11;

void add_misc(py::module& m) {
  using gemmi::SmallStructure;
  py::class_<SmallStructure> small_structure(m, "SmallStructure");
  py::class_<SmallStructure::Site>(small_structure, "Site")
    .def(py::init<>())
    .def(py::init(&gemmi::atom_to_site))
    .def_readonly("label", &SmallStructure::Site::label)
    .def_readonly("type_symbol", &SmallStructure::Site::type_symbol)
    .def_readonly("fract", &SmallStructure::Site::fract)
    .def_readonly("occ", &SmallStructure::Site::occ)
    .def_readonly("u_iso", &SmallStructure::Site::u_iso)
    .def_readonly("element", &SmallStructure::Site::element)
    .def_readonly("charge", &SmallStructure::Site::charge)
    .def_readonly("disorder_group", &SmallStructure::Site::disorder_group)
    .def_readonly("aniso", &SmallStructure::Site::aniso)
    .def("orth", &SmallStructure::Site::orth)
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
    .def("remove_hydrogens", &SmallStructure::remove_hydrogens)
    .def("__repr__", [](const SmallStructure& self) {
        return "<gemmi.SmallStructure: " + std::string(self.name) + ">";
    });
  m.def("mx_to_sx_structure", &gemmi::mx_to_sx_structure,
        py::arg("st"), py::arg("n")=0);

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
  m.attr("hc") = py::float_(gemmi::hc());
}

PYBIND11_MODULE(gemmi, mg) {
  mg.doc() = "Python bindings to GEMMI - a library used in macromolecular\n"
             "crystallography and related fields";
  mg.attr("__version__") = GEMMI_VERSION;
  py::module cif = mg.def_submodule("cif", "CIF file format");
  add_cif(cif);
  add_symmetry(mg);
  add_unitcell(mg);
  add_elem(mg);
  add_meta(mg);
  add_mol(mg);
  add_misc(mg);
  add_grid(mg);
  add_recgrid(mg);
  add_ccp4(mg);
  add_sf(mg);
  add_cif_read(cif);
  add_hkl(mg);
  add_monlib(mg);
  add_alignment(mg);
  add_select(mg);
  add_search(mg);
  add_read_structure(mg);
}
