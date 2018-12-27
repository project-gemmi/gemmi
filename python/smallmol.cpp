// Copyright 2018 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/smcif.hpp"
#include "gemmi/chemcomp.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<ChemComp::Atom>)

void add_smcif(py::module& m) {
  py::bind_vector<std::vector<Restraints::Bond>>(m, "VectorRestraintsBond");
  py::bind_vector<std::vector<ChemComp::Atom>>(m, "VectorChemCompAtom");

  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });

  py::class_<AtomicStructure> atomic_structure(m, "AtomicStructure");
  py::class_<AtomicStructure::Site>(atomic_structure, "Site")
    .def_readonly("label", &AtomicStructure::Site::label)
    .def_readonly("type_symbol", &AtomicStructure::Site::type_symbol)
    .def_readonly("fract", &AtomicStructure::Site::fract)
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
}

void add_chemcomp(py::module& m) {
  py::enum_<BondType>(m, "BondType")
    .value("Unspec", BondType::Unspec)
    .value("Single", BondType::Single)
    .value("Double", BondType::Double)
    .value("Triple", BondType::Triple)
    .value("Aromatic", BondType::Aromatic)
    .value("Deloc", BondType::Deloc)
    .value("Metal", BondType::Metal);

  py::class_<Restraints> restraints(m, "Restraints");
  py::class_<Restraints::AtomId>(restraints, "AtomId")
    .def_readwrite("comp", &Restraints::AtomId::comp)
    .def_readwrite("atom", &Restraints::AtomId::atom)
    .def("get_from",
         (Atom* (Restraints::AtomId::*)(Residue&, Residue*, char) const)
         &Restraints::AtomId::get_from,
         py::arg("res1"), py::arg("res2"), py::arg("altloc"),
         py::return_value_policy::reference)
    ;
  py::class_<Restraints::Bond>(restraints, "Bond")
    .def_readwrite("id1", &Restraints::Bond::id1)
    .def_readwrite("id2", &Restraints::Bond::id2)
    .def_readwrite("type", &Restraints::Bond::type)
    .def_readwrite("aromatic", &Restraints::Bond::aromatic)
    .def_readwrite("value", &Restraints::Bond::value)
    .def_readwrite("esd", &Restraints::Bond::esd)
    .def("lexicographic_str", &Restraints::Bond::lexicographic_str)
    ;
  restraints
    .def_readwrite("bonds", &Restraints::bonds)
    .def("empty", &Restraints::empty)
    .def("get_bond", &Restraints::get_bond,
         py::return_value_policy::reference_internal)
    .def("get_bond", [](Restraints& self,
                        const std::string& a1, const std::string& a2)
                                                      -> Restraints::Bond& {
            auto it = self.find_bond(a1, a2);
            if (it == self.bonds.end())
              fail("Bond restraint not found: " + a1 + "-" + a2);
            return *it;
         }, py::return_value_policy::reference_internal)
    ;

  py::class_<ChemComp> chemcomp(m, "ChemComp");
  py::class_<ChemComp::Atom>(chemcomp, "Atom")
    .def_readonly("id", &ChemComp::Atom::id)
    .def_readonly("el", &ChemComp::Atom::el)
    .def_readonly("charge", &ChemComp::Atom::charge)
    .def_readonly("chem_type", &ChemComp::Atom::chem_type)
    .def("is_hydrogen", &ChemComp::Atom::is_hydrogen)
    ;
  chemcomp
    .def_readonly("name", &ChemComp::name)
    .def_readonly("group", &ChemComp::group)
    .def_readonly("atoms", &ChemComp::atoms)
    .def_readonly("rt", &ChemComp::rt)
    .def("get_atom", &ChemComp::get_atom)
    .def("remove_hydrogens", &ChemComp::remove_hydrogens)
    ;
  m.def("make_chemcomp_from_block", &make_chemcomp_from_block);
}
