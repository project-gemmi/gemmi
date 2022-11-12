// Copyright 2018 Global Phasing Ltd.

#include "gemmi/chemcomp.hpp"    // for ChemComp
#include "gemmi/to_chemcomp.hpp" // for add_chemcomp_to_block

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Angle>)
PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Torsion>)
PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Chirality>)
PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Plane>)
PYBIND11_MAKE_OPAQUE(std::vector<ChemComp::Atom>)

void add_chemcomp(py::module& m) {
  py::class_<ChemComp> chemcomp(m, "ChemComp");
  py::class_<ChemComp::Atom> chemcompatom(chemcomp, "Atom");

  py::class_<Restraints> restraints(m, "Restraints");
  py::class_<Restraints::Bond> restraintsbond(restraints, "Bond");
  py::class_<Restraints::Angle> restraintsangle(restraints, "Angle");
  py::class_<Restraints::Torsion> restraintstorsion(restraints, "Torsion");
  py::class_<Restraints::Chirality> restraintschirality(restraints, "Chirality");
  py::class_<Restraints::Plane> restraintsplane(restraints, "Plane");

  py::bind_vector<std::vector<Restraints::Bond>>(m, "RestraintsBonds");
  py::bind_vector<std::vector<Restraints::Angle>>(m, "RestraintsAngles");
  py::bind_vector<std::vector<Restraints::Torsion>>(m, "RestraintsTorsions");
  py::bind_vector<std::vector<Restraints::Chirality>>(m, "RestraintsChirs");
  py::bind_vector<std::vector<Restraints::Plane>>(m, "RestraintsPlanes");
  py::bind_vector<std::vector<ChemComp::Atom>>(m, "ChemCompAtoms");

  py::enum_<BondType>(m, "BondType")
    .value("Unspec", BondType::Unspec)
    .value("Single", BondType::Single)
    .value("Double", BondType::Double)
    .value("Triple", BondType::Triple)
    .value("Aromatic", BondType::Aromatic)
    .value("Deloc", BondType::Deloc)
    .value("Metal", BondType::Metal);

  py::enum_<ChiralityType>(m, "ChiralityType")
    .value("Positive", ChiralityType::Positive)
    .value("Negative", ChiralityType::Negative)
    .value("Both", ChiralityType::Both);

  py::enum_<Restraints::DistanceOf>(restraints, "DistanceOf")
    .value("ElectronCloud", Restraints::DistanceOf::ElectronCloud)
    .value("Nucleus", Restraints::DistanceOf::Nucleus);


  py::class_<Restraints::AtomId>(restraints, "AtomId")
    .def(py::init([](int comp, const std::string& atom) {
          return new Restraints::AtomId{comp, atom};
    }))
    .def(py::init([](const std::string& atom) {
          return new Restraints::AtomId{1, atom};
    }))
    .def_readwrite("comp", &Restraints::AtomId::comp)
    .def_readwrite("atom", &Restraints::AtomId::atom)
    .def("get_from",
         (Atom* (Restraints::AtomId::*)(Residue&, Residue*, char) const)
         &Restraints::AtomId::get_from,
         py::arg("res1"), py::arg("res2"), py::arg("altloc"),
         py::return_value_policy::reference)
    .def("__repr__", [](const Restraints::AtomId& self) {
        return cat("<gemmi.Restraints.AtomId ", self.comp, ' ', self.atom, '>');
    });
  restraintsbond
    .def(py::init<>())
    .def_readwrite("id1", &Restraints::Bond::id1)
    .def_readwrite("id2", &Restraints::Bond::id2)
    .def_readwrite("type", &Restraints::Bond::type)
    .def_readwrite("aromatic", &Restraints::Bond::aromatic)
    .def_readwrite("value", &Restraints::Bond::value)
    .def_readwrite("esd", &Restraints::Bond::esd)
    .def_readwrite("value_nucleus", &Restraints::Bond::value_nucleus)
    .def_readwrite("esd_nucleus", &Restraints::Bond::esd_nucleus)
    .def("lexicographic_str", &Restraints::Bond::lexicographic_str)
    .def("__repr__", [](const Restraints::Bond& self) {
        return "<gemmi.Restraints.Bond " + self.str() + ">";
    });
  restraintsangle
    .def(py::init<>())
    .def_readwrite("id1", &Restraints::Angle::id1)
    .def_readwrite("id2", &Restraints::Angle::id2)
    .def_readwrite("id3", &Restraints::Angle::id3)
    .def_readwrite("value", &Restraints::Angle::value)
    .def_readwrite("esd", &Restraints::Angle::esd)
    .def("__repr__", [](const Restraints::Angle& self) {
        return "<gemmi.Restraints.Angle " + self.str() + ">";
    });
  restraintstorsion
    .def(py::init<>())
    .def_readwrite("label", &Restraints::Torsion::label)
    .def_readwrite("id1", &Restraints::Torsion::id1)
    .def_readwrite("id2", &Restraints::Torsion::id2)
    .def_readwrite("id3", &Restraints::Torsion::id3)
    .def_readwrite("id4", &Restraints::Torsion::id4)
    .def_readwrite("value", &Restraints::Torsion::value)
    .def_readwrite("esd", &Restraints::Torsion::esd)
    .def_readwrite("period", &Restraints::Torsion::period)
    .def("__repr__", [](const Restraints::Torsion& self) {
        return "<gemmi.Restraints.Torsion " + self.str() + ">";
    });
  restraintschirality
    .def(py::init<>())
    .def_readwrite("id_ctr", &Restraints::Chirality::id_ctr)
    .def_readwrite("id1", &Restraints::Chirality::id1)
    .def_readwrite("id2", &Restraints::Chirality::id2)
    .def_readwrite("id3", &Restraints::Chirality::id3)
    .def_readwrite("sign", &Restraints::Chirality::sign)
    .def("is_wrong", &Restraints::Chirality::is_wrong)
    .def("__repr__", [](const Restraints::Chirality& self) {
        return "<gemmi.Restraints.Chirality " + self.str() + ">";
    });
  restraintsplane
    .def(py::init<>())
    .def_readwrite("label", &Restraints::Plane::label)
    .def_readwrite("ids", &Restraints::Plane::ids)
    .def_readwrite("esd", &Restraints::Plane::esd)
    .def("__repr__", [](const Restraints::Plane& self) {
        return "<gemmi.Restraints.Plane " + self.str() + ">";
    });
  restraints
    .def(py::init<>())
    .def_readwrite("bonds", &Restraints::bonds)
    .def_readwrite("angles", &Restraints::angles)
    .def_readwrite("torsions", &Restraints::torsions)
    .def_readwrite("chirs", &Restraints::chirs)
    .def_readwrite("planes", &Restraints::planes)
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
    .def("find_shortest_path", &Restraints::find_shortest_path)
    .def("chiral_abs_volume", &Restraints::chiral_abs_volume)
    ;

  chemcompatom
    .def_readonly("id", &ChemComp::Atom::id)
    .def_readonly("el", &ChemComp::Atom::el)
    .def_readonly("charge", &ChemComp::Atom::charge)
    .def_readonly("chem_type", &ChemComp::Atom::chem_type)
    .def("is_hydrogen", &ChemComp::Atom::is_hydrogen)
    ;
  py::enum_<ChemComp::Group>(chemcomp, "Group")
      .value("Peptide",    ChemComp::Group::Peptide)
      .value("PPeptide",   ChemComp::Group::PPeptide)
      .value("MPeptide",   ChemComp::Group::MPeptide)
      .value("Dna",        ChemComp::Group::Dna)
      .value("Rna",        ChemComp::Group::Rna)
      .value("DnaRna",     ChemComp::Group::DnaRna)
      .value("Pyranose",   ChemComp::Group::Pyranose)
      .value("Ketopyranose", ChemComp::Group::Ketopyranose)
      .value("Furanose",   ChemComp::Group::Furanose)
      .value("NonPolymer", ChemComp::Group::NonPolymer)
      .value("Null",       ChemComp::Group::Null);
  chemcomp
    .def_readwrite("name", &ChemComp::name)
    .def_readwrite("group", &ChemComp::group)
    .def_readonly("atoms", &ChemComp::atoms)
    .def_readonly("rt", &ChemComp::rt)
    .def_static("group_str", &ChemComp::group_str)
    .def("set_group", &ChemComp::set_group)
    .def("get_atom", &ChemComp::get_atom)
    .def("find_atom", [](ChemComp& self, const std::string& atom_id) {
        auto it = self.find_atom(atom_id);
        return it != self.atoms.end() ? &*it : nullptr;
    }, py::return_value_policy::reference_internal)
    .def("remove_hydrogens", &ChemComp::remove_hydrogens)
    ;
  m.def("make_chemcomp_from_block", &make_chemcomp_from_block);
  m.def("add_chemcomp_to_block", &add_chemcomp_to_block);
}
