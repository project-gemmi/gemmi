// Copyright 2018 Global Phasing Ltd.

#include "gemmi/chemcomp.hpp"    // for ChemComp
#include "gemmi/to_chemcomp.hpp" // for add_chemcomp_to_block

#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for find_shortest_path

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<Restraints::Bond>)
NB_MAKE_OPAQUE(std::vector<Restraints::Angle>)
NB_MAKE_OPAQUE(std::vector<Restraints::Torsion>)
NB_MAKE_OPAQUE(std::vector<Restraints::Chirality>)
NB_MAKE_OPAQUE(std::vector<Restraints::Plane>)
NB_MAKE_OPAQUE(std::vector<ChemComp::Atom>)

void add_chemcomp(nb::module_& m) {
  nb::class_<ChemComp> chemcomp(m, "ChemComp");
  nb::class_<ChemComp::Atom> chemcompatom(chemcomp, "Atom");

  nb::class_<Restraints> restraints(m, "Restraints");
  nb::class_<Restraints::Bond> restraintsbond(restraints, "Bond");
  nb::class_<Restraints::Angle> restraintsangle(restraints, "Angle");
  nb::class_<Restraints::Torsion> restraintstorsion(restraints, "Torsion");
  nb::class_<Restraints::Chirality> restraintschirality(restraints, "Chirality");
  nb::class_<Restraints::Plane> restraintsplane(restraints, "Plane");

  nb::bind_vector<std::vector<Restraints::Bond>, rv_ri>(m, "RestraintsBonds");
  nb::bind_vector<std::vector<Restraints::Angle>, rv_ri>(m, "RestraintsAngles");
  nb::bind_vector<std::vector<Restraints::Torsion>, rv_ri>(m, "RestraintsTorsions");
  nb::bind_vector<std::vector<Restraints::Chirality>, rv_ri>(m, "RestraintsChirs");
  nb::bind_vector<std::vector<Restraints::Plane>, rv_ri>(m, "RestraintsPlanes");
  nb::bind_vector<std::vector<ChemComp::Atom>, rv_ri>(m, "ChemCompAtoms");

  nb::enum_<BondType>(m, "BondType")
    .value("Unspec", BondType::Unspec)
    .value("Single", BondType::Single)
    .value("Double", BondType::Double)
    .value("Triple", BondType::Triple)
    .value("Aromatic", BondType::Aromatic)
    .value("Deloc", BondType::Deloc)
    .value("Metal", BondType::Metal);

  nb::enum_<ChiralityType>(m, "ChiralityType")
    .value("Positive", ChiralityType::Positive)
    .value("Negative", ChiralityType::Negative)
    .value("Both", ChiralityType::Both);

  nb::enum_<Restraints::DistanceOf>(restraints, "DistanceOf")
    .value("ElectronCloud", Restraints::DistanceOf::ElectronCloud)
    .value("Nucleus", Restraints::DistanceOf::Nucleus);


  nb::class_<Restraints::AtomId>(restraints, "AtomId")
    .def("__init__", [](Restraints::AtomId* p, int comp, const std::string& atom) {
        new(p) Restraints::AtomId{comp, atom};
    })
    .def("__init__", [](Restraints::AtomId* p, const std::string& atom) {
        new(p) Restraints::AtomId{1, atom};
    })
    .def_rw("comp", &Restraints::AtomId::comp)
    .def_rw("atom", &Restraints::AtomId::atom)
    .def("get_from",
         (Atom* (Restraints::AtomId::*)(Residue&, Residue*, char, char) const)
         &Restraints::AtomId::get_from,
         nb::arg("res1"), nb::arg("res2"), nb::arg("altloc1"), nb::arg("altloc2"),
         nb::rv_policy::reference)
    .def("__repr__", [](const Restraints::AtomId& self) {
        return cat("<gemmi.Restraints.AtomId ", self.comp, ' ', self.atom, '>');
    });
  restraintsbond
    .def(nb::init<>())
    .def_rw("id1", &Restraints::Bond::id1)
    .def_rw("id2", &Restraints::Bond::id2)
    .def_rw("type", &Restraints::Bond::type)
    .def_rw("aromatic", &Restraints::Bond::aromatic)
    .def_rw("value", &Restraints::Bond::value)
    .def_rw("esd", &Restraints::Bond::esd)
    .def_rw("value_nucleus", &Restraints::Bond::value_nucleus)
    .def_rw("esd_nucleus", &Restraints::Bond::esd_nucleus)
    .def("lexicographic_str", &Restraints::Bond::lexicographic_str)
    .def("__repr__", [](const Restraints::Bond& self) {
        return "<gemmi.Restraints.Bond " + self.str() + ">";
    });
  restraintsangle
    .def(nb::init<>())
    .def_rw("id1", &Restraints::Angle::id1)
    .def_rw("id2", &Restraints::Angle::id2)
    .def_rw("id3", &Restraints::Angle::id3)
    .def_rw("value", &Restraints::Angle::value)
    .def_rw("esd", &Restraints::Angle::esd)
    .def("__repr__", [](const Restraints::Angle& self) {
        return "<gemmi.Restraints.Angle " + self.str() + ">";
    });
  restraintstorsion
    .def(nb::init<>())
    .def_rw("label", &Restraints::Torsion::label)
    .def_rw("id1", &Restraints::Torsion::id1)
    .def_rw("id2", &Restraints::Torsion::id2)
    .def_rw("id3", &Restraints::Torsion::id3)
    .def_rw("id4", &Restraints::Torsion::id4)
    .def_rw("value", &Restraints::Torsion::value)
    .def_rw("esd", &Restraints::Torsion::esd)
    .def_rw("period", &Restraints::Torsion::period)
    .def("__repr__", [](const Restraints::Torsion& self) {
        return "<gemmi.Restraints.Torsion " + self.str() + ">";
    });
  restraintschirality
    .def(nb::init<>())
    .def_rw("id_ctr", &Restraints::Chirality::id_ctr)
    .def_rw("id1", &Restraints::Chirality::id1)
    .def_rw("id2", &Restraints::Chirality::id2)
    .def_rw("id3", &Restraints::Chirality::id3)
    .def_rw("sign", &Restraints::Chirality::sign)
    .def("is_wrong", &Restraints::Chirality::is_wrong)
    .def("__repr__", [](const Restraints::Chirality& self) {
        return "<gemmi.Restraints.Chirality " + self.str() + ">";
    });
  restraintsplane
    .def(nb::init<>())
    .def_rw("label", &Restraints::Plane::label)
    .def_rw("ids", &Restraints::Plane::ids)
    .def_rw("esd", &Restraints::Plane::esd)
    .def("__repr__", [](const Restraints::Plane& self) {
        return "<gemmi.Restraints.Plane " + self.str() + ">";
    });
  restraints
    .def(nb::init<>())
    .def_rw("bonds", &Restraints::bonds)
    .def_rw("angles", &Restraints::angles)
    .def_rw("torsions", &Restraints::torsions)
    .def_rw("chirs", &Restraints::chirs)
    .def_rw("planes", &Restraints::planes)
    .def("empty", &Restraints::empty)
    .def("get_bond", &Restraints::get_bond,
         nb::rv_policy::reference_internal)
    .def("get_bond", [](Restraints& self,
                        const std::string& a1, const std::string& a2)
                                                      -> Restraints::Bond& {
            auto it = self.find_bond(a1, a2);
            if (it == self.bonds.end())
              fail("Bond restraint not found: " + a1 + "-" + a2);
            return *it;
         }, nb::rv_policy::reference_internal)
    .def("find_shortest_path", &Restraints::find_shortest_path)
    .def("chiral_abs_volume", &Restraints::chiral_abs_volume)
    ;

  nb::enum_<ChemComp::Group>(chemcomp, "Group")
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
  chemcompatom
    .def_rw("id", &ChemComp::Atom::id)
    .def_rw("el", &ChemComp::Atom::el)
    .def_rw("charge", &ChemComp::Atom::charge)
    .def_rw("chem_type", &ChemComp::Atom::chem_type)
    .def("is_hydrogen", &ChemComp::Atom::is_hydrogen)
    ;
  nb::class_<ChemComp::Aliasing>(chemcomp, "Aliasing")
    .def_ro("group", &ChemComp::Aliasing::group)
    .def("name_from_alias", &ChemComp::Aliasing::name_from_alias)
    ;
  chemcomp
    .def_rw("name", &ChemComp::name)
    .def_rw("group", &ChemComp::group)
    .def_ro("atoms", &ChemComp::atoms)
    .def_ro("rt", &ChemComp::rt)
    .def_static("read_group", &ChemComp::read_group)
    .def_static("group_str", &ChemComp::group_str)
    .def("set_group", &ChemComp::set_group)
    .def("get_atom", &ChemComp::get_atom)
    .def("find_atom", [](ChemComp& self, const std::string& atom_id) {
        auto it = self.find_atom(atom_id);
        return it != self.atoms.end() ? &*it : nullptr;
    }, nb::rv_policy::reference_internal)
    .def("remove_hydrogens", &ChemComp::remove_hydrogens)
    ;
  m.def("make_chemcomp_from_block", &make_chemcomp_from_block);
  m.def("add_chemcomp_to_block", &add_chemcomp_to_block);
}
