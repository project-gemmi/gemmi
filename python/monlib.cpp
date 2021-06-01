// Copyright 2018 Global Phasing Ltd.

#include "gemmi/chemcomp.hpp"
#include "gemmi/monlib.hpp"
#include "gemmi/tostr.hpp"
#include "gemmi/read_cif.hpp"  // for read_cif_gz

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
using monomers_type = std::map<std::string, ChemComp>;
using links_type = std::map<std::string, ChemLink>;
using modifications_type = std::map<std::string, ChemMod>;
PYBIND11_MAKE_OPAQUE(monomers_type)
PYBIND11_MAKE_OPAQUE(links_type)
PYBIND11_MAKE_OPAQUE(modifications_type)

void add_monlib(py::module& m) {
  py::class_<ChemMod> chemmod(m, "ChemMod");
  py::class_<ChemLink> chemlink(m, "ChemLink");
  py::class_<ChemLink::Side> chemlinkside(chemlink, "Side");
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
  py::bind_map<monomers_type>(m, "ChemCompMap");
  py::bind_map<links_type>(m, "ChemLinkMap");
  py::bind_map<modifications_type>(m, "ChemModMap");

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
        return tostr("<gemmi.Restraints.AtomId ", self.comp, ' ', self.atom, '>');
    });
  restraintsbond
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
    .def_readwrite("id1", &Restraints::Angle::id1)
    .def_readwrite("id2", &Restraints::Angle::id2)
    .def_readwrite("id3", &Restraints::Angle::id3)
    .def_readwrite("value", &Restraints::Angle::value)
    .def_readwrite("esd", &Restraints::Angle::esd)
    .def("__repr__", [](const Restraints::Angle& self) {
        return "<gemmi.Restraints.Angle " + self.str() + ">";
    });
  restraintstorsion
    .def_readwrite("id1", &Restraints::Torsion::id1)
    .def_readwrite("id2", &Restraints::Torsion::id2)
    .def_readwrite("id3", &Restraints::Torsion::id3)
    .def_readwrite("id4", &Restraints::Torsion::id3)
    .def_readwrite("value", &Restraints::Torsion::value)
    .def_readwrite("esd", &Restraints::Torsion::esd)
    .def_readwrite("period", &Restraints::Torsion::period)
    .def("__repr__", [](const Restraints::Torsion& self) {
        return "<gemmi.Restraints.Torsion " + self.str() + ">";
    });
  restraintschirality
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
    .def_readwrite("label", &Restraints::Plane::label)
    .def_readwrite("ids", &Restraints::Plane::ids)
    .def_readwrite("esd", &Restraints::Plane::esd)
    .def("__repr__", [](const Restraints::Plane& self) {
        return "<gemmi.Restraints.Plane " + self.str() + ">";
    });
  restraints
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
    ;

  chemcompatom
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

  chemlink
    .def_readwrite("id", &ChemLink::id)
    .def_readwrite("name", &ChemLink::name)
    .def_readwrite("side1", &ChemLink::side1)
    .def_readwrite("side2", &ChemLink::side2)
    .def_readwrite("rt", &ChemLink::rt)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemLink " + self.id + ">";
    });

  py::enum_<ChemLink::Group>(chemlinkside, "Group")
      .value("Peptide",  ChemLink::Group::Peptide)
      .value("PPeptide", ChemLink::Group::PPeptide)
      .value("MPeptide", ChemLink::Group::MPeptide)
      .value("Pyranose", ChemLink::Group::Pyranose)
      .value("DnaRna",   ChemLink::Group::DnaRna)
      .value("Null",     ChemLink::Group::Null);

  chemlinkside
    .def_readwrite("comp", &ChemLink::Side::comp)
    .def_readwrite("mod", &ChemLink::Side::mod)
    .def_readwrite("group", &ChemLink::Side::group)
    .def("__repr__", [](const ChemLink::Side& self) {
        return "<gemmi.ChemLink.Side " + self.comp + "/" +
               ChemLink::group_str(self.group) + ">";
    });

  chemmod
    .def_readwrite("id", &ChemMod::id)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  py::class_<MonLib>(m, "MonLib")
    .def(py::init<>())
    .def_readonly("monomers", &MonLib::monomers)
    .def_readonly("links", &MonLib::links)
    .def_readonly("modifications", &MonLib::modifications)
    .def("find_link", &MonLib::find_link, py::arg("link_id"),
         py::return_value_policy::reference_internal)
    .def("match_link", &MonLib::match_link,
         py::arg("comp1"), py::arg("atom1"),
         py::arg("comp2"), py::arg("atom2"),
         py::return_value_policy::reference_internal)
    .def("add_monomer_if_present", &MonLib::add_monomer_if_present)
    .def("add_monomers_if_present", &MonLib::add_monomers_if_present)
    .def("insert_chemlinks", [](MonLib &self, const cif::Document &doc) {
        insert_chemlinks(doc, self.links);
    })
    .def("insert_chemmods", [](MonLib &self, const cif::Document &doc) {
        insert_chemmods(doc, self.modifications);
    })
    .def("insert_comp_list", [](MonLib &self, const cif::Document &doc) {
        insert_comp_list(doc, self.residue_infos);
    })
    .def("__repr__", [](const MonLib& self) {
        return "<gemmi.MonLib with " +
               std::to_string(self.monomers.size()) + " monomers, " +
               std::to_string(self.links.size()) + " links, " +
               std::to_string(self.modifications.size()) + " modifications>";
    });

  m.def("read_monomer_lib", [](const std::string& monomer_dir,
                               const std::vector<std::string>& resnames) {
    return read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz);
  });
  m.def("read_monomer_cif", [](const std::string& path) {
    return read_monomer_cif(path, gemmi::read_cif_gz);
  });

  py::class_<BondIndex>(m, "BondIndex")
    .def(py::init<const Model&>(), py::keep_alive<1, 2>())
    .def("add_link", &BondIndex::add_link)
    .def("add_monomer_bonds", &BondIndex::add_monomer_bonds)
    .def("are_linked", &BondIndex::are_linked)
    .def("graph_distance", &BondIndex::graph_distance,
         py::arg("a"), py::arg("b"), py::arg("same_index"),
         py::arg("max_distance")=4)
    ;
}
