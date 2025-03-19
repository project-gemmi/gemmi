// Copyright 2018 Global Phasing Ltd.

#include "gemmi/monlib.hpp"

#include "common.h"
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>  // for MonLib::match_link
#include <nanobind/stl/vector.h>

using namespace gemmi;

using monomers_type = std::map<std::string, ChemComp>;
using links_type = std::map<std::string, ChemLink>;
using modifications_type = std::map<std::string, ChemMod>;
NB_MAKE_OPAQUE(monomers_type)
NB_MAKE_OPAQUE(links_type)
NB_MAKE_OPAQUE(modifications_type)

void add_monlib(nb::module_& m) {
  nb::class_<ChemMod> chemmod(m, "ChemMod");
  nb::class_<ChemLink> chemlink(m, "ChemLink");

  nb::bind_map<monomers_type, rv_ri>(m, "ChemCompMap");
  nb::bind_map<links_type, rv_ri>(m, "ChemLinkMap");
  nb::bind_map<modifications_type, rv_ri>(m, "ChemModMap");

  nb::class_<ChemLink::Side>(chemlink, "Side")
    .def(nb::init<>())
    .def_rw("comp", &ChemLink::Side::comp)
    .def_rw("mod", &ChemLink::Side::mod)
    .def_rw("group", &ChemLink::Side::group)
    .def("__repr__", [](const ChemLink::Side& self) {
        return "<gemmi.ChemLink.Side " + self.comp + "/" +
               ChemComp::group_str(self.group) + ">";
    });

  nb::class_<ChemMod::AtomMod>(chemmod, "AtomMod")
    .def_rw("func", &ChemMod::AtomMod::func)
    .def_rw("old_id", &ChemMod::AtomMod::old_id)
    .def_rw("new_id", &ChemMod::AtomMod::new_id)
    .def_rw("el", &ChemMod::AtomMod::el)
    .def_rw("charge", &ChemMod::AtomMod::charge)
    .def_rw("chem_type", &ChemMod::AtomMod::chem_type)
  ;
  chemlink
    .def(nb::init<>())
    .def_rw("id", &ChemLink::id)
    .def_rw("name", &ChemLink::name)
    .def_rw("side1", &ChemLink::side1)
    .def_rw("side2", &ChemLink::side2)
    .def_rw("rt", &ChemLink::rt)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemLink " + self.id + ">";
    });

  chemmod
    .def(nb::init<>())
    .def_rw("id", &ChemMod::id)
    .def_rw("name", &ChemMod::name)
    .def_rw("comp_id", &ChemMod::comp_id)
    .def_rw("group_id", &ChemMod::group_id)
    .def_rw("atom_mods", &ChemMod::atom_mods)
    .def_rw("rt", &ChemMod::rt)
    .def("__repr__", [](const ChemMod& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  nb::class_<EnerLib>(m, "EnerLib");  // NOLINT(bugprone-unused-raii)
  nb::class_<MonLib>(m, "MonLib")
    .def(nb::init<>())
    .def_prop_rw("monomer_dir",
                  [](const MonLib& self) { return self.monomer_dir; },
                  &MonLib::set_monomer_dir)
    .def_ro("ener_lib", &MonLib::ener_lib)
    .def_ro("monomers", &MonLib::monomers)
    .def_ro("links", &MonLib::links)
    .def_ro("modifications", &MonLib::modifications)
    .def("get_link", &MonLib::get_link, nb::arg("link_id"),
         nb::rv_policy::reference_internal)
    .def("get_mod", &MonLib::get_mod, nb::arg("name"),
         nb::rv_policy::reference_internal)
    .def("match_link", &MonLib::match_link,
         nb::arg("res1"), nb::arg("atom1"), nb::arg("alt1"),
         nb::arg("res2"), nb::arg("atom2"), nb::arg("alt2"),
         nb::arg("min_bond_sq")=0.,
         nb::rv_policy::reference_internal)
    .def("test_link", [](const MonLib& self, const ChemLink& link,
                         const std::string& res1, const std::string& atom1,
                         const std::string& res2, const std::string& atom2)
    -> std::tuple<bool, const ChemComp::Aliasing*, const ChemComp::Aliasing*> {
      const ChemComp::Aliasing* aliasing1 = nullptr;
      const ChemComp::Aliasing* aliasing2 = nullptr;
      bool match = (!link.rt.bonds.empty() &&
                    self.link_side_matches_residue(link.side1, res1, &aliasing1) &&
                    self.link_side_matches_residue(link.side2, res2, &aliasing2) &&
                    atom_match_with_alias(link.rt.bonds[0].id1.atom, atom1, aliasing1) &&
                    atom_match_with_alias(link.rt.bonds[0].id2.atom, atom2, aliasing2));
      return std::make_tuple(match, aliasing1, aliasing2);
    }, nb::arg("link"), nb::arg("res1"), nb::arg("atom1"), nb::arg("res2"), nb::arg("atom2"),
         nb::rv_policy::reference_internal)
    .def("add_monomer_if_present", &MonLib::add_monomer_if_present)
    .def("read_monomer_doc", &MonLib::read_monomer_doc)
    .def("read_monomer_cif", &MonLib::read_monomer_cif)
    .def("read_monomer_lib", &MonLib::read_monomer_lib,
         nb::arg("monomer_dir"), nb::arg("resnames"), nb::arg("logging")=nb::none())
    .def("find_ideal_distance", [](const MonLib& self, CRA &cra1, CRA cra2) {
      return self.find_ideal_distance(cra1, cra2);
    })
    .def("update_old_atom_names", &MonLib::update_old_atom_names,
         nb::arg("st"), nb::arg("logging")=nb::none())
    .def("path", &MonLib::path, nb::arg("code")=std::string())
    .def("__repr__", [](const MonLib& self) {
        return "<gemmi.MonLib with " +
               std::to_string(self.monomers.size()) + " monomers, " +
               std::to_string(self.links.size()) + " links, " +
               std::to_string(self.modifications.size()) + " modifications>";
    })
    .def("clone", [](const MonLib& self) { return new MonLib(self); });

  m.def("read_monomer_lib", &read_monomer_lib,
        nb::arg("monomer_dir"), nb::arg("resnames"),
        nb::arg("libin")=std::string(), nb::arg("ignore_missing")=false);
}
