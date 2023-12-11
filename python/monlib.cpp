// Copyright 2018 Global Phasing Ltd.

#include "gemmi/monlib.hpp"
#include "gemmi/read_cif.hpp"  // for read_cif_gz

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

using monomers_type = std::map<std::string, ChemComp>;
using links_type = std::map<std::string, ChemLink>;
using modifications_type = std::map<std::string, ChemMod>;
PYBIND11_MAKE_OPAQUE(monomers_type)
PYBIND11_MAKE_OPAQUE(links_type)
PYBIND11_MAKE_OPAQUE(modifications_type)

void add_monlib(py::module& m) {
  py::class_<ChemMod> chemmod(m, "ChemMod");
  py::class_<ChemLink> chemlink(m, "ChemLink");

  py::bind_map<monomers_type>(m, "ChemCompMap");
  py::bind_map<links_type>(m, "ChemLinkMap");
  py::bind_map<modifications_type>(m, "ChemModMap");

  py::class_<ChemLink::Side>(chemlink, "Side")
    .def(py::init<>())
    .def_readwrite("comp", &ChemLink::Side::comp)
    .def_readwrite("mod", &ChemLink::Side::mod)
    .def_readwrite("group", &ChemLink::Side::group)
    .def("__repr__", [](const ChemLink::Side& self) {
        return "<gemmi.ChemLink.Side " + self.comp + "/" +
               ChemComp::group_str(self.group) + ">";
    });

  py::class_<ChemMod::AtomMod>(chemmod, "AtomMod")
    .def_readwrite("func", &ChemMod::AtomMod::func)
    .def_readwrite("old_id", &ChemMod::AtomMod::old_id)
    .def_readwrite("new_id", &ChemMod::AtomMod::new_id)
    .def_readwrite("el", &ChemMod::AtomMod::el)
    .def_readwrite("charge", &ChemMod::AtomMod::charge)
    .def_readwrite("chem_type", &ChemMod::AtomMod::chem_type)
  ;
  chemlink
    .def(py::init<>())
    .def_readwrite("id", &ChemLink::id)
    .def_readwrite("name", &ChemLink::name)
    .def_readwrite("side1", &ChemLink::side1)
    .def_readwrite("side2", &ChemLink::side2)
    .def_readwrite("rt", &ChemLink::rt)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemLink " + self.id + ">";
    });

  chemmod
    .def(py::init<>())
    .def_readwrite("id", &ChemMod::id)
    .def_readwrite("name", &ChemMod::name)
    .def_readwrite("comp_id", &ChemMod::comp_id)
    .def_readwrite("group_id", &ChemMod::group_id)
    .def_readwrite("atom_mods", &ChemMod::atom_mods)
    .def_readwrite("rt", &ChemMod::rt)
    .def("__repr__", [](const ChemMod& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  py::class_<EnerLib>(m, "EnerLib");
  py::class_<MonLib>(m, "MonLib")
    .def(py::init<>())
    .def_readonly("monomer_dir", &MonLib::monomer_dir)
    .def_readonly("ener_lib", &MonLib::ener_lib)
    .def_readonly("monomers", &MonLib::monomers)
    .def_readonly("links", &MonLib::links)
    .def_readonly("modifications", &MonLib::modifications)
    .def("get_link", &MonLib::get_link, py::arg("link_id"),
         py::return_value_policy::reference_internal)
    .def("get_mod", &MonLib::get_mod, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("match_link", &MonLib::match_link,
         py::arg("res1"), py::arg("atom1"), py::arg("alt1"),
         py::arg("res2"), py::arg("atom2"), py::arg("alt2"),
         py::arg("min_bond_sq")=0.,
         py::return_value_policy::reference_internal)
    .def("test_link", [](const MonLib& self, const ChemLink& link,
                         const std::string& res1, const std::string& atom1,
                         const std::string& res2, const std::string& atom2) {
      const ChemComp::Aliasing* aliasing1 = nullptr;
      const ChemComp::Aliasing* aliasing2 = nullptr;
      bool match = (!link.rt.bonds.empty() &&
                    self.link_side_matches_residue(link.side1, res1, &aliasing1) &&
                    self.link_side_matches_residue(link.side2, res2, &aliasing2) &&
                    atom_match_with_alias(link.rt.bonds[0].id1.atom, atom1, aliasing1) &&
                    atom_match_with_alias(link.rt.bonds[0].id2.atom, atom2, aliasing2));
      return py::make_tuple(match, aliasing1, aliasing2);
    }, py::arg("link"), py::arg("res1"), py::arg("atom1"), py::arg("res2"), py::arg("atom2"))
    .def("add_monomer_if_present", &MonLib::add_monomer_if_present)
    .def("read_monomer_doc", &MonLib::read_monomer_doc)
    .def("read_monomer_cif", [](MonLib& self, const std::string& path) {
      return self.read_monomer_cif(path, gemmi::read_cif_gz);
    })
    .def("read_monomer_lib", [](MonLib& self, const std::string& monomer_dir,
                                const std::vector<std::string>& resnames) {
      return self.read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz);
    })
    .def("find_ideal_distance", [](const MonLib& self, CRA &cra1, CRA cra2) {
      return self.find_ideal_distance(cra1, cra2);
    })
    .def("update_old_atom_names", &MonLib::update_old_atom_names)
    .def("path", &MonLib::path, py::arg("code")=std::string())
    .def("__repr__", [](const MonLib& self) {
        return "<gemmi.MonLib with " +
               std::to_string(self.monomers.size()) + " monomers, " +
               std::to_string(self.links.size()) + " links, " +
               std::to_string(self.modifications.size()) + " modifications>";
    });

  m.def("read_monomer_lib", [](const std::string& monomer_dir,
                               const std::vector<std::string>& resnames,
                               const std::string& libin,
                               bool ignore_missing) {
    return read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz, libin, ignore_missing);
  }, py::arg("monomer_dir"), py::arg("resnames"), py::arg("libin")=std::string(),
     py::arg("ignore_missing")=false);
}
