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
using residue_infos_type = std::map<std::string, ResidueInfo>;
PYBIND11_MAKE_OPAQUE(monomers_type)
PYBIND11_MAKE_OPAQUE(links_type)
PYBIND11_MAKE_OPAQUE(modifications_type)
PYBIND11_MAKE_OPAQUE(residue_infos_type)

void add_monlib(py::module& m) {
  py::class_<ChemMod> chemmod(m, "ChemMod");
  py::class_<ChemLink> chemlink(m, "ChemLink");
  py::class_<ChemLink::Side> chemlinkside(chemlink, "Side");

  py::bind_map<monomers_type>(m, "ChemCompMap");
  py::bind_map<links_type>(m, "ChemLinkMap");
  py::bind_map<modifications_type>(m, "ChemModMap");
  py::bind_map<residue_infos_type>(m, "ResidueInfoMap");

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
      .value("Ketopyranose", ChemLink::Group::Ketopyranose)
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
    .def_readwrite("name", &ChemMod::name)
    .def_readwrite("comp_id", &ChemMod::comp_id)
    .def_readwrite("group_id", &ChemMod::group_id)
    .def_readwrite("rt", &ChemMod::rt)
    .def("__repr__", [](const ChemMod& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  py::class_<MonLib>(m, "MonLib")
    .def(py::init<>())
    .def_readonly("mon_lib_list", &MonLib::mon_lib_list)
    .def_readonly("monomers", &MonLib::monomers)
    .def_readonly("links", &MonLib::links)
    .def_readonly("modifications", &MonLib::modifications)
    .def_readonly("residue_infos", &MonLib::residue_infos)
    .def("get_link", &MonLib::get_link, py::arg("link_id"),
         py::return_value_policy::reference_internal)
    .def("find_mod", &MonLib::find_mod, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("find_residue_info", &MonLib::find_residue_info, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("match_link", &MonLib::match_link,
         py::arg("res1"), py::arg("atom1"),
         py::arg("res2"), py::arg("atom2"), py::arg("altloc"),
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
    .def("read_monomer_cif", [](MonLib& self, const std::string& path) {
      return self.read_monomer_cif(path, gemmi::read_cif_gz);
    })
    .def("read_monomer_lib", [](MonLib& self, const std::string& monomer_dir,
                                const std::vector<std::string>& resnames) {
      return self.read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz);
    })
    .def("path", &MonLib::path, py::arg("code")=nullptr)
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
