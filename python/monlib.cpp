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
    .def_readwrite("rt", &ChemMod::rt)
    .def("__repr__", [](const ChemMod& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  py::class_<MonLib>(m, "MonLib")
    .def(py::init<>())
    .def_readonly("monomer_dir", &MonLib::monomer_dir)
    .def_readonly("monomers", &MonLib::monomers)
    .def_readonly("links", &MonLib::links)
    .def_readonly("modifications", &MonLib::modifications)
    .def("get_link", &MonLib::get_link, py::arg("link_id"),
         py::return_value_policy::reference_internal)
    .def("get_mod", &MonLib::get_mod, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("match_link", [](const MonLib &self,
                          const Residue& res1, const std::string& atom1,
                          const Residue& res2, const std::string& atom2,
                          char alt, double min_bond_sq) {
      const ChemComp::Aliasing* aliasing1 = nullptr;
      const ChemComp::Aliasing* aliasing2 = nullptr;
      return self.match_link(res1, atom1, res2, atom2, alt, min_bond_sq,
                             &aliasing1, &aliasing2);
    }, py::arg("res1"), py::arg("atom1"),
       py::arg("res2"), py::arg("atom2"), py::arg("altloc"),
       py::arg("min_bond_sq")=0.,
       py::return_value_policy::reference_internal)
    .def("add_monomer_if_present", &MonLib::add_monomer_if_present)
    .def("insert_chemcomps", &MonLib::insert_chemcomps)
    .def("insert_chemlinks", &MonLib::insert_chemlinks)
    .def("insert_chemmods", &MonLib::insert_chemmods)
    // deprecated
    .def("insert_comp_list", [](MonLib &self, const cif::Document &doc) {
        insert_comp_list(doc, self.cc_groups);
    })
    .def("read_monomer_cif", [](MonLib& self, const std::string& path) {
      return self.read_monomer_cif(path, gemmi::read_cif_gz);
    })
    .def("read_monomer_lib", [](MonLib& self, const std::string& monomer_dir,
                                const std::vector<std::string>& resnames) {
      return self.read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz);
    })
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
