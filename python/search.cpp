// Copyright 2018 Global Phasing Ltd.

#include "gemmi/neighbor.hpp"
#include "gemmi/linkhunt.hpp"
#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<NeighborSearch::Mark*>)

void add_search(py::module& m) {
  py::class_<NeighborSearch> neighbor_search(m, "NeighborSearch");
  py::class_<NeighborSearch::Mark>(neighbor_search, "Mark")
    .def_readonly("x", &NeighborSearch::Mark::x)
    .def_readonly("y", &NeighborSearch::Mark::y)
    .def_readonly("z", &NeighborSearch::Mark::z)
    .def_readonly("altloc", &NeighborSearch::Mark::altloc)
    .def_readonly("element", &NeighborSearch::Mark::element)
    .def_readonly("image_idx", &NeighborSearch::Mark::image_idx)
    .def_readonly("chain_idx", &NeighborSearch::Mark::chain_idx)
    .def_readonly("residue_idx", &NeighborSearch::Mark::residue_idx)
    .def_readonly("atom_idx", &NeighborSearch::Mark::atom_idx)
    .def("pos", &NeighborSearch::Mark::pos)
    .def("to_cra", (CRA (NeighborSearch::Mark::*)(Model&) const)
                   &NeighborSearch::Mark::to_cra)
    .def("to_site", (SmallStructure::Site& (NeighborSearch::Mark::*)(SmallStructure&) const)
                    &NeighborSearch::Mark::to_site)
    .def("__repr__", [](const NeighborSearch::Mark& self) {
        return cat("<gemmi.NeighborSearch.Mark ", self.element.name(),
                   " of atom ", self.chain_idx, '/', self.residue_idx, '/',
                   self.atom_idx, '>');
    });
  py::bind_vector<std::vector<NeighborSearch::Mark*>>(m, "VectorMarkPtr");
  neighbor_search
    .def_readonly("radius_specified", &NeighborSearch::radius_specified)
    .def(py::init<Model&, const UnitCell&, double>(),
         py::arg("model"), py::arg("cell"), py::arg("max_radius")/*,
         py::keep_alive<1, 2>()*/)
    .def(py::init([](Structure& st, double max_radius, int model_index) {
      return new NeighborSearch(st.models.at(model_index), st.cell, max_radius);
    }), py::arg("st"), py::arg("max_radius"), py::arg("model_index")=0,
        py::keep_alive<1, 2>())
    .def(py::init<SmallStructure&, double>(),
         py::arg("small_structure"), py::arg("max_radius"),
         py::keep_alive<1, 2>())
    .def("populate", &NeighborSearch::populate, py::arg("include_h")=true,
         "Usually run after constructing NeighborSearch.")
    .def("add_chain", &NeighborSearch::add_chain,
         py::arg("chain"), py::arg("include_h")=true)
    .def("add_atom", &NeighborSearch::add_atom,
         py::arg("atom"), py::arg("n_ch"), py::arg("n_res"), py::arg("n_atom"),
         "Lower-level alternative to populate()")
    .def("find_atoms", &NeighborSearch::find_atoms,
         py::arg("pos"), py::arg("alt")='\0', py::arg("radius")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("find_neighbors", &NeighborSearch::find_neighbors,
         py::arg("atom"), py::arg("min_dist")=0, py::arg("max_dist")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("find_nearest_atom", &NeighborSearch::find_nearest_atom,
         py::return_value_policy::reference_internal)
    .def("find_site_neighbors", &NeighborSearch::find_site_neighbors,
         py::arg("atom"), py::arg("min_dist")=0, py::arg("max_dist")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("dist", &NeighborSearch::dist)
    .def("get_image_transformation", &NeighborSearch::get_image_transformation)
    .def_property_readonly("grid_cell",
        [](const NeighborSearch& self) { return self.grid.unit_cell; })
    .def("__repr__", [](const NeighborSearch& self) {
        return cat("<gemmi.NeighborSearch with grid ",
                   self.grid.nu, ", ", self.grid.nv, ", ", self.grid.nw, '>');
    });
  m.def("merge_atoms_in_expanded_model", &merge_atoms_in_expanded_model,
        py::arg("model"), py::arg("cell"), py::arg("max_dist")=0.2);

  py::class_<ContactSearch> contactsearch(m, "ContactSearch");
  py::enum_<ContactSearch::Ignore> csignore(contactsearch, "Ignore");
  py::class_<ContactSearch::Result> csresult(contactsearch, "Result");

  contactsearch
    .def(py::init<float>())
    .def_readwrite("search_radius", &ContactSearch::search_radius)
    .def_readwrite("ignore", &ContactSearch::ignore)
    .def_readwrite("twice", &ContactSearch::twice)
    .def_readwrite("special_pos_cutoff_sq", &ContactSearch::special_pos_cutoff_sq)
    .def_readwrite("min_occupancy", &ContactSearch::min_occupancy)
    .def("setup_atomic_radii", &ContactSearch::setup_atomic_radii)
    .def("get_radius", [](const ContactSearch& self, Element el) {
        return self.get_radius(el.elem);
    })
    .def("set_radius", [](ContactSearch& self, Element el, float r) {
        self.set_radius(el.elem, r);
    })
    .def("find_contacts", &ContactSearch::find_contacts)
    ;

  csignore
    .value("Nothing", ContactSearch::Ignore::Nothing)
    .value("SameResidue", ContactSearch::Ignore::SameResidue)
    .value("AdjacentResidues", ContactSearch::Ignore::AdjacentResidues)
    .value("SameChain", ContactSearch::Ignore::SameChain)
    .value("SameAsu", ContactSearch::Ignore::SameAsu);

  csresult
    .def_readonly("partner1", &ContactSearch::Result::partner1)
    .def_readonly("partner2", &ContactSearch::Result::partner2)
    .def_readonly("image_idx", &ContactSearch::Result::image_idx)
    .def_property_readonly("dist", [](ContactSearch::Result& self) {
        return std::sqrt(self.dist_sq);
    })
    ;

  py::class_<LinkHunt> linkhunt(m, "LinkHunt");
  py::class_<LinkHunt::Match> linkhuntmatch(linkhunt, "Match");
  linkhunt
    .def(py::init<>())
    .def("index_chem_links", &LinkHunt::index_chem_links,
         py::arg("monlib"), py::arg("use_alias")=true, py::keep_alive<1, 2>())
    .def("find_possible_links", &LinkHunt::find_possible_links,
         py::arg("st"), py::arg("bond_margin"), py::arg("radius_margin"),
         py::arg("ignore")=ContactSearch::Ignore::SameResidue)
    ;

  linkhuntmatch
    .def_readonly("chem_link", &LinkHunt::Match::chem_link)
    .def_readonly("chem_link_count", &LinkHunt::Match::chem_link_count)
    .def_readonly("cra1", &LinkHunt::Match::cra1)
    .def_readonly("cra2", &LinkHunt::Match::cra2)
    .def_readonly("same_image", &LinkHunt::Match::same_image)
    .def_readonly("bond_length", &LinkHunt::Match::bond_length)
    .def_readonly("conn", &LinkHunt::Match::conn)
    ;
}
