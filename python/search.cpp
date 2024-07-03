// Copyright 2018 Global Phasing Ltd.

#include "gemmi/neighbor.hpp"
#include "gemmi/linkhunt.hpp"
#include "gemmi/bond_idx.hpp"
#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<NeighborSearch::Mark*>)

void add_search(nb::module_& m) {
  nb::class_<NeighborSearch> neighbor_search(m, "NeighborSearch");
  nb::class_<NeighborSearch::Mark>(neighbor_search, "Mark")
    .def_ro("pos", &NeighborSearch::Mark::pos)
    .def_ro("altloc", &NeighborSearch::Mark::altloc)
    .def_ro("element", &NeighborSearch::Mark::element)
    .def_ro("image_idx", &NeighborSearch::Mark::image_idx)
    .def_ro("chain_idx", &NeighborSearch::Mark::chain_idx)
    .def_ro("residue_idx", &NeighborSearch::Mark::residue_idx)
    .def_ro("atom_idx", &NeighborSearch::Mark::atom_idx)
    .def("to_cra", (CRA (NeighborSearch::Mark::*)(Model&) const)
                   &NeighborSearch::Mark::to_cra)
    .def("to_site", (SmallStructure::Site& (NeighborSearch::Mark::*)(SmallStructure&) const)
                    &NeighborSearch::Mark::to_site)
    .def("__repr__", [](const NeighborSearch::Mark& self) {
        return cat("<gemmi.NeighborSearch.Mark ", int(self.image_idx), " of atom ",
                   self.chain_idx, '/', self.residue_idx, '/', self.atom_idx,
                   " element ", self.element.name(), ">");
    });
  nb::bind_vector<std::vector<NeighborSearch::Mark*>>(m, "VectorMarkPtr");
  neighbor_search
    .def_ro("radius_specified", &NeighborSearch::radius_specified)
    .def(nb::init<Model&, const UnitCell&, double>(),
         nb::arg("model"), nb::arg("cell"), nb::arg("max_radius")/*,
         nb::keep_alive<1, 2>()*/)
    .def("__init__", [](NeighborSearch* ns, Structure& st,
                        double max_radius, int model_index) {
      new(ns) NeighborSearch(st.models.at(model_index), st.cell, max_radius);
    }, nb::arg("st"), nb::arg("max_radius"), nb::arg("model_index")=0,
        nb::keep_alive<1, 2>())
    .def(nb::init<SmallStructure&, double>(),
         nb::arg("small_structure"), nb::arg("max_radius"),
         nb::keep_alive<1, 2>())
    .def("populate", &NeighborSearch::populate, nb::arg("include_h")=true,
         "Usually run after constructing NeighborSearch.")
    .def("add_chain", &NeighborSearch::add_chain,
         nb::arg("chain"), nb::arg("include_h")=true)
    .def("add_atom", &NeighborSearch::add_atom,
         nb::arg("atom"), nb::arg("n_ch"), nb::arg("n_res"), nb::arg("n_atom"),
         "Lower-level alternative to populate()")
    .def("add_site", &NeighborSearch::add_site,
         nb::arg("site"), nb::arg("n"),
         "Lower-level alternative to populate() for SmallStructure")
    .def("find_atoms", &NeighborSearch::find_atoms,
         nb::arg("pos"), nb::arg("alt")='\0',
         nb::kw_only(), nb::arg("min_dist")=0, nb::arg("radius")=0,
         nb::rv_policy::move, nb::keep_alive<0, 1>())
    .def("find_neighbors", &NeighborSearch::find_neighbors,
         nb::arg("atom"), nb::arg("min_dist")=0, nb::arg("max_dist")=0,
         nb::rv_policy::move, nb::keep_alive<0, 1>())
    .def("find_nearest_atom", &NeighborSearch::find_nearest_atom,
         nb::arg("pos"), nb::arg("radius")=INFINITY,
         nb::rv_policy::reference_internal)
    .def("find_site_neighbors", &NeighborSearch::find_site_neighbors,
         nb::arg("atom"), nb::arg("min_dist")=0, nb::arg("max_dist")=0,
         nb::rv_policy::move, nb::keep_alive<0, 1>())
    .def("dist", &NeighborSearch::dist)
    .def("get_image_transformation", &NeighborSearch::get_image_transformation)
    .def_prop_ro("grid_cell",
        [](const NeighborSearch& self) { return self.grid.unit_cell; })
    .def("__repr__", [](const NeighborSearch& self) {
        return cat("<gemmi.NeighborSearch with grid ",
                   self.grid.nu, ", ", self.grid.nv, ", ", self.grid.nw, '>');
    });

  nb::class_<ContactSearch> contactsearch(m, "ContactSearch");
  nb::enum_<ContactSearch::Ignore> csignore(contactsearch, "Ignore");
  nb::class_<ContactSearch::Result> csresult(contactsearch, "Result");

  contactsearch
    .def(nb::init<double>())
    .def_rw("search_radius", &ContactSearch::search_radius)
    .def_rw("ignore", &ContactSearch::ignore)
    .def_rw("twice", &ContactSearch::twice)
    .def_rw("special_pos_cutoff_sq", &ContactSearch::special_pos_cutoff_sq)
    .def_rw("min_occupancy", &ContactSearch::min_occupancy)
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
    .def_ro("partner1", &ContactSearch::Result::partner1)
    .def_ro("partner2", &ContactSearch::Result::partner2)
    .def_ro("image_idx", &ContactSearch::Result::image_idx)
    .def_prop_ro("dist", [](ContactSearch::Result& self) {
        return std::sqrt(self.dist_sq);
    })
    ;

  nb::class_<LinkHunt> linkhunt(m, "LinkHunt");
  nb::class_<LinkHunt::Match> linkhuntmatch(linkhunt, "Match");
  linkhunt
    .def(nb::init<>())
    .def("index_chem_links", &LinkHunt::index_chem_links,
         nb::arg("monlib"), nb::arg("use_alias")=true, nb::keep_alive<1, 2>())
    .def("find_possible_links", &LinkHunt::find_possible_links,
         nb::arg("st"), nb::arg("bond_margin"), nb::arg("radius_margin"),
         nb::arg("ignore")=ContactSearch::Ignore::SameResidue)
    ;

  linkhuntmatch
    .def_ro("chem_link", &LinkHunt::Match::chem_link)
    .def_ro("chem_link_count", &LinkHunt::Match::chem_link_count)
    .def_ro("cra1", &LinkHunt::Match::cra1)
    .def_ro("cra2", &LinkHunt::Match::cra2)
    .def_ro("same_image", &LinkHunt::Match::same_image)
    .def_ro("bond_length", &LinkHunt::Match::bond_length)
    .def_ro("conn", &LinkHunt::Match::conn)
    ;

  nb::class_<BondIndex>(m, "BondIndex")
    .def(nb::init<const Model&>(), nb::keep_alive<1, 2>())
    .def("add_link", &BondIndex::add_link)
    .def("add_monomer_bonds", &BondIndex::add_monomer_bonds)
    .def("are_linked", &BondIndex::are_linked)
    .def("graph_distance", &BondIndex::graph_distance,
         nb::arg("a"), nb::arg("b"), nb::arg("same_index"),
         nb::arg("max_distance")=4)
    ;
}
