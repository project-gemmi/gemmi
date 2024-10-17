// Copyright 2021 Global Phasing Ltd.

#include "gemmi/topo.hpp"
#include "gemmi/riding_h.hpp"  // for adjust_hydrogen_distances
#include "gemmi/crd.hpp"       // for prepare_refmac_crd, ...

#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/array.h>   // for Topo::Bond::atoms
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for find_missing_atoms
#include <nanobind/stl/unique_ptr.h>  // for prepare_topology

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<Topo::Bond>)
NB_MAKE_OPAQUE(std::vector<Topo::Angle>)
NB_MAKE_OPAQUE(std::vector<Topo::Torsion>)
NB_MAKE_OPAQUE(std::vector<Topo::Chirality>)
NB_MAKE_OPAQUE(std::vector<Topo::Plane>)
NB_MAKE_OPAQUE(std::vector<Topo::Link>)
NB_MAKE_OPAQUE(std::vector<Topo::ChainInfo>)
NB_MAKE_OPAQUE(std::vector<Topo::ResInfo>)
NB_MAKE_OPAQUE(std::vector<Topo::Rule>)
NB_MAKE_OPAQUE(std::vector<Topo::Mod>)
NB_MAKE_OPAQUE(std::vector<Topo::FinalChemComp>)

void add_topo(nb::module_& m) {
  nb::class_<Topo> topo(m, "Topo");

  nb::enum_<HydrogenChange>(m, "HydrogenChange")
    .value("NoChange", HydrogenChange::NoChange)
    .value("Shift", HydrogenChange::Shift)
    .value("Remove", HydrogenChange::Remove)
    .value("ReAdd", HydrogenChange::ReAdd)
    .value("ReAddButWater", HydrogenChange::ReAddButWater)
    .value("ReAddKnown", HydrogenChange::ReAddKnown);

  nb::class_<Topo::Bond>(topo, "Bond")
    .def_ro("restr", &Topo::Bond::restr)
    .def_ro("atoms", &Topo::Bond::atoms)
    .def("calculate", &Topo::Bond::calculate)
    .def("calculate_z", &Topo::Bond::calculate_z)
    ;
  nb::bind_vector<std::vector<Topo::Bond>, rv_ri>(m, "TopoBonds");

  nb::class_<Topo::Angle>(topo, "Angle")
    .def_ro("restr", &Topo::Angle::restr)
    .def_ro("atoms", &Topo::Angle::atoms)
    .def("calculate", &Topo::Angle::calculate)
    .def("calculate_z", &Topo::Angle::calculate_z)
    ;
  nb::bind_vector<std::vector<Topo::Angle>, rv_ri>(m, "TopoAngles");

  nb::class_<Topo::Torsion>(topo, "Torsion")
    .def_ro("restr", &Topo::Torsion::restr)
    .def_ro("atoms", &Topo::Torsion::atoms)
    .def("calculate", &Topo::Torsion::calculate)
    .def("calculate_z", &Topo::Torsion::calculate_z)
    ;
  nb::bind_vector<std::vector<Topo::Torsion>, rv_ri>(m, "TopoTorsions");

  nb::class_<Topo::Chirality>(topo, "Chirality")
    .def_ro("restr", &Topo::Chirality::restr)
    .def_ro("atoms", &Topo::Chirality::atoms)
    .def("calculate", &Topo::Chirality::calculate)
    .def("calculate_z", &Topo::Chirality::calculate_z,
         nb::arg("ideal_abs_vol"), nb::arg("esd"))
    .def("check", &Topo::Chirality::check)
    ;
  nb::bind_vector<std::vector<Topo::Chirality>, rv_ri>(m, "TopoChirs");

  nb::class_<Topo::Plane>(topo, "Plane")
    .def_ro("restr", &Topo::Plane::restr)
    .def_ro("atoms", &Topo::Plane::atoms)
    .def("has", &Topo::Plane::has)
    ;
  nb::bind_vector<std::vector<Topo::Plane>, rv_ri>(m, "TopoPlanes");

  nb::enum_<Topo::RKind>(m, "RKind")
    .value("Bond", Topo::RKind::Bond)
    .value("Angle", Topo::RKind::Angle)
    .value("Torsion", Topo::RKind::Torsion)
    .value("Chirality", Topo::RKind::Chirality)
    .value("Plane", Topo::RKind::Plane)
    ;
  nb::class_<Topo::Rule>(topo, "Rule")
    .def_ro("rkind", &Topo::Rule::rkind)
    .def_ro("index", &Topo::Rule::index)
    ;
  nb::bind_vector<std::vector<Topo::Rule>, rv_ri>(m, "TopoRules");

  nb::class_<Topo::Mod>(topo, "Mod")
    .def_ro("id", &Topo::Mod::id)
    .def_ro("alias", &Topo::Mod::alias)
    .def_ro("altloc", &Topo::Mod::altloc)
    ;
  nb::bind_vector<std::vector<Topo::Mod>, rv_ri>(m, "TopoMods");

  nb::class_<Topo::FinalChemComp>(topo, "FinalChemComp")
    .def_ro("altloc", &Topo::FinalChemComp::altloc)
    .def_ro("cc", &Topo::FinalChemComp::cc)
    ;
  nb::bind_vector<std::vector<Topo::FinalChemComp>, rv_ri>(m, "TopoFinalChemComps");

  nb::class_<Topo::Link>(topo, "Link")
    .def_ro("link_id", &Topo::Link::link_id)
    .def_ro("res1", &Topo::Link::res1)
    .def_ro("res2", &Topo::Link::res2)
    .def_ro("alt1", &Topo::Link::alt1)
    .def_ro("alt2", &Topo::Link::alt2)
    .def_ro("link_rules", &Topo::Link::link_rules)
    ;
  nb::bind_vector<std::vector<Topo::Link>, rv_ri>(m, "TopoLinks");

  nb::class_<Topo::ResInfo>(topo, "ResInfo")
    .def_ro("res", &Topo::ResInfo::res)
    .def_ro("prev", &Topo::ResInfo::prev)
    .def_ro("mods", &Topo::ResInfo::mods)
    .def_ro("chemcomps", &Topo::ResInfo::chemcomps)
    .def_ro("monomer_rules", &Topo::ResInfo::monomer_rules)
    .def("get_final_chemcomp", &Topo::ResInfo::get_final_chemcomp)
    ;
  nb::bind_vector<std::vector<Topo::ResInfo>, rv_ri>(m, "TopoResInfos");

  nb::class_<Topo::ChainInfo>(topo, "ChainInfo")
    .def_prop_ro("chain_ref", [](const Topo::ChainInfo& self)
        { return &self.chain_ref; }, nb::rv_policy::reference_internal)
    .def_ro("subchain_name", &Topo::ChainInfo::subchain_name)
    .def_ro("entity_id", &Topo::ChainInfo::entity_id)
    .def_ro("polymer", &Topo::ChainInfo::polymer)
    .def_ro("polymer_type", &Topo::ChainInfo::polymer_type)
    .def_ro("res_infos", &Topo::ChainInfo::res_infos)
    ;

  nb::class_<std::vector<Topo::ChainInfo>>(m, "TopoChainInfos")
    .def("__len__", &std::vector<Topo::ChainInfo>::size)
    .def("__getitem__", [](const std::vector<Topo::ChainInfo> &self, size_t i) {return self.at(i);})
    ;

  topo
    .def(nb::init<>())
    .def("adjust_hydrogen_distances", &adjust_hydrogen_distances,
         nb::arg("of"), nb::arg("default_scale")=1.)
    .def_ro("bonds", &Topo::bonds)
    .def_ro("angles", &Topo::angles)
    .def_ro("torsions", &Topo::torsions)
    .def_ro("chirs", &Topo::chirs)
    .def_ro("planes", &Topo::planes)
    .def_ro("extras", &Topo::extras)
    .def_ro("chain_infos", &Topo::chain_infos)
    .def("ideal_chiral_abs_volume", &Topo::ideal_chiral_abs_volume)
    .def("links_to_previous", [](Topo& self, Residue* res) {
        if (Topo::ResInfo* ri = self.find_resinfo(res))
          return ri->prev;
        fail("links_to_previous(): Residue not found");
    }, nb::rv_policy::reference_internal)
    .def("first_bond_in_link", &Topo::first_bond_in_link,
         nb::rv_policy::reference_internal)
    .def("set_cispeps_in_structure", &Topo::set_cispeps_in_structure)
    .def("find_missing_atoms", &find_missing_atoms,
         nb::arg("including_hydrogen")=false)
    ;

  m.def("prepare_topology", &prepare_topology,
        nb::arg("st"), nb::arg("monlib"), nb::arg("model_index")=0,
        nb::arg("h_change")=HydrogenChange::NoChange, nb::arg("reorder")=false,
        nb::arg("warnings")=nb::none(), nb::arg("ignore_unknown_links")=false,
        nb::arg("use_cispeps")=false);

  // crd.hpp
  m.def("setup_for_crd", &setup_for_crd);
  m.def("prepare_refmac_crd", &prepare_refmac_crd);
  m.def("add_automatic_links", &add_automatic_links);
  m.def("add_dictionary_blocks", &add_dictionary_blocks);
}
