// Copyright 2021 Global Phasing Ltd.

#include "gemmi/topo.hpp"
#include "gemmi/placeh.hpp"  // for adjust_hydrogen_distances

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Topo::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Angle>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Torsion>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Chirality>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Plane>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::ExtraLink>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::ChainInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::ResInfo::Prev>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::ResInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Force>)

void add_topo(py::module& m) {
  py::class_<Topo> topo(m, "Topo");

  py::enum_<HydrogenChange>(m, "HydrogenChange")
    .value("None", HydrogenChange::None)
    .value("Shift", HydrogenChange::Shift)
    .value("Remove", HydrogenChange::Remove)
    .value("ReAdd", HydrogenChange::ReAdd)
    .value("ReAddButWater", HydrogenChange::ReAddButWater);

  py::class_<Topo::Bond>(topo, "Bond")
    .def_readonly("restr", &Topo::Bond::restr)
    .def_readonly("atoms", &Topo::Bond::atoms)
    .def("calculate", &Topo::Bond::calculate)
    .def("calculate_z", &Topo::Bond::calculate_z)
    ;
  py::class_<Topo::Angle>(topo, "Angle")
    .def_readonly("restr", &Topo::Angle::restr)
    .def_readonly("atoms", &Topo::Angle::atoms)
    .def("calculate", &Topo::Angle::calculate)
    .def("calculate_z", &Topo::Angle::calculate_z)
    ;
  py::class_<Topo::Torsion>(topo, "Torsion")
    .def_readonly("restr", &Topo::Torsion::restr)
    .def_readonly("atoms", &Topo::Torsion::atoms)
    .def("calculate", &Topo::Torsion::calculate)
    .def("calculate_z", &Topo::Torsion::calculate_z)
    ;
  py::class_<Topo::Chirality>(topo, "Chirality")
    .def_readonly("restr", &Topo::Chirality::restr)
    .def_readonly("atoms", &Topo::Chirality::atoms)
    .def("calculate", &Topo::Chirality::calculate)
    .def("check", &Topo::Chirality::check)
    ;
  py::class_<Topo::Plane>(topo, "Plane")
    .def_readonly("restr", &Topo::Plane::restr)
    .def_readonly("atoms", &Topo::Plane::atoms)
    .def("has", &Topo::Plane::has)
    ;
  py::enum_<Topo::Provenance>(m, "Provenance")
    .value("None", Topo::Provenance::None)
    .value("PrevLink", Topo::Provenance::PrevLink)
    .value("Monomer", Topo::Provenance::Monomer)
    .value("NextLink", Topo::Provenance::NextLink)
    .value("ExtraLink", Topo::Provenance::ExtraLink)
    ;
  py::enum_<Topo::RKind>(m, "RKind")
    .value("Bond", Topo::RKind::Bond)
    .value("Angle", Topo::RKind::Angle)
    .value("Torsion", Topo::RKind::Torsion)
    .value("Chirality", Topo::RKind::Chirality)
    .value("Plane", Topo::RKind::Plane)
    ;
  py::class_<Topo::Force>(topo, "Force")
    .def_readonly("provenance", &Topo::Force::provenance)
    .def_readonly("rkind", &Topo::Force::rkind)
    .def_readonly("index", &Topo::Force::index)
    ;
  py::class_<Topo::ResInfo> resinfo(topo, "ResInfo");
  py::class_<Topo::ResInfo::Prev>(resinfo, "Prev")
    .def_readonly("link", &Topo::ResInfo::Prev::link)
    .def_readonly("idx", &Topo::ResInfo::Prev::idx)
    .def("get", (Topo::ResInfo* (Topo::ResInfo::Prev::*)(Topo::ResInfo*)const) &Topo::ResInfo::Prev::get, py::return_value_policy::reference_internal)
    ;
  resinfo
    .def_readonly("res", &Topo::ResInfo::res)
    .def_readonly("prev", &Topo::ResInfo::prev)
    .def_readonly("mods", &Topo::ResInfo::mods)
    .def_readonly("chemcomp", &Topo::ResInfo::chemcomp)
    .def_readonly("forces", &Topo::ResInfo::forces)
    ;
  py::class_<Topo::ChainInfo>(topo, "ChainInfo")
    .def_readonly("name", &Topo::ChainInfo::name)
    .def_readonly("entity_id", &Topo::ChainInfo::entity_id)
    .def_readonly("polymer", &Topo::ChainInfo::polymer)
    .def_readonly("polymer_type", &Topo::ChainInfo::polymer_type)
    .def_readonly("res_infos", &Topo::ChainInfo::res_infos)
    ;
  py::class_<Topo::ExtraLink>(topo, "ExtraLink")
    .def_readonly("res1", &Topo::ExtraLink::res1)
    .def_readonly("res2", &Topo::ExtraLink::res2)
    .def_readonly("alt1", &Topo::ExtraLink::alt1)
    .def_readonly("alt2", &Topo::ExtraLink::alt2)
    .def_readonly("link_id", &Topo::ExtraLink::link_id)
    ;

  py::bind_vector<std::vector<Topo::Bond>>(m, "TopoBonds");
  py::bind_vector<std::vector<Topo::Angle>>(m, "TopoAngles");
  py::bind_vector<std::vector<Topo::Torsion>>(m, "TopoTorsions");
  py::bind_vector<std::vector<Topo::Chirality>>(m, "TopoChirs");
  py::bind_vector<std::vector<Topo::Plane>>(m, "TopoPlanes");
  py::bind_vector<std::vector<Topo::ChainInfo>>(m, "TopoChainInfos");
  py::bind_vector<std::vector<Topo::Force>>(m, "TopoForces");
  py::bind_vector<std::vector<Topo::ResInfo::Prev>>(m, "TopoResInfoPrevs");
  py::bind_vector<std::vector<Topo::ResInfo>>(m, "TopoResInfos");
  py::bind_vector<std::vector<Topo::ExtraLink>>(m, "TopoExtraLinks");

  topo
    .def(py::init<>())
    .def("adjust_hydrogen_distances", &adjust_hydrogen_distances,
         py::arg("of"), py::arg("default_scale")=1.)
    .def_readonly("bonds", &Topo::bonds)
    .def_readonly("angles", &Topo::angles)
    .def_readonly("torsions", &Topo::torsions)
    .def_readonly("chirs", &Topo::chirs)
    .def_readonly("planes", &Topo::planes)
    .def_readonly("extras", &Topo::extras)
    .def_readonly("chain_infos", &Topo::chain_infos)
    ;

  m.def("prepare_topology", &prepare_topology,
        py::arg("st"), py::arg("monlib"), py::arg("model_index")=0,
        py::arg("h_change")=HydrogenChange::None, py::arg("reorder")=false,
        py::arg("raise_errors")=false);
}
