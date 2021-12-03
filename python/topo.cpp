// Copyright 2021 Global Phasing Ltd.

#include "gemmi/topo.hpp"
#include "gemmi/placeh.hpp"  // for adjust_hydrogen_distances

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>  // for detail::pythonbuf

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Topo::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Angle>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Torsion>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Chirality>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Plane>)

void add_topo(py::module& m) {
  py::class_<Topo> topo(m, "Topo");

  py::enum_<HydrogenChange>(m, "HydrogenChange")
    .value("NoChange", HydrogenChange::NoChange)
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
    .def("calculate_z", &Topo::Chirality::calculate_z,
         py::arg("ideal_abs_vol"), py::arg("esd"))
    .def("check", &Topo::Chirality::check)
    ;
  py::class_<Topo::Plane>(topo, "Plane")
    .def_readonly("restr", &Topo::Plane::restr)
    .def_readonly("atoms", &Topo::Plane::atoms)
    .def("has", &Topo::Plane::has)
    ;
  py::bind_vector<std::vector<Topo::Bond>>(m, "TopoBonds");
  py::bind_vector<std::vector<Topo::Angle>>(m, "TopoAngles");
  py::bind_vector<std::vector<Topo::Torsion>>(m, "TopoTorsions");
  py::bind_vector<std::vector<Topo::Chirality>>(m, "TopoChirs");
  py::bind_vector<std::vector<Topo::Plane>>(m, "TopoPlanes");

  topo
    .def(py::init<>())
    .def("adjust_hydrogen_distances", &adjust_hydrogen_distances,
         py::arg("of"), py::arg("default_scale")=1.)
    .def_readonly("bonds", &Topo::bonds)
    .def_readonly("angles", &Topo::angles)
    .def_readonly("torsions", &Topo::torsions)
    .def_readonly("chirs", &Topo::chirs)
    .def_readonly("planes", &Topo::planes)
    .def("ideal_chiral_abs_volume", &Topo::ideal_chiral_abs_volume)
    ;

  m.def("prepare_topology",
    [](Structure& st, MonLib& monlib, size_t model_index,
       HydrogenChange h_change, bool reorder,
       const py::object& pywarnings, bool ignore_unknown_links) {
      std::ostream* warnings = nullptr;
      std::ostream os(nullptr);
      std::unique_ptr<py::detail::pythonbuf> buffer;
      if (!pywarnings.is_none()) {
        buffer.reset(new py::detail::pythonbuf(pywarnings));
        os.rdbuf(buffer.get());
        warnings = &os;
      }
      prepare_topology(st, monlib, model_index, h_change, reorder,
                       warnings, ignore_unknown_links);
    }, py::arg("st"), py::arg("monlib"), py::arg("model_index")=0,
       py::arg("h_change")=HydrogenChange::NoChange, py::arg("reorder")=false,
       py::arg("warnings")=py::none(), py::arg("ignore_unknown_links")=false);
}
