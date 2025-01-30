// Copyright 2018 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/util.hpp"  // for cat
#include "common.h"
#include <nanobind/stl/string.h>
#include <nanobind/stl/array.h>  // for Ccp4Base::axis_positions

using namespace gemmi;

namespace {

template<typename T>
auto add_ccp4_common(nb::module_& m, const char* name) {
  using Map = Ccp4<T>;
  return nb::class_<Map, Ccp4Base>(m, name)
    .def(nb::init<>())
    .def_rw("grid", &Map::grid)
    .def("setup", &Map::setup,
         nb::arg("default_value"), nb::arg("mode")=MapSetup::Full)
    .def("update_ccp4_header", &Map::update_ccp4_header,
         nb::arg("mode")=-1, nb::arg("update_stats")=true)
    .def("full_cell", &Map::full_cell)
    .def("write_ccp4_map", &Map::write_ccp4_map, nb::arg("filename"))
    .def("set_extent", &Map::set_extent)
    .def("__repr__", [=](const Map& self) {
        const SpaceGroup* sg = self.grid.spacegroup;
        return cat("<gemmi.", name, " with grid ",
                   self.grid.nu, 'x', self.grid.nv, 'x', self.grid.nw,
                   " in SG #", sg ? std::to_string(sg->ccp4) : "?", '>');
    });
}

}  // anonymous namespace

void add_ccp4(nb::module_& m) {
  nb::enum_<MapSetup>(m, "MapSetup")
    .value("Full", MapSetup::Full)
    .value("NoSymmetry", MapSetup::NoSymmetry)
    .value("ReorderOnly", MapSetup::ReorderOnly);

  nb::class_<Ccp4Base>(m, "Ccp4Base")
    .def("header_i32", &Ccp4Base::header_i32)
    .def("header_float", &Ccp4Base::header_float)
    .def("header_str", &Ccp4Base::header_str, nb::arg("w"), nb::arg("len")=80)
    .def("set_header_i32", &Ccp4Base::set_header_i32)
    .def("set_header_float", &Ccp4Base::set_header_float)
    .def("set_header_str", &Ccp4Base::set_header_str)
    .def("axis_positions", &Ccp4Base::axis_positions)
    .def("get_extent", &Ccp4Base::get_extent)
    .def("has_skew_transformation", &Ccp4Base::has_skew_transformation)
    .def("get_skew_transformation", &Ccp4Base::get_skew_transformation)
    ;

  add_ccp4_common<float>(m, "Ccp4Map");
  add_ccp4_common<int8_t>(m, "Ccp4Mask");
  m.def("read_ccp4_map", &read_ccp4_map,
        nb::arg("path"), nb::arg("setup")=false, nb::rv_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", &read_ccp4_mask,
        nb::arg("path"), nb::arg("setup")=false, nb::rv_policy::move,
        "Reads a CCP4 file, mode 0 (int8_t data, usually 0/1 masks).");
  m.def("read_ccp4_header", &read_ccp4_header,
        nb::arg("path"), nb::rv_policy::move);
}
