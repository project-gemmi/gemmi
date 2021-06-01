// Copyright 2018 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped
#include "gemmi/tostr.hpp"
#include "common.h"
#include <pybind11/stl.h>

#define GEMMI_READ_MAP_IMPLEMENTATION
#include "gemmi/read_map.hpp"  // defines read_ccp4_map, read_ccp4_mask

namespace py = pybind11;
using namespace gemmi;

template<typename T>
py::class_<T> add_ccp4_common(py::module& m, const char* name) {
  using Map = Ccp4<T>;
  return py::class_<Map, Ccp4Base>(m, name)
    .def(py::init<>())
    .def_readwrite("grid", &Map::grid)
    .def("update_ccp4_header", &Map::update_ccp4_header,
         py::arg("mode")=-1, py::arg("update_stats")=true)
    .def("full_cell", &Map::full_cell)
    .def("write_ccp4_map", &Map::write_ccp4_map, py::arg("filename"))
    .def("set_extent", &Map::set_extent)
    .def("__repr__", [=](const Map& self) {
        const SpaceGroup* sg = self.grid.spacegroup;
        return tostr("<gemmi.", name, " with grid ",
                     self.grid.nu, 'x', self.grid.nv, 'x', self.grid.nw,
                     " in SG #", sg ? std::to_string(sg->ccp4) : "?", '>');
    });
}

void add_ccp4(py::module& m) {
  py::class_<Ccp4Base>(m, "Ccp4Base")
    .def("header_i32", &Ccp4Base::header_i32)
    .def("header_float", &Ccp4Base::header_float)
    .def("header_str", &Ccp4Base::header_str)
    .def("set_header_i32", &Ccp4Base::set_header_i32)
    .def("set_header_float", &Ccp4Base::set_header_float)
    .def("set_header_str", &Ccp4Base::set_header_str)
    .def("axis_positions", &Ccp4Base::axis_positions)
    .def("get_extent", &Ccp4Base::get_extent)
    ;

  add_ccp4_common<float>(m, "Ccp4Map")
    .def("setup", [](Ccp4<float>& self, float default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=NAN);
  add_ccp4_common<int8_t>(m, "Ccp4Mask")
    .def("setup", [](Ccp4<int8_t>& self, int8_t default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=-1);
  m.def("read_ccp4_map", &read_ccp4_map,
        py::arg("path"), py::arg("setup")=false, py::return_value_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", &read_ccp4_mask,
        py::arg("path"), py::arg("setup")=false, py::return_value_policy::move,
        "Reads a CCP4 file, mode 0 (int8_t data, usually 0/1 masks).");
}
