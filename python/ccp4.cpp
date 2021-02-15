// Copyright 2018 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped
#include "gemmi/tostr.hpp"
#include "common.h"

namespace py = pybind11;
using namespace gemmi;

template<typename T>
py::class_<T> add_ccp4_common(py::module& m, const char* name) {
  using Map = Ccp4<T>;
  return py::class_<Map>(m, name)
    .def(py::init<>())
    .def_readwrite("grid", &Map::grid)
    .def("header_i32", &Map::header_i32)
    .def("header_float", &Map::header_float)
    .def("header_str", &Map::header_str)
    .def("set_header_i32", &Map::set_header_i32)
    .def("set_header_float", &Map::set_header_float)
    .def("set_header_str", &Map::set_header_str)
    .def("update_ccp4_header", &Map::update_ccp4_header,
         py::arg("mode")=-1, py::arg("update_stats")=true)
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
  add_ccp4_common<float>(m, "Ccp4Map")
    .def("setup", [](Ccp4<float>& self, float default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=NAN);
  add_ccp4_common<int8_t>(m, "Ccp4Mask")
    .def("setup", [](Ccp4<int8_t>& self, int8_t default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=-1);
  m.def("read_ccp4_map", [](const std::string& path) {
          Ccp4<float> ccp4;
          ccp4.read_ccp4(MaybeGzipped(path));
          return ccp4;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", [](const std::string& path) {
          Ccp4<int8_t> ccp4;
          ccp4.read_ccp4(MaybeGzipped(path));
          return ccp4;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 0 (int8_t data, usually 0/1 masks).");
}
