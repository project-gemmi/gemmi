// Copyright 2018 Global Phasing Ltd.

#include "gemmi/grid.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

template<typename Gr>
void add_grid_methods(py::class_<Gr>& cl) {
  cl.def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      Gr grid;
      grid.set_size(nx, ny, nz);
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def_readonly("nu", &Gr::nu, "size in the first (fastest-changing) dim")
    .def_readonly("nv", &Gr::nv, "size in the second dimension")
    .def_readonly("nw", &Gr::nw, "size in the third (slowest-changing) dim")
    .def("header_i32", &Gr::header_i32)
    .def("header_float", &Gr::header_float)
    .def("header_str", &Gr::header_str)
    .def("set_header_i32", &Gr::set_header_i32)
    .def("set_header_float", &Gr::set_header_float)
    .def("set_header_str", &Gr::set_header_str)
    .def("get_value", &Gr::get_value_s)
    .def("set_value", &Gr::set_value_s)
    .def_readwrite("space_group", &Gr::space_group)
    .def_readonly("unit_cell", &Gr::unit_cell)
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def("set_points_around", &Gr::set_points_around)
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("fill", &Gr::fill, py::arg("value"))
    .def("update_ccp4_header", &Gr::update_ccp4_header,
         py::arg("mode"), py::arg("update_stats"))
    .def("write_ccp4_map", &Gr::write_ccp4_map, py::arg("filename"))
    .def("__iter__", [](const Gr& self) {
        return py::make_iterator(self.data);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Gr& self) {
        return "<gemmi.FloatGrid(" + std::to_string(self.nu) + ", "
                                   + std::to_string(self.nv) + ", "
                                   + std::to_string(self.nw) + ")>";
    });
}

void add_grid(py::module& m) {
  py::class_<Grid<float>> float_grid(m, "FloatGrid");
  add_grid_methods(float_grid);
  float_grid
    .def("setup", [](Grid<float>& self) { self.setup(GridSetup::Full, NAN); });

  py::class_<Grid<int8_t>> mask_grid(m, "MaskGrid");
  add_grid_methods(mask_grid);
  mask_grid
    .def("setup", [](Grid<int8_t>& self) { self.setup(GridSetup::Full, -1); });

  m.def("read_ccp4_map", [](const std::string& path) {
          Grid<float> grid;
          grid.read_ccp4_map(path);
          return grid;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", [](const std::string& path) {
          Grid<int8_t> grid;
          grid.read_ccp4_map(path);
          return grid;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 0 (int8_t data, usually 0/1 masks).");
}
