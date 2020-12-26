// Copyright 2018 Global Phasing Ltd.

#include <complex>

// for symmetrize_min and symmetrize_max
bool operator<(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) < std::norm(b);
}
bool operator>(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) > std::norm(b);
}

#include "gemmi/grid.hpp"
#include "gemmi/dencalc.hpp"  // for mask_points_in_constant_radius
#include "gemmi/tostr.hpp"

#include "common.h"
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

template<typename T>
void add_grid(py::module& m, const std::string& name) {
  using GrBase = GridBase<T>;
  using Gr = Grid<T>;
  using Masked = MaskedGrid<T>;
  using GrPoint = typename Gr::Point;

  py::class_<GrBase> pyGridBase(m, (name + "Base").c_str(), py::buffer_protocol());
  py::class_<Gr, GrBase> gr(m, name.c_str());
  py::class_<Masked> pyMaskedGrid (m, ("Masked" + name).c_str());
  py::class_<GrPoint> pyGrPoint(gr, "Point");

  pyGridBase
    .def_buffer([](GrBase &g) {
      return py::buffer_info(g.data.data(),
                             {g.nu, g.nv, g.nw},       // dimensions
                             {sizeof(T),               // strides
                              sizeof(T) * g.nu,
                              sizeof(T) * g.nu * g.nv});
    })
    .def_readonly("nu", &GrBase::nu, "size in the first (fastest-changing) dim")
    .def_readonly("nv", &GrBase::nv, "size in the second dimension")
    .def_readonly("nw", &GrBase::nw, "size in the third (slowest-changing) dim")
    .def_property_readonly("shape", [](const GrBase& self){
      return py::make_tuple(self.nu, self.nv, self.nw);
    })
    .def_readwrite("spacegroup", &GrBase::spacegroup)
    .def_readwrite("unit_cell", &GrBase::unit_cell)
    .def_readonly("axis_order", &GrBase::axis_order)
    .def_property_readonly("point_count", &GrBase::point_count)
    .def("fill", &GrBase::fill, py::arg("value"))
    .def("sum", &GrBase::sum)
    .def("__iter__", [](GrBase& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;

  gr
    .def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      Gr* grid = new Gr();
      grid->set_size(nx, ny, nz);
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def(py::init([](py::array_t<T> arr, const UnitCell *cell, const SpaceGroup* sg) {
      auto r = arr.template unchecked<3>();
      Gr* grid = new Gr();
      grid->set_size((int)r.shape(0), (int)r.shape(1), (int)r.shape(2));
      for (int k = 0; k < r.shape(2); ++k)
        for (int j = 0; j < r.shape(1); ++j)
          for (int i = 0; i < r.shape(0); ++i)
            grid->data[grid->index_q(i, j, k)] = r(i, j, k);
      if (cell)
        grid->set_unit_cell(*cell);
      if (sg)
        grid->spacegroup = sg;
      return grid;
    }), py::arg().noconvert(), py::arg("cell")=nullptr, py::arg("spacegroup")=nullptr)
    .def_property_readonly("spacing", [](const Gr& self){
        return py::make_tuple(self.spacing[0], self.spacing[1], self.spacing[2]);
    })
    .def("get_value", &Gr::get_value)
    .def("set_value", &Gr::set_value)
    .def("get_point", &Gr::get_point)
    .def("point_to_fractional", &Gr::point_to_fractional)
    .def("point_to_position", &Gr::point_to_position)
    .def("point_to_position", [](const Gr& self, int u, int v, int w) {
      GrPoint p = {u, v, w, nullptr};
      return self.point_to_position(p);
    }, py::arg("u"), py::arg("v"), py::arg("w"))
    .def("point_to_index", &Gr::point_to_index)
    .def("interpolate_value",
         (T (Gr::*)(const Fractional&) const) &Gr::interpolate_value)
    .def("interpolate_value",
         (T (Gr::*)(const Position&) const) &Gr::interpolate_value)
    .def("interpolate_values",
         [](const Gr& self, py::array_t<T> arr, const Transform& tr) {
        auto r = arr.template mutable_unchecked<3>();
        for (int i = 0; i < r.shape(0); ++i)
          for (int j = 0; j < r.shape(1); ++j)
            for (int k = 0; k < r.shape(2); ++k) {
              Position pos(tr.apply(Vec3(i, j, k)));
              r(i, j, k) = self.interpolate_value(pos);
            }
    }, py::arg().noconvert(), py::arg())
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def("set_points_around", &Gr::set_points_around,
         py::arg("position"), py::arg("radius"), py::arg("value"))
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("asu", &Gr::asu)
    .def("mask_points_in_constant_radius", &mask_points_in_constant_radius<T>,
         py::arg("model"), py::arg("radius"), py::arg("value"))
    .def("__repr__", [=](const Gr& self) {
        return tostr("<gemmi.", name, '(', self.nu, ", ", self.nv, ", ", self.nw, ")>");
    });

  pyGrPoint
    .def_readonly("u", &GrPoint::u)
    .def_readonly("v", &GrPoint::v)
    .def_readonly("w", &GrPoint::w)
    .def_property("value",
                  [](const GrPoint& self) { return *self.value; },
                  [](GrPoint& self, T x) { *self.value = x; })
    .def("__repr__", [=](const GrPoint& self) {
        return tostr("<gemmi.", name, "Point (", self.u, ", ", self.v, ", ",
                     self.w, ") -> ", +*self.value, '>');
    });
    ;

  pyMaskedGrid
    .def_readonly("grid", &Masked::grid, py::return_value_policy::reference)
    .def_readonly("mask", &Masked::mask)
    .def("__iter__", [](Masked& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;
}

void add_grid(py::module& m) {
  py::enum_<AxisOrder>(m, "AxisOrder")
    .value("XYZ", AxisOrder::XYZ)
    .value("ZYX", AxisOrder::ZYX);
  add_grid<int8_t>(m, "Int8Grid");
  add_grid<float>(m, "FloatGrid");
  add_grid<std::complex<float>>(m, "ComplexGrid");
}
