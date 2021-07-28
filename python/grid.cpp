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
#include "gemmi/solmask.hpp"  // for SolventMasker, mask_points_in_constant_radius
#include "gemmi/blob.hpp"     // for Blob, find_blobs_by_flood_fill
#include "gemmi/tostr.hpp"

#include "common.h"
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

template<typename T>
void add_grid_base(py::module& m, const char* name) {
  using GrBase = GridBase<T>;
  using GrPoint = typename GridBase<T>::Point;

  py::class_<GrBase, GridMeta> grid_base(m, name, py::buffer_protocol());
  py::class_<GrPoint>(grid_base, "Point")
    .def_readonly("u", &GrPoint::u)
    .def_readonly("v", &GrPoint::v)
    .def_readonly("w", &GrPoint::w)
    .def_property("value",
                  [](const GrPoint& self) { return *self.value; },
                  [](GrPoint& self, T x) { *self.value = x; })
    .def("__repr__", [=](const GrPoint& self) {
        return tostr("<gemmi.", name, ".Point (", self.u, ", ", self.v, ", ",
                     self.w, ") -> ", +*self.value, '>');
    });

  grid_base
    .def_buffer([](GrBase &g) {
      return py::buffer_info(g.data.data(),
                             {g.nu, g.nv, g.nw},       // dimensions
                             {sizeof(T),               // strides
                              sizeof(T) * g.nu,
                              sizeof(T) * g.nu * g.nv});
    })
    .def("point_to_index", &GrBase::point_to_index)
    .def("index_to_point", &GrBase::index_to_point)
    .def("fill", &GrBase::fill, py::arg("value"))
    .def("sum", &GrBase::sum)
    .def("__iter__", [](GrBase& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;
}

template<typename T>
void add_grid(py::module& m, const std::string& name) {
  using Gr = Grid<T>;
  using GrPoint = typename GridBase<T>::Point;
  using Masked = MaskedGrid<T>;
  py::class_<Masked> masked_grid (m, ("Masked" + name).c_str());

  py::class_<Gr, GridBase<T>>(m, name.c_str())
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
    .def("set_size", &Gr::set_size)
    .def("get_value", &Gr::get_value)
    .def("set_value", &Gr::set_value)
    .def("get_point", &Gr::get_point)
    .def("get_position", &Gr::get_position)
    .def("get_fractional", &Gr::get_fractional)
    .def("get_nearest_point", (GrPoint (Gr::*)(const Position&)) &Gr::get_nearest_point)
    .def("point_to_fractional", &Gr::point_to_fractional)
    .def("point_to_position", &Gr::point_to_position)
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
    .def("tricubic_interpolation",
         (double (Gr::*)(const Fractional&) const) &Gr::tricubic_interpolation)
    .def("tricubic_interpolation",
         (double (Gr::*)(const Position&) const) &Gr::tricubic_interpolation)
    .def("tricubic_interpolation_der",
         (std::array<double,4> (Gr::*)(const Fractional&) const)
         &Gr::tricubic_interpolation_der)
    .def("copy_metadata_from", &Gr::copy_metadata_from)
    .def("setup_from", &Gr::template setup_from<Structure>,
         py::arg("st"), py::arg("spacing"))
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def("set_points_around", &Gr::set_points_around,
         py::arg("position"), py::arg("radius"), py::arg("value"), py::arg("use_pbc")=true)
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("symmetrize_abs_max", &Gr::symmetrize_abs_max)
    .def("symmetrize_sum", &Gr::symmetrize_sum)
    .def("asu", &Gr::asu)
    .def("mask_points_in_constant_radius", &mask_points_in_constant_radius<T>,
         py::arg("model"), py::arg("radius"), py::arg("value"))
    .def("get_subarray",
         [](const Gr& self, int u0, int v0, int w0, int new_nu, int new_nv, int new_nw) {
        py::array_t<T> arr({new_nu, new_nv, new_nw},
                           {sizeof(T), sizeof(T)*new_nu, sizeof(T)*new_nu*new_nv});
        py::buffer_info buf = arr.request();
        T* ptr = (T*) buf.ptr;
        int idx = 0;
        for (int w = 0; w < new_nw; w++)
          for (int v = 0; v < new_nv; v++)
            for (int u = 0; u < new_nu; u++)
              ptr[idx++] = self.get_value(u0 + u, v0 + v, w0 + w);
        return arr;
    })
    .def("__repr__", [=](const Gr& self) {
        return tostr("<gemmi.", name, '(', self.nu, ", ", self.nv, ", ", self.nw, ")>");
    });

  masked_grid
    .def_readonly("grid", &Masked::grid, py::return_value_policy::reference)
    .def_readonly("mask", &Masked::mask)
    .def("__iter__", [](Masked& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;
}

// for pandda-2
Grid<float> interpolate_positions(const Grid<float>& moving_map,
                                  Grid<float> interpolated_map,
                                  const std::vector<std::array<int,3>>& points,
                                  const std::vector<Fractional>& positions) {
  py::gil_scoped_release release;  // Release gil for threading support
  if (points.size() != positions.size())
    fail("interpolate_positions(): list sizes differ");
  for (std::size_t i = 0; i < positions.size(); i++)
    interpolated_map.set_value(points[i][0],
                               points[i][1],
                               points[i][2],
                               moving_map.interpolate_value(positions[i]));
  return interpolated_map;
}

void add_grid(py::module& m) {
  py::enum_<AxisOrder>(m, "AxisOrder")
    .value("XYZ", AxisOrder::XYZ)
    .value("ZYX", AxisOrder::ZYX);

  py::class_<GridMeta>(m, "GridMeta")
    .def_readwrite("spacegroup", &GridMeta::spacegroup)
    .def_readwrite("unit_cell", &GridMeta::unit_cell)
    .def_readonly("nu", &GridMeta::nu, "size in the first (fastest-changing) dim")
    .def_readonly("nv", &GridMeta::nv, "size in the second dimension")
    .def_readonly("nw", &GridMeta::nw, "size in the third (slowest-changing) dim")
    .def_readonly("axis_order", &GridMeta::axis_order)
    .def_property_readonly("point_count", &GridMeta::point_count)
    .def_property_readonly("shape", [](const GridMeta& self){
      return py::make_tuple(self.nu, self.nv, self.nw);
    });

  add_grid_base<int8_t>(m, "Int8GridBase");
  add_grid<int8_t>(m, "Int8Grid");
  add_grid_base<float>(m, "FloatGridBase");
  add_grid<float>(m, "FloatGrid");
  add_grid_base<std::complex<float>>(m, "ComplexGridBase");
  m.def("interpolate_positions", &interpolate_positions);

  // from solmask.hpp
  py::enum_<AtomicRadiiSet>(m, "AtomicRadiiSet")
    .value("VanDerWaals", AtomicRadiiSet::VanDerWaals)
    .value("Cctbx", AtomicRadiiSet::Cctbx)
    .value("Refmac", AtomicRadiiSet::Refmac)
    .value("Constant", AtomicRadiiSet::Constant);
  py::class_<SolventMasker>(m, "SolventMasker")
    .def(py::init<AtomicRadiiSet, double>(),
         py::arg("choice"), py::arg("constant_r")=0.)
    .def_readwrite("atomic_radii_set", &SolventMasker::atomic_radii_set)
    .def_readwrite("rprobe", &SolventMasker::rprobe)
    .def_readwrite("rshrink", &SolventMasker::rshrink)
    .def_readwrite("island_min_volume", &SolventMasker::island_min_volume)
    .def_readwrite("constant_r", &SolventMasker::constant_r)
    .def("set_radii", &SolventMasker::set_radii,
         py::arg("choice"), py::arg("constant_r")=0.)
    .def("put_mask_on_int8_grid", &SolventMasker::put_mask_on_grid<int8_t>)
    .def("put_mask_on_float_grid", &SolventMasker::put_mask_on_grid<float>)
    .def("set_to_zero", &SolventMasker::set_to_zero)
    ;

  py::class_<Blob>(m, "Blob")
    .def_readonly("volume", &Blob::volume)
    .def_readonly("score", &Blob::score)
    .def_readonly("peak_value", &Blob::peak_value)
    .def_readonly("centroid", &Blob::centroid)
    .def_readonly("peak_pos", &Blob::peak_pos)
    ;
  m.def("find_blobs_by_flood_fill",
        [](const Grid<float>& grid, double cutoff,
           double min_volume, double min_score, double min_peak) {
       BlobCriteria crit;
       crit.cutoff = cutoff;
       crit.min_volume = min_volume;
       crit.min_score = min_score;
       crit.min_peak = min_peak;
       return find_blobs_by_flood_fill(grid, crit);
    }, py::arg("grid"), py::arg("cutoff")=true,
       py::arg("min_volume")=10., py::arg("min_score")=15, py::arg("min_peak")=0.);
}
