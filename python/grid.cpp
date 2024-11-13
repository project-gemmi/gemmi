// Copyright 2018 Global Phasing Ltd.

#include <complex>

// for symmetrize_min and symmetrize_max
bool operator<(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) < std::norm(b);
}
bool operator>(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) > std::norm(b);
}

#include "common.h"
#include "array.h"
#include "make_iterator.h"
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for find_blobs_by_flood_fill, ...

#include "gemmi/grid.hpp"
#include "gemmi/floodfill.hpp"  // for flood_fill_above
#include "gemmi/solmask.hpp"  // for SolventMasker, mask_points_in_constant_radius
#include "gemmi/blob.hpp"     // for Blob, find_blobs_by_flood_fill
#include "gemmi/asumask.hpp"  // for MaskedGrid

using namespace gemmi;

namespace {

template<typename T>
auto grid_to_array(GridMeta& g, std::vector<T>& data) {
  // should we take AxisOrder into account here?
  return nb::ndarray<nb::numpy, T>(data.data(),
                                   {(size_t)g.nu, (size_t)g.nv, (size_t)g.nw},
                                   nb::handle(),
                                   {1, g.nu, g.nu * g.nv});
}

template<typename T>
nb::class_<GridBase<T>, GridMeta> add_grid_base(nb::module_& m, const char* name) {
  using GrBase = GridBase<T>;
  using GrPoint = typename GridBase<T>::Point;

  nb::class_<GrBase, GridMeta> grid_base(m, name);
  nb::class_<GrPoint>(grid_base, "Point")
    .def_ro("u", &GrPoint::u)
    .def_ro("v", &GrPoint::v)
    .def_ro("w", &GrPoint::w)
    .def_prop_rw("value",
                  [](const GrPoint& self) { return *self.value; },
                  [](GrPoint& self, T x) { *self.value = x; })
    .def("__repr__", [=](const GrPoint& self) {
        return nb::str("<gemmi.{}.Point ({}, {}, {}) -> {}>")
            .format(name, self.u, self.v, self.w, self.value);
    });

  auto to_array = [](GrBase& gr) { return grid_to_array(gr, gr.data); };
  grid_base
    .def_prop_ro("array", to_array, nb::rv_policy::reference_internal)
    .def("__array__", [](nb::handle_t<GrBase>& h, nb::handle dtype, nb::handle copy) {
        return handle_numpy_array_args(h.attr("array"), dtype, copy);
    }, nb::arg("dtype")=nb::none(), nb::arg("copy")=nb::none())
    .def("point_to_index", &GrBase::point_to_index)
    .def("index_to_point", &GrBase::index_to_point)
    .def("fill", &GrBase::fill, nb::arg("value"))
    .def("sum", &GrBase::sum)
    .def("__iter__", [](GrBase& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>())
    ;
  return grid_base;
}

template<typename T>
nb::class_<Grid<T>, GridBase<T>> add_grid_common(nb::module_& m, const std::string& name) {
  using Gr = Grid<T>;
  using GrPoint = typename GridBase<T>::Point;
  using Masked = MaskedGrid<T>;
  nb::class_<Gr, GridBase<T>> grid(m, name.c_str());
  nb::class_<Masked> masked_grid (m, ("Masked" + name).c_str());

  grid
    .def(nb::init<>())
    .def("__init__", [](Gr* grid, int nx, int ny, int nz) {
      new(grid) Gr();
      grid->set_size(nx, ny, nz);
    }, nb::arg("nx"), nb::arg("ny"), nb::arg("nz"))
    .def("__init__", [](Gr* grid, const nb::ndarray<nb::numpy, T, nb::ndim<3>>& arr,
                        const UnitCell *cell, const SpaceGroup* sg) {
      new(grid) Gr();
      auto r = arr.view();
      grid->set_size((int)r.shape(0), (int)r.shape(1), (int)r.shape(2));
      for (size_t k = 0; k < r.shape(2); ++k)
        for (size_t j = 0; j < r.shape(1); ++j)
          for (size_t i = 0; i < r.shape(0); ++i)
            grid->data[grid->index_q(i, j, k)] = r(i, j, k);
      if (cell)
        grid->set_unit_cell(*cell);
      if (sg)
        grid->spacegroup = sg;
    }, nb::arg().noconvert(), nb::arg("cell")=nb::none(), nb::arg("spacegroup")=nb::none())
    .def_prop_ro("spacing", [](const Gr& self) {
        return nb::make_tuple(self.spacing[0], self.spacing[1], self.spacing[2]);
    })
    .def("set_size", &Gr::set_size)
    .def("set_size_from_spacing", &Gr::set_size_from_spacing,
         nb::arg("spacing"), nb::arg("rounding"))
    .def("get_value", &Gr::get_value)
    .def("set_value", &Gr::set_value)
    .def("get_point", &Gr::get_point)
    .def("get_nearest_point", (GrPoint (Gr::*)(const Position&)) &Gr::get_nearest_point)
    .def("point_to_fractional", &Gr::point_to_fractional)
    .def("point_to_position", &Gr::point_to_position)
    .def("change_values", &Gr::change_values, nb::arg("old_value"), nb::arg("new_value"))
    .def("copy_metadata_from", &Gr::copy_metadata_from)
    .def("setup_from", &Gr::template setup_from<Structure>,
         nb::arg("st"), nb::arg("spacing")=0.)
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def("set_points_around", &Gr::set_points_around,
         nb::arg("position"), nb::arg("radius"), nb::arg("value"), nb::arg("use_pbc")=true)
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("symmetrize_abs_max", &Gr::symmetrize_abs_max)
    .def("symmetrize_sum", &Gr::symmetrize_sum)
    .def("masked_asu", &masked_asu<T>, nb::keep_alive<0, 1>())
    .def("mask_points_in_constant_radius", &mask_points_in_constant_radius<T>,
         nb::arg("model"), nb::arg("radius"), nb::arg("value"),
         nb::arg("ignore_hydrogen")=false, nb::arg("ignore_zero_occupancy_atoms")=false)
    .def("get_subarray",
         [](const Gr& self, std::array<int,3> start, std::array<int,3> shape) {
        auto arr = make_numpy_array<T>(
            {(size_t)shape[0], (size_t)shape[1], (size_t)shape[2]},
            {1, int64_t(shape[0]), int64_t(shape[0]*shape[1])});
        self.get_subarray(arr.data(), start, shape);
        return arr;
    }, nb::arg("start"), nb::arg("shape"))
    .def("set_subarray",
         [](Gr& self,
            const nb::ndarray<T, nb::ndim<3>, nb::f_contig, nb::device::cpu>& arr,
            std::array<int,3> start) {
        self.set_subarray(arr.data(), start,
                          {(int)arr.shape(0), (int)arr.shape(1), (int)arr.shape(2)});
    }, nb::arg("arr"), nb::arg("start"))
    .def("clone", [](const Gr& self) { return new Gr(self); })
    .def("__repr__", [=](const Gr& self) {
        return cat("<gemmi.", name, '(', self.nu, ", ", self.nv, ", ", self.nw, ")>");
    });

  masked_grid
    .def_ro("grid", &Masked::grid, nb::rv_policy::reference)
    .def_prop_ro("mask_array", [](Masked& self) {
        return grid_to_array(*self.grid, self.mask);
    }, nb::rv_policy::reference_internal)
    .def("__iter__", [](Masked& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>())
    ;
  return grid;
}

template<typename T>
void add_grid_interpolation(nb::class_<Grid<T>, GridBase<T>>& grid) {
  using Gr = Grid<T>;
  grid
    .def("interpolate_value",
         (T (Gr::*)(const Fractional&, int) const) &Gr::interpolate_value,
         nb::arg(), nb::arg("order")=1)
    .def("interpolate_value",
         (T (Gr::*)(const Position&, int) const) &Gr::interpolate_value,
         nb::arg(), nb::arg("order")=1)
    // deprecated, use interpolate_value(..., order=3)
    .def("tricubic_interpolation",
         (double (Gr::*)(const Fractional&) const) &Gr::tricubic_interpolation)
    // deprecated, use interpolate_value(..., order=3)
    .def("tricubic_interpolation",
         (double (Gr::*)(const Position&) const) &Gr::tricubic_interpolation)
    .def("tricubic_interpolation_der",
         (std::array<double,4> (Gr::*)(const Fractional&) const)
         &Gr::tricubic_interpolation_der)
    .def("interpolate_position_array",
         [](const Gr& self, const nb::ndarray<double, nb::shape<-1,3>, nb::device::cpu>& xyz,
            int order, const Transform* to_frac) {
        auto xyz_view = xyz.view();
        size_t len = xyz_view.shape(0);
        auto values = make_numpy_array<T>({len});
        T* data = values.data();
        const Transform& frac = to_frac ? *to_frac : self.unit_cell.frac;
        for (size_t i = 0; i < len; ++i) {
          Position pos(xyz_view(i, 0), xyz_view(i, 1), xyz_view(i, 2));
          Fractional fpos = Fractional(frac.apply(pos));
          data[i] = self.interpolate_value(fpos, order);
        }
        return values;
    }, nb::arg("xyz"), nb::arg("order")=1, nb::arg("to_frac")=nb::none())
    // The name of this function is not very descriptive, but since it's used
    // in a few external projects, renaming it isn't worth the hassle.
    // cf. interpolate_grid
    .def("interpolate_values",
         [](const Gr& self, nb::ndarray<nb::numpy, T, nb::ndim<3>>& arr,
            const Transform& tr, int order) {
        auto r = arr.view();
        for (size_t i = 0; i < r.shape(0); ++i)
          for (size_t j = 0; j < r.shape(1); ++j)
            for (size_t k = 0; k < r.shape(2); ++k) {
              Position pos(tr.apply(Vec3(i, j, k)));
              Fractional fpos = self.unit_cell.fractionalize(pos);
              r(i, j, k) = self.interpolate_value(fpos, order);
            }
    }, nb::arg().noconvert(), nb::arg(), nb::arg("order")=1)
    ;
}

}  // anonymous namespace

void add_grid(nb::module_& m) {
  nb::enum_<AxisOrder>(m, "AxisOrder")
    .value("Unknown", AxisOrder::Unknown)
    .value("XYZ", AxisOrder::XYZ)
    .value("ZYX", AxisOrder::ZYX);

  nb::enum_<GridSizeRounding>(m, "GridSizeRounding")
    .value("Nearest", GridSizeRounding::Nearest)
    .value("Up", GridSizeRounding::Up)
    .value("Down", GridSizeRounding::Down);

  nb::class_<GridMeta>(m, "GridMeta")
    .def_rw("spacegroup", &GridMeta::spacegroup)
    .def_rw("unit_cell", &GridMeta::unit_cell)
    .def_ro("nu", &GridMeta::nu, "size in the first (fastest-changing) dim")
    .def_ro("nv", &GridMeta::nv, "size in the second dimension")
    .def_ro("nw", &GridMeta::nw, "size in the third (slowest-changing) dim")
    .def_ro("axis_order", &GridMeta::axis_order)
    .def_prop_ro("point_count", &GridMeta::point_count)
    .def("get_position", &GridMeta::get_position)
    .def("get_fractional", &GridMeta::get_fractional)
    .def_prop_ro("shape", [](const GridMeta& self) {
      return nb::make_tuple(self.nu, self.nv, self.nw);
    });

  add_grid_base<int8_t>(m, "Int8GridBase")
    .def("get_nonzero_extent", &get_nonzero_extent<int8_t>);
  add_grid_common<int8_t>(m, "Int8Grid");

  add_grid_base<float>(m, "FloatGridBase")
    .def("calculate_correlation", &calculate_correlation<float>)
    .def("get_nonzero_extent", &get_nonzero_extent<float>)
    ;
  auto grid_float = add_grid_common<float>(m, "FloatGrid");
  add_grid_interpolation<float>(grid_float);
  grid_float.def("symmetrize_avg", &Grid<float>::symmetrize_avg);
  grid_float.def("normalize", &Grid<float>::normalize);
  grid_float.def("add_soft_edge_to_mask", &add_soft_edge_to_mask<float>);

  add_grid_base<std::complex<float>>(m, "ComplexGridBase");

  // from solmask.hpp
  nb::enum_<AtomicRadiiSet>(m, "AtomicRadiiSet")
    .value("VanDerWaals", AtomicRadiiSet::VanDerWaals)
    .value("Cctbx", AtomicRadiiSet::Cctbx)
    .value("Refmac", AtomicRadiiSet::Refmac)
    .value("Constant", AtomicRadiiSet::Constant);
  nb::class_<SolventMasker>(m, "SolventMasker")
    .def(nb::init<AtomicRadiiSet, double>(),
         nb::arg("choice"), nb::arg("constant_r")=0.)
    .def_rw("atomic_radii_set", &SolventMasker::atomic_radii_set)
    .def_rw("rprobe", &SolventMasker::rprobe)
    .def_rw("rshrink", &SolventMasker::rshrink)
    .def_rw("island_min_volume", &SolventMasker::island_min_volume)
    .def_rw("constant_r", &SolventMasker::constant_r)
    .def_rw("ignore_hydrogen", &SolventMasker::ignore_hydrogen)
    .def_rw("ignore_zero_occupancy_atoms", &SolventMasker::ignore_zero_occupancy_atoms)
    .def("set_radii", &SolventMasker::set_radii,
         nb::arg("choice"), nb::arg("constant_r")=0.)
    .def("put_mask_on_int8_grid", &SolventMasker::put_mask_on_grid<int8_t>)
    .def("put_mask_on_float_grid", &SolventMasker::put_mask_on_grid<float>)
    .def("set_to_zero", &SolventMasker::set_to_zero)
    ;
  m.def("interpolate_grid", &interpolate_grid<float>,
        nb::arg("dest"), nb::arg("src"), nb::arg("tr"), nb::arg("order")=1);
  m.def("interpolate_grid_of_aligned_model2", &interpolate_grid_of_aligned_model2<float>,
        nb::arg("dest"), nb::arg("src"), nb::arg("tr"),
        nb::arg("dest_model"), nb::arg("radius"), nb::arg("order")=1);


  // from blob.hpp
  nb::class_<Blob>(m, "Blob")
    .def_ro("volume", &Blob::volume)
    .def_ro("score", &Blob::score)
    .def_ro("peak_value", &Blob::peak_value)
    .def_ro("centroid", &Blob::centroid)
    .def_ro("peak_pos", &Blob::peak_pos)
    ;
  m.def("find_blobs_by_flood_fill",
        [](const Grid<float>& grid, double cutoff, double min_volume,
           double min_score, double min_peak, bool negate) {
       BlobCriteria crit;
       crit.cutoff = cutoff;
       crit.min_volume = min_volume;
       crit.min_score = min_score;
       crit.min_peak = min_peak;
       return find_blobs_by_flood_fill(grid, crit, negate);
    }, nb::arg("grid"), nb::arg("cutoff"), nb::arg("min_volume")=10.,
       nb::arg("min_score")=15., nb::arg("min_peak")=0., nb::arg("negate")=false);

  // from floodfill.hpp
  m.def("flood_fill_above", &flood_fill_above,
        nb::arg("grid"), nb::arg("seeds"), nb::arg("threshold"), nb::arg("negate")=false);

  // from asumask.hpp
  nb::class_<AsuBrick>(m, "AsuBrick")
    .def_ro("size", &AsuBrick::size)
    .def_ro("incl", &AsuBrick::incl)
    .def("get_extent", &AsuBrick::get_extent)
    .def("str", &AsuBrick::str)
    ;
  m.def("find_asu_brick", &find_asu_brick);
}
