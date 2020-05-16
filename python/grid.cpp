// Copyright 2018 Global Phasing Ltd.

#include <complex>

// for symmetrize_min and symmetrize_max
bool operator<(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) < std::norm(b);
}
bool operator>(const std::complex<float>& a, const std::complex<float>& b) {
    return std::norm(a) > std::norm(b);
}

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped
#include "gemmi/neighbor.hpp"
#include "gemmi/tostr.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<NeighborSearch::Mark*>)

template<typename T>
std::string grid_dim_str(const GridBase<T>& g) {
  return std::to_string(g.nu) + ", " + std::to_string(g.nv) + ", " +
         std::to_string(g.nw);
}

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
    .def("get_value", &Gr::get_value)
    .def("set_value", &Gr::set_value)
    .def("get_point", &Gr::get_point)
    .def("point_to_fractional", &Gr::point_to_fractional)
    .def("point_to_position", &Gr::point_to_position)
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
    .def("__repr__", [=](const Gr& self) {
        return tostr("<gemmi.", name, '(', grid_dim_str(self), ")>");
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
                     self.w, ") -> ", *self.value, '>');
    });
    ;

  pyMaskedGrid
    .def_readonly("grid", &Masked::grid, py::return_value_policy::reference)
    .def_readonly("mask", &Masked::mask)
    .def("__iter__", [](Masked& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;

  using ReGr = ReciprocalGrid<T>;
  py::class_<ReGr, GrBase>(m, ("Reciprocal" + name).c_str())
    .def_readonly("half_l", &ReGr::half_l)
    .def(py::init<>())
    .def("get_value", &ReGr::get_value)
    .def("get_value_or_zero", &ReGr::get_value_or_zero)
    .def("set_value", &ReGr::set_value)
    .def("to_hkl", &ReGr::to_hkl)
    .def("__repr__", [=](const ReGr& self) {
        return tostr("<gemmi.Reciprocal", name, '(', grid_dim_str(self), ")>");
    });
}

template<typename T>
py::class_<T> add_ccp4(py::module& m, const char* name) {
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
         py::arg("mode"), py::arg("update_stats"))
    .def("write_ccp4_map", &Map::write_ccp4_map, py::arg("filename"))
    .def("__repr__", [=](const Map& self) {
        const SpaceGroup* sg = self.grid.spacegroup;
        return tostr("<gemmi.", name, " with grid (", grid_dim_str(self.grid),
                     ") in SG #", sg ? std::to_string(sg->ccp4) : "?", '>');
    });
}

void add_grid(py::module& m) {
  py::enum_<AxisOrder> pyAxisOrder(m, "AxisOrder");

  add_grid<int8_t>(m, "Int8Grid");
  add_grid<float>(m, "FloatGrid");
  add_grid<std::complex<float>>(m, "ComplexGrid");
  add_ccp4<float>(m, "Ccp4Map")
    .def("setup", [](Ccp4<float>& self, float default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=NAN);
  add_ccp4<int8_t>(m, "Ccp4Mask")
    .def("setup", [](Ccp4<int8_t>& self, int8_t default_value) {
            self.setup(GridSetup::Full, default_value);
         }, py::arg("default_value")=-1);
  m.def("read_ccp4_map", [](const std::string& path) {
          Ccp4<float> grid;
          grid.read_ccp4(MaybeGzipped(path));
          return grid;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", [](const std::string& path) {
          Ccp4<int8_t> grid;
          grid.read_ccp4(MaybeGzipped(path));
          return grid;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 0 (int8_t data, usually 0/1 masks).");

  pyAxisOrder
    .value("XYZ", AxisOrder::XYZ)
    .value("ZYX", AxisOrder::ZYX);

  py::class_<NeighborSearch> neighbor_search(m, "NeighborSearch");
  py::class_<NeighborSearch::Mark>(neighbor_search, "Mark")
    .def_readonly("x", &NeighborSearch::Mark::x)
    .def_readonly("y", &NeighborSearch::Mark::y)
    .def_readonly("z", &NeighborSearch::Mark::z)
    .def_readonly("altloc", &NeighborSearch::Mark::altloc)
    .def_readonly("element", &NeighborSearch::Mark::element)
    .def_readonly("image_idx", &NeighborSearch::Mark::image_idx)
    .def_readonly("chain_idx", &NeighborSearch::Mark::chain_idx)
    .def_readonly("residue_idx", &NeighborSearch::Mark::residue_idx)
    .def_readonly("atom_idx", &NeighborSearch::Mark::atom_idx)
    .def("pos", &NeighborSearch::Mark::pos)
    .def("to_cra", (CRA (NeighborSearch::Mark::*)(Model&) const)
                   &NeighborSearch::Mark::to_cra)
    .def("__repr__", [](const NeighborSearch::Mark& self) {
        return tostr("<gemmi.NeighborSearch.Mark ", self.element.name(),
                     " of atom ", self.chain_idx, '/', self.residue_idx, '/',
                     self.atom_idx, '>');
    });
  py::bind_vector<std::vector<NeighborSearch::Mark*>>(m, "VectorMarkPtr");
  neighbor_search
    .def(py::init<Model&, const UnitCell&, double>(),
         py::arg("model"), py::arg("cell"), py::arg("max_radius")/*,
         py::keep_alive<1, 2>()*/)
    .def(py::init([](Structure& st, double max_radius, int model_index) {
      return new NeighborSearch(st.models.at(model_index), st.cell, max_radius);
    }), py::arg("st"), py::arg("max_radius"), py::arg("model_index")=0,
        py::keep_alive<1, 2>())
    .def("populate", &NeighborSearch::populate, py::arg("include_h")=true,
         "Usually run after constructing NeighborSearch.")
    .def("add_atom", &NeighborSearch::add_atom,
         py::arg("atom"), py::arg("n_ch"), py::arg("n_res"), py::arg("n_atom"),
         "Lower-level alternative to populate()")
    .def("find_atoms", &NeighborSearch::find_atoms,
         py::arg("pos"), py::arg("alt")='\0', py::arg("radius")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("find_neighbors", &NeighborSearch::find_neighbors,
         py::arg("atom"), py::arg("min_dist")=0, py::arg("max_dist")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("dist", &NeighborSearch::dist)
    .def("__repr__", [](const NeighborSearch& self) {
        return tostr("<gemmi.NeighborSearch with grid ",
                     grid_dim_str(self.grid), '>');
    });
}
