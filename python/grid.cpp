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
#include "gemmi/subcells.hpp"
#include "gemmi/tostr.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<SubCells::Mark*>)

template<typename T>
std::string grid_dim_str(const Grid<T>& g) {
  return std::to_string(g.nu) + ", " + std::to_string(g.nv) + ", " +
         std::to_string(g.nw);
}

template<typename T>
void add_grid(py::module& m, const char* name) {
  using Gr = Grid<T>;
  py::class_<Gr> gr(m, name, py::buffer_protocol());
  gr
    .def_buffer([](Gr &g) {
      return py::buffer_info(g.data.data(),
                             {g.nu, g.nv, g.nw},       // dimensions
                             {sizeof(T),               // strides
                              sizeof(T) * g.nu,
                              sizeof(T) * g.nu * g.nv});
    })
    .def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      Gr grid;
      grid.set_size(nx, ny, nz);
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def_readonly("nu", &Gr::nu, "size in the first (fastest-changing) dim")
    .def_readonly("nv", &Gr::nv, "size in the second dimension")
    .def_readonly("nw", &Gr::nw, "size in the third (slowest-changing) dim")
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
    .def_readwrite("spacegroup", &Gr::spacegroup)
    .def_readwrite("unit_cell", &Gr::unit_cell)
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def_readonly("full_canonical", &Gr::full_canonical)
    .def_property_readonly("point_count", &Gr::point_count)
    .def("set_points_around", &Gr::set_points_around,
         py::arg("position"), py::arg("radius"), py::arg("value"))
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("fill", &Gr::fill, py::arg("value"))
    .def("sum", &Gr::sum)
    .def("asu", &Gr::asu)
    .def("__iter__", [](Gr& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    .def("__repr__", [=](const Gr& self) {
        return tostr("<gemmi.", name, '(', grid_dim_str(self), ")>");
    });

  using GrPoint = typename Gr::Point;
  py::class_<GrPoint>(gr, "Point")
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

  using Masked = MaskedGrid<T>;
  py::class_<Masked>(m, ("Masked" + std::string(name)).c_str())
    .def_readonly("grid", &Masked::grid, py::return_value_policy::reference)
    .def_readonly("mask", &Masked::mask)
    .def("__iter__", [](Masked& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;
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
  add_grid<float>(m, "FloatGrid");
  add_grid<int8_t>(m, "Int8Grid");
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

  py::class_<SubCells> subcells(m, "SubCells");
  py::class_<SubCells::Mark>(subcells, "Mark")
    .def_readonly("x", &SubCells::Mark::x)
    .def_readonly("y", &SubCells::Mark::y)
    .def_readonly("z", &SubCells::Mark::z)
    .def_readonly("altloc", &SubCells::Mark::altloc)
    .def_readonly("element", &SubCells::Mark::element)
    .def_readonly("image_idx", &SubCells::Mark::image_idx)
    .def_readonly("chain_idx", &SubCells::Mark::chain_idx)
    .def_readonly("residue_idx", &SubCells::Mark::residue_idx)
    .def_readonly("atom_idx", &SubCells::Mark::atom_idx)
    .def("pos", &SubCells::Mark::pos)
    .def("to_cra",
         (CRA (SubCells::Mark::*)(Model&) const) &SubCells::Mark::to_cra)
    .def("__repr__", [](const SubCells::Mark& self) {
        return tostr("<gemmi.SubCells.Mark ", element_name(self.element),
                     " of atom ", self.chain_idx, '/', self.residue_idx, '/',
                     self.atom_idx, '>');
    });
  py::bind_vector<std::vector<SubCells::Mark*>>(m, "VectorSubCellsMarkPtr");
  subcells
    .def(py::init<const Model&, const UnitCell&, double>(),
         py::arg("model"), py::arg("cell"), py::arg("max_radius"),
         py::keep_alive<1, 2>())
    .def("populate", &SubCells::populate, py::arg("include_h")=true,
         "Usually run after constructing SubCells.")
    .def("add_atom", &SubCells::add_atom,
         py::arg("atom"), py::arg("n_ch"), py::arg("n_res"), py::arg("n_atom"),
         "Lower-level alternative to populate()")
    .def("find_atoms", &SubCells::find_atoms,
         py::arg("pos"), py::arg("alt"), py::arg("radius"),
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("find_neighbors", &SubCells::find_neighbors,
         py::arg("atom"), py::arg("min_dist"), py::arg("max_dist"),
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("dist", &SubCells::dist)
    .def("__repr__", [](const SubCells& self) {
        return tostr("<gemmi.SubCells with grid ",grid_dim_str(self.grid),'>');
    });
}
