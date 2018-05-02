// Copyright 2018 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/subcells.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

template<typename T>
std::string grid_dim_str(const Grid<T>& g) {
  return std::to_string(g.nu) + ", " + std::to_string(g.nv) + ", " +
         std::to_string(g.nw);
}

template<typename T>
void add_grid(py::module& m, const char* name) {
  using Gr = Grid<T>;
  py::class_<Gr>(m, name)
    .def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      Gr grid;
      grid.set_size(nx, ny, nz);
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def_readonly("nu", &Gr::nu, "size in the first (fastest-changing) dim")
    .def_readonly("nv", &Gr::nv, "size in the second dimension")
    .def_readonly("nw", &Gr::nw, "size in the third (slowest-changing) dim")
    .def("get_value", &Gr::get_value_s)
    .def("set_value", &Gr::set_value_s)
    .def_readwrite("space_group", &Gr::space_group)
    .def_readonly("unit_cell", &Gr::unit_cell)
    .def("set_unit_cell", (void (Gr::*)(const UnitCell&)) &Gr::set_unit_cell)
    .def("set_points_around", &Gr::set_points_around)
    .def("symmetrize_min", &Gr::symmetrize_min)
    .def("symmetrize_max", &Gr::symmetrize_max)
    .def("fill", &Gr::fill, py::arg("value"))
    .def("__iter__", [](const Gr& self) {
        return py::make_iterator(self.data);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [=](const Gr& self) {
        return "<gemmi." + std::string(name) + "(" + grid_dim_str(self) + ")>";
    });
}

template<typename T>
py::class_<T> add_ccp4(py::module& m, const char* name) {
  using Map = Ccp4<T>;
  return py::class_<Map>(m, name)
    .def(py::init<>())
    .def_readonly("grid", &Map::grid)
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
        const SpaceGroup* sg = self.grid.space_group;
        std::string sg_str = sg ?  std::to_string(sg->ccp4) : "?";
        return "<gemmi." + std::string(name) + " with grid (" +
               grid_dim_str(self.grid) + ") in SG #" + sg_str + ">";
    });
}

void add_grid(py::module& m) {
  add_grid<float>(m, "FloatGrid");
  add_grid<int8_t>(m, "Int8Grid");
  add_ccp4<float>(m, "Ccp4Map")
    .def("setup", [](Ccp4<float>& self) { self.setup(GridSetup::Full, NAN); });
  add_ccp4<int8_t>(m, "Ccp4Mask")
    .def("setup", [](Ccp4<int8_t>& self) { self.setup(GridSetup::Full, -1); });

  m.def("read_ccp4_map", [](const std::string& path) {
          Ccp4<float> grid;
          grid.read_ccp4_map(path);
          return grid;
        }, py::arg("path"), py::return_value_policy::move,
        "Reads a CCP4 file, mode 2 (floating-point data).");
  m.def("read_ccp4_mask", [](const std::string& path) {
          Ccp4<int8_t> grid;
          grid.read_ccp4_map(path);
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
        return "<gemmi.SubCells.Mark " +
               std::string(element_name(self.element)) + " of atom " +
               std::to_string(self.chain_idx) + "/" +
               std::to_string(self.residue_idx) + "/" +
               std::to_string(self.atom_idx) + ">";
    });
  subcells
    .def(py::init<const Model&, const UnitCell&, double>(),
         py::arg("model"), py::arg("cell"), py::arg("max_radius"))
    .def("find_atoms", &SubCells::find_atoms,
         py::arg("pos"), py::arg("alt"), py::arg("radius"))
    .def("dist", &SubCells::dist)
    .def("__repr__", [](const SubCells& self) {
        return "<gemmi.SubCells with grid " + grid_dim_str(self.grid) + ">";
    });
}
