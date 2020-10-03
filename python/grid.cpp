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
#include "gemmi/fourier.hpp"  // for get_f_phi_on_grid

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "common.h"  // for normalize_index
namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<NeighborSearch::Mark*>)

template<typename T>
std::string grid_dim_str(const GridBase<T>& g) {
  return std::to_string(g.nu) + ", " + std::to_string(g.nv) + ", " +
         std::to_string(g.nw);
}


template<typename T> void add_to_asu_data(T&) {}
template<> void add_to_asu_data(py::class_<FPhiGrid<float>::AsuData>& cl) {
  using AsuData = typename FPhiGrid<float>::AsuData;
  cl.def("get_size_for_hkl", &get_size_for_hkl<AsuData>,
         py::arg("min_size")=std::array<int,3>{{0,0,0}}, py::arg("sample_rate")=0.);
  cl.def("data_fits_into", &data_fits_into<AsuData>, py::arg("size"));
  cl.def("get_f_phi_on_grid", get_f_phi_on_grid<float, AsuData>,
         py::arg("size"), py::arg("half_l")=false, py::arg("order")=AxisOrder::XYZ);
  cl.def("transform_f_phi_to_map", &transform_f_phi_to_map2<float, AsuData>,
         py::arg("min_size")=std::array<int,3>{{0,0,0}},
         py::arg("sample_rate")=0.,
         py::arg("exact_size")=std::array<int,3>{{0,0,0}},
         py::arg("order")=AxisOrder::XYZ);
}

template<typename AsuData, typename F>
py::array_t<float> make_new_column(const AsuData& asu_data, F f) {
  if (!asu_data.unit_cell().is_crystal())
    throw std::runtime_error("AsuData: unknown unit cell parameters");
  py::array_t<float> arr(asu_data.size());
  py::buffer_info buf = arr.request();
  float* ptr = (float*) buf.ptr;
  for (size_t i = 0; i < asu_data.size(); ++i)
    ptr[i] = static_cast<float>(f(asu_data.unit_cell(), asu_data.get_hkl(i)));
  return arr;
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
                     self.w, ") -> ", +*self.value, '>');
    });
    ;

  pyMaskedGrid
    .def_readonly("grid", &Masked::grid, py::return_value_policy::reference)
    .def_readonly("mask", &Masked::mask)
    .def("__iter__", [](Masked& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>())
    ;

  using ReGr = ReciprocalGrid<T>;
  py::class_<ReGr, GrBase> regr(m, ("Reciprocal" + name).c_str());

  py::class_<typename ReGr::HklValue>(regr, "HklValue")
    .def_readonly("hkl", &ReGr::HklValue::hkl)
    .def_readonly("value", &ReGr::HklValue::value)
    .def("__repr__", [name](const typename ReGr::HklValue& self) {
        return tostr("<gemmi.Reciprocal", name, ".HklValue (",
                     self.hkl[0], ',', self.hkl[1], ',', self.hkl[2], ") ",
                     self.value, '>');
    });

  using AsuData = typename ReGr::AsuData;
  py::class_<AsuData> asu_data(regr, "AsuData");
  asu_data
    .def(py::init([](const UnitCell& unit_cell, const SpaceGroup* sg,
                     py::array_t<int> hkl, py::array_t<T> values) {
      auto h = hkl.unchecked<2>();
      if (h.shape(1) != 3)
        throw std::domain_error("error: the size of the second dimension != 3");
      auto v = values.template unchecked<1>();
      if (h.shape(0) != v.shape(0))
        throw std::domain_error("error: arrays have different lengths");
      AsuData* ret = new AsuData;
      ret->spacegroup_ = sg;
      ret->unit_cell_ = unit_cell;
      ret->unit_cell_.set_cell_images_from_spacegroup(ret->spacegroup_);
      ret->v.reserve(h.shape(0));
      for (ssize_t i = 0; i < h.shape(0); ++i)
        ret->v.push_back({{{h(i, 0), h(i, 1), h(i, 2)}}, v(i)});
      return ret;
    }), py::arg("cell"), py::arg("sg").none(false),
        py::arg("miller_array"), py::arg("value_array"))
    .def("__iter__", [](AsuData& self) { return py::make_iterator(self.v); },
         py::keep_alive<0, 1>())
    .def("__len__", [](const AsuData& self) { return self.v.size(); })
    .def("__getitem__", [](AsuData& self, int index) -> typename ReGr::HklValue& {
        return self.v.at(normalize_index(index, self.v));
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def_readwrite("spacegroup", &AsuData::spacegroup_)
    .def_readwrite("unit_cell", &AsuData::unit_cell_)
    .def_property_readonly("miller_array", [](const AsuData& self) {
      const typename ReGr::HklValue* data = self.v.data();
      py::array::ShapeContainer shape({(ssize_t)self.v.size(), 3});
      py::array::StridesContainer strides({(const char*)(data+1) - (const char*)data,
                                           sizeof(int)});
      return py::array_t<int>(shape, strides, &data->hkl[0], py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("value_array", [](const AsuData& self) {
      const typename ReGr::HklValue* data = self.v.data();
      ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<T>({(ssize_t)self.v.size()}, {stride},
                            &data->value, py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def("make_1_d2_array", [](const AsuData& asu_data) {
      return make_new_column(asu_data, [](const UnitCell& cell, Miller hkl) {
        return cell.calculate_1_d2(hkl);
      });
    })
    .def("make_d_array", [](const AsuData& asu_data) {
      return make_new_column(asu_data, [](const UnitCell& cell, Miller hkl) {
        return cell.calculate_d(hkl);
      });
    })
    .def("__repr__", [name](const AsuData& self) {
        return tostr("<gemmi.Reciprocal", name, ".AsuData with ", self.v.size(), " values>");
    });
    add_to_asu_data(asu_data);

  regr
    .def_readonly("half_l", &ReGr::half_l)
    .def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      ReGr* grid = new ReGr();
      grid->set_size_without_checking(nx, ny, nz);
      grid->axis_order = AxisOrder::XYZ;
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def(py::init([](py::array_t<T> arr, const UnitCell *cell, const SpaceGroup* sg) {
      auto r = arr.template unchecked<3>();
      ReGr* grid = new ReGr();
      grid->set_size_without_checking((int)r.shape(0), (int)r.shape(1), (int)r.shape(2));
      grid->axis_order = AxisOrder::XYZ;
      for (int k = 0; k < r.shape(2); ++k)
        for (int j = 0; j < r.shape(1); ++j)
          for (int i = 0; i < r.shape(0); ++i)
            grid->data[grid->index_q(i, j, k)] = r(i, j, k);
      if (cell)
        grid->unit_cell = *cell;
      if (sg)
        grid->spacegroup = sg;
      return grid;
    }), py::arg().noconvert(), py::arg("cell")=nullptr, py::arg("spacegroup")=nullptr)
    .def("get_value", &ReGr::get_value)
    .def("get_value_or_zero", &ReGr::get_value_or_zero)
    .def("set_value", &ReGr::set_value)
    .def("to_hkl", &ReGr::to_hkl)
    .def("prepare_asu_data", &ReGr::prepare_asu_data,
         py::arg("dmin")=0., py::arg("unblur")=0.,
         py::arg("with_000")=false, py::arg("with_sys_abs")=false)
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
    .def("to_site", (SmallStructure::Site& (NeighborSearch::Mark::*)(SmallStructure&) const)
                    &NeighborSearch::Mark::to_site)
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
    .def(py::init<SmallStructure&, double>(),
         py::arg("small_structure"), py::arg("max_radius"),
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
    .def("find_site_neighbors", &NeighborSearch::find_site_neighbors,
         py::arg("atom"), py::arg("min_dist")=0, py::arg("max_dist")=0,
         py::return_value_policy::move, py::keep_alive<0, 1>())
    .def("dist", &NeighborSearch::dist)
    .def("__repr__", [](const NeighborSearch& self) {
        return tostr("<gemmi.NeighborSearch with grid ",
                     grid_dim_str(self.grid), '>');
    });
}
