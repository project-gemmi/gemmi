// Copyright 2019 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/refln.hpp"
#include "gemmi/fourier.hpp"  // for get_size_for_hkl, get_f_phi_on_grid, ...
#include "gemmi/tostr.hpp"
#include "gemmi/fprime.hpp"
#include "gemmi/reciproc.hpp" // for count_reflections, make_miller_vector
#include "gemmi/cif2mtz.hpp"  // for CifToMtz
#include "gemmi/merge.hpp"  // for Intensities

#include "common.h"
#include "arrvec.h"  // py_array_from_vector
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<ReflnBlock>)

namespace gemmi {
  // operator<< is used by stl_bind for vector's __repr__
  inline std::ostream& operator<< (std::ostream& os, const ReflnBlock& rb) {
    os << "<gemmi.ReflnBlock " << rb.block.name << " with ";
    if (rb.default_loop)
      os << rb.default_loop->width() << " x " << rb.default_loop->length();
    else
      os << " no ";
    os << " loop>";
    return os;
  }
}

void add_hkl(py::module& m) {
  py::class_<ReflnBlock> pyReflnBlock(m, "ReflnBlock");
  py::bind_vector<std::vector<ReflnBlock>>(m, "ReflnBlocks");

  pyReflnBlock
    .def_readonly("block", &ReflnBlock::block)
    .def_readonly("entry_id", &ReflnBlock::entry_id)
    .def_readonly("cell", &ReflnBlock::cell)
    .def_readonly("spacegroup", &ReflnBlock::spacegroup)
    .def_readonly("wavelength", &ReflnBlock::wavelength)
    .def("column_labels", &ReflnBlock::column_labels)
    .def("make_int_array",
         [](ReflnBlock& self, const std::string& tag, int null) {
           return py_array_from_vector(self.make_vector(tag, null));
    }, py::arg("tag"), py::arg("null"))
    .def("make_float_array",
         [](ReflnBlock& self, const std::string& tag, double null) {
           return py_array_from_vector(self.make_vector(tag, null));
    }, py::arg("tag"), py::arg("null")=NAN)
    .def("make_float_array", &ReflnBlock::make_vector<double>,
         py::arg("tag"), py::arg("null")=NAN)
    .def("make_miller_array", [](ReflnBlock& self) {
        return py::array_t<int>(py::cast(self.make_miller_vector()));
    })
    .def("make_1_d2_array", [](ReflnBlock& self) {
        return py_array_from_vector(self.make_1_d2_vector());
    })
    .def("make_d_array", [](ReflnBlock& self) {
        return py_array_from_vector(self.make_d_vector());
    })
    .def("get_size_for_hkl",
         [](const ReflnBlock& self,
            std::array<int,3> min_size, double sample_rate) {
          return get_size_for_hkl(ReflnDataProxy(self), min_size, sample_rate);
    }, py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.)
    .def("data_fits_into", [](const ReflnBlock& self, std::array<int,3> size) {
        return data_fits_into(ReflnDataProxy(self), size);
    }, py::arg("size"))
    .def("get_f_phi_on_grid", [](const ReflnBlock& self,
                                 const std::string& f_col,
                                 const std::string& phi_col,
                                 std::array<int, 3> size,
                                 bool half_l, AxisOrder order) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        FPhiProxy<ReflnDataProxy> fphi(ReflnDataProxy{self}, f_idx, phi_idx);
        return get_f_phi_on_grid<float>(fphi, size, half_l, order);
    }, py::arg("f"), py::arg("phi"), py::arg("size"),
       py::arg("half_l")=false, py::arg("order")=AxisOrder::XYZ)
    .def("get_value_on_grid", [](const ReflnBlock& self,
                                 const std::string& column,
                                 std::array<int, 3> size,
                                 bool half_l, AxisOrder order) {
        size_t col_idx = self.get_column_index(column);
        return get_value_on_grid<float>(ReflnDataProxy(self), col_idx,
                                        size, half_l, order);
    }, py::arg("column"), py::arg("size"), py::arg("half_l")=false,
       py::arg("order")=AxisOrder::XYZ)
    .def("transform_f_phi_to_map", [](const ReflnBlock& self,
                                      const std::string& f_col,
                                      const std::string& phi_col,
                                      std::array<int, 3> min_size,
                                      std::array<int, 3> exact_size,
                                      double sample_rate,
                                      AxisOrder order) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        FPhiProxy<ReflnDataProxy> fphi(ReflnDataProxy{self}, f_idx, phi_idx);
        return transform_f_phi_to_map2<float>(fphi, min_size, sample_rate,
                                              exact_size, order);
    }, py::arg("f"), py::arg("phi"),
       py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("exact_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.,
       py::arg("order")=AxisOrder::XYZ)
    .def("get_float", &make_asu_data<float, ReflnBlock>,
         py::arg("col"), py::arg("as_is")=false)
    .def("get_int", &make_asu_data<int, ReflnBlock>,
         py::arg("col"), py::arg("as_is")=false)
    .def("get_f_phi", [](const ReflnBlock& self, const std::string& f_col,
                                                 const std::string& phi_col,
                                                 bool as_is) {
        return make_asu_data<std::complex<float>, 2>(self, {f_col, phi_col}, as_is);
    }, py::arg("f"), py::arg("phi"), py::arg("as_is")=false)
    .def("get_value_sigma", [](const ReflnBlock& self, const std::string& f_col,
                                                       const std::string& sigma_col,
                                                       bool as_is) {
        return make_asu_data<ValueSigma<float>, 2>(self, {f_col, sigma_col}, as_is);
    }, py::arg("f"), py::arg("sigma"), py::arg("as_is")=false)
    .def("is_unmerged", &ReflnBlock::is_unmerged)
    .def("use_unmerged", &ReflnBlock::use_unmerged)
    .def("__bool__", [](const ReflnBlock& self) { return self.ok(); })
    .def("__repr__", [](const ReflnBlock& self) { return tostr(self); })
    ;
  m.def("as_refln_blocks",
        [](cif::Document& d) { return as_refln_blocks(std::move(d.blocks)); });
  m.def("hkl_cif_as_refln_block", &hkl_cif_as_refln_block, py::arg("block"));
  m.def("transform_f_phi_grid_to_map", [](FPhiGrid<float> grid) {
          return transform_f_phi_grid_to_map<float>(std::move(grid));
        }, py::arg("grid"));
  m.def("transform_map_to_f_phi", &transform_map_to_f_phi<float>,
        py::arg("map"), py::arg("half_l")=false, py::arg("use_scale")=true);
  m.def("cromer_liberman", [](int z, double energy) {
          std::pair<double, double> r;
          r.first = cromer_liberman(z, energy, &r.second);
          return r;
        }, py::arg("z"), py::arg("energy"));
  m.def("count_reflections", &count_reflections,
        py::arg("cell"), py::arg("spacegroup"), py::arg("dmin"),
        py::arg("dmax")=0., py::arg("unique")=true);
  m.def("make_miller_array", [](const UnitCell& cell, const SpaceGroup* sg,
                                double dmin, double dmax, bool unique) {
          return py::array_t<int>(py::cast(
                      gemmi::make_miller_vector(cell, sg, dmin, dmax, unique)));
        }, py::arg("cell"), py::arg("spacegroup"), py::arg("dmin"),
           py::arg("dmax")=0., py::arg("unique")=true);

  py::class_<CifToMtz>(m, "CifToMtz")
    .def(py::init<>())
    .def_readwrite("title", &CifToMtz::title)
    .def_readwrite("history", &CifToMtz::history)
    .def_readwrite("spec_lines", &CifToMtz::spec_lines)
    .def("convert_block_to_mtz", [](const CifToMtz& self, const ReflnBlock& rb) {
        std::ostringstream out;
        return new Mtz(self.convert_block_to_mtz(rb, out));
    })
    ;

  py::class_<Intensities> intensities(m, "Intensities");

  py::enum_<Intensities::Type>(intensities, "Type")
    .value("None", Intensities::Type::None)
    .value("Unmerged", Intensities::Type::Unmerged)
    .value("Mean", Intensities::Type::Mean)
    .value("Anomalous", Intensities::Type::Anomalous)
    ;

  intensities
    .def(py::init<>())
    .def_readwrite("spacegroup", &Intensities::spacegroup)
    .def_readwrite("unit_cell", &Intensities::unit_cell)
    .def_readwrite("type", &Intensities::type)
    .def("resolution_range", &Intensities::resolution_range)
    .def("remove_systematic_absences", &Intensities::remove_systematic_absences)
    .def("merge_in_place", &Intensities::merge_in_place, py::arg("itype"))
    .def("read_mtz", &Intensities::read_mtz, py::arg("mtz"), py::arg("type"))
    .def_property_readonly("miller_array", [](const Intensities& self) {
      const Intensities::Refl* data = self.data.data();
      py::array::ShapeContainer shape({(py::ssize_t)self.data.size(), 3});
      py::array::StridesContainer strides({(const char*)(data+1) - (const char*)data,
                                           sizeof(int)});
      return py::array_t<int>(shape, strides, &data->hkl[0], py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("value_array", [](const Intensities& self) {
      const Intensities::Refl* data = self.data.data();
      py::ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<double>({(py::ssize_t)self.data.size()}, {stride},
                                 &data->value, py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("sigma_array", [](const Intensities& self) {
      const Intensities::Refl* data = self.data.data();
      py::ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<double>({(py::ssize_t)self.data.size()}, {stride},
                                 &data->sigma, py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("nobs_array", [](const Intensities& self) {
      const Intensities::Refl* data = self.data.data();
      py::ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<short>({(py::ssize_t)self.data.size()}, {stride},
                                &data->nobs, py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("isign_array", [](const Intensities& self) {
      const Intensities::Refl* data = self.data.data();
      py::ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<short>({(py::ssize_t)self.data.size()}, {stride},
                                &data->isign, py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def("set_data", [](Intensities& self, const UnitCell& unit_cell, 
                        const SpaceGroup* sg, py::array_t<int> hkl, 
                        py::array_t<double> values, py::array_t<double> sigmas) {
      auto h = hkl.unchecked<2>();
      if (h.shape(1) != 3)
        throw std::domain_error("error: the size of the second dimension != 3");
      auto v = values.template unchecked<1>();
      auto s = sigmas.template unchecked<1>();
      if (h.shape(0) != v.shape(0) || h.shape(0) != s.shape(0))
        throw std::domain_error("error: arrays have different lengths");

      self.unit_cell = unit_cell;
      self.spacegroup = sg;
      self.data.reserve(h.shape(0));
      for (py::ssize_t i = 0; i < h.shape(0); ++i)
        if (!std::isnan(v(i)) && s(i) > 0)
          self.data.push_back({{{h(i, 0), h(i, 1), h(i, 2)}}, 1, 0, v(i), s(i)});

      self.switch_to_asu_indices();
      self.type = Intensities::Type::Unmerged;
    }, py::arg("cell"), py::arg("sg").none(false),
       py::arg("miller_array"), py::arg("value_array"), py::arg("sigma_array"))
    ;
}
