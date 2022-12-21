// Copyright 2020 Global Phasing Ltd.

#include "gemmi/recgrid.hpp"
#include "gemmi/fourier.hpp"  // for get_size_for_hkl, get_f_phi_on_grid
#include "tostr.hpp"

#include "common.h"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

namespace gemmi {
  std::ostream& operator<< (std::ostream& os, const ValueSigma<float>& vs) {
    os << vs.value << " +/- " << vs.sigma;
    return os;
  }
}

template<typename T, typename F>
py::array_t<float> make_new_column(const AsuData<T>& asu_data, F f) {
  if (!asu_data.unit_cell().is_crystal())
    throw std::runtime_error("AsuData: unknown unit cell parameters");
  py::array_t<float> arr(asu_data.size());
  py::buffer_info buf = arr.request();
  float* ptr = (float*) buf.ptr;
  for (size_t i = 0; i < asu_data.size(); ++i)
    ptr[i] = static_cast<float>(f(asu_data.unit_cell(), asu_data.get_hkl(i)));
  return arr;
}

template<typename T> void add_to_asu_data(T&) {}
template<> void add_to_asu_data(py::class_<AsuData<std::complex<float>>>& cl) {
  using AsuData = AsuData<std::complex<float>>;
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
  cl.def("calculate_correlation", [](const AsuData& self, const AsuData& other) {
      return calculate_hkl_complex_correlation(self.v, other.v);
  });
}
template<> void add_to_asu_data(py::class_<AsuData<float>>& cl) {
  using AsuData = AsuData<float>;
  cl.def("calculate_correlation", [](const AsuData& self, const AsuData& other) {
      return calculate_hkl_value_correlation(self.v, other.v);
  });
}

template<typename T>
void add_asudata(py::module& m, const std::string& prefix) {
  py::class_<HklValue<T>>(m, (prefix + "HklValue").c_str())
    .def_readonly("hkl", &HklValue<T>::hkl)
    .def_readwrite("value", &HklValue<T>::value)
    .def("__repr__", [prefix](const HklValue<T>& self) {
        return tostr("<gemmi.", prefix, "HklValue (",
                     self.hkl[0], ',', self.hkl[1], ',', self.hkl[2], ") ",
                     self.value, '>');
    });

  using AsuData = AsuData<T>;
  py::class_<AsuData> asu_data(m, (prefix + "AsuData").c_str());
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
      for (py::ssize_t i = 0; i < h.shape(0); ++i)
        ret->v.push_back({{{h(i, 0), h(i, 1), h(i, 2)}}, v(i)});
      return ret;
    }), py::arg("cell"), py::arg("sg").none(false),
        py::arg("miller_array"), py::arg("value_array"))
    .def("__iter__", [](AsuData& self) { return py::make_iterator(self.v); },
         py::keep_alive<0, 1>())
    .def("__len__", [](const AsuData& self) { return self.v.size(); })
    .def("__getitem__", [](AsuData& self, int index) -> HklValue<T>& {
        return self.v.at(normalize_index(index, self.v));
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def_readwrite("spacegroup", &AsuData::spacegroup_)
    .def_readwrite("unit_cell", &AsuData::unit_cell_)
    .def_property_readonly("miller_array", [](const AsuData& self) {
      const HklValue<T>* data = self.v.data();
      py::array::ShapeContainer shape({(py::ssize_t)self.v.size(), 3});
      py::array::StridesContainer strides({(const char*)(data+1) - (const char*)data,
                                           sizeof(int)});
      return py::array_t<int>(shape, strides, &data->hkl[0], py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("value_array", [](const AsuData& self) {
      const HklValue<T>* data = self.v.data();
      py::ssize_t stride = (const char*)(data+1) - (const char*)data;
      return py::array_t<T>({(py::ssize_t)self.v.size()}, {stride},
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
    .def("count_equal_values", [](const AsuData& self, const AsuData& other) {
      return count_equal_values(self.v, other.v);
    })
    .def("ensure_sorted", &AsuData::ensure_sorted)
    .def("ensure_asu", &AsuData::ensure_asu, py::arg("tnt_asu")=false)
    .def("copy", [](const AsuData& self) {
      return new AsuData(self);
    })
    .def("__repr__", [prefix](const AsuData& self) {
        return tostr("<gemmi.", prefix, "AsuData with ", self.v.size(), " values>");
    });
    add_to_asu_data(asu_data);
}

template<typename TA, typename TG=TA>
void add_asudata_and_recgrid(py::module& m,
                             const std::string& prefix_asu,
                             const std::string& rgrid_name) {
  using RecGr = ReciprocalGrid<TG>;
  py::class_<RecGr, GridBase<TG>> recgrid(m, rgrid_name.c_str());

  add_asudata<TA>(m, prefix_asu);

  recgrid
    .def_readonly("half_l", &RecGr::half_l)
    .def(py::init<>())
    .def(py::init([](int nx, int ny, int nz) {
      RecGr* grid = new RecGr();
      grid->set_size_without_checking(nx, ny, nz);
      grid->axis_order = AxisOrder::XYZ;
      return grid;
    }), py::arg("nx"), py::arg("ny"), py::arg("nz"))
    .def(py::init([](py::array_t<TG> arr, const UnitCell *cell, const SpaceGroup* sg) {
      auto r = arr.template unchecked<3>();
      RecGr* grid = new RecGr();
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
    .def("get_value", &RecGr::get_value)
    .def("get_value_or_zero", &RecGr::get_value_or_zero)
    .def("set_value", &RecGr::set_value)
    .def("to_hkl", &RecGr::to_hkl)
    .def("calculate_1_d2", &RecGr::calculate_1_d2)
    .def("calculate_d", &RecGr::calculate_d)
    .def("get_value_by_hkl", [](RecGr &self, py::array_t<int> hkl, double unblur,
                                bool mott_bethe, TA mott_bethe_000) {
      auto h = hkl.unchecked<2>();
      if (h.shape(1) != 3)
        throw std::domain_error("error: the size of the second dimension != 3");
      py::array_t<TA> vals(h.shape(0));
      TA* ptr = (TA*) vals.request().ptr;
      for (py::ssize_t i = 0; i < h.shape(0); ++i) {
        Miller mi{{h(i, 0), h(i, 1), h(i, 2)}};
        if (mott_bethe && mi[0] == 0 && mi[1] == 0 && mi[2] == 0)
          ptr[i] = mott_bethe_000;
        else
          ptr[i] = self.get_value_by_hkl(mi, unblur, mott_bethe);
      }
      return vals;
    }, py::arg("hkl"), py::arg("unblur")=0, py::arg("mott_bethe")=false,
       py::arg("mott_bethe_000")=0)

    .def("prepare_asu_data", &RecGr::template prepare_asu_data<TA>,
         py::arg("dmin")=0., py::arg("unblur")=0.,
         py::arg("with_000")=false, py::arg("with_sys_abs")=false,
         py::arg("mott_bethe")=false)
    .def("__repr__", [=](const RecGr& self) {
        return tostr("<gemmi.", rgrid_name, '(', self.nu, ", ", self.nv, ", ", self.nw, ")>");
    });
}

void add_recgrid(py::module& m) {
  using VS = ValueSigma<float>;
  try {
    PYBIND11_NUMPY_DTYPE(VS, value, sigma);
  } catch (py::error_already_set &) {
    // numpy can't be imported, that's ok
  }
  py::class_<VS>(m, "ValueSigma")
    .def_readwrite("value", &VS::value)
    .def_readwrite("sigma", &VS::sigma)
    .def("__repr__", [](const VS& self) {
        return tostr("<gemmi.ValueSigma(", self.value, ", ", self.sigma, ")>");
    });
  py::class_<ComplexCorrelation>(m, "ComplexCorrelation")
    .def_readonly("n", &ComplexCorrelation::n)
    .def("coefficient", &ComplexCorrelation::coefficient)
    .def("mean_ratio", &ComplexCorrelation::mean_ratio)
    ;
  add_asudata_and_recgrid<int, int8_t>(m, "Int", "ReciprocalInt8Grid");
  add_asudata_and_recgrid<float>(m, "Float", "ReciprocalFloatGrid");
  add_asudata_and_recgrid<std::complex<float>>(m, "Complex", "ReciprocalComplexGrid");
  add_asudata<VS>(m, "ValueSigma");
}
