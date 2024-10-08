// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include "array.h"
#include "make_iterator.h"
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>

#include "gemmi/recgrid.hpp"
#include "gemmi/fourier.hpp"  // for get_size_for_hkl, get_f_phi_on_grid
#include "gemmi/sprintf.hpp"

using namespace gemmi;

namespace {
template<typename T, typename Func>
auto make_new_column(const AsuData<T>& asu_data, Func f) {
  if (!asu_data.unit_cell().is_crystal())
    throw std::runtime_error("AsuData: unknown unit cell parameters");
  auto numpy_arr = make_numpy_array<float>({asu_data.size()});
  float* arr = numpy_arr.data();
  for (size_t i = 0; i < asu_data.size(); ++i)
    arr[i] = static_cast<float>(f(asu_data.unit_cell(), asu_data.get_hkl(i)));
  return numpy_arr;
}

template<typename T> void add_to_asu_data(T&) {}
template<> void add_to_asu_data(nb::class_<AsuData<std::complex<float>>>& cl) {
  using AsuData = AsuData<std::complex<float>>;
  cl.def("get_size_for_hkl", &get_size_for_hkl<AsuData>,
         nb::arg("min_size")=std::array<int,3>{{0,0,0}}, nb::arg("sample_rate")=0.);
  cl.def("data_fits_into", &data_fits_into<AsuData>, nb::arg("size"));
  cl.def("get_f_phi_on_grid", get_f_phi_on_grid<float, AsuData>,
         nb::arg("size"), nb::arg("half_l")=false, nb::arg("order")=AxisOrder::XYZ);
  cl.def("transform_f_phi_to_map", &transform_f_phi_to_map2<float, AsuData>,
         nb::arg("min_size")=std::array<int,3>{{0,0,0}},
         nb::arg("sample_rate")=0.,
         nb::arg("exact_size")=std::array<int,3>{{0,0,0}},
         nb::arg("order")=AxisOrder::XYZ);
  cl.def("calculate_correlation", [](const AsuData& self, const AsuData& other) {
      return calculate_hkl_complex_correlation(self.v, other.v);
  });
}
template<> void add_to_asu_data(nb::class_<AsuData<float>>& cl) {
  using AsuData = AsuData<float>;
  cl.def("calculate_correlation", [](const AsuData& self, const AsuData& other) {
      return calculate_hkl_value_correlation(self.v, other.v);
  });
}
template<> void add_to_asu_data(nb::class_<AsuData<ValueSigma<float>>>& cl) {
  cl.def("discard_by_sigma_ratio", &discard_by_sigma_ratio<float>);
}

template<typename T> struct array_for {
  using type = cpu_c_array<T>;
};
template<> struct array_for<ValueSigma<float>> {
  using type = nb::ndarray<float, nb::shape<-1, 2>, nb::device::cpu, nb::c_contig>;
};

template<typename T>
void add_asudata(nb::module_& m, const std::string& prefix) {
  nb::class_<HklValue<T>>(m, (prefix + "HklValue").c_str())
    .def_ro("hkl", &HklValue<T>::hkl)
    .def_rw("value", &HklValue<T>::value)
    .def("__repr__", [prefix](const HklValue<T>& self) {
        return nb::str("<gemmi.{}HklValue ({},{},{}) {}>")
                .format(prefix, self.hkl[0], self.hkl[1], self.hkl[2], self.value);
    });

  using AsuData = AsuData<T>;
  nb::class_<AsuData> asu_data(m, (prefix + "AsuData").c_str());
  asu_data
    .def("__init__", [](AsuData* p, const UnitCell& unit_cell, const SpaceGroup* sg,
                        const cpu_miller_array& hkl, const typename array_for<T>::type& values) {
      new(p) AsuData;
      auto h = hkl.view();
      auto v = values.view();
      if (h.shape(0) != v.shape(0))
        throw std::domain_error("error: arrays have different lengths");
      p->spacegroup_ = sg;
      p->unit_cell_ = unit_cell;
      p->unit_cell_.set_cell_images_from_spacegroup(p->spacegroup_);
      p->v.reserve(h.shape(0));
      for (size_t i = 0; i < h.shape(0); ++i) {
        Miller hkl {h(i, 0), h(i, 1), h(i, 2)};
        if constexpr (std::is_same<T, ValueSigma<float>>::value)
          p->v.push_back({hkl, T{v(i, 0), v(i, 1)}});
        else
          p->v.push_back({hkl, v(i)});
      }
    }, nb::arg("cell"), nb::arg("sg")/*.none(false)*/,
       nb::arg("miller_array"), nb::arg("value_array"))
    .def("__iter__", [](AsuData& self) {
        return usual_iterator(self, self.v);
    }, nb::keep_alive<0, 1>())
    .def("__len__", [](const AsuData& self) { return self.v.size(); })
    .def("__getitem__", [](AsuData& self, int index) -> HklValue<T>& {
        return self.v.at(normalize_index(index, self.v));
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def_rw("spacegroup", &AsuData::spacegroup_)
    .def_rw("unit_cell", &AsuData::unit_cell_)
    .def_prop_ro("miller_array", [](AsuData& self) {
      // cf. vector_member_array
      int64_t stride = int64_t(sizeof(HklValue<T>) / sizeof(int));
      return nb::ndarray<nb::numpy, int, nb::shape<-1,3>>(
          &self.v[0].hkl[0], {self.v.size(), 3}, nb::handle(), {stride, 1});
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("value_array", [](AsuData& self) {
      if constexpr (std::is_same<T, ValueSigma<float>>::value) {
        auto stride = sizeof(ValueSigma<HklValue<float>>) / sizeof(float);
        return nb::ndarray<nb::numpy, float, nb::shape<-1, 2>>(
                &(self.v.data()->value.value),
                {self.v.size(), 2},
                nb::handle(),
                {(int64_t)stride, 1});
      } else {
        return vector_member_array(self.v, &HklValue<T>::value);
      }
    }, nb::rv_policy::reference_internal)
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
    .def("ensure_asu", &AsuData::ensure_asu, nb::arg("tnt_asu")=false)
    .def("copy", [](const AsuData& self) {
      return new AsuData(self);
    })
    .def("__repr__", [prefix](const AsuData& self) {
        return cat("<gemmi.", prefix, "AsuData with ", self.v.size(), " values>");
    });
    add_to_asu_data(asu_data);
}

template<typename TA, typename TG=TA>
void add_asudata_and_recgrid(nb::module_& m,
                             const std::string& prefix_asu,
                             const std::string& rgrid_name) {
  using RecGr = ReciprocalGrid<TG>;
  // needed before add_asudata for ComplexAsuData.get_f_phi_on_grid
  nb::class_<RecGr, GridBase<TG>> recgrid(m, rgrid_name.c_str());

  add_asudata<TA>(m, prefix_asu);

  recgrid
    .def_ro("half_l", &RecGr::half_l)
    .def(nb::init<>())
    .def("__init__", [](RecGr* grid, int nx, int ny, int nz) {
      new(grid) RecGr();
      grid->set_size_without_checking(nx, ny, nz);
      grid->axis_order = AxisOrder::XYZ;
    }, nb::arg("nx"), nb::arg("ny"), nb::arg("nz"))
    .def("__init__", [](RecGr* grid,
                        const nb::ndarray<TG, nb::ndim<3>, nb::device::cpu>& arr,
                        const UnitCell *cell, const SpaceGroup* sg) {
      auto r = arr.view();
      new(grid) RecGr();
      grid->set_size_without_checking((int)r.shape(0), (int)r.shape(1), (int)r.shape(2));
      grid->axis_order = AxisOrder::XYZ;
      for (size_t k = 0; k < r.shape(2); ++k)
        for (size_t j = 0; j < r.shape(1); ++j)
          for (size_t i = 0; i < r.shape(0); ++i)
            grid->data[grid->index_q(i, j, k)] = r(i, j, k);
      if (cell)
        grid->unit_cell = *cell;
      if (sg)
        grid->spacegroup = sg;
    }, nb::arg().noconvert(), nb::arg("cell")=nb::none(), nb::arg("spacegroup")=nb::none())
    .def("get_value", &RecGr::get_value)
    .def("get_value_or_zero", &RecGr::get_value_or_zero)
    .def("set_value", &RecGr::set_value)
    .def("to_hkl", &RecGr::to_hkl)
    .def("calculate_1_d2", &RecGr::calculate_1_d2)
    .def("calculate_d", &RecGr::calculate_d)
    .def("get_value_by_hkl", [](RecGr &self, const cpu_miller_array& hkl, double unblur,
                                bool mott_bethe, TA mott_bethe_000) {
      auto h = hkl.view();
      auto vals = make_numpy_array<TA>({h.shape(0)});
      TA* ptr = vals.data();
      for (size_t i = 0; i < h.shape(0); ++i) {
        Miller mi{h(i, 0), h(i, 1), h(i, 2)};
        if (mott_bethe && mi[0] == 0 && mi[1] == 0 && mi[2] == 0)
          ptr[i] = mott_bethe_000;
        else
          ptr[i] = self.get_value_by_hkl(mi, unblur, mott_bethe);
      }
      return vals;
    }, nb::arg("hkl"), nb::arg("unblur")=0, nb::arg("mott_bethe")=false,
       nb::arg("mott_bethe_000")=0)

    .def("prepare_asu_data", &RecGr::template prepare_asu_data<TA>,
         nb::arg("dmin")=0., nb::arg("unblur")=0.,
         nb::arg("with_000")=false, nb::arg("with_sys_abs")=false,
         nb::arg("mott_bethe")=false)
    .def("__repr__", [=](const RecGr& self) {
        return cat("<gemmi.", rgrid_name, '(', self.nu, ", ", self.nv, ", ", self.nw, ")>");
    });
}

}  // anonymous namespace

void add_recgrid(nb::module_& m) {
  using VS = ValueSigma<float>;
  nb::class_<VS>(m, "ValueSigma")
    .def_rw("value", &VS::value)
    .def_rw("sigma", &VS::sigma)
    .def("__repr__", [](const VS& self) {
        char buf[64];
        snprintf_z(buf, 64, "<gemmi.ValueSigma(%g, %g)>", self.value, self.sigma);
        return std::string(buf);
    });
  nb::class_<ComplexCorrelation>(m, "ComplexCorrelation")
    .def_ro("n", &ComplexCorrelation::n)
    .def("coefficient", &ComplexCorrelation::coefficient)
    .def("mean_ratio", &ComplexCorrelation::mean_ratio)
    ;
  add_asudata_and_recgrid<int, int8_t>(m, "Int", "ReciprocalInt8Grid");
  add_asudata_and_recgrid<float>(m, "Float", "ReciprocalFloatGrid");
  add_asudata_and_recgrid<std::complex<float>>(m, "Complex", "ReciprocalComplexGrid");
  add_asudata<VS>(m, "ValueSigma");
}
