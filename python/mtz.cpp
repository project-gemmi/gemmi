// Copyright 2019 Global Phasing Ltd.

#include "common.h"
#include "array.h"
#include "make_iterator.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "gemmi/mtz.hpp"
#include "gemmi/fourier.hpp"

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<Mtz::Dataset>)
NB_MAKE_OPAQUE(std::vector<Mtz::Column>)
NB_MAKE_OPAQUE(std::vector<Mtz::Batch>)

namespace {

template<typename T>
struct VectorRef {
  std::vector<T>& v;
  VectorRef(std::vector<T>& vec) : v(vec) {}
};

// Minimal std::vector bindings, for Batch::ints and Batch::floats.
// Don't allow the user to resize the vector, only to get and set values.
template<typename T>
void bind_wrapped_vector(nb::handle scope, const char* name) {
  using SizeType = typename std::vector<T>::size_type;
  using DiffType = typename std::vector<T>::difference_type;
  auto wrap_idx = [&](DiffType i, SizeType length) -> SizeType {
    SizeType idx = (i >= 0 ? (SizeType)i : (SizeType)i + length);
    if (idx >= length)
      throw nb::index_error();
    return idx;
  };
  using VR = VectorRef<T>;
  nb::class_<VR>(scope, name)
    .def("__getitem__", [&](VR& r, DiffType i) { return r.v[wrap_idx(i, r.v.size())]; })
    .def("__setitem__", [&](VR& r, DiffType i, T x) { r.v[wrap_idx(i, r.v.size())] = x; })
    .def("__len__", [](const VR& r) { return r.v.size(); })
    ;
}

template<typename Func>
auto make_new_column(const Mtz& mtz, int dataset, Func func) {
  if (!mtz.has_data())
    throw std::runtime_error("MTZ: the data must be read first");
  const UnitCell& cell = mtz.get_cell(dataset);
  if (!cell.is_crystal())
    throw std::runtime_error("MTZ: unknown unit cell parameters");
  size_t stride = mtz.columns.size();
  size_t n = (size_t) mtz.nreflections;
  auto numpy_arr = make_numpy_array<float>({n});
  float* arr = numpy_arr.data();
  const float* h = mtz.data.data();
  for (size_t i = 0; i < n; ++i, h += stride)
    arr[i] = func(cell, h[0], h[1], h[2]);
  return numpy_arr;
}

auto mtz_to_array(Mtz& self) {
  size_t nrow = self.has_data() ? (size_t) self.nreflections : 0;
  size_t ncol = self.columns.size();
  return nb::ndarray<nb::numpy, float, nb::ndim<2>>(self.data.data(), {nrow, ncol}, nb::handle());
}

auto column_to_array(Mtz::Column& self) {
  return nb::ndarray<nb::numpy, float, nb::ndim<1>>(self.parent->data.data() + self.idx,
                                                    {(size_t)self.size()},
                                                    nb::handle(),
                                                    {(int64_t) self.stride()});
}

}  // anonymous namespace

void add_mtz(nb::module_& m) {
  nb::class_<Mtz> mtz(m, "Mtz");
  nb::class_<Mtz::Dataset> pyMtzDataset(mtz, "Dataset");
  nb::class_<Mtz::Column> pyMtzColumn(mtz, "Column");
  nb::class_<Mtz::Batch> pyMtzBatch(mtz, "Batch");

  nb::bind_vector<std::vector<Mtz::Dataset>, rv_ri>(m, "MtzDatasets");
  nb::bind_vector<std::vector<Mtz::Column>, rv_ri>(m, "MtzColumns");
  nb::bind_vector<std::vector<Mtz::Batch>, rv_ri>(m, "MtzBatches");
  bind_wrapped_vector<int>(m, "BatchInts");
  bind_wrapped_vector<float>(m, "BatchFloats");

  mtz
    .def(nb::init<bool>(), nb::arg("with_base")=false)
    .def_prop_ro("array", &mtz_to_array, nb::rv_policy::reference_internal)
    .def("__array__", [](nb::handle_t<Mtz>& h, nb::handle dtype, nb::handle copy) {
        return handle_numpy_array_args(h.attr("array"), dtype, copy);
    }, nb::arg("dtype")=nb::none(), nb::arg("copy")=nb::none())
    .def_rw("title", &Mtz::title)
    .def_rw("nreflections", &Mtz::nreflections)
    .def_rw("sort_order", &Mtz::sort_order)
    .def_rw("min_1_d2", &Mtz::min_1_d2)
    .def_rw("max_1_d2", &Mtz::max_1_d2)
    .def_rw("valm", &Mtz::valm)
    .def_ro("nsymop", &Mtz::nsymop)
    .def_rw("cell", &Mtz::cell)
    .def_rw("spacegroup", &Mtz::spacegroup, nb::arg().none())
    .def_ro("spacegroup_name", &Mtz::spacegroup_name)
    .def_ro("spacegroup_number", &Mtz::spacegroup_number)
    .def_rw("datasets", &Mtz::datasets)
    .def_rw("columns", &Mtz::columns)
    .def_rw("batches", &Mtz::batches)
    .def_rw("history", &Mtz::history)
    .def_rw("appended_text", &Mtz::appended_text)
    .def("set_logging", [](Mtz& self, Logger&& logger) {
        self.logger = std::move(logger);
    }, nb::arg().none())
    .def("resolution_high", &Mtz::resolution_high)
    .def("resolution_low", &Mtz::resolution_low)
    .def("dataset", (Mtz::Dataset& (Mtz::*)(int)) &Mtz::dataset,
         nb::arg("id"))
    .def("count", &Mtz::count, nb::arg("label"))
    .def("column_with_label",
         (Mtz::Column* (Mtz::*)(const std::string&, const Mtz::Dataset*, char))
                                                     &Mtz::column_with_label,
         nb::arg("label"), nb::arg("dataset")=nb::none(), nb::arg("type")='*',
         nb::rv_policy::reference_internal)
    .def("rfree_column", (Mtz::Column* (Mtz::*)()) &Mtz::rfree_column,
          nb::rv_policy::reference_internal)
    .def("columns_with_type", &Mtz::columns_with_type,
         nb::arg("type"), nb::rv_policy::reference_internal)
    .def("column_labels", [](const Mtz& self) {
        std::vector<std::string> labels;
        labels.reserve(self.columns.size());
        for (const Mtz::Column& c : self.columns)
          labels.push_back(c.label);
        return labels;
    })
    .def("get_cell", (UnitCell& (Mtz::*)(int)) &Mtz::get_cell,
         nb::arg("dataset")=-1)
    .def("set_cell_for_all", &Mtz::set_cell_for_all)
    .def("make_miller_array", [](const Mtz& self) {
        size_t n = (size_t) self.nreflections;
        auto numpy_arr = make_numpy_array<int>({n, 3});
        int* arr = numpy_arr.data();
        for (size_t i = 0; i < n; ++i)
          for (size_t j = 0; j != 3; ++j)
            *arr++ = (int) self.data[self.columns.size() * i + j];
        return numpy_arr;
    })
    .def("make_1_d2_array", [](const Mtz& mtz, int dataset) {
        return make_new_column(mtz, dataset,
            [](const UnitCell& cell, float h, float k, float l) noexcept {
              return (float) cell.calculate_1_d2_double(h, k, l);
            });
    }, nb::arg("dataset")=-1)
    .def("make_d_array", [](const Mtz& mtz, int dataset) {
        return make_new_column(mtz, dataset,
            [](const UnitCell& cell, float h, float k, float l) noexcept {
              double _1_d2 = cell.calculate_1_d2_double(h, k, l);
              return (float) (1.0 / std::sqrt(_1_d2));
        });
    }, nb::arg("dataset")=-1)
    .def("get_size_for_hkl",
         [](const Mtz& self, std::array<int,3> min_size, double sample_rate) {
          return get_size_for_hkl(MtzDataProxy{self}, min_size, sample_rate);
    }, nb::arg("min_size")=std::array<int,3>{{0,0,0}},
       nb::arg("sample_rate")=0.)
    .def("data_fits_into", [](const Mtz& self, std::array<int,3> size) {
        return data_fits_into(MtzDataProxy{self}, size);
    }, nb::arg("size"))
    .def("get_f_phi_on_grid", [](const Mtz& self,
                                 const std::string& f_col,
                                 const std::string& phi_col,
                                 std::array<int, 3> size,
                                 bool half_l,
                                 AxisOrder order) {
        const Mtz::Column& f = self.get_column_with_label(f_col);
        const Mtz::Column& phi = self.get_column_with_label(phi_col);
        FPhiProxy<MtzDataProxy> fphi(MtzDataProxy{self}, f.idx, phi.idx);
        return get_f_phi_on_grid<float>(fphi, size, half_l, order);
    }, nb::arg("f"), nb::arg("phi"), nb::arg("size"),
       nb::arg("half_l")=false, nb::arg("order")=AxisOrder::XYZ)
    .def("get_value_on_grid", [](const Mtz& self,
                                 const std::string& label,
                                 std::array<int, 3> size,
                                 bool half_l,
                                 AxisOrder order) {
        const Mtz::Column& col = self.get_column_with_label(label);
        return get_value_on_grid<float>(MtzDataProxy{self}, col.idx,
                                        size, half_l, order);
    }, nb::arg("label"), nb::arg("size"), nb::arg("half_l")=false,
       nb::arg("order")=AxisOrder::XYZ)
    .def("transform_f_phi_to_map", [](const Mtz& self,
                                      const std::string& f_col,
                                      const std::string& phi_col,
                                      std::array<int, 3> min_size,
                                      std::array<int, 3> exact_size,
                                      double sample_rate,
                                      AxisOrder order) {
        const Mtz::Column& f = self.get_column_with_label(f_col);
        const Mtz::Column& phi = self.get_column_with_label(phi_col);
        FPhiProxy<MtzDataProxy> fphi(MtzDataProxy{self}, f.idx, phi.idx);
        return transform_f_phi_to_map2<float>(fphi, min_size, sample_rate,
                                              exact_size, order);
    }, nb::arg("f"), nb::arg("phi"),
       nb::arg("min_size")=std::array<int,3>{{0,0,0}},
       nb::arg("exact_size")=std::array<int,3>{{0,0,0}},
       nb::arg("sample_rate")=0.,
       nb::arg("order")=AxisOrder::XYZ)
    .def("get_float", &make_asu_data<float, Mtz>,
         nb::arg("col"), nb::arg("as_is")=false)
    .def("get_int", &make_asu_data<int, Mtz>,
         nb::arg("col"), nb::arg("as_is")=false)
    .def("get_f_phi", [](const Mtz& self, const std::string& f_col,
                                          const std::string& phi_col,
                                          bool as_is) {
        return make_asu_data<std::complex<float>, 2>(self, {f_col, phi_col}, as_is);
    }, nb::arg("f"), nb::arg("phi"), nb::arg("as_is")=false)
    .def("get_value_sigma", [](const Mtz& self, const std::string& f_col,
                                                const std::string& sigma_col,
                                                bool as_is) {
        return make_asu_data<ValueSigma<float>, 2>(self, {f_col, sigma_col}, as_is);
    }, nb::arg("f"), nb::arg("sigma"), nb::arg("as_is")=false)
    .def("add_dataset", &Mtz::add_dataset, nb::arg("name"),
         nb::rv_policy::reference_internal)
    .def("add_column", &Mtz::add_column, nb::arg("label"), nb::arg("type"),
         nb::arg("dataset_id")=-1, nb::arg("pos")=-1, nb::arg("expand_data")=true,
         nb::rv_policy::reference_internal)
    .def("replace_column", &Mtz::replace_column,
         nb::arg("dest_idx"), nb::arg("src_col"),
         nb::arg("trailing_cols")=std::vector<std::string>(),
         nb::rv_policy::reference_internal)
    .def("copy_column", &Mtz::copy_column,
         nb::arg("dest_idx"), nb::arg("src_col"),
         nb::arg("trailing_cols")=std::vector<std::string>(),
         nb::rv_policy::reference_internal)
    .def("remove_column", &Mtz::remove_column, nb::arg("index"))
    .def("set_data", [](Mtz& self, const AsuData<std::complex<float>>& asu_data) {
         if (self.columns.size() != 5)
           fail("Mtz.set_data(): Mtz must have 5 columns to put H,K,L,F,Phi.");
         self.nreflections = (int) asu_data.v.size();
         self.data.clear();
         add_asu_f_phi_to_float_vector(self.data, asu_data);
    }, nb::arg("asu_data"))
    .def("set_data", [](Mtz& self, const AsuData<float>& asu_data) {
         if (self.columns.size() != 4)
           fail("Mtz.set_data(): Mtz must have 4 columns.");
         self.nreflections = (int) asu_data.v.size();
         self.data.clear();
         self.data.reserve(self.data.size() + asu_data.v.size() * 4);
         for (const auto& item : asu_data.v) {
           for (int i = 0; i != 3; ++i)
             self.data.push_back((float) item.hkl[i]);
           self.data.push_back(item.value);
         }
    }, nb::arg("asu_data"))
    .def("set_data", [](Mtz& self, const nb::ndarray<float, nb::ndim<2>>& arr) {
         size_t nrow = arr.shape(0);
         size_t ncol = arr.shape(1);
         if (ncol != self.columns.size())
           fail("Mtz.set_data(): expected " +
                std::to_string(self.columns.size()) + " columns.");
         self.nreflections = (int) nrow;
         self.data.resize(nrow * ncol);
         auto r = arr.view();
         for (size_t row = 0; row < nrow; row++)
           for (size_t col = 0; col < ncol; col++)
             self.data[row*ncol+col] = r(row, col);
    }, nb::arg("array"))
    .def("filtered", [](Mtz& self, const cpu_array<bool>& selection) {
        if (!self.has_data())
          throw std::runtime_error("Mtz.filtered(): no data, read it first");
        auto v = selection.view();
        if (v.shape(0) != (size_t) self.nreflections)
          throw nb::value_error("boolean array must match the number of reflections");
        std::vector<float> saved_data;
        saved_data.swap(self.data); // avoid copying data
        Mtz ret(self);
        saved_data.swap(self.data);
        ret.nreflections = 0;
        for (size_t i = 0; i < v.shape(0); ++i)
          ret.nreflections += v(i) ? 1 : 0;
        size_t ncol = ret.columns.size();
        ret.data.reserve(ncol * ret.nreflections);
        for (size_t i = 0, pos = 0; i < v.shape(0); ++i, pos += ncol)
          if (v(i))
            ret.data.insert(ret.data.end(), &self.data[pos], &self.data[pos + ncol]);
        return ret;
    })
    .def("update_reso", &Mtz::update_reso)
    .def("sort", &Mtz::sort, nb::arg("use_first")=3)
    .def("ensure_asu", &Mtz::ensure_asu, nb::arg("tnt_asu")=false)
    .def("switch_to_original_hkl", &Mtz::switch_to_original_hkl)
    .def("switch_to_asu_hkl", &Mtz::switch_to_asu_hkl)
    .def("write_to_file", &Mtz::write_to_file, nb::arg("path"))
    .def("reindex", &Mtz::reindex, nb::arg("op"))
    .def("expand_to_p1", &Mtz::expand_to_p1)
    // handy for testing, but slow and can't handle duplicated column names
    .def("row_as_dict", [](const Mtz& self, const Miller& hkl) {
        size_t offset = self.find_offset_of_hkl(hkl);
        nb::dict data;
        if (offset != (size_t)-1)
          for (const Mtz::Column& column : self.columns)
            data[column.label.c_str()] = self.data[offset++];
        return data;
    }, nb::arg("hkl"))
    .def("__repr__", [](const Mtz& self) {
        return cat("<gemmi.Mtz with ", self.columns.size(), " columns, ",
                   self.nreflections, " reflections>");
    });

  pyMtzDataset
    .def_rw("id", &Mtz::Dataset::id)
    .def_rw("project_name", &Mtz::Dataset::project_name)
    .def_rw("crystal_name", &Mtz::Dataset::crystal_name)
    .def_rw("dataset_name", &Mtz::Dataset::dataset_name)
    .def_rw("cell", &Mtz::Dataset::cell)
    .def_rw("wavelength", &Mtz::Dataset::wavelength)
    .def("__repr__", [](const Mtz::Dataset& self) {
      return cat("<gemmi.Mtz.Dataset ", self.id, ' ', self.project_name,
                 '/', self.crystal_name, '/', self.dataset_name, '>');
    })
    ;
  pyMtzColumn
    .def_prop_ro("array", &column_to_array, nb::rv_policy::reference_internal)
    .def("__array__", [](nb::handle_t<Mtz::Column>& h, nb::handle dtype, nb::handle copy) {
        return handle_numpy_array_args(h.attr("array"), dtype, copy);
    }, nb::arg("dtype")=nb::none(), nb::arg("copy")=nb::none())
    .def_prop_ro("dataset",
            (Mtz::Dataset& (Mtz::Column::*)()) &Mtz::Column::dataset)
    .def_rw("dataset_id", &Mtz::Column::dataset_id)
    .def_rw("type", &Mtz::Column::type)
    .def_rw("label", &Mtz::Column::label)
    .def_rw("min_value", &Mtz::Column::min_value)
    .def_rw("max_value", &Mtz::Column::max_value)
    .def_rw("source", &Mtz::Column::source)
    .def_rw("idx", &Mtz::Column::idx)
    .def("is_integer", &Mtz::Column::is_integer)
    .def("__len__", &Mtz::Column::size)
    .def("__getitem__", [](const Mtz::Column& self, int index) -> float {
        return self.at(index >= 0 ? index : index + self.size());
    }, nb::arg("index"))
    .def("__iter__", [](Mtz::Column& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>())
    .def("__repr__", [](const Mtz::Column& self) {
        return cat("<gemmi.Mtz.Column ", self.label, " type ", self.type, '>');
    })
    ;
  pyMtzBatch
    .def(nb::init<>())
    .def_rw("number", &Mtz::Batch::number)
    .def_rw("title", &Mtz::Batch::title)
    .def_prop_rw("ints",
        [](Mtz::Batch& self) { return VectorRef<int>(self.ints); },
        [](Mtz::Batch& self, const VectorRef<int>& x) { self.ints = x.v; },
        nb::rv_policy::reference_internal)
    .def_prop_rw("floats",
        [](Mtz::Batch& self) { return VectorRef<float>(self.floats); },
        [](Mtz::Batch& self, const VectorRef<float>& x) { self.floats = x.v; },
        nb::rv_policy::reference_internal)
    .def_rw("axes", &Mtz::Batch::axes)
    .def_prop_rw("cell", &Mtz::Batch::get_cell, &Mtz::Batch::set_cell)
    .def_prop_rw("dataset_id", &Mtz::Batch::dataset_id, &Mtz::Batch::set_dataset_id)
    .def_prop_rw("wavelength", &Mtz::Batch::wavelength, &Mtz::Batch::set_wavelength)
    .def("clone", [](const Mtz::Batch& self) { return new Mtz::Batch(self); })
    ;

  m.def("read_mtz_file", [](const std::string& path, Logger&& logging, bool with_data) {
    std::unique_ptr<Mtz> mtz(new Mtz);
    mtz->logger = std::move(logging);
    mtz->read_file_gz(path, with_data);
    return mtz.release();
  }, nb::arg("path"), nb::arg("logging")=nb::none(), nb::arg("with_data")=true);
}
