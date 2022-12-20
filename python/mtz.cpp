// Copyright 2019 Global Phasing Ltd.

#include "gemmi/mtz.hpp"
#include "gemmi/reindex.hpp"  // for reindex_mtz
#include "gemmi/fourier.hpp"
#include "gemmi/gz.hpp"
#include "tostr.hpp"

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Dataset>)
PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Column>)
PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Batch>)
PYBIND11_MAKE_OPAQUE(std::vector<const Mtz::Column*>)

namespace gemmi {
  // operator<< is used by stl_bind for vector's __repr__
  inline std::ostream& operator<< (std::ostream& os, const Mtz::Dataset& ds) {
    os << "<gemmi.Mtz.Dataset " << ds.id << ' ' << ds.project_name
       << '/' << ds.crystal_name << '/' << ds.dataset_name << '>';
    return os;
  }

  inline std::ostream& operator<< (std::ostream& os, const Mtz::Column* col) {
    os << "<gemmi.Mtz.Column " << col->label << " type " << col->type << '>';
    return os;
  }
}

template<typename F>
py::array_t<float> make_new_column(const Mtz& mtz, int dataset, F f) {
  if (!mtz.has_data())
    throw std::runtime_error("MTZ: the data must be read first");
  const UnitCell& cell = mtz.get_cell(dataset);
  if (!cell.is_crystal())
    throw std::runtime_error("MTZ: unknown unit cell parameters");
  py::array_t<float> arr(mtz.nreflections);
  py::buffer_info buf = arr.request();
  float* ptr = (float*) buf.ptr;
  for (int i = 0; i < mtz.nreflections; ++i) {
    int hidx = (int) mtz.columns.size() * i;
    ptr[i] = f(cell, mtz.data[hidx], mtz.data[hidx+1], mtz.data[hidx+2]);
  }
  return arr;
}

static py::array_t<float> make_1_d2_array(const Mtz& mtz, int dataset) {
  return make_new_column(mtz, dataset,
                         [](const UnitCell& cell, float h, float k, float l) {
                           return (float) cell.calculate_1_d2_double(h, k, l);
                         });
}
static py::array_t<float> make_d_array(const Mtz& mtz, int dataset) {
  return make_new_column(mtz, dataset,
                         [](const UnitCell& cell, float h, float k, float l) {
                           double _1_d2 = cell.calculate_1_d2_double(h, k, l);
                           return (float) (1.0 / std::sqrt(_1_d2));
                         });
}

void add_mtz(py::module& m) {
  py::class_<Mtz> mtz(m, "Mtz", py::buffer_protocol());
  py::class_<Mtz::Dataset> pyMtzDataset(mtz, "Dataset");
  py::class_<Mtz::Column> pyMtzColumn(mtz, "Column", py::buffer_protocol());
  py::class_<Mtz::Batch> pyMtzBatch(mtz, "Batch");

  py::bind_vector<std::vector<Mtz::Dataset>>(m, "MtzDatasets");
  py::bind_vector<std::vector<Mtz::Column>>(m, "MtzColumns");
  py::bind_vector<std::vector<Mtz::Batch>>(m, "MtzBatches");
  py::bind_vector<std::vector<const Mtz::Column*>>(m, "MtzColumnRefs");

  mtz
    .def(py::init<bool>(), py::arg("with_base")=false)
    .def_buffer([](Mtz &self) {
      int nrow = self.has_data() ? self.nreflections : 0;
      int ncol = (int) self.columns.size();
      return py::buffer_info(self.data.data(),
                             {nrow, ncol}, // dimensions
                             {4 * ncol, 4});  // strides
    })
    .def_property_readonly("array", [](const Mtz& self) {
      int nrow = self.has_data() ? self.nreflections : 0;
      int ncol = (int) self.columns.size();
      return py::array_t<float>({nrow, ncol}, {4 * ncol, 4},
                                self.data.data(), py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_readwrite("title", &Mtz::title)
    .def_readwrite("nreflections", &Mtz::nreflections)
    .def_readwrite("sort_order", &Mtz::sort_order)
    .def_readwrite("min_1_d2", &Mtz::min_1_d2)
    .def_readwrite("max_1_d2", &Mtz::max_1_d2)
    .def_readwrite("valm", &Mtz::valm)
    .def_readonly("nsymop", &Mtz::nsymop)
    .def_readwrite("cell", &Mtz::cell)
    .def_readwrite("spacegroup", &Mtz::spacegroup)
    .def_readonly("spacegroup_name", &Mtz::spacegroup_name)
    .def_readonly("spacegroup_number", &Mtz::spacegroup_number)
    .def_readwrite("datasets", &Mtz::datasets)
    .def_readwrite("columns", &Mtz::columns)
    .def_readwrite("batches", &Mtz::batches)
    .def_readwrite("history", &Mtz::history)
    .def_readwrite("appended_text", &Mtz::appended_text)
    .def("resolution_high", &Mtz::resolution_high)
    .def("resolution_low", &Mtz::resolution_low)
    .def("dataset", (Mtz::Dataset& (Mtz::*)(int)) &Mtz::dataset,
         py::arg("id"))
    .def("count", &Mtz::count, py::arg("label"))
    .def("column_with_label",
         (Mtz::Column* (Mtz::*)(const std::string&, const Mtz::Dataset*))
                                                     &Mtz::column_with_label,
         py::arg("label"), py::arg("dataset")=nullptr,
         py::return_value_policy::reference_internal)
    .def("rfree_column", (Mtz::Column* (Mtz::*)()) &Mtz::rfree_column,
          py::return_value_policy::reference_internal)
    .def("columns_with_type", &Mtz::columns_with_type,
         py::arg("type"), py::keep_alive<0, 1>())
    .def("column_labels", [](const Mtz& self) {
        std::vector<std::string> labels;
        labels.reserve(self.columns.size());
        for (const Mtz::Column& c : self.columns)
          labels.push_back(c.label);
        return labels;
    })
    .def("get_cell", (UnitCell& (Mtz::*)(int)) &Mtz::get_cell,
         py::arg("dataset")=-1)
    .def("set_cell_for_all", &Mtz::set_cell_for_all)
    .def("make_miller_array", [](const Mtz& self) {
        py::array_t<int> arr({self.nreflections, 3});
        py::buffer_info buf = arr.request();
        int* ptr = (int*) buf.ptr;
        for (int i = 0; i < self.nreflections; ++i)
          for (int j = 0; j != 3; ++j)
            ptr[3*i + j] = (int) self.data[self.columns.size() * i + j];
        return arr;
    })
    .def("make_1_d2_array", &make_1_d2_array, py::arg("dataset")=-1)
    .def("make_d_array", &make_d_array, py::arg("dataset")=-1)
    .def("get_size_for_hkl",
         [](const Mtz& self, std::array<int,3> min_size, double sample_rate) {
          return get_size_for_hkl(MtzDataProxy{self}, min_size, sample_rate);
    }, py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.)
    .def("data_fits_into", [](const Mtz& self, std::array<int,3> size) {
        return data_fits_into(MtzDataProxy{self}, size);
    }, py::arg("size"))
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
    }, py::arg("f"), py::arg("phi"), py::arg("size"),
       py::arg("half_l")=false, py::arg("order")=AxisOrder::XYZ)
    .def("get_value_on_grid", [](const Mtz& self,
                                 const std::string& label,
                                 std::array<int, 3> size,
                                 bool half_l,
                                 AxisOrder order) {
        const Mtz::Column& col = self.get_column_with_label(label);
        return get_value_on_grid<float>(MtzDataProxy{self}, col.idx,
                                        size, half_l, order);
    }, py::arg("label"), py::arg("size"), py::arg("half_l")=false,
       py::arg("order")=AxisOrder::XYZ)
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
    }, py::arg("f"), py::arg("phi"),
       py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("exact_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.,
       py::arg("order")=AxisOrder::XYZ)
    .def("get_float", &make_asu_data<float, Mtz>,
         py::arg("col"), py::arg("as_is")=false)
    .def("get_int", &make_asu_data<int, Mtz>,
         py::arg("col"), py::arg("as_is")=false)
    .def("get_f_phi", [](const Mtz& self, const std::string& f_col,
                                          const std::string& phi_col,
                                          bool as_is) {
        return make_asu_data<std::complex<float>, 2>(self, {f_col, phi_col}, as_is);
    }, py::arg("f"), py::arg("phi"), py::arg("as_is")=false)
    .def("get_value_sigma", [](const Mtz& self, const std::string& f_col,
                                                const std::string& sigma_col,
                                                bool as_is) {
        return make_asu_data<ValueSigma<float>, 2>(self, {f_col, sigma_col}, as_is);
    }, py::arg("f"), py::arg("sigma"), py::arg("as_is")=false)
    .def("add_dataset", &Mtz::add_dataset, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("add_column", &Mtz::add_column, py::arg("label"), py::arg("type"),
         py::arg("dataset_id")=-1, py::arg("pos")=-1, py::arg("expand_data")=true,
         py::return_value_policy::reference_internal)
    .def("replace_column", &Mtz::replace_column,
         py::arg("dest_idx"), py::arg("src_col"),
         py::arg("trailing_cols")=std::vector<std::string>(),
         py::return_value_policy::reference_internal)
    .def("copy_column", &Mtz::copy_column,
         py::arg("dest_idx"), py::arg("src_col"),
         py::arg("trailing_cols")=std::vector<std::string>(),
         py::return_value_policy::reference_internal)
    .def("remove_column", &Mtz::remove_column, py::arg("index"))
    .def("set_data", [](Mtz& self, const AsuData<std::complex<float>>& asu_data) {
         if (self.columns.size() != 5)
           fail("Mtz.set_data(): Mtz must have 5 columns to put H,K,L,F,Phi.");
         self.nreflections = (int) asu_data.v.size();
         self.data.clear();
         add_asu_f_phi_to_float_vector(self.data, asu_data);
    }, py::arg("asu_data"))
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
    }, py::arg("asu_data"))
    .def("set_data", [](Mtz& self, py::array_t<float> arr) {
         if (arr.ndim() != 2)
           fail("Mtz.set_data(): expected 2D array.");
         int nrow = (int) arr.shape(0);
         int ncol = (int) arr.shape(1);
         if (ncol != (int) self.columns.size())
           fail("Mtz.set_data(): expected " +
                std::to_string(self.columns.size()) + " columns.");
         self.nreflections = nrow;
         self.data.resize(nrow * ncol);
         auto r = arr.unchecked<2>();
         for (py::ssize_t row = 0; row < nrow; row++)
           for (py::ssize_t col = 0; col < ncol; col++)
             self.data[row*ncol+col] = r(row, col);
    }, py::arg("array"))
    .def("update_reso", &Mtz::update_reso)
    .def("sort", &Mtz::sort, py::arg("use_first")=3)
    .def("ensure_asu", &Mtz::ensure_asu, py::arg("tnt_asu")=false)
    .def("switch_to_original_hkl", &Mtz::switch_to_original_hkl)
    .def("switch_to_asu_hkl", &Mtz::switch_to_asu_hkl)
    .def("write_to_file", &Mtz::write_to_file, py::arg("path"))
    .def("reindex", [](Mtz& self, const Op& op) {
        std::ostringstream out;
        reindex_mtz(self, op, &out);
        return out.str();
    }, py::arg("op"))
    .def("__repr__", [](const Mtz& self) {
        return tostr("<gemmi.Mtz with ", self.columns.size(), " columns, ",
                     self.nreflections, " reflections>");
    });

  pyMtzDataset
    .def_readwrite("id", &Mtz::Dataset::id)
    .def_readwrite("project_name", &Mtz::Dataset::project_name)
    .def_readwrite("crystal_name", &Mtz::Dataset::crystal_name)
    .def_readwrite("dataset_name", &Mtz::Dataset::dataset_name)
    .def_readwrite("cell", &Mtz::Dataset::cell)
    .def_readwrite("wavelength", &Mtz::Dataset::wavelength)
    .def("__repr__", [](const Mtz::Dataset& self) { return tostr(self); })
    ;
  pyMtzColumn
    .def_buffer([](Mtz::Column& self) {
      return py::buffer_info(self.parent->data.data() + self.idx,
                             std::vector<py::ssize_t>(1, self.size()), // dimensions
                             {4 * self.stride()});  // strides
    })
    .def_property_readonly("array", [](const Mtz::Column& self) {
      return py::array_t<float>({self.size()}, {4 * self.stride()},
                                self.parent->data.data() + self.idx,
                                py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("dataset",
            (Mtz::Dataset& (Mtz::Column::*)()) &Mtz::Column::dataset)
    .def_readwrite("dataset_id", &Mtz::Column::dataset_id)
    .def_readwrite("type", &Mtz::Column::type)
    .def_readwrite("label", &Mtz::Column::label)
    .def_readwrite("min_value", &Mtz::Column::min_value)
    .def_readwrite("max_value", &Mtz::Column::max_value)
    .def_readwrite("source", &Mtz::Column::source)
    .def_readwrite("idx", &Mtz::Column::idx)
    .def("is_integer", &Mtz::Column::is_integer)
    .def("__len__", &Mtz::Column::size)
    .def("__getitem__", [](const Mtz::Column& self, int index) -> float {
        return self.at(index >= 0 ? index : index + self.size());
    }, py::arg("index"))
    .def("__iter__", [](Mtz::Column& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Mtz::Column& self) { return tostr(&self); })
    ;
  pyMtzBatch
    .def_readwrite("number", &Mtz::Batch::number)
    .def_property_readonly("dataset_id", &Mtz::Batch::dataset_id)
    .def_readwrite("title", &Mtz::Batch::title)
    .def_readonly("ints", &Mtz::Batch::ints)
    .def_readonly("floats", &Mtz::Batch::floats)
    .def_readonly("axes", &Mtz::Batch::axes)
    ;

  m.def("read_mtz_file", [](const std::string& path) {
      return read_mtz(MaybeGzipped(path), true);
  }, py::arg("path"), py::return_value_policy::move);
}
