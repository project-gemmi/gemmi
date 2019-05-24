// Copyright 2019 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/mtz.hpp"
#include "gemmi/refln.hpp"
#include "gemmi/fourier.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Column>)
PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Dataset>)
PYBIND11_MAKE_OPAQUE(std::vector<ReflnBlock>)
PYBIND11_MAKE_OPAQUE(std::vector<const Mtz::Column*>)

namespace gemmi {
  // operator<< is used by stl_bind for vector's __repr__
  inline std::ostream& operator<< (std::ostream& os, const Mtz::Dataset& ds) {
    os << "<gemmi.Mtz.Dataset " << ds.id << ' ' << ds.project_name
       << '/' << ds.crystal_name << '/' << ds.dataset_name << '>';
    return os;
  }

  inline std::ostream& operator<< (std::ostream& os, const ReflnBlock& rb) {
    os << "<gemmi.ReflnBlock " + rb.block.name + " with ";
    if (rb.default_loop)
      os << rb.default_loop->width() << " x " << rb.default_loop->length();
    else
      os << " no ";
    os << " loop>";
    return os;
  }

  inline std::ostream& operator<< (std::ostream& os, const Mtz::Column* col) {
    os << "<gemmi.Mtz.Column " << col->label << " type " << col->type << '>';
    return os;
  }
}

template <typename T> std::string tostr(const T& s) {
  std::ostringstream stream;
  stream << s;
  return stream.str();
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
    int hidx = mtz.columns.size() * i;
    ptr[i] = f(cell, mtz.data[hidx], mtz.data[hidx+1], mtz.data[hidx+2]);
  }
  return arr;
}

static py::array_t<float> make_1_d2_array(const Mtz& mtz, int dataset) {
  return make_new_column(mtz, dataset,
                         [](const UnitCell& cell, float h, float k, float l) {
                           return (float) cell.calculate_1_d2(h, k, l);
                         });
}
static py::array_t<float> make_d_array(const Mtz& mtz, int dataset) {
  return make_new_column(mtz, dataset,
                         [](const UnitCell& cell, float h, float k, float l) {
                           return (float) cell.calculate_d(h, k, l);
                         });
}

static
Grid<std::complex<float>> mtz_get_f_phi_on_grid(const Mtz& self,
                                                const std::string& f_col,
                                                const std::string& phi_col,
                                                bool half_l,
                                                std::array<int, 3> min_size,
                                                double sample_rate) {
  const Mtz::Column* f = self.column_with_label(f_col);
  const Mtz::Column* phi = self.column_with_label(phi_col);
  if (!f || !phi)
    fail("Column labels not found.");
  return get_f_phi_on_grid<float>(MtzDataProxy{self}, f->idx, phi->idx,
                                  half_l, min_size, sample_rate);
}

template<typename T>
py::array_t<T> py_array_from_vector(std::vector<T>&& original_vec) {
  auto v = new std::vector<T>(std::move(original_vec));
  py::capsule cap(v, [](void* p) { delete (std::vector<T>*) p; });
  return py::array_t<T>(v->size(), v->data(), cap);
}

void add_hkl(py::module& m) {
  py::bind_vector<std::vector<Mtz::Column>>(m, "MtzColumns");
  py::bind_vector<std::vector<Mtz::Dataset>>(m, "MtzDatasets");
  py::bind_vector<std::vector<ReflnBlock>>(m, "ReflnBlocks");
  py::bind_vector<std::vector<const Mtz::Column*>>(m, "MtzColumnRefs");

  py::class_<Mtz> mtz(m, "Mtz", py::buffer_protocol());
  mtz.def(py::init<>())
    .def_buffer([](Mtz &self) {
      int nrow = self.has_data() ? self.nreflections : 0;
      int ncol = (int) self.columns.size();
      return py::buffer_info(self.data.data(),
                             {nrow, ncol}, // dimensions
                             {4 * ncol, 4});  // strides
    })
    .def_readwrite("title", &Mtz::title)
    .def_readwrite("nreflections", &Mtz::nreflections)
    .def_readwrite("nbatches", &Mtz::nbatches)
    .def_readwrite("min_1_d2", &Mtz::min_1_d2)
    .def_readwrite("max_1_d2", &Mtz::max_1_d2)
    .def_readwrite("valm", &Mtz::valm)
    .def_readwrite("nsymop", &Mtz::nsymop)
    .def_readwrite("cell", &Mtz::cell)
    .def_readwrite("spacegroup", &Mtz::spacegroup)
    .def_readwrite("datasets", &Mtz::datasets)
    .def_readwrite("columns", &Mtz::columns)
    .def_readwrite("history", &Mtz::history)
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
    .def("make_1_d2_array", &make_1_d2_array, py::arg("dataset")=-1)
    .def("make_d_array", &make_d_array, py::arg("dataset")=-1)
    .def("get_f_phi_on_grid", &mtz_get_f_phi_on_grid,
         py::arg("f"), py::arg("phi"), py::arg("half_l")=false,
         py::arg("size")=std::array<int,3>{{0,0,0}},
         py::arg("sample_rate")=0.)
    .def("transform_f_phi_to_map", [](const Mtz& self,
                                      const std::string& f_col,
                                      const std::string& phi_col,
                                      std::array<int, 3> min_size,
                                      double sample_rate) {
        return transform_f_phi_grid_to_map(
                          mtz_get_f_phi_on_grid(self, f_col, phi_col, true,
                                                min_size, sample_rate));
    }, py::arg("f"), py::arg("phi"),
       py::arg("size")=std::array<int,3>{{0,0,0}}, py::arg("sample_rate")=0.)
    .def("add_dataset", &Mtz::add_dataset, py::arg("name"),
         py::return_value_policy::reference_internal)
    .def("add_column", &Mtz::add_column, py::arg("label"), py::arg("type"),
         py::arg("dataset_id")=-1, py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
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
         for (ssize_t row = 0; row < nrow; row++)
           for (ssize_t col = 0; col < ncol; col++)
             self.data[row*ncol+col] = r(row, col);
    }, py::arg("array"))
    .def("write_to_file", &Mtz::write_to_file, py::arg("path"))
    .def("__repr__", [](const Mtz& self) {
        return "<gemmi.Mtz with " +
               std::to_string(self.columns.size()) + " columns, " +
               std::to_string(self.nreflections) + " reflections>";
    });
  py::class_<Mtz::Dataset>(mtz, "Dataset")
    .def_readwrite("id", &Mtz::Dataset::id)
    .def_readwrite("project_name", &Mtz::Dataset::project_name)
    .def_readwrite("crystal_name", &Mtz::Dataset::crystal_name)
    .def_readwrite("dataset_name", &Mtz::Dataset::dataset_name)
    .def_readwrite("cell", &Mtz::Dataset::cell)
    .def_readwrite("wavelength", &Mtz::Dataset::wavelength)
    .def("__repr__", [](const Mtz::Dataset& self) {
      return tostr(self);
    });
  py::class_<Mtz::Column>(mtz, "Column", py::buffer_protocol())
    .def_buffer([](Mtz::Column& self) {
      return py::buffer_info(self.parent->data.data() + self.idx,
                             {self.size()},         // dimensions
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
    .def("__len__", &Mtz::Column::size)
    .def("__getitem__", [](const Mtz::Column& self, int index) -> float {
        return self.at(index >= 0 ? index : index + self.size());
    }, py::arg("index"))
    .def("__iter__", [](Mtz::Column& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Mtz::Column& self) { return tostr(&self); });

  m.def("read_mtz_file", &read_mtz_file);

  py::class_<ReflnBlock>(m, "ReflnBlock")
    .def_readonly("block", &ReflnBlock::block)
    .def_readonly("entry_id", &ReflnBlock::entry_id)
    .def_readonly("cell", &ReflnBlock::cell)
    .def_readonly("spacegroup", &ReflnBlock::spacegroup)
    .def_readonly("wavelength", &ReflnBlock::wavelength)
    .def("column_labels", &ReflnBlock::column_labels)
    .def("make_array_int",
         [](ReflnBlock& self, const std::string& tag, int null) {
           return py_array_from_vector(self.make_vector(tag, null));
    }, py::arg("tag"), py::arg("null"))
    .def("make_array_float",
         [](ReflnBlock& self, const std::string& tag, double null) {
           return py_array_from_vector(self.make_vector(tag, null));
    }, py::arg("tag"), py::arg("null")=NAN)
    .def("make_array_float", &ReflnBlock::make_vector<double>,
         py::arg("tag"), py::arg("null")=NAN)
    //.def("make_hkl_array", &ReflnBlock::make_hkl_vector)
    .def("make_1_d2_array", [](ReflnBlock& self) {
        return py_array_from_vector(self.make_1_d2_vector());
    })
    .def("get_f_phi_on_grid", [](const ReflnBlock& self,
                                 const std::string& f_col,
                                 const std::string& phi_col,
                                 bool half_l,
                                 std::array<int, 3> min_size,
                                 double sample_rate) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        return get_f_phi_on_grid<float>(ReflnDataProxy{self}, f_idx, phi_idx,
                                        half_l, min_size, sample_rate);
    }, py::arg("f"), py::arg("phi"), py::arg("half_l")=false,
       py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.)
    .def("transform_f_phi_to_map", [](const ReflnBlock& self,
                                      const std::string& f_col,
                                      const std::string& phi_col,
                                      std::array<int, 3> min_size,
                                      double sample_rate) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        return transform_f_phi_to_map<float>(ReflnDataProxy{self},
                                             f_idx, phi_idx,
                                             min_size, sample_rate);
    }, py::arg("f"), py::arg("phi"),
       py::arg("min_size")=std::array<int,3>{{0,0,0}},
       py::arg("sample_rate")=0.)
    .def("is_unmerged", &ReflnBlock::is_unmerged)
    .def("use_unmerged", &ReflnBlock::use_unmerged)
    .def("__bool__", [](const ReflnBlock& self) { return self.ok(); })
    .def("__repr__", [](const ReflnBlock& self) {
      return tostr(self);
    });
  m.def("as_refln_blocks",
        [](cif::Document& d) { return as_refln_blocks(std::move(d.blocks)); });
  m.def("transform_f_phi_grid_to_map", [](Grid<std::complex<float>> grid) {
          return transform_f_phi_grid_to_map<float>(std::move(grid));
        }, py::arg("grid"));
  m.def("transform_map_to_f_phi", &transform_map_to_f_phi<float>,
        py::arg("map"), py::arg("half_l")=false);
}
