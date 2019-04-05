// Copyright 2019 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/mtz.hpp"
#include "gemmi/refln.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Column>)
PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Dataset>)
PYBIND11_MAKE_OPAQUE(std::vector<ReflnBlock>)

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
    int hidx = mtz.ncol * i;
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

  py::class_<Mtz> mtz(m, "Mtz", py::buffer_protocol());
  mtz.def(py::init<>())
    .def_buffer([](Mtz &self) {
      int nrow = self.has_data() ? self.nreflections : 0;
      return py::buffer_info(self.data.data(),
                             {nrow, self.ncol}, // dimensions
                             {4 * self.ncol, 4});  // strides
    })
    .def_readwrite("title", &Mtz::title)
    .def_readwrite("ncol", &Mtz::ncol)
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
         (Mtz::Column* (Mtz::*)(const std::string&)) &Mtz::column_with_label,
         py::arg("label"), py::return_value_policy::reference_internal)
    .def("column_with_type",
         (Mtz::Column* (Mtz::*)(char)) &Mtz::column_with_type,
         py::arg("type"), py::return_value_policy::reference_internal)
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
    .def("get_map_coef_as_grid", &Mtz::get_map_coef_as_grid<float>,
         py::arg("f"), py::arg("phi"), py::arg("size")=std::array<int,3>{0,0,0})
    .def("write_to_file", &Mtz::write_to_file, py::arg("path"))
    .def("__repr__", [](const Mtz& self) {
        return "<gemmi.Mtz with " + std::to_string(self.ncol) + " columns, " +
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
    .def("__len__", &Mtz::Column::size)
    .def("__getitem__", [](const Mtz::Column& self, int index) -> float {
        return self.at(index >= 0 ? index : index + self.size());
    }, py::arg("index"))
    .def("__iter__", [](Mtz::Column& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>())
    .def("__repr__", [](const Mtz::Column& self) {
        return "<gemmi.Mtz.Column " + self.label +
               " type " + std::string(1, self.type) + ">";
    });

  m.def("read_mtz_file", &read_mtz_file);

  py::class_<ReflnBlock>(m, "ReflnBlock")
    .def_readonly("block", &ReflnBlock::block)
    .def_readonly("entry_id", &ReflnBlock::entry_id)
    .def_readonly("cell", &ReflnBlock::cell)
    .def_readonly("spacegroup", &ReflnBlock::spacegroup)
    .def_readonly("wavelength", &ReflnBlock::wavelength)
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
    .def("__repr__", [](const ReflnBlock& self) {
      return tostr(self);
    });
  m.def("as_refln_blocks",
        [](cif::Document& d) { return as_refln_blocks(std::move(d.blocks)); });
}
