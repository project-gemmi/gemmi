// Copyright 2017 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/mtz.hpp"

#include <cstdio>
#include <array>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Column>)
PYBIND11_MAKE_OPAQUE(std::vector<Mtz::Dataset>)

static std::string triple(double x, double y, double z) {
  using namespace std;  // VS2015/17 doesn't like std::snprintf
  char buf[128];
  snprintf(buf, 128, "%g, %g, %g", x, y, z);
  return std::string(buf);
}

void add_unitcell(py::module& m) {
  py::class_<Vec3>(m, "Vec3")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Vec3::x)
    .def_readwrite("y", &Vec3::y)
    .def_readwrite("z", &Vec3::z)
    .def("__getitem__", (double (Vec3::*)(int) const) &Vec3::at)
    .def("__repr__", [](const Vec3& self) {
        return "<gemmi.Vec3(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Mat33>(m, "Mat33")
    .def(py::init<>())
    .def("multiply", (Mat33 (Mat33::*)(const Mat33&) const) &Mat33::multiply)
    .def("multiply", (Vec3 (Mat33::*)(const Vec3&) const) &Mat33::multiply)
    .def("transpose", &Mat33::transpose)
    .def("approx", &Mat33::approx, py::arg("other"), py::arg("epsilon"))
    .def("determinant", &Mat33::determinant)
    .def("inverse", &Mat33::inverse)
    .def("as_list", [](const Mat33& m) -> std::array<std::array<double,3>,3> {
        return {{{{m[0][0], m[0][1], m[0][2]}},
                 {{m[1][0], m[1][1], m[1][2]}},
                 {{m[2][0], m[2][1], m[2][2]}}}};
    })
    .def("__repr__", [](const Mat33& self) {
        const auto& a = self.a;
        return "<gemmi.Mat33 [" + triple(a[0][0], a[0][1], a[0][2]) + "]\n"
               "             [" + triple(a[1][0], a[1][1], a[1][2]) + "]\n"
               "             [" + triple(a[2][0], a[2][1], a[2][2]) + "]>";
    });
  py::class_<Transform>(m, "Transform")
    .def_readonly("mat", &Transform::mat)
    .def_readonly("vec", &Transform::vec)
    .def("inverse", &Transform::inverse)
    .def("apply", &Transform::apply);

  py::class_<Position, Vec3>(m, "Position")
    .def(py::init<double,double,double>())
    .def("dist", [](const Position& self, const Position& other) {
        return self.dist(other);
    })
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Fractional, Vec3>(m, "Fractional")
    .def(py::init<double,double,double>())
    .def("__getitem__", (double (Fractional::*)(int) const) &Fractional::at)
    .def("__repr__", [](const Fractional& self) {
        return "<gemmi.Fractional(" + triple(self.x, self.y, self.z) + ")>";
    });

  py::class_<FTransform, Transform>(m, "FTransform")
    .def("apply", &FTransform::apply);


  py::class_<SymImage> symimage(m, "SymImage");
  symimage
    .def("dist", &SymImage::dist)
    .def("__repr__", [](const SymImage& self) {
        return "<gemmi.SymImage cell:[" +
          triple(self.box[0], self.box[1], self.box[2]) +
          "] sym:" + std::to_string(self.sym_id) + ">";
    });
  py::enum_<SameAsu>(symimage, "Asu")
    .value("Same", SameAsu::Yes)
    .value("Different", SameAsu::No)
    .value("Any", SameAsu::Any);

  py::class_<UnitCell>(m, "UnitCell")
    .def(py::init<>())
    .def(py::init([](double a, double b, double c,
                     double alpha, double beta, double gamma) {
           UnitCell cell;
           cell.set(a, b, c, alpha, beta, gamma);
           return cell;
         }), py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("alpha"), py::arg("beta"), py::arg("gamma"))
    .def_readonly("a", &UnitCell::a)
    .def_readonly("b", &UnitCell::b)
    .def_readonly("c", &UnitCell::c)
    .def_readonly("alpha", &UnitCell::alpha)
    .def_readonly("beta", &UnitCell::beta)
    .def_readonly("gamma", &UnitCell::gamma)
    .def_readonly("volume", &UnitCell::volume)
    .def_readonly("images", &UnitCell::images)
    .def("set", &UnitCell::set)
    .def("is_crystal", &UnitCell::is_crystal)
    .def("fractionalize", &UnitCell::fractionalize)
    .def("orthogonalize", &UnitCell::orthogonalize)
    .def("volume_per_image", &UnitCell::volume_per_image)
    .def("find_nearest_image", &UnitCell::find_nearest_image,
        py::arg("ref"), py::arg("pos"), py::arg("asu")=SameAsu::Any)
    .def("is_special_position", &UnitCell::is_special_position,
         py::arg("pos"), py::arg("max_dist")=0.8)
    .def("calculate_1_d2", &UnitCell::calculate_1_d2,
         py::arg("hk"), py::arg("k"), py::arg("l"))
    .def("calculate_d", &UnitCell::calculate_d,
         py::arg("hk"), py::arg("k"), py::arg("l"))
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });
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

void add_mtz(py::module& m) {
  py::bind_vector<std::vector<Mtz::Column>>(m, "VectorMtzColumn");
  py::bind_vector<std::vector<Mtz::Dataset>>(m, "VectorMtzDataset");

  py::class_<Mtz> mtz(m, "Mtz", py::buffer_protocol());
  mtz.def(py::init<>())
    .def_buffer([](Mtz &self) {
      int nrow = self.has_data() ? self.nreflections : 0;
      return py::buffer_info(self.data.data(),
                             4, py::format_descriptor<float>::format(),
                             2, {nrow, self.ncol}, // dimensions
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
         py::arg("number"))
    .def("count", &Mtz::count, py::arg("label"))
    .def("column_with_label", &Mtz::column_with_label, py::arg("label"),
         py::return_value_policy::reference_internal)
    .def("column_with_type", &Mtz::column_with_type, py::arg("type"),
         py::return_value_policy::reference_internal)
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
    .def("__repr__", [](const Mtz& self) {
        return "<gemmi.Mtz with " + std::to_string(self.ncol) + " columns, " +
               std::to_string(self.nreflections) + " reflections>";
    });
  py::class_<Mtz::Dataset>(mtz, "Dataset")
    .def_readwrite("number", &Mtz::Dataset::number)
    .def_readwrite("project_name", &Mtz::Dataset::project_name)
    .def_readwrite("crystal_name", &Mtz::Dataset::crystal_name)
    .def_readwrite("dataset_name", &Mtz::Dataset::dataset_name)
    .def_readwrite("cell", &Mtz::Dataset::cell)
    .def_readwrite("wavelength", &Mtz::Dataset::wavelength)
    .def("__repr__", [](const Mtz::Dataset& self) {
        return "<gemmi.Mtz.Dataset " + std::to_string(self.number) + " " +
               self.project_name + "/" + self.crystal_name + "/" +
               self.dataset_name + ">";
    });
  py::class_<Mtz::Column>(mtz, "Column", py::buffer_protocol())
    .def_buffer([](Mtz::Column& self) {
      return py::buffer_info(self.parent->data.data() + self.idx,
                             4, py::format_descriptor<float>::format(),
                             1, {self.size()},      // dimensions
                             {4 * self.stride()});  // strides
    })
    .def_property_readonly("array", [](const Mtz::Column& self) {
      return py::array_t<float>({self.size()}, {4 * self.stride()},
                                self.parent->data.data() + self.idx,
                                py::cast(self));
    }, py::return_value_policy::reference_internal)
    .def_readwrite("dataset_number", &Mtz::Column::dataset_number)
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
        return "<gemmi.Mtz.Column " + self.label + " type " + std::string(1, self.type) +
               ">";
    });

  m.def("read_mtz_file", &read_mtz_file);
}
