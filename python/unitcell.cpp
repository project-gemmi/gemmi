// Copyright 2017 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"

#include <cstdio>  // for snprintf
#include <array>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
//#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

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


  py::class_<SymImage>(m, "SymImage")
    .def("dist", &SymImage::dist)
    .def("__repr__", [](const SymImage& self) {
        return "<gemmi.SymImage box:[" +
          triple(self.box[0], self.box[1], self.box[2]) +
          "] sym:" + std::to_string(self.sym_id) + ">";
    });

  py::enum_<Asu>(m, "Asu")
    .value("Same", Asu::Same)
    .value("Different", Asu::Different)
    .value("Any", Asu::Any);

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
         py::arg("ref"), py::arg("pos"), py::arg("asu")=Asu::Any)
    .def("is_special_position",
         (int (UnitCell::*)(const Position&, double) const)
           &UnitCell::is_special_position,
         py::arg("pos"), py::arg("max_dist")=0.8)
    .def("is_special_position",
         (int (UnitCell::*)(const Fractional&, double) const)
           &UnitCell::is_special_position,
         py::arg("fpos"), py::arg("max_dist"))
    .def("calculate_1_d2",
         (double (UnitCell::*)(const Miller&)const) &UnitCell::calculate_1_d2,
         py::arg("hkl"))
    .def("calculate_d", &UnitCell::calculate_d, py::arg("hkl"))
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });
}

