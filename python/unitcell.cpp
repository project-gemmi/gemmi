// Copyright 2017 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"

#include <cstdio>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace gemmi;

static std::string triple(double x, double y, double z) {
  using namespace std;  // VS2015/17 doesn't like std::snprintf
  char buf[128];
  snprintf(buf, 128, "%g, %g, %g", x, y, z);
  return std::string(buf);
}

void add_unitcell(py::module& m) {
  py::class_<Position>(m, "Position")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Position::x)
    .def_readwrite("y", &Position::y)
    .def_readwrite("z", &Position::z)
    .def("dist", [](const Position& self, const Position& other) {
        return self.dist(other);
    })
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Fractional>(m, "Fractional")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Fractional::x)
    .def_readwrite("y", &Fractional::y)
    .def_readwrite("z", &Fractional::z)
    .def("__repr__", [](const Fractional& self) {
        return "<gemmi.Fractional(" + triple(self.x, self.y, self.z) + ")>";
    });

  py::enum_<SymmetryImage>(m, "SymmetryImage")
    .value("Same", SymmetryImage::Same)
    .value("Different", SymmetryImage::Different)
    .value("Unspecified", SymmetryImage::Unspecified);

  py::class_<NearbyImage>(m, "NearbyImage")
    .def("dist", &NearbyImage::dist)
    .def("__repr__", [](const NearbyImage& self) {
        return "<gemmi.NearbyImage cell:[" +
          triple(self.box[0], self.box[1], self.box[2]) +
          "] sym:" + std::to_string(self.sym_id) + ">";
    });

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
    .def("set", &UnitCell::set)
    .def("fractionalize", &UnitCell::fractionalize)
    .def("orthogonalize", &UnitCell::orthogonalize)
    .def("volume_per_image", &UnitCell::volume_per_image)
    .def("find_nearest_image", &UnitCell::find_nearest_image,
        py::arg("ref"), py::arg("pos"),
        py::arg("sym_image")=SymmetryImage::Unspecified)
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });
}
