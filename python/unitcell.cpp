// Copyright 2017 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/seqid.hpp"

#include <cstdio>  // for snprintf
#include <array>
#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include "miller_a.h"

namespace py = pybind11;
using namespace gemmi;

static std::string triple(double x, double y, double z) {
  using namespace std;  // VS2015/17 doesn't like std::snprintf
  char buf[128];
  snprintf(buf, 128, "%g, %g, %g", x, y, z);
  return std::string(buf);
}

static std::vector<std::string>
expand_protein_one_letter_string(const std::string& s) {
  std::vector<std::string> r;
  r.reserve(s.size());
  for (char c : s)
    r.push_back(expand_protein_one_letter(c));
  return r;
}
void mat33_from_list(Mat33& self, std::array<std::array<double,3>,3>& m) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      self.a[i][j] = m[i][j];
}

template<typename T> void add_smat33(py::module& m, const char* name) {
  using M = SMat33<T>;
  py::class_<M>(m, name)
    .def(py::init([](T u11, T u22, T u33, T u12, T u13, T u23) {
           return M{u11, u22, u33, u12, u13, u23};
         }), py::arg("u11"), py::arg("u22"), py::arg("u33"),
             py::arg("u12"), py::arg("u13"), py::arg("u23"))
    .def_readwrite("u11", &M::u11)
    .def_readwrite("u22", &M::u22)
    .def_readwrite("u33", &M::u33)
    .def_readwrite("u12", &M::u12)
    .def_readwrite("u13", &M::u13)
    .def_readwrite("u23", &M::u23)
    .def("elements", &M::elements)
    .def("as_mat33", &M::as_mat33)
    .def("trace", &M::trace)
    .def("nonzero", &M::nonzero)
    .def("determinant", &M::determinant)
    .def("inverse", &M::inverse)
    .def("r_u_r", &M::r_u_r)
    .def("r_u_r", [](const M& self, py::array_t<int> arr) {
        int nrow = (int) arr.shape(0);
        int ncol = (int) arr.shape(1);
        if (ncol != 3)
           fail("SMat33::r_u_r(): expected 3 columns.");
        std::vector<T> v;
        v.reserve(nrow);
        auto r = arr.unchecked<2>();
        for (ssize_t row = 0; row < nrow; ++row)
           v.push_back((T)self.r_u_r(Vec3(r(row, 0), r(row, 1), r(row, 2))));
        return py_array_from_vector(std::move(v));
    }, py::arg().noconvert())
    .def("transformed_by", &M::template transformed_by<double>)
    .def("calculate_eigenvalues", &M::calculate_eigenvalues)
    ;
}

template<typename T> void add_box(py::module& m, const char* name) {
  using M = Box<T>;
  py::class_<M>(m, name)
    .def(py::init<>())
    .def_readwrite("minimum", &M::minimum)
    .def_readwrite("maximum", &M::maximum)
    .def("get_size", &M::get_size)
    .def("extend", &M::extend)
    .def("add_margin", &M::add_margin)
    ;
}

void add_unitcell(py::module& m) {
  py::class_<Vec3>(m, "Vec3")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Vec3::x)
    .def_readwrite("y", &Vec3::y)
    .def_readwrite("z", &Vec3::z)
    .def("dot", &Vec3::dot)
    .def("cross", &Vec3::cross)
    .def("approx", &Vec3::approx, py::arg("other"), py::arg("epsilon"))
    .def("tolist", [](const Vec3& self) {
        return std::array<double,3>{{self.x, self.y, self.z}};
    })
    .def("fromlist", [](Vec3& self, std::array<double,3>& v) {
        self.x = v[0];
        self.y = v[1];
        self.z = v[2];
    })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    //.def(py::self -= py::self)  // Clang warning -Wself-assign-overloaded
    .def(operator-=(py::self, py::self))
    .def(py::self * float())
    .def(py::self *= float())
    .def(float() * py::self)
    .def(py::self / float())
    .def(py::self /= float())
    .def(-py::self)
    .def("__getitem__", (double (Vec3::*)(int) const) &Vec3::at)
    .def("__setitem__", [](Vec3& self, int idx, double value) {
        self.at(idx) = value;
    })
    .def("__repr__", [](const Vec3& self) {
        return "<gemmi.Vec3(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Mat33>(m, "Mat33", py::buffer_protocol())
    .def(py::init<>())
    .def(py::init([](std::array<std::array<double,3>,3>& m) {
       Mat33 *mat = new Mat33();
       mat33_from_list(*mat, m);
       return mat;
     }))
    .def_buffer([](Mat33 &self) {
      return py::buffer_info(&self.a[0][0],
                             {3, 3}, // dimensions
                             {sizeof(double)*3, sizeof(double)});  // strides
    })
    .def("multiply", (Mat33 (Mat33::*)(const Mat33&) const) &Mat33::multiply)
    .def("multiply", (Vec3 (Mat33::*)(const Vec3&) const) &Mat33::multiply)
    .def("left_multiply", &Mat33::left_multiply)
    .def("transpose", &Mat33::transpose)
    .def("approx", &Mat33::approx, py::arg("other"), py::arg("epsilon"))
    .def("determinant", &Mat33::determinant)
    .def("inverse", &Mat33::inverse)
    .def("is_identity", &Mat33::is_identity)
    .def("tolist", [](const Mat33& m) -> std::array<std::array<double,3>,3> {
        return {{{{m[0][0], m[0][1], m[0][2]}},
                 {{m[1][0], m[1][1], m[1][2]}},
                 {{m[2][0], m[2][1], m[2][2]}}}};
    })
    .def("fromlist", mat33_from_list)
    .def("__repr__", [](const Mat33& self) {
        const auto& a = self.a;
        return "<gemmi.Mat33 [" + triple(a[0][0], a[0][1], a[0][2]) + "]\n"
               "             [" + triple(a[1][0], a[1][1], a[1][2]) + "]\n"
               "             [" + triple(a[2][0], a[2][1], a[2][2]) + "]>";
    });

  add_smat33<double>(m, "SMat33d");
  add_smat33<float>(m, "SMat33f");

  py::class_<Transform>(m, "Transform")
    .def(py::init<>())
    .def(py::init([](const Mat33& m, const Vec3& v) {
      Transform* tr = new Transform();
      tr->mat = m;
      tr->vec = v;
      return tr;
    }), py::arg("mat33"), py::arg("vec3"))
    .def_readonly("mat", &Transform::mat)
    .def_readonly("vec", &Transform::vec)
    .def("inverse", &Transform::inverse)
    .def("apply", &Transform::apply)
    .def("is_identity", &Transform::is_identity)
    .def("approx", &Transform::approx);

  py::class_<Position, Vec3>(m, "Position")
    .def(py::init<double,double,double>())
    .def(py::init<const Vec3&>())
    .def("dist", [](const Position& self, const Position& other) {
        return self.dist(other);
    })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    .def(operator-=(py::self, py::self))
    .def(py::self * float())
    .def(py::self *= float())
    .def(float() * py::self)
    .def(py::self / float())
    .def(py::self /= float())
    .def(-py::self)
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Fractional, Vec3>(m, "Fractional")
    .def(py::init<double,double,double>())
    .def("__getitem__", (double (Fractional::*)(int) const) &Fractional::at)
    .def("__repr__", [](const Fractional& self) {
        return "<gemmi.Fractional(" + triple(self.x, self.y, self.z) + ")>";
    });

  add_box<Position>(m, "PositionBox");
  add_box<Fractional>(m, "FractionalBox");

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
    .def(py::init<double,double,double,double,double,double>(),
         py::arg("a"), py::arg("b"), py::arg("c"),
         py::arg("alpha"), py::arg("beta"), py::arg("gamma"))
    .def_readonly("a", &UnitCell::a)
    .def_readonly("b", &UnitCell::b)
    .def_readonly("c", &UnitCell::c)
    .def_readonly("alpha", &UnitCell::alpha)
    .def_readonly("beta", &UnitCell::beta)
    .def_readonly("gamma", &UnitCell::gamma)
    .def_readonly("volume", &UnitCell::volume)
    .def_readonly("images", &UnitCell::images)
    .def_property_readonly("parameters", [](const UnitCell& self) {
        return py::make_tuple(self.a, self.b, self.c, self.alpha, self.beta, self.gamma);
    })
    .def_property_readonly("fractionalization_matrix",
                           [](const UnitCell& self) { return self.frac.mat; })
    .def_property_readonly("orthogonalization_matrix",
                           [](const UnitCell& self) { return self.orth.mat; })
    .def("set", &UnitCell::set)
    .def("is_crystal", &UnitCell::is_crystal)
    .def("calculate_u_eq", &UnitCell::calculate_u_eq)
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
    .def("calculate_1_d2", &UnitCell::calculate_1_d2, py::arg("hkl"))
    .def("calculate_1_d2_array", [](const UnitCell& u, py::array_t<int> hkl) {
        return miller_function<double>(u, &UnitCell::calculate_1_d2, hkl);
    })
    .def("calculate_d", &UnitCell::calculate_d, py::arg("hkl"))
    .def("calculate_d_array", [](const UnitCell& u, py::array_t<int> hkl) {
        return miller_function<double>(u, &UnitCell::calculate_d, hkl);
    })
    .def("metric_tensor", &UnitCell::metric_tensor)
    .def("reciprocal_metric_tensor", &UnitCell::reciprocal_metric_tensor)
    .def("reciprocal", &UnitCell::reciprocal)
    .def("get_hkl_limits", &UnitCell::get_hkl_limits, py::arg("dmin"))
    .def(py::self == py::self)
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });


  // resinfo.hpp
  py::class_<ResidueInfo>(m, "ResidueInfo")
    .def_readonly("one_letter_code", &ResidueInfo::one_letter_code)
    .def_readonly("hydrogen_count", &ResidueInfo::hydrogen_count)
    .def_readonly("weight", &ResidueInfo::weight)
    .def("found", &ResidueInfo::found)
    .def("is_standard", &ResidueInfo::is_standard)
    .def("is_water", &ResidueInfo::is_water)
    .def("is_nucleic_acid", &ResidueInfo::is_nucleic_acid)
    .def("is_amino_acid", &ResidueInfo::is_amino_acid);

  m.def("find_tabulated_residue", &find_tabulated_residue, py::arg("name"),
        "Find chemical component information in the internal table.");
  m.def("expand_protein_one_letter", &expand_protein_one_letter);
  m.def("expand_protein_one_letter_string", &expand_protein_one_letter_string);
}
