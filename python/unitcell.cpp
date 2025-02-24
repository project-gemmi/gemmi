// Copyright 2017 Global Phasing Ltd.

#include "gemmi/unitcell.hpp"
#include "gemmi/cellred.hpp"  // GruberVector
#include "gemmi/twin.hpp"  // find_lattice_2fold_ops
#include "gemmi/sprintf.hpp"  // snprintf_z

#include <array>
#include "common.h"
#include "array.h"  // numpy_array_from_vector, miller_function
#include "serial.h"  // for getstate, setstate
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>    // for find_lattice_2fold_ops
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for UnitCell::images

using namespace gemmi;

namespace {

std::string triple(double x, double y, double z) {
  char buf[128];
  auto r = [](double d) { return std::fabs(d) < 1e-15 ? 0 : d; };
  snprintf_z(buf, 128, "%g, %g, %g", r(x), r(y), r(z));
  return std::string(buf);
}

void mat33_from_list(Mat33& self, std::array<std::array<double,3>,3>& m) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      self.a[i][j] = m[i][j];
}

nb::tuple make_six_tuple(const std::array<double,6>& v) {
  return nb::make_tuple(v[0], v[1], v[2], v[3], v[4], v[5]);
}

auto mat33_to_array(Mat33& self) {
  return nb::ndarray<nb::numpy, double, nb::shape<3,3>, nb::c_contig>(
          &self.a[0][0], {3, 3}, nb::handle());
}

template<typename T> void add_smat33(nb::module_& m, const char* name) {
  using M = SMat33<T>;
  nb::class_<M>(m, name)
    .def("__init__", [](M* p, T u11, T u22, T u33, T u12, T u13, T u23) {
           new(p) M{u11, u22, u33, u12, u13, u23};
         }, nb::arg("u11"), nb::arg("u22"), nb::arg("u33"),
            nb::arg("u12"), nb::arg("u13"), nb::arg("u23"))
    .def_rw("u11", &M::u11)
    .def_rw("u22", &M::u22)
    .def_rw("u33", &M::u33)
    .def_rw("u12", &M::u12)
    .def_rw("u13", &M::u13)
    .def_rw("u23", &M::u23)
    .def("elements_pdb", &M::elements_pdb)
    .def("elements_voigt", &M::elements_voigt)
    .def("as_mat33", &M::as_mat33)
    .def("trace", &M::trace)
    .def("nonzero", &M::nonzero)
    .def("determinant", &M::determinant)
    .def("inverse", &M::inverse)
    .def("scaled", &M::template scaled<T>)
    .def("added_kI", &M::added_kI)
    .def("r_u_r", (double (M::*)(const Vec3&) const) &M::r_u_r)
    .def("r_u_r", [](const M& self, const cpu_miller_array& arr) {
        std::vector<T> v;
        size_t len = arr.shape(0);
        v.reserve(len);
        auto r = arr.view();
        for (size_t row = 0; row < len; ++row)
           v.push_back((T)self.r_u_r({r(row, 0), r(row, 1), r(row, 2)}));
        return numpy_array_from_vector(std::move(v));
    }, nb::arg().noconvert())
    .def("multiply", &M::multiply)
    .def(nb::self + nb::self)
    .def(nb::self - nb::self)  // NOLINT(misc-redundant-expression)
    .def("transformed_by", &M::template transformed_by<T>)
    .def("calculate_eigenvalues", &M::calculate_eigenvalues)
    .def("__repr__", [name](const M& m) {
        char buf[128];
        snprintf_z(buf, 128, "<gemmi.%s(%g, %g, %g, %g, %g, %g)>",
                   name, m.u11, m.u22, m.u33, m.u12, m.u13, m.u23);
        return std::string(buf);
    });
}

template<typename T> void add_box(nb::module_& m, const char* name) {
  using M = Box<T>;
  nb::class_<M>(m, name)
    .def(nb::init<>())
    .def_rw("minimum", &M::minimum)
    .def_rw("maximum", &M::maximum)
    .def("get_size", &M::get_size)
    .def("extend", &M::extend)
    .def("add_margin", &M::add_margin)
    ;
}

}  // anonymous namespace

void add_unitcell(nb::module_& m) {
  nb::class_<Vec3>(m, "Vec3")
    .def(nb::init<double,double,double>())
    .def_rw("x", &Vec3::x)
    .def_rw("y", &Vec3::y)
    .def_rw("z", &Vec3::z)
    .def("normalized", &Vec3::normalized)
    .def("dot", &Vec3::dot)
    .def("cross", &Vec3::cross)
    .def("length", &Vec3::length)
    .def("approx", &Vec3::approx, nb::arg("other"), nb::arg("epsilon"))
    .def("tolist", [](const Vec3& self) {
        return std::array<double,3>{{self.x, self.y, self.z}};
    })
    .def("fromlist", [](Vec3& self, std::array<double,3>& v) {
        self.x = v[0];
        self.y = v[1];
        self.z = v[2];
    })
    .def(nb::self + nb::self)
    .def(nb::self - nb::self)  // NOLINT(misc-redundant-expression)
    .def(nb::self += nb::self, nb::rv_policy::none)
    //.def(nb::self -= nb::self)  // Clang warning -Wself-assign-overloaded
    .def(operator-=(nb::self, nb::self), nb::rv_policy::none)
    .def(nb::self * float())
    .def(nb::self *= float(), nb::rv_policy::none)
    .def(float() * nb::self)
    .def(nb::self / float())
    .def(nb::self /= float(), nb::rv_policy::none)
    .def(-nb::self)
    .def("__getitem__", (double (Vec3::*)(int) const) &Vec3::at)
    .def("__setitem__", [](Vec3& self, int idx, double value) {
        self.at(idx) = value;
    })
    .def("__repr__", [](const Vec3& self) {
        return "<gemmi.Vec3(" + triple(self.x, self.y, self.z) + ")>";
    });
  nb::class_<Mat33> mat33(m, "Mat33");
  mat33
    .def(nb::init<>())
    .def("__init__", [](Mat33* mat, std::array<std::array<double,3>,3>& arr) {
        new(mat) Mat33();
        mat33_from_list(*mat, arr);
    })
    .def_prop_ro("array", &mat33_to_array, nb::rv_policy::reference_internal)
    .def("__array__", [](nb::handle_t<Mat33>& h, nb::handle dtype, nb::handle copy) {
        return handle_numpy_array_args(h.attr("array"), dtype, copy);
    }, nb::arg("dtype")=nb::none(), nb::arg("copy")=nb::none())

    .def("row_copy", &Mat33::row_copy)
    .def("column_copy", &Mat33::column_copy)
    .def(nb::self + nb::self)
    .def(nb::self - nb::self)  // NOLINT(misc-redundant-expression)
    .def("multiply", (Mat33 (Mat33::*)(const Mat33&) const) &Mat33::multiply)
    .def("multiply", (Vec3 (Mat33::*)(const Vec3&) const) &Mat33::multiply)
    .def("__matmul__", (Mat33 (Mat33::*)(const Mat33&) const) &Mat33::multiply, nb::is_operator())
    .def("__matmul__", (Vec3 (Mat33::*)(const Vec3&) const) &Mat33::multiply, nb::is_operator())
    .def("left_multiply", &Mat33::left_multiply)
    .def("multiply_by_diagonal", &Mat33::multiply_by_diagonal)
    .def("transpose", &Mat33::transpose)
    .def("trace", &Mat33::trace)
    .def("approx", &Mat33::approx, nb::arg("other"), nb::arg("epsilon"))
    .def("determinant", &Mat33::determinant)
    .def("inverse", &Mat33::inverse)
    .def("is_identity", &Mat33::is_identity)
    .def("tolist", [](const Mat33& mat) -> std::array<std::array<double,3>,3> {
        return {{{{mat[0][0], mat[0][1], mat[0][2]}},
                 {{mat[1][0], mat[1][1], mat[1][2]}},
                 {{mat[2][0], mat[2][1], mat[2][2]}}}};
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

  nb::class_<Transform> transform(m, "Transform");
  transform
    .def(nb::init<>())
    .def("__init__", [](Transform* tr, const Mat33& m, const Vec3& v) {
        new(tr) Transform();
        tr->mat = m;
        tr->vec = v;
    }, nb::arg("mat33"), nb::arg("vec3"))
    .def_ro("mat", &Transform::mat)
    .def_ro("vec", &Transform::vec)
    .def("inverse", &Transform::inverse)
    .def("apply", &Transform::apply)
    .def("combine", &Transform::combine)
    .def("__matmul__", &Transform::combine, nb::is_operator())
    .def("is_identity", &Transform::is_identity)
    .def("approx", &Transform::approx, nb::arg("other"), nb::arg("epsilon"));

  nb::class_<Position, Vec3>(m, "Position")
    .def(nb::init<double,double,double>())
    .def(nb::init<const Vec3&>())
    .def("dist", [](const Position& self, const Position& other) {
        return self.dist(other);
    })
    .def(nb::self + nb::self)
    .def(nb::self - nb::self)  // NOLINT(misc-redundant-expression)
    .def(nb::self += nb::self, nb::rv_policy::none)
    .def(operator-=(nb::self, nb::self), nb::rv_policy::none)
    .def(nb::self * float())
    .def(nb::self *= float(), nb::rv_policy::none)
    .def(float() * nb::self)
    .def(nb::self / float())
    .def(nb::self /= float(), nb::rv_policy::none)
    .def(-nb::self)
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + triple(self.x, self.y, self.z) + ")>";
    });
  nb::class_<Fractional, Vec3>(m, "Fractional")
    .def(nb::init<double,double,double>())
    .def(nb::init<const Vec3&>())
    .def("wrap_to_unit", &Fractional::wrap_to_unit)
    .def("wrap_to_zero", &Fractional::wrap_to_zero)
    .def("__getitem__", (double (Fractional::*)(int) const) &Fractional::at)
    .def(nb::self + nb::self)
    .def(nb::self - nb::self)  // NOLINT(misc-redundant-expression)
    .def("__repr__", [](const Fractional& self) {
        return "<gemmi.Fractional(" + triple(self.x, self.y, self.z) + ")>";
    });

  add_box<Position>(m, "PositionBox");
  add_box<Fractional>(m, "FractionalBox");

  nb::class_<FTransform, Transform>(m, "FTransform")
    .def(nb::init<>())
    .def("apply", &FTransform::apply);

  nb::class_<NearestImage>(m, "NearestImage")
    .def("dist", &NearestImage::dist)
    .def("same_asu", &NearestImage::same_asu)
    .def("symmetry_code", &NearestImage::symmetry_code, nb::arg("underscore")=true)
    .def_ro("sym_idx", &NearestImage::sym_idx)
    .def_prop_ro("pbc_shift", [](const NearestImage& self) {
        return nb::make_tuple(self.pbc_shift[0], self.pbc_shift[1], self.pbc_shift[2]);
    })
    .def("__repr__", [](const NearestImage& self) {
        char buf[64];
        snprintf_z(buf, 64, "<gemmi.NearestImage %s in distance %.2f>",
                   self.symmetry_code(true).c_str(), self.dist());
        return std::string(buf);
    });

  nb::enum_<Asu>(m, "Asu")
    .value("Same", Asu::Same)
    .value("Different", Asu::Different)
    .value("Any", Asu::Any);

  nb::class_<UnitCell>(m, "UnitCell")
    .def(nb::init<>())
    .def(nb::init<double,double,double,double,double,double>(),
         nb::arg("a"), nb::arg("b"), nb::arg("c"),
         nb::arg("alpha"), nb::arg("beta"), nb::arg("gamma"))
    .def_ro("a", &UnitCell::a)
    .def_ro("b", &UnitCell::b)
    .def_ro("c", &UnitCell::c)
    .def_ro("alpha", &UnitCell::alpha)
    .def_ro("beta", &UnitCell::beta)
    .def_ro("gamma", &UnitCell::gamma)
    .def_ro("volume", &UnitCell::volume)
    .def_ro("explicit_matrices", &UnitCell::explicit_matrices)
    .def_ro("images", &UnitCell::images)
    .def_prop_ro("parameters", [](const UnitCell& u) {
        return nb::make_tuple(u.a, u.b, u.c, u.alpha, u.beta, u.gamma);
    })

    .def_ro("frac", &UnitCell::frac)
    .def_ro("orth", &UnitCell::orth)
    .def("set", &UnitCell::set)
    .def("changed_basis_forward", &UnitCell::changed_basis_forward,
         nb::arg("op"), nb::arg("set_images"))
    .def("changed_basis_backward", &UnitCell::changed_basis_backward,
         nb::arg("op"), nb::arg("set_images"))
    .def("is_compatible_with_spacegroup", &UnitCell::is_compatible_with_spacegroup,
         nb::arg("sg"), nb::arg("eps")=1e-3)
    .def("is_crystal", &UnitCell::is_crystal)
    .def("approx", [](const UnitCell& self, const UnitCell& o, double epsilon) {
        return self.approx(o, epsilon);
    }, nb::arg("other"), nb::arg("epsilon"))
    .def("is_similar", &UnitCell::is_similar,
         nb::arg("other"), nb::arg("rel"), nb::arg("deg"))
    .def("calculate_u_eq", &UnitCell::calculate_u_eq)
    .def("fractionalize", &UnitCell::fractionalize)
    .def("orthogonalize", &UnitCell::orthogonalize)
    .def("orthogonalize_box", &UnitCell::orthogonalize_box)
    .def("op_as_transform", &UnitCell::op_as_transform)
    .def("volume_per_image", &UnitCell::volume_per_image)
    .def("find_nearest_image", &UnitCell::find_nearest_image,
         nb::arg("ref"), nb::arg("pos"), nb::arg("asu")=Asu::Any)
    .def("find_nearest_pbc_image",
         (NearestImage (UnitCell::*)(const Fractional&, Fractional, int) const)
         &UnitCell::find_nearest_pbc_image,
         nb::arg("fref"), nb::arg("fpos"), nb::arg("image_idx")=0)
    .def("find_nearest_pbc_image",
         (NearestImage (UnitCell::*)(const Position&, const Position&, int) const)
         &UnitCell::find_nearest_pbc_image,
         nb::arg("ref"), nb::arg("pos"), nb::arg("image_idx")=0)
    .def("find_nearest_pbc_position", &UnitCell::find_nearest_pbc_position,
         nb::arg("ref"), nb::arg("pos"), nb::arg("image_idx"), nb::arg("inverse")=false)
    .def("is_special_position",
         (int (UnitCell::*)(const Position&, double) const)
           &UnitCell::is_special_position,
         nb::arg("pos"), nb::arg("max_dist")=0.8)
    .def("is_special_position",
         (int (UnitCell::*)(const Fractional&, double) const)
           &UnitCell::is_special_position,
         nb::arg("fpos"), nb::arg("max_dist"))
    .def("calculate_1_d2", &UnitCell::calculate_1_d2, nb::arg("hkl"))
    .def("calculate_1_d2_array", [](const UnitCell& u, const cpu_miller_array& hkl) {
        return miller_function<double>(u, &UnitCell::calculate_1_d2, hkl);
    })
    .def("calculate_d", &UnitCell::calculate_d, nb::arg("hkl"))
    .def("calculate_d_array", [](const UnitCell& u, const cpu_miller_array& hkl) {
        return miller_function<double>(u, &UnitCell::calculate_d, hkl);
    })
    .def("metric_tensor", &UnitCell::metric_tensor)
    .def("reciprocal_metric_tensor", &UnitCell::reciprocal_metric_tensor)
    .def("reciprocal", &UnitCell::reciprocal)
    .def("get_hkl_limits", &UnitCell::get_hkl_limits, nb::arg("dmin"))
    .def("primitive_orth_matrix", &UnitCell::primitive_orth_matrix, nb::arg("centring_type"))
    // NOLINTNEXTLINE(misc-redundant-expression)
    .def(nb::self == nb::self, nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def("__getstate__", &getstate<UnitCell>)
    .def("__setstate__", &setstate<UnitCell>)
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });

  nb::class_<SellingVector> selling_vector(m, "SellingVector");

  nb::class_<GruberVector>(m, "GruberVector")
    .def(nb::init<const std::array<double,6>&>())
    .def(nb::init<const UnitCell&, const SpaceGroup*, bool>(),
         nb::arg("cell"), nb::arg("sg").none(), nb::arg("track_change_of_basis")=false)
    .def(nb::init<const UnitCell&, char, bool>(),
         nb::arg("cell"), nb::arg("centring"), nb::arg("track_change_of_basis")=false)
    .def_prop_ro("parameters", [](const GruberVector& g) {
      return make_six_tuple(g.parameters());
    })
    .def("cell_parameters", [](const GruberVector& self) {
      return make_six_tuple(self.cell_parameters());
    })
    .def("get_cell",
         [](const GruberVector& self) { return new UnitCell(self.get_cell()); })
    .def_prop_ro("change_of_basis", [](const GruberVector& self) {
        return self.change_of_basis.get();
    }, nb::rv_policy::reference_internal)
    .def("selling", &GruberVector::selling)
    .def("is_normalized", &GruberVector::is_normalized)
    .def("is_buerger", &GruberVector::is_buerger, nb::arg("epsilon")=1e-9)
    .def("normalize", &GruberVector::normalize, nb::arg("epsilon")=1e-9)
    .def("buerger_reduce", &GruberVector::buerger_reduce)
    .def("niggli_step", &GruberVector::niggli_step, nb::arg("epsilon"))
    .def("niggli_reduce", &GruberVector::niggli_reduce,
         nb::arg("epsilon")=1e-9, nb::arg("iteration_limit")=100)
    .def("is_niggli", &GruberVector::is_niggli, nb::arg("epsilon")=1e-9)
    .def("__repr__", [](const GruberVector& self) {
        char buf[256];
        snprintf_z(buf, 256, "<gemmi.GruberVector((%.2f, %.2f, %.2f, %.2f, %.2f, %.2f))>",
                   self.A, self.B, self.C, self.xi, self.eta, self.zeta);
        return std::string(buf);
    });

  selling_vector
    .def(nb::init<const std::array<double,6>&>())
    .def(nb::init<const UnitCell&, const SpaceGroup*>())
    .def_prop_ro("parameters", [](const SellingVector& self) {
      return make_six_tuple(self.s);
    })
    .def("cell_parameters", [](const SellingVector& self) {
      return make_six_tuple(self.cell_parameters());
    })
    .def("get_cell",
         [](const SellingVector& self) { return new UnitCell(self.get_cell()); })
    .def("sum_b_squared", &SellingVector::sum_b_squared)
    .def("gruber", &SellingVector::gruber)
    .def("is_reduced", &SellingVector::is_reduced, nb::arg("epsilon")=1e-9)
    .def("reduce_step", &SellingVector::reduce_step, nb::arg("epsilon")=1e-9)
    .def("reduce", &SellingVector::reduce,
         nb::arg("epsilon")=1e-9, nb::arg("iteration_limit")=100)
    .def("sort", &SellingVector::sort, nb::arg("epsilon")=1e-9)
    .def("__repr__", [](const SellingVector& self) {
        char buf[256];
        snprintf_z(buf, 256, "<gemmi.SellingVector((%.2f, %.2f, %.2f, %.2f, %.2f, %.2f))>",
                   self.s[0], self.s[1], self.s[2], self.s[3], self.s[4], self.s[5]);
        return std::string(buf);
    });

  // twin.hpp
  m.def("find_lattice_2fold_ops", &find_lattice_2fold_ops,
        nb::arg("reduced_cell"), nb::arg("max_obliq"));
  m.def("find_lattice_symmetry_r", &find_lattice_symmetry_r);
  m.def("find_lattice_symmetry", &find_lattice_symmetry,
        nb::arg("cell"), nb::arg("centring"), nb::arg("max_obliq"));
  m.def("find_twin_laws", &find_twin_laws,
        nb::arg("cell"), nb::arg("sg"), nb::arg("max_obliq"), nb::arg("all_ops"));
}
