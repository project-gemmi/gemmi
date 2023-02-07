// Copyright 2023 MRC Laboratory of Molecular Biology

#include "gemmi/refine/geom.hpp"    // for Geometry
#include "gemmi/refine/ll.hpp"    // for LL
#include "gemmi/it92.hpp"

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>  // for detail::pythonbuf

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Bond::Value>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Angle>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Angle::Value>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Torsion>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Torsion::Value>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Chirality>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Plane>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Interval>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Stacking>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Harmonic>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Special>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Vdw>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::bond_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::angle_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::torsion_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::chiral_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::plane_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::stacking_reporting_t>)
PYBIND11_MAKE_OPAQUE(std::vector<Geometry::Reporting::vdw_reporting_t>)

py::tuple for_coo_matrix(GeomTarget &self) {
  const size_t n = self.am.size();
  py::array_t<int> rowarr(n);
  py::array_t<int> colarr(n);
  int* row = (int*) rowarr.request().ptr;
  int* col = (int*) colarr.request().ptr;
  self.get_am_col_row(row, col);
  return py::make_tuple(&self.am, py::make_tuple(rowarr, colarr));
}

template<typename Table>
py::tuple for_coo_matrix(LL<Table> &self) {
  const auto am = self.fisher_diag_from_table();
  const size_t n = am.size();
  py::array_t<int> rowarr(n);
  py::array_t<int> colarr(n);
  int* row = (int*) rowarr.request().ptr;
  int* col = (int*) colarr.request().ptr;
  self.get_am_col_row(row, col);
  return py::make_tuple(am, py::make_tuple(rowarr, colarr));
}

py::tuple precondition_eigen_coo(py::array_t<double> am, py::array_t<int> rows,
                                 py::array_t<int> cols, int N, double cutoff) {
  int* colp = (int*) cols.request().ptr;
  int* rowp = (int*) rows.request().ptr;
  double* amp = (double*) am.request().ptr;
  auto len = cols.shape(0);

  //std::vector<gemmi::SMat33<double>> blocks(N);
  std::vector<double> blocks(2*N);
  for(int i = 0; i < len; ++i) {
    const int c = colp[i], r = rowp[i];
    const int b = c % 3, j = c / 3;
    int k;
    if (r < c - b || r > c) continue;
    if (c == r) k = b;
    else if (b == 1 && r == c - 1) k = 3;
    else if (b == 2 && r == c - 2) k = 4;
    else k = 5; //if (b == 2 && r == c - 1) k = 5;
    blocks[j*6+k] = amp[i];
  }

  std::vector<double> ret(N * 3);
  std::vector<int> retrow(N * 3), retcol(N * 3);
  for (int i = 0; i < N / 3; ++i) {
    const gemmi::SMat33<double> m{blocks[i*6], blocks[i*6+1], blocks[i*6+2], blocks[i*6+3], blocks[i*6+4], blocks[i*6+5]};
    const gemmi::Mat33 pinv = gemmi::eigen_decomp_inv(m, cutoff, true);
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        ret[i*9+3*j+k] = pinv[j][k];
        retrow[i*9+3*j+k] = 3 * i + j;
        retcol[i*9+3*j+k] = 3 * i + k;
      }
  }

  return py::make_tuple(ret, py::make_tuple(retrow, retcol));
}

template<typename Table>
void add_ll(py::module& m, const char* name) {
  using T = gemmi::LL<Table>;
  py::class_<T>(m, name)
    .def(py::init<gemmi::UnitCell, gemmi::SpaceGroup*, std::vector<gemmi::Atom*>, bool, bool, int, bool>(),
         py::arg("cell"), py::arg("sg"), py::arg("atoms"), py::arg("mott_bethe"),
         py::arg("refine_xyz"), py::arg("adp_mode"), py::arg("refine_h"))
    .def("set_ncs", &T::set_ncs)
    .def("calc_grad", &T::calc_grad)
    .def("make_fisher_table_diag_fast", &T::make_fisher_table_diag_fast)
    .def("fisher_diag_from_table", &T::fisher_diag_from_table)
    .def("fisher_for_coo", [](T &self) {return for_coo_matrix(self);}, py::return_value_policy::reference_internal)
    .def_readonly("table_bs", &T::table_bs)
    .def_readonly("pp1", &T::pp1)
    .def_readonly("bb", &T::bb)
    .def_readonly("aa", &T::aa)
  ;
}

void add_refine(py::module& m) {
  py::class_<GeomTarget> geomtarget(m, "GeomTarget");
  py::class_<Geometry> geom(m, "Geometry");

  py::class_<Geometry::Reporting>(geom, "Reporting")
    .def_readonly("bonds", &Geometry::Reporting::bonds)
    .def_readonly("angles", &Geometry::Reporting::angles)
    .def_readonly("torsions", &Geometry::Reporting::torsions)
    .def_readonly("chirs", &Geometry::Reporting::chirs)
    .def_readonly("planes", &Geometry::Reporting::planes)
    .def_readonly("stackings", &Geometry::Reporting::stackings)
    .def_readonly("vdws", &Geometry::Reporting::vdws)
  ;
  py::class_<Geometry::Bond> bond(geom, "Bond");
  py::class_<Geometry::Angle> angle(geom, "Angle");
  py::class_<Geometry::Torsion> torsion(geom, "Torsion");
  py::class_<Geometry::Chirality> chirality(geom, "Chirality");
  py::class_<Geometry::Plane> plane(geom, "Plane");
  py::class_<Geometry::Vdw> vdw(geom, "Vdw");
  py::class_<Geometry::Bond::Value>(bond, "Value")
    .def(py::init<double,double,double,double>())
    .def_readwrite("value", &Geometry::Bond::Value::value)
    .def_readwrite("sigma", &Geometry::Bond::Value::sigma)
    .def_readwrite("value_nucleus", &Geometry::Bond::Value::value_nucleus)
    .def_readwrite("sigma_nucleus", &Geometry::Bond::Value::sigma_nucleus)
    ;
  py::class_<Geometry::Angle::Value>(angle, "Value")
    .def(py::init<double,double>())
    .def_readwrite("value", &Geometry::Angle::Value::value)
    .def_readwrite("sigma", &Geometry::Angle::Value::sigma)
    ;
  py::class_<Geometry::Torsion::Value>(torsion, "Value")
    .def(py::init<double,double,int>())
    .def_readwrite("value", &Geometry::Torsion::Value::value)
    .def_readwrite("sigma", &Geometry::Torsion::Value::sigma)
    .def_readwrite("period", &Geometry::Torsion::Value::period)
    ;
  bond
    .def(py::init<Atom*,Atom*>())
    .def("set_image", &Geometry::Bond::set_image)
    .def_readwrite("type", &Geometry::Bond::type)
    .def_readwrite("alpha", &Geometry::Bond::alpha)
    .def_readwrite("sym_idx", &Geometry::Bond::sym_idx)
    .def_readwrite("pbc_shift", &Geometry::Bond::pbc_shift)
    .def_readwrite("atoms", &Geometry::Bond::atoms)
    .def_readwrite("values", &Geometry::Bond::values)
    ;
  angle
    .def(py::init<Atom*,Atom*,Atom*>())
    .def_readwrite("atoms", &Geometry::Angle::atoms)
    .def_readwrite("values", &Geometry::Angle::values)
    ;
  torsion
    .def(py::init<Atom*,Atom*,Atom*,Atom*>())
    .def_readwrite("atoms", &Geometry::Torsion::atoms)
    .def_readwrite("values", &Geometry::Torsion::values)
    ;
  chirality
    .def(py::init<Atom*,Atom*,Atom*,Atom*>())
    .def_readwrite("value", &Geometry::Chirality::value)
    .def_readwrite("sigma", &Geometry::Chirality::sigma)
    .def_readwrite("sign", &Geometry::Chirality::sign)
    .def_readwrite("atoms", &Geometry::Chirality::atoms)
    ;
  plane
    .def(py::init<std::vector<Atom*>>())
    .def_readwrite("sigma", &Geometry::Plane::sigma)
    .def_readwrite("atoms", &Geometry::Plane::atoms)
    ;
  py::class_<Geometry::Interval>(geom, "Interval")
    .def(py::init<Atom*,Atom*>())
    .def_readwrite("dmin", &Geometry::Interval::dmin)
    .def_readwrite("dmax", &Geometry::Interval::dmax)
    .def_readwrite("smin", &Geometry::Interval::smin)
    .def_readwrite("smax", &Geometry::Interval::smax)
    .def_readwrite("atoms", &Geometry::Interval::atoms)
    ;
  py::class_<Geometry::Harmonic>(geom, "Harmonic")
    .def(py::init<Atom*>())
    .def_readwrite("sigma", &Geometry::Harmonic::sigma)
    .def_readwrite("atom", &Geometry::Harmonic::atom)
    ;
  py::class_<Geometry::Special>(geom, "Special")
    .def(py::init<Atom*>())
    .def_readwrite("sigma_t", &Geometry::Special::sigma_t)
    .def_readwrite("sigma_u", &Geometry::Special::sigma_u)
    .def_readwrite("u_val_incl", &Geometry::Special::u_val_incl)
    .def_readwrite("trans_t", &Geometry::Special::trans_t)
    .def_readwrite("mat_u", &Geometry::Special::mat_u)
    .def_readwrite("atom", &Geometry::Special::atom)
    ;
  py::class_<Geometry::Stacking>(geom, "Stacking")
    .def(py::init<std::vector<Atom*>,std::vector<Atom*>>())
    .def_readwrite("dist", &Geometry::Stacking::dist)
    .def_readwrite("sd_dist", &Geometry::Stacking::sd_dist)
    .def_readwrite("angle", &Geometry::Stacking::angle)
    .def_readwrite("sd_angle", &Geometry::Stacking::sd_angle)
    .def_readwrite("planes", &Geometry::Stacking::planes)
    ;
  vdw
    .def(py::init<Atom*,Atom*>())
    .def("set_image", &Geometry::Vdw::set_image)
    .def_readwrite("type", &Geometry::Vdw::type)
    .def_readwrite("value", &Geometry::Vdw::value)
    .def_readwrite("sigma", &Geometry::Vdw::sigma)
    .def_readwrite("sym_idx", &Geometry::Vdw::sym_idx)
    .def_readwrite("pbc_shift", &Geometry::Vdw::pbc_shift)
    .def_readwrite("atoms", &Geometry::Vdw::atoms)
    ;

  py::bind_vector<std::vector<Geometry::Reporting::bond_reporting_t>>(geom, "ReportingBonds");
  py::bind_vector<std::vector<Geometry::Reporting::angle_reporting_t>>(geom, "ReportingAngles");
  py::bind_vector<std::vector<Geometry::Reporting::torsion_reporting_t>>(geom, "ReportingTorsions");
  py::bind_vector<std::vector<Geometry::Reporting::chiral_reporting_t>>(geom, "ReportingChirals");
  py::bind_vector<std::vector<Geometry::Reporting::plane_reporting_t>>(geom, "ReportingPlanes");
  py::bind_vector<std::vector<Geometry::Reporting::stacking_reporting_t>>(geom, "ReportingStackings");
  py::bind_vector<std::vector<Geometry::Reporting::vdw_reporting_t>>(geom, "ReportingVdws");
  py::bind_vector<std::vector<Geometry::Bond>>(geom, "Bonds");
  py::bind_vector<std::vector<Geometry::Angle>>(geom, "Angles");
  py::bind_vector<std::vector<Geometry::Chirality>>(geom, "Chiralitys");
  py::bind_vector<std::vector<Geometry::Torsion>>(geom, "Torsions");
  py::bind_vector<std::vector<Geometry::Plane>>(geom, "Planes");
  py::bind_vector<std::vector<Geometry::Interval>>(geom, "Intervals");
  py::bind_vector<std::vector<Geometry::Stacking>>(geom, "Stackings");
  py::bind_vector<std::vector<Geometry::Harmonic>>(geom, "Harmonics");
  py::bind_vector<std::vector<Geometry::Special>>(geom, "Specials");
  py::bind_vector<std::vector<Geometry::Vdw>>(geom, "Vdws");
  py::bind_vector<std::vector<Geometry::Bond::Value>>(bond, "Values");
  py::bind_vector<std::vector<Geometry::Angle::Value>>(angle, "Values");
  py::bind_vector<std::vector<Geometry::Torsion::Value>>(torsion, "Values");

  geomtarget
    .def_readonly("target", &GeomTarget::target)
    .def_readonly("vn", &GeomTarget::vn)
    .def_readonly("am", &GeomTarget::am)
    .def("am_for_coo", [](GeomTarget &self) {return for_coo_matrix(self);}, py::return_value_policy::reference_internal)
  ;
  geom
    .def(py::init<Structure&>(), py::arg("st"))
    .def_readonly("bonds", &Geometry::bonds)
    .def_readonly("angles", &Geometry::angles)
    .def_readonly("chirs", &Geometry::chirs)
    .def_readonly("torsions", &Geometry::torsions)
    .def_readonly("planes", &Geometry::planes)
    .def_readonly("intervals", &Geometry::intervals)
    .def_readonly("stackings", &Geometry::stackings)
    .def_readonly("harmonics", &Geometry::harmonics)
    .def_readonly("specials", &Geometry::specials)
    .def_readonly("vdws", &Geometry::vdws)
    .def_readonly("target", &Geometry::target)
    .def_readonly("reporting", &Geometry::reporting)
    .def("load_topo", &Geometry::load_topo)
    .def("finalize_restraints", &Geometry::finalize_restraints)
    .def("setup_target", &Geometry::setup_target)
    .def("clear_target", &Geometry::clear_target)
    .def("setup_vdw", &Geometry::setup_vdw)
    .def("calc", &Geometry::calc)
    .def("calc_adp_restraint", &Geometry::calc_adp_restraint)
  ;

  py::class_<TableS3>(m, "TableS3")
    .def(py::init<double, double>(), py::arg("d_min"), py::arg("d_max"))
    .def_readonly("s3_values", &TableS3::s3_values)
    .def_readonly("y_values", &TableS3::y_values)
    .def("make_table",&TableS3::make_table)
    .def("get_value", &TableS3::get_value)
    ;
  add_ll<gemmi::IT92<double>>(m, "LLX");
  m.def("precondition_eigen_coo", &precondition_eigen_coo);
}
