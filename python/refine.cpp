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
using namespace pybind11::literals; // to bring in the `_a` literal
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
    .def("get_summary_table", [](const Geometry::Reporting& self, bool use_nucleus) {
      std::vector<std::string> keys;
      std::vector<int> nrest;
      std::vector<double> rmsd, rmsz;
      auto append = [&](const std::string& k, const std::vector<double>& delsq,
                        const std::vector<double>& zsq) {
        keys.emplace_back(k);
        nrest.push_back(delsq.size());
        rmsd.push_back(std::sqrt(std::accumulate(delsq.begin(), delsq.end(), 0.) / nrest.back()));
        rmsz.push_back(std::sqrt(std::accumulate(zsq.begin(), zsq.end(), 0.) / nrest.back()));
      };
      // Bond
      std::map<int, std::vector<double>> delsq, zsq;
      for (const auto& b : self.bonds) {
        const auto& restr = std::get<0>(b);
        const auto& val = std::get<1>(b);
        const double sigma = use_nucleus ? val->sigma_nucleus : val->sigma;
        const double d2 = sq(std::get<2>(b)), z2 = sq(std::get<2>(b) / sigma);
        const int k = (restr->type == 2 ? 2 :
                       (restr->atoms[0]->is_hydrogen() || restr->atoms[1]->is_hydrogen()) ? 1 : 0);
        delsq[k].push_back(d2);
        zsq[k].push_back(z2);
      }
      for (const auto& p : delsq)
        if (!p.second.empty())
          append(p.first == 2 ? "External distances" :
                 p.first == 1 ? "Bond distances, H" :
                 "Bond distances, non H", p.second, zsq[p.first]);

      // Angle
      delsq.clear(); zsq.clear();
      for (const auto& a : self.angles) {
        const auto& restr = std::get<0>(a);
        const auto& val = std::get<1>(a);
        const double d2 = sq(std::get<2>(a)), z2 = sq(std::get<2>(a) / val->sigma);
        const int k = (restr->atoms[0]->is_hydrogen() || restr->atoms[1]->is_hydrogen()) ? 1 : 0;
        delsq[k].push_back(d2);
        zsq[k].push_back(z2);
      }
      for (const auto& p : delsq)
        if (!p.second.empty())
          append(p.first == 1 ? "Bond angles, H" : "Bond angles, non H", p.second, zsq[p.first]);

      // Torsion
      delsq.clear(); zsq.clear();
      for (const auto& t : self.torsions) {
        const auto& val = std::get<1>(t);
        const double d2 = sq(std::get<2>(t)), z2 = sq(std::get<2>(t) / val->sigma);
        delsq[val->period].push_back(d2);
        zsq[val->period].push_back(z2);
      }
      for (const auto& p : delsq)
        if (!p.second.empty())
          append("Torsion angles, period " + std::to_string(p.first), p.second, zsq[p.first]);

      // Chiral
      delsq.clear(); zsq.clear();
      for (const auto& c : self.chirs) {
        const auto& val = std::get<0>(c);
        const double d2 = sq(std::get<1>(c)), z2 = sq(std::get<1>(c) / val->sigma);
        delsq[0].push_back(d2);
        zsq[0].push_back(z2);
      }
      if (!delsq[0].empty())
        append("Chiral centres", delsq[0], zsq[0]);

      // Plane
      delsq.clear(); zsq.clear();
      for (const auto& p : self.planes) {
        const auto& val = std::get<0>(p);
        for (double d : std::get<1>(p)) {
          const double d2 = sq(d), z2 = sq(d / val->sigma);
          delsq[0].push_back(d2);
          zsq[0].push_back(z2);
        }
      }
      if (!delsq[0].empty())
        append("Planar groups", delsq[0], zsq[0]);

      // Stack
      delsq.clear(); zsq.clear();
      for (const auto& s : self.stackings) {
        const auto& restr = std::get<0>(s);
        const double da2 = sq(std::get<1>(s)), za2 = sq(std::get<1>(s) / restr->sd_angle);
        const double dd2 = 0.5 * (sq(std::get<2>(s)) + sq(std::get<3>(s)));
        const double zd2 = dd2 / sq(restr->sd_dist);
        delsq[0].push_back(da2);
        zsq[0].push_back(za2);
        delsq[1].push_back(dd2);
        zsq[1].push_back(zd2);
      }
      if (!delsq[0].empty())
        append("Stacking angles", delsq[0], zsq[0]);
      if (!delsq[1].empty())
        append("Stacking distances", delsq[1], zsq[1]);

      // VDW
      delsq.clear(); zsq.clear();
      for (const auto& v : self.vdws) {
        const auto& restr = std::get<0>(v);
        const double d2 = sq(std::get<1>(v)), z2 = sq(std::get<1>(v) / restr->sigma);
        delsq[restr->type].push_back(d2);
        zsq[restr->type].push_back(z2);
      }
      for (const auto& p : delsq)
        if (!p.second.empty()) {
          const int i = p.first > 6 ? p.first - 6 : p.first;
          append((i == 1 ? "VDW nonbonded" :
                  i == 2 ? "VDW torsion" :
                  i == 3 ? "VDW hbond" :
                  i == 4 ? "VDW metal" :
                  i == 5 ? "VDW dummy" :
                  "VDW dummy-dummy") + std::string(p.first > 6 ? ", symmetry" : ""),
                 p.second, zsq[p.first]);
        }
      return py::dict("Restraint type"_a=keys, "N restraints"_a=nrest,
                      "r.m.s.d."_a=rmsd, "r.m.s.Z"_a=rmsz);
    })
    .def("get_bond_outliers", [](const Geometry::Reporting& self, bool use_nucleus, double min_z) {
      std::vector<const Atom*> atom1, atom2;
      std::vector<double> values, ideals, zs, alphas;
      std::vector<int> types;
      for (const auto& b : self.bonds) {
        const auto& restr = std::get<0>(b);
        const auto& val = std::get<1>(b);
        const double ideal = use_nucleus ? val->value_nucleus : val->value;
        const double sigma = use_nucleus ? val->sigma_nucleus : val->sigma;
        const double z = std::get<2>(b) / sigma; // value - ideal
        if (std::abs(z) >= min_z) {
          atom1.push_back(restr->atoms[0]);
          atom2.push_back(restr->atoms[1]);
          values.push_back(std::get<2>(b) + ideal);
          ideals.push_back(ideal);
          zs.push_back(z);
          types.push_back(restr->type);
          alphas.push_back(restr->alpha);
        }
      }
      return py::dict("atom1"_a=atom1, "atom2"_a=atom2, "value"_a=values,
                      "ideal"_a=ideals, "z"_a=zs, "type"_a=types, "alpha"_a=alphas);
    }, py::arg("use_nucleus"), py::arg("min_z"))
    .def("get_angle_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2, atom3;
      std::vector<double> values, ideals, zs;
      for (const auto& t : self.angles) {
        const auto& restr = std::get<0>(t);
        const auto& val = std::get<1>(t);
        const double z = std::get<2>(t) / val->sigma; // value - ideal
        if (std::abs(z) >= min_z) {
          atom1.push_back(restr->atoms[0]);
          atom2.push_back(restr->atoms[1]);
          atom3.push_back(restr->atoms[2]);
          values.push_back(std::get<2>(t) + val->value);
          ideals.push_back(val->value);
          zs.push_back(z);
        }
      }
      return py::dict("atom1"_a=atom1, "atom2"_a=atom2, "atom3"_a=atom3,
                      "value"_a=values, "ideal"_a=ideals, "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_torsion_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2, atom3, atom4;
      std::vector<double> values, ideals, zs;
      std::vector<int> pers;
      std::vector<std::string> labels;
      for (const auto& t : self.torsions) {
        const auto& restr = std::get<0>(t);
        const auto& val = std::get<1>(t);
        const double z = std::get<2>(t) / val->sigma; // value - ideal
        if (std::abs(z) >= min_z) {
          atom1.push_back(restr->atoms[0]);
          atom2.push_back(restr->atoms[1]);
          atom3.push_back(restr->atoms[2]);
          atom4.push_back(restr->atoms[3]);
          labels.push_back(val->label);
          values.push_back(std::get<2>(t) + val->value);
          ideals.push_back(val->value);
          pers.push_back(val->period);
          zs.push_back(z);
        }
      }
      return py::dict("label"_a=labels, "atom1"_a=atom1, "atom2"_a=atom2, "atom3"_a=atom3, "atom4"_a=atom4,
                      "value"_a=values, "ideal"_a=ideals, "per"_a=pers, "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_chiral_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2, atom3, atom4;
      std::vector<double> values, ideals, zs;
      std::vector<bool> signs;
      for (const auto& t : self.chirs) {
        const auto& restr = std::get<0>(t);
        const double z = std::get<1>(t) / restr->sigma; // value - ideal
        if (std::abs(z) >= min_z) {
          atom1.push_back(restr->atoms[0]);
          atom2.push_back(restr->atoms[1]);
          atom3.push_back(restr->atoms[2]);
          atom4.push_back(restr->atoms[3]);
          values.push_back(std::get<1>(t) + std::get<2>(t));
          ideals.push_back(std::get<2>(t));
          signs.push_back(restr->sign == ChiralityType::Both);
          zs.push_back(z);
        }
      }
      return py::dict("atomc"_a=atom1, "atom1"_a=atom2, "atom2"_a=atom3, "atom3"_a=atom4,
                      "value"_a=values, "ideal"_a=ideals, "both"_a=signs, "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_plane_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atoms;
      std::vector<double> values, zs;
      std::vector<std::string> labels;
      for (const auto& t : self.planes) {
        const auto& restr = std::get<0>(t);
        for (size_t i = 0; i < restr->atoms.size(); ++i) {
          const double z = std::get<1>(t)[i] / restr->sigma;
          if (std::abs(z) >= min_z) {
            atoms.push_back(restr->atoms[i]);
            labels.push_back(restr->label);
            values.push_back(std::get<1>(t)[i]);
            zs.push_back(z);
          }
        }
      }
      return py::dict("label"_a=labels, "atom"_a=atoms, "dev"_a=values, "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_stacking_angle_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2;
      std::vector<double> values, ideals, zs;
      for (const auto& t : self.stackings) {
        const auto& restr = std::get<0>(t);
        const double za = std::get<1>(t) / restr->sd_angle;
        if (std::abs(za) >= min_z) {
          atom1.push_back(restr->planes[0][0]); // report only first atom
          atom2.push_back(restr->planes[1][0]);
          values.push_back(std::get<1>(t) + restr->angle);
          ideals.push_back(restr->angle);
          zs.push_back(za);
        }
      }
      return py::dict("plane1"_a=atom1, "plane2"_a=atom2, "value"_a=values, "ideal"_a=ideals, "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_stacking_dist_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2;
      std::vector<double> values, ideals, zs;
      for (const auto& t : self.stackings) {
        const auto& restr = std::get<0>(t);
        const double zd1 = std::get<2>(t) / restr->sd_dist;
        const double zd2 = std::get<3>(t) / restr->sd_dist;
        if (std::min(std::abs(zd1), std::abs(zd2)) >= min_z) {
          const double zd = std::abs(zd1) > std::abs(zd2) ? zd1 : zd2;
          atom1.push_back(restr->planes[0][0]); // report only first atom
          atom2.push_back(restr->planes[1][0]);
          values.push_back(zd * restr->sd_dist + restr->dist);
          ideals.push_back(restr->dist);
          zs.push_back(zd);
        }
      }
      return py::dict("plane1"_a=atom1, "plane2"_a=atom2, "value"_a=values, "ideal"_a=ideals,  "z"_a=zs);
    }, py::arg("min_z"))
    .def("get_vdw_outliers", [](const Geometry::Reporting& self, double min_z) {
      std::vector<const Atom*> atom1, atom2;
      std::vector<double> values, ideals, zs;
      std::vector<int> types;
      for (const auto& t : self.vdws) {
        const auto& restr = std::get<0>(t);
        const double z = std::get<1>(t) / restr->sigma;
        if (std::abs(z) >= min_z) {
          atom1.push_back(restr->atoms[0]);
          atom2.push_back(restr->atoms[1]);
          values.push_back(std::get<1>(t) + restr->value);
          ideals.push_back(restr->value);
          zs.push_back(z);
          types.push_back(restr->type);
        }
      }
      return py::dict("atom1"_a=atom1, "atom2"_a=atom2, "value"_a=values,
                      "ideal"_a=ideals, "z"_a=zs, "type"_a=types);
    }, py::arg("min_z"))
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
    .def_readwrite("label", &Geometry::Torsion::Value::label)
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
    .def_readwrite("label", &Geometry::Plane::label)
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
    .def("calc", &Geometry::calc, py::arg("use_nucleus"), py::arg("check_only"),
         py::arg("wbond")=1, py::arg("wangle")=1, py::arg("wtors")=1,
         py::arg("wchir")=1, py::arg("wplane")=1, py::arg("wstack")=1, py::arg("wvdw")=1)
    .def("calc_adp_restraint", &Geometry::calc_adp_restraint)
    // vdw parameters
    .def_readwrite("vdw_sdi_vdw", &Geometry::vdw_sdi_vdw)
    .def_readwrite("vdw_sdi_torsion", &Geometry::vdw_sdi_torsion)
    .def_readwrite("vdw_sdi_hbond", &Geometry::vdw_sdi_hbond)
    .def_readwrite("vdw_sdi_metal", &Geometry::vdw_sdi_metal)
    .def_readwrite("hbond_dinc_ad", &Geometry::hbond_dinc_ad)
    .def_readwrite("hbond_dinc_ah", &Geometry::hbond_dinc_ah)
    .def_readwrite("dinc_torsion_o", &Geometry::dinc_torsion_o)
    .def_readwrite("dinc_torsion_n", &Geometry::dinc_torsion_n)
    .def_readwrite("dinc_torsion_c", &Geometry::dinc_torsion_c)
    .def_readwrite("dinc_torsion_all", &Geometry::dinc_torsion_all)
    .def_readwrite("dinc_dummy", &Geometry::dinc_dummy)
    .def_readwrite("vdw_sdi_dummy", &Geometry::vdw_sdi_dummy)
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
