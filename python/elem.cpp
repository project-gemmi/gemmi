// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include "gemmi/elem.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/it92.hpp"
#include "gemmi/c4322.hpp"
#include "gemmi/neutron92.hpp"
#include "gemmi/util.hpp"  // for cat
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace gemmi;

// Round coefficients.
// Decimal number from tables, such as 11.7695, after -> 32-bit float -> double
// becomes 11.769499778747559. Let's round it, it's for printing only anyway.
double roc(float x) {
  float ax = std::abs(x);
  double n = ax < 16.f ? 1e6 : ax < 128.f ? 1e5 : 1e4;
  return std::round(x * n) / n;
}

void add_elem(py::module& m) {
  // it92.hpp
  using IT92 = gemmi::IT92<float>;
  py::class_<IT92::Coef>(m, "IT92Coef")
    .def_property_readonly("a", [](IT92::Coef& c) -> std::array<double,4> {
        return {{ roc(c.a(0)), roc(c.a(1)), roc(c.a(2)), roc(c.a(3)) }};
    })
    .def_property_readonly("b", [](IT92::Coef& c) -> std::array<double,4> {
        return {{ roc(c.b(0)), roc(c.b(1)), roc(c.b(2)), roc(c.b(3)) }};
    })
    .def_property_readonly("c", [](IT92::Coef& c) { return roc(c.c()); })
    .def("get_coefs", [](const IT92::Coef &self) { return self.coefs; })
    .def("set_coefs", &IT92::Coef::set_coefs)
    .def("calculate_sf", py::vectorize(&IT92::Coef::calculate_sf), py::arg("stol2"))
    .def("calculate_density_iso",
         [](const IT92::Coef &self, py::array_t<float> r2, float B) {
             return py::vectorize([&self,B](float r2) {
                 return self.calculate_density_iso(r2, B);
             })(r2);
    }, py::arg("r2"), py::arg("B"))
    ;
  m.def("IT92_normalize", &IT92::normalize);
  // can't define property for py::module_, and we don't expose IT92 as class
  m.def("IT92_get_ignore_charge", []() { return IT92::ignore_charge; });
  m.def("IT92_set_ignore_charge", [](bool v) { IT92::ignore_charge = v; });

  // c4322.hpp
  using C4322 = gemmi::C4322<float>;
  py::class_<C4322::Coef>(m, "C4322Coef")
    .def_property_readonly("a", [](C4322::Coef& c) -> std::array<double,5> {
        return {{ roc(c.a(0)), roc(c.a(1)), roc(c.a(2)), roc(c.a(3)), roc(c.a(4)) }};
    })
    .def_property_readonly("b", [](C4322::Coef& c) -> std::array<double,5> {
        return {{ roc(c.b(0)), roc(c.b(1)), roc(c.b(2)), roc(c.b(3)), roc(c.b(4)) }};
    })
    .def("get_coefs", [](const C4322::Coef &self) { return self.coefs; })
    .def("set_coefs", &C4322::Coef::set_coefs)
    .def("calculate_sf", &C4322::Coef::calculate_sf, py::arg("stol2"))
    .def("calculate_density_iso", &C4322::Coef::calculate_density_iso,
         py::arg("r2"), py::arg("B"))
    ;

  // neutron92.hpp
  using Neutron92 = gemmi::Neutron92<double>;
  py::class_<Neutron92::Coef>(m, "Neutron92")
    .def("get_coefs", [](const Neutron92::Coef &self) { return self.coefs; })
    // the same arguments as above - for consistency
    .def("calculate_sf", &Neutron92::Coef::calculate_sf, py::arg("stol2"))
    .def("calculate_density_iso", &Neutron92::Coef::calculate_density_iso,
         py::arg("r2"), py::arg("B"))
    ;

  // elem.hpp (w/ properties from it92.hpp, ...)
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def("__eq__",
         [](const Element &a, const Element &b) { return a.elem == b.elem; },
         py::is_operator())
#if PY_MAJOR_VERSION < 3  // in Py3 != is inferred from ==
    .def("__ne__",
         [](const Element &a, const Element &b) { return a.elem != b.elem; },
         py::is_operator())
#endif
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("covalent_r", &Element::covalent_r)
    .def_property_readonly("vdw_r", &Element::vdw_r)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def_property_readonly("is_hydrogen", &Element::is_hydrogen)
    .def_property_readonly("is_metal", &Element::is_metal)
    .def_property_readonly("it92", [](const Element& self) {
        return IT92::has(self.elem) ? &IT92::get(self.elem, 0) : nullptr;
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("c4322", [](const Element& self) {
        return C4322::has(self.elem) ? &C4322::get(self.elem) : nullptr;
    }, py::return_value_policy::reference_internal)
    .def_property_readonly("neutron92", [](const Element& self) {
        return Neutron92::get(self.elem);  // a copy is created
    })
    .def("__hash__", [](const Element &self) { return self.ordinal(); })
    .def("__repr__", [](const Element& self) {
        return gemmi::cat("gemmi.Element('", self.name(), "')");
    });

  m.def("IT92_get_exact", [](gemmi::Element el, signed char charge) {
      return IT92::get_exact(el.elem, charge);
  }, py::return_value_policy::reference_internal);

  // resinfo.hpp
  py::enum_<ResidueKind>(m, "ResidueKind")
    .value("UNKNOWN", ResidueKind::UNKNOWN)
    .value("AA", ResidueKind::AA)
    .value("AAD", ResidueKind::AAD)
    .value("PAA", ResidueKind::PAA)
    .value("MAA", ResidueKind::MAA)
    .value("RNA", ResidueKind::RNA)
    .value("DNA", ResidueKind::DNA)
    .value("BUF", ResidueKind::BUF)
    .value("HOH", ResidueKind::HOH)
    .value("PYR", ResidueKind::PYR)
    .value("KET", ResidueKind::KET)
    .value("ELS", ResidueKind::ELS);

  py::class_<ResidueInfo>(m, "ResidueInfo")
    .def_readonly("kind", &ResidueInfo::kind)
    .def_readonly("one_letter_code", &ResidueInfo::one_letter_code)
    .def_readonly("hydrogen_count", &ResidueInfo::hydrogen_count)
    .def_readonly("weight", &ResidueInfo::weight)
    .def("found", &ResidueInfo::found)
    .def("is_standard", &ResidueInfo::is_standard)
    .def("fasta_code", &ResidueInfo::fasta_code)
    .def("is_water", &ResidueInfo::is_water)
    .def("is_nucleic_acid", &ResidueInfo::is_nucleic_acid)
    .def("is_amino_acid", &ResidueInfo::is_amino_acid);

  m.def("find_tabulated_residue", &find_tabulated_residue, py::arg("name"),
        "Find chemical component information in the internal table.");
  m.def("expand_one_letter", &expand_one_letter);
  m.def("expand_one_letter_sequence", &expand_one_letter_sequence);
}
