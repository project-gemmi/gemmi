// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for expand_one_letter_sequence

#include "gemmi/elem.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/it92.hpp"
#include "gemmi/c4322.hpp"
#include "gemmi/neutron92.hpp"
#include "gemmi/util.hpp"  // for cat

#include "gemmi/xds_ascii.hpp"  // for XdsAscii
#include "gemmi/xds2mtz.hpp"  // for xds_to_mtz
#include "array.h"  // for vector_member_array

using namespace gemmi;

// Round coefficients.
// Decimal number from tables, such as 11.7695, after -> 32-bit float -> double
// becomes 11.769499778747559. Let's round it, it's for printing only anyway.
static double roc(float x) {
  float ax = std::abs(x);
  double n = ax < 16.f ? 1e6 : ax < 128.f ? 1e5 : 1e4;
  return std::round(x * n) / n;
}

void add_elem(nb::module_& m) {
  // it92.hpp
  using IT92 = gemmi::IT92<float>;
  nb::class_<IT92::Coef>(m, "IT92Coef")
    .def_prop_ro("a", [](IT92::Coef& c) -> std::array<double,4> {
        return {{ roc(c.a(0)), roc(c.a(1)), roc(c.a(2)), roc(c.a(3)) }};
    })
    .def_prop_ro("b", [](IT92::Coef& c) -> std::array<double,4> {
        return {{ roc(c.b(0)), roc(c.b(1)), roc(c.b(2)), roc(c.b(3)) }};
    })
    .def_prop_ro("c", [](IT92::Coef& c) { return roc(c.c()); })
    .def("get_coefs", [](const IT92::Coef &self) { return self.coefs; })
    .def("set_coefs", &IT92::Coef::set_coefs)
    .def("calculate_sf", &IT92::Coef::calculate_sf, nb::arg("stol2"))
    .def("calculate_density_iso", &IT92::Coef::calculate_density_iso,
         nb::arg("r2"), nb::arg("B"))
    ;
  m.def("IT92_normalize", &IT92::normalize);
  // can't define property for nb::module_, and we don't expose IT92 as class
  m.def("IT92_get_ignore_charge", []() { return IT92::ignore_charge; });
  m.def("IT92_set_ignore_charge", [](bool v) { IT92::ignore_charge = v; });

  // c4322.hpp
  using C4322 = gemmi::C4322<float>;
  nb::class_<C4322::Coef>(m, "C4322Coef")
    .def_prop_ro("a", [](C4322::Coef& c) -> std::array<double,5> {
        return {{ roc(c.a(0)), roc(c.a(1)), roc(c.a(2)), roc(c.a(3)), roc(c.a(4)) }};
    })
    .def_prop_ro("b", [](C4322::Coef& c) -> std::array<double,5> {
        return {{ roc(c.b(0)), roc(c.b(1)), roc(c.b(2)), roc(c.b(3)), roc(c.b(4)) }};
    })
    .def("get_coefs", [](const C4322::Coef &self) { return self.coefs; })
    .def("set_coefs", &C4322::Coef::set_coefs)
    .def("calculate_sf", &C4322::Coef::calculate_sf, nb::arg("stol2"))
    .def("calculate_density_iso", &C4322::Coef::calculate_density_iso,
         nb::arg("r2"), nb::arg("B"))
    ;

  // neutron92.hpp
  using Neutron92 = gemmi::Neutron92<double>;
  nb::class_<Neutron92::Coef>(m, "Neutron92")
    .def("get_coefs", [](const Neutron92::Coef &self) { return self.coefs; })
    // the same arguments as above - for consistency
    .def("calculate_sf", &Neutron92::Coef::calculate_sf, nb::arg("stol2"))
    .def("calculate_density_iso", &Neutron92::Coef::calculate_density_iso,
         nb::arg("r2"), nb::arg("B"))
    ;

  // elem.hpp (w/ properties from it92.hpp, ...)
  nb::class_<Element>(m, "Element")
    .def(nb::init<const std::string &>())
    .def(nb::init<int>())
    // NOLINTNEXTLINE(misc-redundant-expression)
    .def(nb::self == nb::self, nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def_prop_ro("name", &Element::name)
    .def_prop_ro("weight", &Element::weight)
    .def_prop_ro("covalent_r", &Element::covalent_r)
    .def_prop_ro("vdw_r", &Element::vdw_r)
    .def_prop_ro("atomic_number", &Element::atomic_number)
    .def_prop_ro("is_hydrogen", &Element::is_hydrogen)
    .def_prop_ro("is_metal", &Element::is_metal)
    .def_prop_ro("it92", [](const Element& self) {
        return IT92::has(self.elem) ? &IT92::get(self.elem, 0) : nullptr;
    }, nb::rv_policy::reference)
    .def_prop_ro("c4322", [](const Element& self) {
        return C4322::has(self.elem) ? &C4322::get(self.elem) : nullptr;
    }, nb::rv_policy::reference)
    .def_prop_ro("neutron92", [](const Element& self) {
        return Neutron92::get(self.elem);  // a copy is created
    })
    .def("__hash__", [](const Element &self) { return self.ordinal(); })
    .def("__repr__", [](const Element& self) {
        return gemmi::cat("gemmi.Element('", self.name(), "')");
    });

  m.def("IT92_get_exact", [](gemmi::Element el, signed char charge) {
      return IT92::get_exact(el.elem, charge);
  }, nb::rv_policy::reference);

  m.def("set_is_metal", [](const char* el, bool value) {
      set_is_metal(find_element(el), value);
  });

  // resinfo.hpp
  nb::enum_<ResidueKind>(m, "ResidueKind")
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

  nb::class_<ResidueInfo>(m, "ResidueInfo")
    .def_ro("kind", &ResidueInfo::kind)
    .def_ro("one_letter_code", &ResidueInfo::one_letter_code)
    .def_ro("hydrogen_count", &ResidueInfo::hydrogen_count)
    .def_ro("weight", &ResidueInfo::weight)
    .def("found", &ResidueInfo::found)
    .def("is_standard", &ResidueInfo::is_standard)
    .def("fasta_code", &ResidueInfo::fasta_code)
    .def("is_water", &ResidueInfo::is_water)
    .def("is_nucleic_acid", &ResidueInfo::is_nucleic_acid)
    .def("is_amino_acid", &ResidueInfo::is_amino_acid);

  m.def("find_tabulated_residue", &find_tabulated_residue, nb::arg("name"),
        "Find chemical component information in the internal table.");
  m.def("expand_one_letter", &expand_one_letter);
  m.def("expand_one_letter_sequence", &expand_one_letter_sequence);
}

void add_xds(nb::module_& m) {
  nb::class_<XdsAscii>(m, "XdsAscii")
    .def_ro("source_path", &XdsAscii::source_path)
    .def_ro("read_columns", &XdsAscii::read_columns)
    .def_ro("spacegroup_number", &XdsAscii::spacegroup_number)
    .def_ro("wavelength", &XdsAscii::wavelength)
    .def_ro("cell_constants", &XdsAscii::cell_constants)
    .def_ro("rotation_axis", &XdsAscii::rotation_axis)
    .def_ro("friedels_law", &XdsAscii::friedels_law)
    .def_ro("generated_by", &XdsAscii::generated_by)
    .def_prop_ro("data_size", [](XdsAscii& self) { return self.data.size(); })
    .def_prop_ro("miller_array", [](XdsAscii& self) {
        constexpr int64_t stride = int64_t(sizeof(XdsAscii::Refl) / sizeof(int));
        return nb::ndarray<nb::numpy, int, nb::shape<-1,3>>(
          &self.data[0].hkl[0], {self.data.size(), 3}, nb::handle(), {stride, 1});
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("iset_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::iset);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("iobs_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::iobs);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("sigma_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::sigma);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("xd_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::xd);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("yd_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::yd);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("zd_array", [](XdsAscii& self) {
        return vector_member_array(self.data, &XdsAscii::Refl::zd);
    }, nb::rv_policy::reference_internal)
    .def("filtered", [](XdsAscii& self, const cpu_array<bool>& selection) {
        auto v = selection.view();
        size_t n = self.data.size();
        if (v.shape(0) != n)
          throw nb::value_error("boolean array must match the number of reflections");
        XdsAscii ret((const XdsAsciiMetadata&) self);
        size_t count = 0;
        for (size_t i = 0; i < v.shape(0); ++i)
          count += (size_t) v(i);
        ret.data.reserve(count);
        for (size_t i = 0; i < v.shape(0); ++i)
          if (v(i))
            ret.data.push_back(self.data[i]);
        return ret;
    })
    .def("apply_polarization_correction", &XdsAscii::apply_polarization_correction,
         nb::arg("p"), nb::arg("normal"))
    .def("to_mtz", &gemmi::xds_to_mtz)
    ;
  m.def("read_xds_ascii", &read_xds_ascii);
}
