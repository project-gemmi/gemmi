// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"

#include "common.h"
#include "array.h"  // for miller_function
#include <nanobind/make_iterator.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

using namespace gemmi;

void add_symmetry(nb::module_& m) {
  nb::enum_<CrystalSystem>(m, "CrystalSystem")
    .value("Triclinic", CrystalSystem::Triclinic)
    .value("Monoclinic", CrystalSystem::Monoclinic)
    .value("Orthorhombic", CrystalSystem::Orthorhombic)
    .value("Tetragonal", CrystalSystem::Tetragonal)
    .value("Trigonal", CrystalSystem::Trigonal)
    .value("Hexagonal", CrystalSystem::Hexagonal)
    .value("Cubic", CrystalSystem::Cubic)
    ;
  nb::class_<Op>(m, "Op")
    .def("__init__", [](Op* op) {
        new(op) Op(Op::identity());
    })
    .def("__init__", [](Op* op, const std::string& s) {
        new(op) Op(parse_triplet(s));
    })
    .def_prop_ro_static("DEN", [](const nb::object&) { return Op::DEN; },
                           "Denominator (integer) for the translation vector.")
    .def_rw("rot", &Op::rot, "3x3 integer matrix.")
    .def_rw("tran", &Op::tran,
       "Numerators (integers) of the translation vector. Denominator DEN=24.")
    .def("triplet", &Op::triplet, nb::arg("style")='x')
    .def("inverse", &Op::inverse, "Returns inverted operator.")
    .def("wrap", &Op::wrap, "Wrap the translation part to [0,1)")
    .def("translated", &Op::translated, nb::arg("a"), "Adds a to tran")
    .def("transposed_rot", &Op::transposed_rot)
    .def("det_rot", &Op::det_rot, "Determinant of the 3x3 matrix.")
    .def("rot_type", &Op::rot_type)
    .def("combine", &Op::combine, nb::arg("b"),
         "Combine two symmetry operations.")
    .def("seitz", [](const Op& self) {
            auto arr = self.int_seitz();
            auto mat = nb::list();
            nb::object fr = nb::module_::import_("fractions").attr("Fraction");
            for (int i = 0; i < 4; ++i) {
              auto row = nb::list();
              for (int j = 0; j < 4; ++j) {
                auto v = arr[i][j];
                if (i == 3 || v == 0)
                  row.append(v);
                else if (std::abs(v) == Op::DEN)
                  row.append(v / Op::DEN);
                else
                  row.append(fr(v, Op::DEN + 0));  // +0 to avoid linker error
              }
              mat.append(row);
            }
            return mat;
         }, "Returns Seitz matrix (fractions)")
    .def("float_seitz", &Op::float_seitz, "Returns Seitz matrix (floats)")
    .def("apply_to_xyz", &Op::apply_to_xyz, nb::arg("xyz"))
    .def("apply_to_hkl", &Op::apply_to_hkl, nb::arg("hkl"))
    .def("phase_shift", &Op::phase_shift, nb::arg("hkl"))
    .def(nb::self * nb::self)
    .def("__mul__", [](const Op &a, const std::string& b) { return a * parse_triplet(b); },
         nb::is_operator())
    .def("__rmul__", [](const Op& a, const std::string& b) { return parse_triplet(b) * a; },
         nb::is_operator())
    .def(nb::self == nb::self)  // NOLINT(misc-redundant-expression)
    .def("__eq__", [](const Op& a, const std::string& b) { return a == parse_triplet(b); },
         nb::is_operator(), nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def("__copy__", [](const Op& self) { return Op(self); })
    .def("__deepcopy__", [](const Op& self, const nb::dict&) { return Op(self); }, nb::arg("memo"))
    .def("__hash__", [](const Op& self) { return std::hash<Op>()(self); })
    .def("__repr__", [](const Op& self) {
        return "<gemmi.Op(\"" + self.triplet() + "\")>";
    });

  m.def("parse_triplet", &parse_triplet, nb::arg("triplet"),
        "Parse coordinate triplet into gemmi.Op.");
  m.def("parse_triplet_part", [](const std::string& s) { return parse_triplet_part(s); },
        "Parse one of the three parts of a triplet.");
  m.def("make_triplet_part", &make_triplet_part,
        nb::arg("xyz"), nb::arg("w"), nb::arg("style")='x',
        "Make one of the three parts of a triplet.");

  nb::class_<GroupOps>(m, "GroupOps")
    .def("__init__", [](GroupOps* p, const std::vector<Op>& ops) {
        new(p) GroupOps(split_centering_vectors(ops));
    })
    .def("__iter__", [](const GroupOps& self) {
        return nb::make_iterator(nb::type<GroupOps>(), "iterator", self);
    }, nb::keep_alive<0, 1>())
    .def("__eq__", [](const GroupOps &a, const GroupOps &b) { return a.is_same_as(b); },
         nb::is_operator(), nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def("__len__", [](const GroupOps& g) { return g.order(); })
    .def("__deepcopy__", [](const GroupOps& g, const nb::dict&) { return GroupOps(g); },
         nb::arg("memo"))
    .def_rw("sym_ops", &GroupOps::sym_ops,
               "Symmetry operations (to be combined with centering vectors).")
    .def_rw("cen_ops", &GroupOps::cen_ops, "Centering vectors.")
    .def("add_missing_elements", &GroupOps::add_missing_elements)
    .def("find_centering", &GroupOps::find_centering)
    .def("has_same_centring", &GroupOps::has_same_centring)
    .def("has_same_rotations", &GroupOps::has_same_rotations)
    .def("is_centrosymmetric", &GroupOps::is_centrosymmetric)
    .def("is_reflection_centric", &GroupOps::is_reflection_centric)
    .def("centric_flag_array", [](const GroupOps& g, const cpu_miller_array& hkl) {
        return miller_function<bool>(g, &GroupOps::is_reflection_centric, hkl);
    })
    .def("epsilon_factor_without_centering", &GroupOps::epsilon_factor_without_centering)
    .def("epsilon_factor", &GroupOps::epsilon_factor)
    .def("epsilon_factor_array", [](const GroupOps& g, const cpu_miller_array& hkl) {
        return miller_function<int>(g, &GroupOps::epsilon_factor, hkl);
    })
    .def("epsilon_factor_without_centering_array", [](const GroupOps& g,
                                                      const cpu_miller_array& hkl) {
        return miller_function<int>(g, &GroupOps::epsilon_factor_without_centering, hkl);
    })
    .def("is_systematically_absent", &GroupOps::is_systematically_absent)
    .def("systematic_absences", [](const GroupOps& g, const cpu_miller_array& hkl) {
        return miller_function<bool>(g, &GroupOps::is_systematically_absent, hkl);
    })
    .def("find_grid_factors", &GroupOps::find_grid_factors,
         "Minimal multiplicity for real-space grid (e.g. 1,1,6 for P61).")
    .def("change_basis_forward", &GroupOps::change_basis_forward, nb::arg("cob"),
         "Applies the change-of-basis operator (in place).")
    .def("change_basis_backward", &GroupOps::change_basis_backward, nb::arg("cob"),
         "Applies inverse of the change-of-basis operator (in place).")
    .def("derive_symmorphic", &GroupOps::derive_symmorphic)
    .def("add_inversion", &GroupOps::add_inversion)
    ;

  nb::class_<SpaceGroup>(m, "SpaceGroup")
    .def(nb::new_([](int n) {
          int main_table_length = int(sizeof(spacegroup_tables::main) / sizeof(SpaceGroup));
          // a second meaning of SpaceGroup(n), for internal use in __reduce__
          if (n < INT_MIN + main_table_length)
            return &spacegroup_tables::main[n - INT_MIN];
          return &get_spacegroup_by_number(n);
    }), nb::arg("ccp4"), nb::rv_policy::reference)
    .def(nb::new_([](const std::string& s) {
          return &get_spacegroup_by_name(s);
    }), nb::arg("hm"), nb::rv_policy::reference)
    .def_ro("number", &SpaceGroup::number, "number 1-230.")
    .def_ro("ccp4", &SpaceGroup::ccp4, "ccp4 number")
    // Intel Compiler would not compile .def_ro("hm", ...),
    // the difference with hm, qualifier and hall is that they are char[N]
    .def_prop_ro("hm", [](const SpaceGroup &self) -> const char* {
        return self.hm;
    }, "Hermann-Mauguin name")
    .def_ro("ext", &SpaceGroup::ext, "Extension (1, 2, H, R or none)")
    .def_prop_ro("qualifier", [](const SpaceGroup &s) -> const char* {
        return s.qualifier;
    }, "e.g. 'cab'")
    .def_prop_ro("hall", [](const SpaceGroup &self) -> const char* {
        return self.hall;
    }, "Hall symbol")
    .def_prop_ro("basisop", &SpaceGroup::basisop)
    .def("xhm", &SpaceGroup::xhm, "extended Hermann-Mauguin name")
    .def("centring_type", &SpaceGroup::centring_type)
    .def("short_name", &SpaceGroup::short_name,
         "H-M name w/o spaces and with 1's removed in '1 ... 1'.")
    .def("is_enantiomorphic", &SpaceGroup::is_enantiomorphic)
    .def("is_sohncke", &SpaceGroup::is_sohncke)
    .def("is_symmorphic", &SpaceGroup::is_symmorphic)
    .def("point_group_hm", &SpaceGroup::point_group_hm,
         "Returns H-M name of the point group.")
    .def("laue_str", &SpaceGroup::laue_str,
         "Returns name of the Laue class (for centrosymmetric groups "
         "the same as point_group_hm).")
    .def("crystal_system", &SpaceGroup::crystal_system)
    .def("crystal_system_str", &SpaceGroup::crystal_system_str,
         "Returns lower-case name of the crystal system.")
    .def("is_centrosymmetric", &SpaceGroup::is_centrosymmetric)
    .def("monoclinic_unique_axis", &SpaceGroup::monoclinic_unique_axis)
    .def("is_reference_setting", &SpaceGroup::is_reference_setting)
    .def("centred_to_primitive", &SpaceGroup::centred_to_primitive)
    .def("change_of_hand_op", &SpaceGroup::change_of_hand_op)
    .def("operations", &SpaceGroup::operations, "Group of operations")
    .def("switch_to_asu", [](const SpaceGroup& sg,
                             const nb::ndarray<int, nb::shape<-1,-1>, nb::device::cpu>& hkl) {
        auto h = hkl.view();
        if (h.shape(1) < 3)
          throw std::domain_error("error: the size of the second dimension < 3");
        GroupOps gops = sg.operations();
        ReciprocalAsu asu(&sg);
        for (size_t i = 0; i < h.shape(0); ++i) {
          Op::Miller hkl = asu.to_asu({{h(i, 0), h(i, 1), h(i, 2)}}, gops).first;
          for (size_t j = 0; j != 3; ++j)
            h(i, j) = hkl[j];
        }
    }, nb::arg("miller_array").noconvert())
    // Check equality by comparing Hall symbol strings.
    // In Python, SpaceGroup always points to an entry in spacegroup_tables::main,
    // so Hall symbols are consistent. The same settings with different H-M names
    // (see space group 65 in the ITB list) are considered equal.
    .def("__eq__", [](const SpaceGroup& a, const SpaceGroup& b) {
        return strcmp(a.hall, b.hall) == 0;
    }, nb::is_operator(), nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def("__reduce__", [](const SpaceGroup& self) {
        // faster than serializing self.xhm()
        std::ptrdiff_t pos = &self - spacegroup_tables::main;
        int main_table_length = int(sizeof(spacegroup_tables::main) / sizeof(SpaceGroup));
        (void) main_table_length;
        assert(pos >= 0 && pos < main_table_length);
        return nb::make_tuple(nb::type<SpaceGroup>(), nb::make_tuple(INT_MIN + (int)pos));
    })
    .def("__repr__", [](const SpaceGroup &self) {
        return "<gemmi.SpaceGroup(\"" + self.xhm() + "\")>";
    });

  nb::class_<ReciprocalAsu>(m, "ReciprocalAsu")
    .def(nb::init<const SpaceGroup*, bool>(), nb::arg(), nb::arg("tnt")=false)
    .def("is_in", &ReciprocalAsu::is_in, nb::arg("hkl"))
    .def("condition_str", &ReciprocalAsu::condition_str)
    .def("to_asu",
         nb::overload_cast<const Op::Miller&, const GroupOps&>(&ReciprocalAsu::to_asu, nb::const_),
         nb::arg("hkl"), nb::arg("group_ops"))
    ;

  nb::handle mod = m;
  m.def("spacegroup_table", [mod]() {
      return nb::make_iterator<nb::rv_policy::reference>(mod, "spacegroup_iterator",
                                                         spacegroup_tables::main);
  });
  m.def("spacegroup_table_itb", [mod]() {
      return nb::make_iterator<nb::rv_policy::reference>(mod, "itb_spacegroup_iterator",
                                                         spacegroup_tables::main + 0,
                                                         spacegroup_tables::main + 530);
  });
  m.def("generators_from_hall", &generators_from_hall, nb::arg("hall"),
        "Parse Hall notation.");
  m.def("symops_from_hall", &symops_from_hall, nb::arg("hall"),
        "Parse Hall notation.");
  m.def("find_spacegroup_by_number", &find_spacegroup_by_number,
        nb::arg("ccp4"), nb::rv_policy::reference,
        "Returns space-group of given number.");
  m.def("find_spacegroup_by_name", &find_spacegroup_by_name,
        nb::arg("hm"), nb::arg("alpha")=0., nb::arg("gamma")=0., nb::arg("prefer")="",
        nb::rv_policy::reference,
        "Returns space-group with given name.");
  m.def("get_spacegroup_reference_setting", &get_spacegroup_reference_setting,
        nb::arg("number"), nb::rv_policy::reference);
  m.def("find_spacegroup_by_ops", &find_spacegroup_by_ops,
        nb::arg("group_ops"), nb::rv_policy::reference,
        "Returns space-group with identical operations.");
  m.def("seitz_to_op", &seitz_to_op);
}
