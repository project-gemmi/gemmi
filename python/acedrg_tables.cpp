// Copyright 2025 Global Phasing Ltd.

#include "gemmi/acedrg_tables.hpp"

#include "common.h"
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

using namespace gemmi;

void add_acedrg_tables(nb::module_& m) {
  nb::enum_<Hybridization>(m, "Hybridization")
    .value("SP1", Hybridization::SP1)
    .value("SP2", Hybridization::SP2)
    .value("SP3", Hybridization::SP3)
    .value("SPD5", Hybridization::SPD5)
    .value("SPD6", Hybridization::SPD6)
    .value("SPD7", Hybridization::SPD7)
    .value("SPD8", Hybridization::SPD8)
    .value("SP_NON", Hybridization::SP_NON);

  nb::enum_<CoordGeometry>(m, "CoordGeometry")
    .value("LINEAR", CoordGeometry::LINEAR)
    .value("TRIGONAL_PLANAR", CoordGeometry::TRIGONAL_PLANAR)
    .value("T_SHAPED", CoordGeometry::T_SHAPED)
    .value("TETRAHEDRAL", CoordGeometry::TETRAHEDRAL)
    .value("SQUARE_PLANAR", CoordGeometry::SQUARE_PLANAR)
    .value("TRIGONAL_BIPYRAMIDAL", CoordGeometry::TRIGONAL_BIPYRAMIDAL)
    .value("SQUARE_PYRAMIDAL", CoordGeometry::SQUARE_PYRAMIDAL)
    .value("OCTAHEDRAL", CoordGeometry::OCTAHEDRAL)
    .value("TRIGONAL_PRISM", CoordGeometry::TRIGONAL_PRISM)
    .value("PENTAGONAL_BIPYRAMIDAL", CoordGeometry::PENTAGONAL_BIPYRAMIDAL)
    .value("CAPPED_OCTAHEDRAL", CoordGeometry::CAPPED_OCTAHEDRAL)
    .value("SQUARE_ANTIPRISM", CoordGeometry::SQUARE_ANTIPRISM)
    .value("UNKNOWN", CoordGeometry::UNKNOWN);

  nb::class_<CodAtomInfo>(m, "CodAtomInfo")
    .def(nb::init<>())
    .def_ro("index", &CodAtomInfo::index)
    .def_ro("hashing_value", &CodAtomInfo::hashing_value)
    .def_ro("el", &CodAtomInfo::el)
    .def_ro("hybrid", &CodAtomInfo::hybrid)
    .def_ro("cod_class", &CodAtomInfo::cod_class)
    .def_ro("cod_main", &CodAtomInfo::cod_main)
    .def_ro("nb_symb", &CodAtomInfo::nb_symb)
    .def_ro("nb2_symb", &CodAtomInfo::nb2_symb)
    .def_ro("connectivity", &CodAtomInfo::connectivity)
    .def_ro("min_ring_size", &CodAtomInfo::min_ring_size)
    .def_ro("is_aromatic", &CodAtomInfo::is_aromatic)
    .def_ro("is_metal", &CodAtomInfo::is_metal)
    .def("__repr__", [](const CodAtomInfo& self) {
        return "<gemmi.CodAtomInfo " + self.cod_class +
               " hash=" + std::to_string(self.hashing_value) + ">";
    });

  nb::class_<ValueStats>(m, "ValueStats")
    .def(nb::init<>())
    .def(nb::init<double, double, int>())
    .def_ro("value", &ValueStats::value)
    .def_ro("sigma", &ValueStats::sigma)
    .def_ro("count", &ValueStats::count)
    .def("__repr__", [](const ValueStats& self) {
        return "<gemmi.ValueStats value=" + std::to_string(self.value) +
               " sigma=" + std::to_string(self.sigma) +
               " count=" + std::to_string(self.count) + ">";
    });

  nb::class_<AcedrgTables>(m, "AcedrgTables")
    .def(nb::init<>())
    .def("load_tables", &AcedrgTables::load_tables,
         nb::arg("tables_dir"),
         "Load COD/CSD statistical tables from the given directory.\n"
         "This directory should contain allOrgBondsHRS.table, "
         "allOrgAnglesHRS.table, etc.")
    .def("tables_loaded", [](const AcedrgTables& self) { return self.tables_loaded_; },
         "Check if tables have been loaded.")
    .def("tables_dir", [](const AcedrgTables& self) { return self.tables_dir_; },
         "Get the directory from which tables were loaded.")
    .def("classify_atoms", &AcedrgTables::classify_atoms,
         nb::arg("chemcomp"),
         "Classify all atoms in a ChemComp and return CodAtomInfo for each.")
    .def("fill_restraints", &AcedrgTables::fill_restraints,
         nb::arg("chemcomp"),
         "Fill missing bond and angle values in a ChemComp using COD/CSD statistics.")
    .def("fill_bond",
         [](const AcedrgTables& self, const ChemComp& cc, Restraints::Bond& bond) {
           auto atom_info = self.classify_atoms(cc);
           self.fill_bond(cc, atom_info, bond);
         },
         nb::arg("chemcomp"), nb::arg("bond"),
         "Fill a single bond value using COD/CSD statistics.")
    .def("fill_angle",
         [](const AcedrgTables& self, const ChemComp& cc, Restraints::Angle& angle) {
           auto atom_info = self.classify_atoms(cc);
           self.fill_angle(cc, atom_info, angle);
         },
         nb::arg("chemcomp"), nb::arg("angle"),
         "Fill a single angle value using COD/CSD statistics.")
    .def("compute_acedrg_types", &AcedrgTables::compute_acedrg_types,
         nb::arg("chemcomp"),
         "Compute acedrg_type strings for all atoms in a ChemComp.\n"
         "Format: CentralElement(Neighbor1_desc)(Neighbor2_desc)...\n"
         "Example: C(CH3)(COO)(NHH)(H) for alanine CA.")
    .def_rw("upper_bond_sigma", &AcedrgTables::upper_bond_sigma,
            "Maximum sigma for bond restraints (default 0.2)")
    .def_rw("lower_bond_sigma", &AcedrgTables::lower_bond_sigma,
            "Minimum sigma for bond restraints (default 0.02)")
    .def_rw("upper_angle_sigma", &AcedrgTables::upper_angle_sigma,
            "Maximum sigma for angle restraints (default 3.0)")
    .def_rw("lower_angle_sigma", &AcedrgTables::lower_angle_sigma,
            "Minimum sigma for angle restraints (default 1.5)")
    .def_rw("min_observations_bond", &AcedrgTables::min_observations_bond,
            "Minimum observation count for bonds (default 4)")
    .def_rw("min_observations_angle", &AcedrgTables::min_observations_angle,
            "Minimum observation count for angles (default 3)")
    .def_rw("min_observations_angle_fallback", &AcedrgTables::min_observations_angle_fallback,
            "Minimum observation count for angle fallback levels (default 3)")
    ;

  // Utility functions
  m.def("hybridization_to_string", &hybridization_to_string,
        nb::arg("h"),
        "Convert Hybridization enum to string.");
  m.def("hybridization_from_string", &hybridization_from_string,
        nb::arg("s"),
        "Convert string to Hybridization enum.");
}
