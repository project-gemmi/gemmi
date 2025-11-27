// Copyright Global Phasing Ltd.

#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "array.h"
#include "gemmi/read_cif.hpp"      // for read_cif_gz
#include "gemmi/smcif.hpp"         // for make_small_structure_from_block
#include "gemmi/small.hpp"         // for SmallStructure
#include "gemmi/gz.hpp"            // for estimate_uncompressed_size
#include "gemmi/interop.hpp"       // for atom_to_site, mx_to_sx_structure
#include "gemmi/flat.hpp"          // for FlatStructure, FlatAtom
#include "gemmi/pymol_select.hpp"  // for select_atoms

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<SmallStructure::Site>)
NB_MAKE_OPAQUE(std::vector<FlatAtom>)

void add_small(nb::module_& m) {
  // from smcif.hpp
  m.def("read_small_structure", [](const std::string& path) {
          cif::Block block = read_cif_gz(path).sole_block();
          return new SmallStructure(make_small_structure_from_block(block));
        }, nb::arg("path"), "Reads a small molecule CIF file.");
  m.def("make_small_structure_from_block", &make_small_structure_from_block,
        nb::arg("block"), "Takes CIF block and returns SmallStructure.");

  // and an unrelated function from gz.hpp
  m.def("estimate_uncompressed_size", &estimate_uncompressed_size,
        nb::arg("path"),
        "Returns uncompressed size of a .gz file (not always reliable)");

  using gemmi::SmallStructure;
  nb::class_<SmallStructure> small_structure(m, "SmallStructure");
  nb::class_<SmallStructure::Site>(small_structure, "Site")
    .def(nb::init<>())
    .def("__init__", [](SmallStructure::Site* site, const Atom& atom, const UnitCell& cell) {
        new(site) SmallStructure::Site(gemmi::atom_to_site(atom, cell));
    })
    .def_rw("label", &SmallStructure::Site::label)
    .def_rw("type_symbol", &SmallStructure::Site::type_symbol)
    .def_rw("fract", &SmallStructure::Site::fract)
    .def_rw("occ", &SmallStructure::Site::occ)
    .def_rw("u_iso", &SmallStructure::Site::u_iso)
    .def_rw("element", &SmallStructure::Site::element)
    .def_rw("charge", &SmallStructure::Site::charge)
    .def_rw("disorder_group", &SmallStructure::Site::disorder_group)
    .def_rw("aniso", &SmallStructure::Site::aniso)
    .def("orth", &SmallStructure::Site::orth)
    .def("clone", [](const SmallStructure::Site& self) { return new SmallStructure::Site(self); })
    .def("__repr__", [](const SmallStructure::Site& self) {
        return "<gemmi.SmallStructure.Site " + self.label + ">";
    });
  nb::bind_vector<std::vector<SmallStructure::Site>, rv_ri>(small_structure, "SiteList");

  using AtomType = SmallStructure::AtomType;
  nb::class_<AtomType>(small_structure, "AtomType")
    .def_ro("symbol", &AtomType::symbol)
    .def_ro("element", &AtomType::element)
    .def_rw("dispersion_real", &AtomType::dispersion_real)
    .def_rw("dispersion_imag", &AtomType::dispersion_imag)
    .def("__repr__", [](const AtomType& self) {
        return "<gemmi.SmallStructure.AtomType " + self.symbol + ">";
    });

  small_structure
    .def(nb::init<>())
    .def_rw("name", &SmallStructure::name)
    .def_rw("cell", &SmallStructure::cell)
    .def_ro("spacegroup", &SmallStructure::spacegroup,
                  nb::rv_policy::reference_internal)
    .def_rw("spacegroup_hm", &SmallStructure::spacegroup_hm)
    .def_rw("spacegroup_hall", &SmallStructure::spacegroup_hall)
    .def_rw("spacegroup_number", &SmallStructure::spacegroup_number)
    .def_rw("symops", &SmallStructure::symops)
    .def_rw("sites", &SmallStructure::sites)
    .def_ro("atom_types", &SmallStructure::atom_types)
    .def_rw("wavelength", &SmallStructure::wavelength)
    .def("add_site", [](SmallStructure& self, const SmallStructure::Site& site) {
        self.sites.push_back(site);
    })
    .def("determine_and_set_spacegroup", &SmallStructure::determine_and_set_spacegroup,
         nb::arg("order"))
    .def("check_spacegroup", &SmallStructure::check_spacegroup)
    .def("get_atom_type", &SmallStructure::get_atom_type)
    .def("get_all_unit_cell_sites", &SmallStructure::get_all_unit_cell_sites)
    .def("remove_hydrogens", &SmallStructure::remove_hydrogens)
    .def("change_occupancies_to_crystallographic",
         &SmallStructure::change_occupancies_to_crystallographic,
         nb::arg("max_dist")=0.4)
    .def("make_cif_block", &make_cif_block_from_small_structure)
    .def("__repr__", [](const SmallStructure& self) {
        return "<gemmi.SmallStructure: " + std::string(self.name) + ">";
    });
  m.def("mx_to_sx_structure", &gemmi::mx_to_sx_structure,
        nb::arg("st"), nb::arg("n")=0);

  // FlatStructure bindings
  nb::class_<FlatStructure>(m, "FlatStructure")
    .def(nb::init<const Structure&>(), nb::arg("structure"),
         "Create a flat representation of a Structure")
    .def("generate_structure", &FlatStructure::generate_structure,
         "Reconstructs a Structure from the flat table of atoms")
    .def("pymol_select", &remove_not_selected)
    .def("__len__", [](const FlatStructure& self) { return self.table.size(); })
    .def("__repr__", [](const FlatStructure& self) {
        return "<gemmi.FlatStructure with " + std::to_string(self.table.size()) + " atoms>";
    })
    // NumPy-like array properties for atomic data
    .def_prop_ro("b_iso", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::b_iso);
    }, nb::rv_policy::reference_internal, "B-factors as numpy array")
    .def_prop_ro("occ", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::occ);
    }, nb::rv_policy::reference_internal, "Occupancies as numpy array")
    .def_prop_ro("pos", [](FlatStructure& self) {
        // Create a view of positions as (N, 3) array
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom) / sizeof(double));
        return nb::ndarray<nb::numpy, double, nb::shape<-1, 3>>(
            &(self.table.data()->pos.x),
            {self.table.size(), 3},
            nb::handle(),
            {stride, 1});
    }, nb::rv_policy::reference_internal, "Positions as (N, 3) numpy array")
    .def_prop_ro("charge", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::charge);
    }, nb::rv_policy::reference_internal, "Charges as numpy array")
    .def_prop_ro("model_num", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::model_num);
    }, nb::rv_policy::reference_internal, "Model numbers as numpy array")
    .def_prop_ro("selected", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::selected);
    }, nb::rv_policy::reference_internal, "Selection flags as numpy array")
    .def_prop_ro("serials", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::serial);
    }, nb::rv_policy::reference_internal, "Serial numbers as numpy array")
    .def_prop_ro("altlocs", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::altloc);
    }, nb::rv_policy::reference_internal, "Alternate location indicators as numpy array")
    .def_prop_ro("het_flags", [](FlatStructure& self) {
        return vector_member_array(self.table, &FlatAtom::het_flag);
    }, nb::rv_policy::reference_internal, "Het flags as numpy array ('A'=ATOM, 'H'=HETATM)")
    .def_prop_ro("elements", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, uint8_t, nb::shape<-1>>(
            reinterpret_cast<uint8_t*>(&self.table.data()->element),
            {self.table.size()},
            nb::handle(),
            {stride});
    }, nb::rv_policy::reference_internal, "Element types as numpy array")
    .def_prop_ro("entity_type", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, uint8_t, nb::shape<-1>>(
            reinterpret_cast<uint8_t*>(&self.table.data()->entity_type),
            {self.table.size()},
            nb::handle(),
            {stride});
    }, nb::rv_policy::reference_internal, "Entity types as numpy array")
    .def_prop_ro("resnums", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom) / sizeof(int));
        return nb::ndarray<nb::numpy, int, nb::shape<-1>>(
            &self.table.data()->seq_id.num.value,
            {self.table.size()},
            nb::handle(),
            {stride});
    }, nb::rv_policy::reference_internal, "Residue sequence numbers as numpy array")
    .def_prop_ro("icodes", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, char, nb::shape<-1>>(
            &self.table.data()->seq_id.icode,
            {self.table.size()},
            nb::handle(),
            {stride});
    }, nb::rv_policy::reference_internal, "Insertion codes as numpy array")
    .def_prop_ro("element_weights", [](FlatStructure& self) {
        std::vector<double> weights;
        weights.reserve(self.table.size());
        for (const auto& atom : self.table)
            weights.push_back(molecular_weight(atom.element));
        return numpy_array_from_vector(std::move(weights));
    }, nb::rv_policy::move, "Element molecular weights as numpy array")
    // String fields as S8 (8-byte fixed-width string) numpy arrays
    .def_prop_ro("atom_names", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        auto raw = nb::ndarray<nb::numpy, char, nb::shape<-1, 8>>(
            self.table.data()->atom_name,
            {self.table.size(), 8},
            nb::handle(),
            {stride, 1});
        return nb::cast(raw).attr("view")("S8").attr("ravel")();
    }, nb::rv_policy::reference_internal, "Atom names as (N, 8) char array")
    .def_prop_ro("residue_names", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, char, nb::shape<-1, 8>>(
            self.table.data()->residue_name,
            {self.table.size(), 8},
            nb::handle(),
            {stride, 1});
    }, nb::rv_policy::reference_internal, "Residue names as (N, 8) char array")
    .def_prop_ro("chain_ids", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, char, nb::shape<-1, 8>>(
            self.table.data()->chain_id,
            {self.table.size(), 8},
            nb::handle(),
            {stride, 1});
    }, nb::rv_policy::reference_internal, "Chain IDs as (N, 8) char array")
    .def_prop_ro("subchains", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, char, nb::shape<-1, 8>>(
            self.table.data()->subchain,
            {self.table.size(), 8},
            nb::handle(),
            {stride, 1});
    }, nb::rv_policy::reference_internal, "Subchain IDs as (N, 8) char array")
    .def_prop_ro("entity_ids", [](FlatStructure& self) {
        constexpr int64_t stride = static_cast<int64_t>(sizeof(FlatAtom));
        return nb::ndarray<nb::numpy, char, nb::shape<-1, 8>>(
            self.table.data()->entity_id,
            {self.table.size(), 8},
            nb::handle(),
            {stride, 1});
    }, nb::rv_policy::reference_internal, "Entity IDs as (N, 8) char array");
}

