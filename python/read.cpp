// Copyright 2017 Global Phasing Ltd.

#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/mmcif.hpp"         // for make_structure_from_block, ...
#include "gemmi/pdb.hpp"           // for read_pdb_string
#include "gemmi/gz.hpp"            // for estimate_uncompressed_size
#include "gemmi/smcif.hpp"         // for make_small_structure_from_block
#include "gemmi/small.hpp"         // for SmallStructure
#include "gemmi/interop.hpp"       // for atom_to_site, mx_to_sx_structure
#include "gemmi/read_cif.hpp"      // for read_cif_gz, read_mmjson_gz
#include "gemmi/mmread_gz.hpp"     // for read_structure_gz
#include "gemmi/json.hpp"          // for read_mmjson_insitu


using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<SmallStructure::Site>)

void add_cif_read(nb::module_& cif) {
  cif.def("read_file", &read_cif_gz, nb::arg("filename"),
          "Reads a CIF file copying data into Document.");
  cif.def("read", &read_cif_or_mmjson_gz,
          nb::arg("filename"), "Reads normal or gzipped CIF file.");
  cif.def("read_string", [](const std::string& str) {
            return read_cif_from_memory(str.c_str(), str.size(), "string");
          }, nb::arg("string"),
          "Reads a string as a CIF file.");
  cif.def("read_string", [](const nb::bytes& data) {
            return read_cif_from_memory(data.c_str(), data.size(), "data");
          }, nb::arg("data"), "Reads bytes as a CIF file.");
  cif.def("read_mmjson", &read_mmjson_gz,
          nb::arg("filename"), "Reads normal or gzipped mmJSON file.");
  cif.def("read_mmjson_string", [](std::string data) {
      return cif::read_mmjson_insitu(data.data(), data.size());
  });
  cif.def("read_mmjson_string", [](const nb::bytes& data) {
      std::string str(data.c_str(), data.size());
      return cif::read_mmjson_insitu(str.data(), str.size());
  });

  cif.def("as_string", (std::string (*)(const std::string&)) &cif::as_string,
          nb::arg("value"), "Get string content (no quotes) from raw string.");
  cif.def("as_number", &cif::as_number,
          nb::arg("value"), nb::arg("default")=NAN,
          "Returns float number from string");
  cif.def("as_int", (int (*)(const std::string&)) &cif::as_int,
          nb::arg("value"), "Returns int number from string value.");
  cif.def("as_int", (int (*)(const std::string&, int)) &cif::as_int,
          nb::arg("value"), nb::arg("default"),
          "Returns int number from string value or the second arg if null.");
  cif.def("is_null", &cif::is_null, nb::arg("value"));
}

void add_read_structure(nb::module_& m) {
  nb::enum_<ChemCompModel>(m, "ChemCompModel", nb::is_flag(), nb::is_arithmetic())
    .value("Xyz", ChemCompModel::Xyz)
    .value("Example", ChemCompModel::Example)
    .value("Ideal", ChemCompModel::Ideal);

  m.def("read_structure", [](const std::string& path, bool merge,
                             CoorFormat format, cif::Document* save_doc) {
          Structure* st = new Structure(read_structure_gz(path, format, save_doc));
          if (merge)
            st->merge_chain_parts();
          return st;
        }, nb::arg("path"), nb::arg("merge_chain_parts")=true,
           nb::arg("format")=CoorFormat::Unknown,
           nb::arg("save_doc")=nb::none(),
        "Reads a coordinate file into Structure.");
  m.def("make_structure_from_block", &make_structure_from_block,
        nb::arg("block"), "Takes mmCIF block and returns Structure.");
  m.def("make_structure_from_chemcomp_block", &make_structure_from_chemcomp_block,
        nb::arg("block"), nb::arg("which")=7,
        "CIF block from CCD or monomer library -> single-residue Model(s).");

  m.def("read_pdb_string", [](const std::string& s, int max_line_length,
                              bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, split_chain_on_ter, false};
          return new Structure(read_pdb_string(s, "string", options));
        }, nb::arg("s"), nb::arg("max_line_length")=0,
           nb::arg("split_chain_on_ter")=false, "Reads a string as PDB file.");
  m.def("read_pdb_string", [](const nb::bytes& s, int max_line_length,
                              bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, split_chain_on_ter, false};
          return new Structure(read_pdb_from_memory(s.c_str(), s.size(), "string", options));
        }, nb::arg("s"), nb::arg("max_line_length")=0,
           nb::arg("split_chain_on_ter")=false, "Reads a string as PDB file.");
  m.def("read_pdb", [](const std::string& path, int max_line_length,
                       bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, split_chain_on_ter, false};
          return new Structure(read_pdb_gz(path, options));
        }, nb::arg("filename"), nb::arg("max_line_length")=0,
           nb::arg("split_chain_on_ter")=false);

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
}

void add_small(nb::module_& m) {
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
}

