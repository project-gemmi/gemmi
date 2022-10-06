// Copyright 2017 Global Phasing Ltd.

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/json.hpp"
#include "gemmi/smcif.hpp"         // for make_small_structure_from_block
#include "gemmi/small.hpp"         // for SmallStructure
#include "gemmi/interop.hpp"       // for atom_to_site, mx_to_sx_structure
#include "gemmi/chemcomp_xyz.hpp"

#define GEMMI_READ_CIF_IMPLEMENTATION
#include "gemmi/read_cif.hpp" // for read_cif_gz, read_mmjson_gz
#define GEMMI_READ_COOR_IMPLEMENTATION
#include "gemmi/read_coor.hpp"  // for read_structure_gz

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<SmallStructure::Site>)

void add_cif_read(py::module& cif) {
  cif.def("read_file", &cif::read_file, py::arg("filename"),
          "Reads a CIF file copying data into Document.");
  cif.def("read", &read_cif_or_mmjson_gz,
          py::arg("filename"), "Reads normal or gzipped CIF file.");
  cif.def("read_mmjson", &read_mmjson_gz,
          py::arg("filename"), "Reads normal or gzipped mmJSON file.");
  cif.def("read_string", &cif::read_string, py::arg("data"),
          "Reads a string as a CIF file.");

  cif.def("as_string", (std::string (*)(const std::string&)) &cif::as_string,
          py::arg("value"), "Get string content (no quotes) from raw string.");
  cif.def("as_number", &cif::as_number,
          py::arg("value"), py::arg("default")=NAN,
          "Returns float number from string");
  cif.def("as_int", (int (*)(const std::string&)) &cif::as_int,
          py::arg("value"), "Returns int number from string value.");
  cif.def("as_int", (int (*)(const std::string&, int)) &cif::as_int,
          py::arg("value"), py::arg("default"),
          "Returns int number from string value or the second arg if null.");
  cif.def("is_null", &cif::is_null, py::arg("value"));
}

void add_read_structure(py::module& m) {
  m.def("read_structure", [](const std::string& path, bool merge,
                             CoorFormat format, cif::Document* save_doc) {
          Structure* st = new Structure(read_structure_gz(path, format, save_doc));
          if (merge)
            st->merge_chain_parts();
          return st;
        }, py::arg("path"), py::arg("merge_chain_parts")=true,
           py::arg("format")=CoorFormat::Unknown,
           py::arg("save_doc")=nullptr,
        "Reads a coordinate file into Structure.");
  m.def("make_structure_from_block", &make_structure_from_block,
        py::arg("block"), "Takes mmCIF block and returns Structure.");
  m.def("read_pdb_string", [](const std::string& s, int max_line_length,
                              bool split_chain_on_ter) {
          PdbReadOptions options;
          options.max_line_length = max_line_length;
          options.split_chain_on_ter = split_chain_on_ter;
          return new Structure(read_pdb_string(s, "string", options));
        }, py::arg("s"), py::arg("max_line_length")=0,
           py::arg("split_chain_on_ter")=false, "Reads a string as PDB file.");
  m.def("read_pdb", [](const std::string& path, int max_line_length,
                       bool split_chain_on_ter) {
          PdbReadOptions options;
          options.max_line_length = max_line_length;
          options.split_chain_on_ter = split_chain_on_ter;
          return new Structure(read_pdb_gz(path, options));
        }, py::arg("filename"), py::arg("max_line_length")=0,
           py::arg("split_chain_on_ter")=false);

  // from smcif.hpp
  m.def("read_small_structure", [](const std::string& path) {
          cif::Block block = cif::read_file(path).sole_block();
          return new SmallStructure(make_small_structure_from_block(block));
        }, py::arg("path"), "Reads a small molecule CIF file.");
  m.def("make_small_structure_from_block", &make_small_structure_from_block,
        py::arg("block"), "Takes CIF block and returns SmallStructure.");

  // from chemcomp_xyz.hpp
  m.def("make_structure_from_chemcomp_block",
        &make_structure_from_chemcomp_block, py::arg("block"),
        "CIF block from CCD or monomer library -> single-residue Structure.");


  // and an unrelated function from gz.hpp
  m.def("estimate_uncompressed_size", &estimate_uncompressed_size,
        py::arg("path"),
        "Returns uncompressed size of a .gz file (not always reliable)");
}

void add_small(py::module& m) {
  using gemmi::SmallStructure;
  py::class_<SmallStructure> small_structure(m, "SmallStructure");
  py::class_<SmallStructure::Site>(small_structure, "Site")
    .def(py::init<>())
    .def(py::init(&gemmi::atom_to_site))
    .def_readwrite("label", &SmallStructure::Site::label)
    .def_readwrite("type_symbol", &SmallStructure::Site::type_symbol)
    .def_readwrite("fract", &SmallStructure::Site::fract)
    .def_readwrite("occ", &SmallStructure::Site::occ)
    .def_readwrite("u_iso", &SmallStructure::Site::u_iso)
    .def_readwrite("element", &SmallStructure::Site::element)
    .def_readwrite("charge", &SmallStructure::Site::charge)
    .def_readwrite("disorder_group", &SmallStructure::Site::disorder_group)
    .def_readwrite("aniso", &SmallStructure::Site::aniso)
    .def("orth", &SmallStructure::Site::orth)
    .def("__repr__", [](const SmallStructure::Site& self) {
        return "<gemmi.SmallStructure.Site " + self.label + ">";
    });
  py::bind_vector<std::vector<SmallStructure::Site>>(small_structure, "SiteList");

  using AtomType = SmallStructure::AtomType;
  py::class_<AtomType>(small_structure, "AtomType")
    .def_readonly("symbol", &AtomType::symbol)
    .def_readonly("element", &AtomType::element)
    .def_readwrite("dispersion_real", &AtomType::dispersion_real)
    .def_readwrite("dispersion_imag", &AtomType::dispersion_imag)
    .def("__repr__", [](const AtomType& self) {
        return "<gemmi.SmallStructure.AtomType " + self.symbol + ">";
    });

  small_structure
    .def(py::init<>())
    .def_readwrite("name", &SmallStructure::name)
    .def_readwrite("cell", &SmallStructure::cell)
    .def_readwrite("spacegroup_hm", &SmallStructure::spacegroup_hm)
    .def_readonly("sites", &SmallStructure::sites)
    .def_readonly("atom_types", &SmallStructure::atom_types)
    .def_readwrite("wavelength", &SmallStructure::wavelength)
    .def("add_site", [](SmallStructure& self, const SmallStructure::Site& site) {
        self.sites.push_back(site);
    })
    .def("find_spacegroup", &SmallStructure::find_spacegroup)
    .def("get_atom_type", &SmallStructure::get_atom_type)
    .def("get_all_unit_cell_sites", &SmallStructure::get_all_unit_cell_sites)
    .def("remove_hydrogens", &SmallStructure::remove_hydrogens)
    .def("change_occupancies_to_crystallographic",
         &SmallStructure::change_occupancies_to_crystallographic,
         py::arg("max_dist")=0.4)
    .def("setup_cell_images", &SmallStructure::setup_cell_images)
    .def("make_cif_block", &make_cif_block_from_small_structure)
    .def("__repr__", [](const SmallStructure& self) {
        return "<gemmi.SmallStructure: " + std::string(self.name) + ">";
    });
  m.def("mx_to_sx_structure", &gemmi::mx_to_sx_structure,
        py::arg("st"), py::arg("n")=0);
}

// used in cif.cpp
void cif_parse_string(cif::Document& doc, const std::string& data) {
  tao::pegtl::memory_input<> in(data, "string");
  cif::parse_input(doc, in);
}
void cif_parse_file(cif::Document& doc, const std::string& filename) {
  GEMMI_CIF_FILE_INPUT(in, filename);
  cif::parse_input(doc, in);
}
