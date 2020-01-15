// Copyright 2017 Global Phasing Ltd.

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/json.hpp"
#include "gemmi/gzread.hpp"
#include "gemmi/gzread_impl.hpp"
#include "gemmi/smcif.hpp"         // for make_small_structure_from_block
#include "gemmi/chemcomp_xyz.hpp"
#include "gemmi/remarks.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace gemmi;

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
  m.def("read_structure", [](const std::string& path, bool merge) {
          Structure* st = new Structure(read_structure_gz(path));
          if (!st->raw_remarks.empty())
            read_metadata_from_remarks(*st);
          if (merge)
            st->merge_chain_parts();
          return st;
        }, py::arg("path"), py::arg("merge_chain_parts")=true,
        "Reads a coordinate file into Structure.");
  m.def("make_structure_from_block", &make_structure_from_block,
        py::arg("block"), "Takes mmCIF block and returns Structure.");
  m.def("read_pdb_string", [](const std::string& s) {
          Structure* st = new Structure(read_pdb_string(s, "string"));
          read_metadata_from_remarks(*st);
          return st;
        }, py::arg("s"), "Reads a string as PDB file.");

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
