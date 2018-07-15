// Copyright 2017 Global Phasing Ltd.

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/json.hpp"
#define GEMMI_GZREAD_IMPLEMENTATION
#include "gemmi/gzread.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace gemmi;

void add_cif_read(py::module& cif) {
  cif.def("read_file", &cif::read_file, py::arg("filename"),
          "Reads a CIF file copying data into Document.");
  cif.def("read", [](const std::string& s) {
              return read_cif_or_mmjson_gz(s);
          }, py::arg("filename"), "Reads normal or gzipped CIF file.");
  cif.def("read_mmjson", [](const std::string& s) {
              return read_mmjson_gz(s);
          }, py::arg("filename"), "Reads normal or gzipped mmJSON file.");
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
}

void add_read_structure(py::module& m) {
  m.def("read_structure", [](const std::string& path, bool merge) {
          auto st = new Structure(read_structure_gz(path));
          if (merge)
            st->merge_same_name_chains();
          return st;
        }, py::arg("path"), py::arg("merge_same_name_chains")=true,
        "Reads a coordinate file into Structure.");
  m.def("make_structure_from_block", &make_structure_from_block,
        py::arg("block"), "Takes mmCIF block and returns Structure.");
  m.def("read_pdb_string", [](const std::string& s) {
          return new Structure(read_pdb_string(s, "string"));
        }, py::arg("s"), "Reads a string as PDB file.");

  // and an unrelated function from gz.hpp
  m.def("estimate_uncompressed_size", &estimate_uncompressed_size,
        py::arg("path"),
        "Returns uncompressed size of a .gz file (not always reliable)");
}
