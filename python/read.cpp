// Copyright 2017 Global Phasing Ltd.

#include "common.h"
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "gemmi/numb.hpp"
#include "gemmi/cifdoc.hpp"
#include "gemmi/mmcif.hpp"         // for make_structure_from_block, ...
#include "gemmi/pdb.hpp"           // for read_pdb_string
#include "gemmi/read_cif.hpp"      // for read_cif_gz, read_mmjson_gz
#include "gemmi/mmread_gz.hpp"     // for read_structure_gz
#include "gemmi/mmread.hpp"        // for read_structure_from_memory
#include "gemmi/json.hpp"          // for read_mmjson_insitu


using namespace gemmi;

void add_cif_read(nb::module_& cif) {
  cif.def("read_file", &read_cif_gz, nb::arg("filename"), nb::arg("check_level")=1,
          "Reads a CIF file copying data into Document.");
  cif.def("read", &read_cif_or_mmjson_gz,
          nb::arg("filename"), "Reads normal or gzipped CIF file.");
  cif.def("read_string", [](const std::string& str, int check_level) {
            return read_cif_from_memory(str.c_str(), str.size(), "string", check_level);
          }, nb::arg("string"), nb::arg("check_level")=1,
          "Reads a string as a CIF file.");
  cif.def("read_string", [](const nb::bytes& data, int check_level) {
            return read_cif_from_memory(data.c_str(), data.size(), "data", check_level);
          }, nb::arg("data"), nb::arg("check_level")=1,
          "Reads bytes as a CIF file.");
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

  m.def("read_structure_string", [](std::string& s, bool merge,
                                    CoorFormat format, cif::Document* save_doc) {
          Structure* st = new Structure(read_structure_from_memory(s.data(), s.size(), "string",
                                                                   format, save_doc));
          if (merge)
            st->merge_chain_parts();
          return st;
        }, nb::arg("path"), nb::arg("merge_chain_parts")=true,
           nb::arg("format")=CoorFormat::Unknown,
           nb::arg("save_doc")=nb::none(),
        "Reads a coordinate file into Structure.");
  m.def("read_structure_string", [](nb::bytes& s, bool merge,
                                    CoorFormat format, cif::Document* save_doc) {
          Structure* st = new Structure(read_structure_from_memory((char*)s.c_str(), s.size(), "string",
                                                                   format, save_doc));
          if (merge)
            st->merge_chain_parts();
          return st;
        }, nb::arg("path"), nb::arg("merge_chain_parts")=true,
           nb::arg("format")=CoorFormat::Unknown,
           nb::arg("save_doc")=nb::none(),
        "Reads a coordinate file into Structure.");
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
  m.def("populate_structure_from_block", &populate_structure_from_block,
        nb::arg("block"), nb::arg("st"));
  m.def("make_structure_from_chemcomp_block", &make_structure_from_chemcomp_block,
        nb::arg("block"), nb::arg("which")=7,
        "CIF block from CCD or monomer library -> single-residue Model(s).");

  m.def("populate_structure_from_pdb_stream", &populate_structure_from_pdb_stream,
        nb::arg("line_reader"), nb::arg("source"), nb::arg("st"), nb::arg("options"));
  m.def("read_pdb_string", [](const std::string& s, int max_line_length,
                              bool ignore_ter, bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, ignore_ter, split_chain_on_ter, false};
          return new Structure(read_pdb_string(s, "string", options));
        }, nb::arg("s"), nb::arg("max_line_length")=0,
           nb::arg("ignore_ter")=false, nb::arg("split_chain_on_ter")=false, "Reads a string as PDB file.");
  m.def("read_pdb_string", [](const nb::bytes& s, int max_line_length,
                              bool ignore_ter, bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, ignore_ter, split_chain_on_ter, false};
          return new Structure(read_pdb_from_memory(s.c_str(), s.size(), "string", options));
        }, nb::arg("s"), nb::arg("max_line_length")=0,
           nb::arg("ignore_ter")=false, nb::arg("split_chain_on_ter")=false, "Reads a string as PDB file.");
  m.def("read_pdb", [](const std::string& path, int max_line_length,
                       bool ignore_ter, bool split_chain_on_ter) {
          PdbReadOptions options{max_line_length, ignore_ter, split_chain_on_ter, false};
          return new Structure(read_pdb_gz(path, options));
        }, nb::arg("filename"), nb::arg("max_line_length")=0,
           nb::arg("ignore_ter")=false, nb::arg("split_chain_on_ter")=false);
}
