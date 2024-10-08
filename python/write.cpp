// Copyright 2018 Global Phasing Ltd.

#include <cassert>
#include <sstream>  // for ostringstream
#include "gemmi/to_mmcif.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/fstream.hpp"

#include "common.h"
#include <nanobind/stl/string.h>

using namespace gemmi;

void add_write(nb::module_& m, nb::class_<Structure>& structure) {
  nb::class_<MmcifOutputGroups>(m, "MmcifOutputGroups")
    .def(nb::init<bool>())
    .def("__init__", [](MmcifOutputGroups* p, bool all, const nb::kwargs& kwargs) {
      nb::object obj = nb::type<MmcifOutputGroups>()(all);
      for (auto [key, value] : kwargs)
        obj.attr(key) = nb::cast<bool>(value);
      new(p) MmcifOutputGroups(nb::cast<const MmcifOutputGroups&>(obj));
    }, nb::arg("all"), nb::arg("kwargs"))
#define DEF_BIT_PROPERTY(name) \
  .def_prop_rw(#name, [](MmcifOutputGroups g) { return g.name; }, \
                      [](MmcifOutputGroups& g, bool v) { g.name = v; })
    DEF_BIT_PROPERTY(atoms)
    DEF_BIT_PROPERTY(block_name)
    DEF_BIT_PROPERTY(entry)
    DEF_BIT_PROPERTY(database_status)
    DEF_BIT_PROPERTY(author)
    DEF_BIT_PROPERTY(cell)
    DEF_BIT_PROPERTY(symmetry)
    DEF_BIT_PROPERTY(entity)
    DEF_BIT_PROPERTY(entity_poly)
    DEF_BIT_PROPERTY(struct_ref)
    DEF_BIT_PROPERTY(chem_comp)
    DEF_BIT_PROPERTY(exptl)
    DEF_BIT_PROPERTY(diffrn)
    DEF_BIT_PROPERTY(reflns)
    DEF_BIT_PROPERTY(refine)
    DEF_BIT_PROPERTY(title_keywords)
    DEF_BIT_PROPERTY(ncs)
    DEF_BIT_PROPERTY(struct_asym)
    DEF_BIT_PROPERTY(origx)
    DEF_BIT_PROPERTY(struct_conf)
    DEF_BIT_PROPERTY(struct_sheet)
    DEF_BIT_PROPERTY(struct_biol)
    DEF_BIT_PROPERTY(assembly)
    DEF_BIT_PROPERTY(conn)
    DEF_BIT_PROPERTY(cis)
    DEF_BIT_PROPERTY(scale)
    DEF_BIT_PROPERTY(atom_type)
    DEF_BIT_PROPERTY(entity_poly_seq)
    DEF_BIT_PROPERTY(tls)
    DEF_BIT_PROPERTY(software)
    DEF_BIT_PROPERTY(group_pdb)
    DEF_BIT_PROPERTY(auth_all)
    ;
#undef DEF_BIT_PROPERTY

  nb::class_<PdbWriteOptions>(m, "PdbWriteOptions")
    .def("__init__", [](PdbWriteOptions* opt, bool minimal, bool headers_only) {
      new(opt) PdbWriteOptions;
      if (minimal)
        *opt = PdbWriteOptions::minimal();
      else if (headers_only)
        *opt = PdbWriteOptions::headers_only();
    }, nb::arg("minimal")=false, nb::arg("headers_only")=false)
    .def("__init__", [](PdbWriteOptions* p, bool minimal, bool headers_only,
                        const nb::kwargs& kwargs) {
      nb::object obj = nb::type<PdbWriteOptions>()(minimal, headers_only);
      for (auto [key, value] : kwargs)
        obj.attr(key) = nb::cast<bool>(value);
      new(p) PdbWriteOptions(nb::cast<const PdbWriteOptions&>(obj));
    }, nb::arg("minimal")=false, nb::arg("headers_only")=false, nb::arg("kwargs"))
#define DEF_PROPERTY(name) \
  .def_prop_rw(#name, [](PdbWriteOptions g) { return g.name; }, \
                       [](PdbWriteOptions& g, bool v) { g.name = v; })
    DEF_PROPERTY(minimal_file)
    DEF_PROPERTY(atom_records)
    DEF_PROPERTY(seqres_records)
    DEF_PROPERTY(ssbond_records)
    DEF_PROPERTY(link_records)
    DEF_PROPERTY(cispep_records)
    DEF_PROPERTY(cryst1_record)
    DEF_PROPERTY(ter_records)
    DEF_PROPERTY(conect_records)
    DEF_PROPERTY(end_record)
    DEF_PROPERTY(numbered_ter)
    DEF_PROPERTY(ter_ignores_type)
    DEF_PROPERTY(use_linkr)
    DEF_PROPERTY(preserve_serial)
    ;
#undef DEF_PROPERTY

  structure
    .def("make_pdb_string", &make_pdb_string,
         nb::arg("options").sig("PdbWriteOptions()")=PdbWriteOptions())
    .def("write_pdb", [](const Structure& st, const std::string& path, PdbWriteOptions options) {
        Ofstream f(path);
        write_pdb(st, f.ref(), options);
    })
    // deprecated - kept for compatibility
    .def("write_pdb", [](const Structure& st, const std::string& path, const nb::kwargs& kwargs) {
        Ofstream f(path);
        nb::object options = nb::type<PdbWriteOptions>()(**kwargs);
        write_pdb(st, f.ref(), nb::cast<const PdbWriteOptions&>(options));
    }, nb::arg("path"), nb::arg("kwargs"))
    // deprecated
    .def("make_pdb_headers", &make_pdb_headers)
    // deprecated
    .def("write_minimal_pdb", [](const Structure& st, const std::string& path) {
       Ofstream f(path);
       write_minimal_pdb(st, f.ref());
    }, nb::arg("path"))
    // deprecated
    .def("make_minimal_pdb", [](const Structure& st) {
       std::ostringstream os;
       write_minimal_pdb(st, os);
       return os.str();
    })
    .def("make_mmcif_document", &make_mmcif_document,
         nb::arg("groups").sig("MmcifOutputGroups(True)")=MmcifOutputGroups(true))
    .def("make_mmcif_block", &make_mmcif_block,
         nb::arg("groups").sig("MmcifOutputGroups(True)")=MmcifOutputGroups(true))
    .def("update_mmcif_block", &update_mmcif_block, nb::arg("block"),
         nb::arg("groups").sig("MmcifOutputGroups(True)")=MmcifOutputGroups(true))
    .def("make_mmcif_headers", &make_mmcif_headers)
    ;
}
