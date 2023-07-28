// Copyright 2018 Global Phasing Ltd.

#include <sstream>  // for ostringstream
#include "gemmi/to_mmcif.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/fstream.hpp"

#include "common.h"

namespace py = pybind11;
using namespace gemmi;

void add_write(py::module& m, py::class_<Structure>& structure) {
  py::class_<MmcifOutputGroups>(m, "MmcifOutputGroups")
    .def(py::init([](bool all, const py::kwargs& kwargs) {
      MmcifOutputGroups g(all);
      if (kwargs) {
        py::object pyg = py::cast(&g);
        for (auto kwarg : kwargs)
          pyg.attr(kwarg.first) = kwarg.second.cast<bool>();
      }
      return g;
    }), py::arg("all"))
#define DEF_BIT_PROPERTY(name) \
  .def_property(#name, [](MmcifOutputGroups g) { return g.name; }, \
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

  py::class_<PdbWriteOptions>(m, "PdbWriteOptions")
    .def(py::init([](bool minimal, bool headers_only, const py::kwargs& kwargs) {
      PdbWriteOptions opt;
      if (minimal)
        opt = PdbWriteOptions::minimal();
      else if (headers_only)
        opt = PdbWriteOptions::headers_only();
      if (kwargs) {
        py::object py_opt = py::cast(&opt);
        for (auto kwarg : kwargs)
          py_opt.attr(kwarg.first) = kwarg.second.cast<bool>();
      }
      return opt;
    }), py::arg("minimal")=false, py::arg("headers_only")=false)
#define DEF_PROPERTY(name) \
  .def_property(#name, [](PdbWriteOptions g) { return g.name; }, \
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
         py::arg_v("options", PdbWriteOptions(), "PdbWriteOptions()"))
    .def("write_pdb", [](const Structure& st, const std::string& path,
                         PdbWriteOptions options) {
        Ofstream f(path);
        write_pdb(st, f.ref(), options);
    })

    // deprecated - kept for compatibility
    .def("write_pdb", [](const Structure& st, const std::string& path,
                         const py::kwargs& kwargs) {
        PdbWriteOptions opt;
        if (kwargs) {
          py::object py_opt = py::cast(&opt);
          for (auto kwarg : kwargs)
            py_opt.attr(kwarg.first) = kwarg.second.cast<bool>();
        }
        Ofstream f(path);
        write_pdb(st, f.ref(), opt);
    }, py::arg("path"))
    // deprecated
    .def("make_pdb_headers", &make_pdb_headers)
    // deprecated
    .def("write_minimal_pdb", [](const Structure& st, const std::string& path) {
       Ofstream f(path);
       write_minimal_pdb(st, f.ref());
    }, py::arg("path"))
    // deprecated
    .def("make_minimal_pdb", [](const Structure& st) {
       std::ostringstream os;
       write_minimal_pdb(st, os);
       return os.str();
    })

    .def("make_mmcif_document", &make_mmcif_document,
         py::arg_v("groups", MmcifOutputGroups(true), "MmcifOutputGroups(True)"))
    .def("make_mmcif_block", &make_mmcif_block,
         py::arg_v("groups", MmcifOutputGroups(true), "MmcifOutputGroups(True)"))
    .def("update_mmcif_block", &update_mmcif_block, py::arg("block"),
         py::arg_v("groups", MmcifOutputGroups(true), "MmcifOutputGroups(True)"))
    .def("make_mmcif_headers", &make_mmcif_headers)
    ;
}
