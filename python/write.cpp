// Copyright 2018 Global Phasing Ltd.

#define GEMMI_WRITE_IMPLEMENTATION
#include "gemmi/sprintf.hpp"
#include "gemmi/to_mmcif.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/mtz.hpp"
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
    ;

  structure
    .def("make_pdb_headers", &make_pdb_headers)
    .def("write_pdb", [](const Structure& st, const std::string& path,
                         bool seqres_records, bool ssbond_records,
                         bool link_records, bool cispep_records,
                         bool ter_records, bool numbered_ter,
                         bool ter_ignores_type, bool use_linkr) {
       PdbWriteOptions options;
       options.seqres_records = seqres_records;
       options.ssbond_records = ssbond_records;
       options.link_records = link_records;
       options.cispep_records = cispep_records;
       options.ter_records = ter_records;
       options.numbered_ter = numbered_ter;
       options.ter_ignores_type = ter_ignores_type;
       options.use_linkr = use_linkr;
       Ofstream f(path);
       write_pdb(st, f.ref(), options);
    }, py::arg("path"),
       py::arg("seqres_records")=true, py::arg("ssbond_records")=true,
       py::arg("link_records")=true, py::arg("cispep_records")=true,
       py::arg("ter_records")=true, py::arg("numbered_ter")=true,
       py::arg("ter_ignores_type")=false, py::arg("use_linkr")=false)
    .def("write_minimal_pdb",
         [](const Structure& st, const std::string& path) {
       Ofstream f(path);
       write_minimal_pdb(st, f.ref());
    }, py::arg("path"))
    .def("make_minimal_pdb", [](const Structure& st) -> std::string {
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
