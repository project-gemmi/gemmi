// Copyright 2017 Global Phasing Ltd.

#include "gemmi/version.hpp"   // for GEMMI_VERSION
#include "gemmi/math.hpp"      // for hc
#include "gemmi/dirwalk.hpp"   // for CifWalk, CoorFileWalk
#include "gemmi/fileutil.hpp"  // for expand_if_pdb_code
#include "gemmi/bessel.hpp"    // for bessel_i1_over_i0
#include "gemmi/third_party/tao/pegtl/parse_error.hpp" // for parse_error

#include "common.h"
#include <pybind11/numpy.h>    // for vectorize

namespace py = pybind11;

void add_misc(py::module& m) {
  py::class_<gemmi::CifWalk>(m, "CifWalk")
    .def(py::init<const char*>())
    .def("__iter__", [](gemmi::CifWalk& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
  py::class_<gemmi::CoorFileWalk>(m, "CoorFileWalk")
    .def(py::init<const char*>())
    .def("__iter__", [](gemmi::CoorFileWalk& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
  m.def("is_pdb_code", &gemmi::is_pdb_code);
  m.def("expand_pdb_code_to_path", &gemmi::expand_pdb_code_to_path);
  m.def("expand_if_pdb_code", &gemmi::expand_if_pdb_code,
        py::arg("code"), py::arg("filetype")='M');
  m.attr("hc") = py::float_(gemmi::hc());
  m.def("bessel_i1_over_i0", py::vectorize(gemmi::bessel_i1_over_i0));
}

PYBIND11_MODULE(gemmi, mg) {
  mg.doc() = "Python bindings to GEMMI - a library used in macromolecular\n"
             "crystallography and related fields";
  mg.attr("__version__") = GEMMI_VERSION;

  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p)
        std::rethrow_exception(p);
    } catch (const std::system_error &e) {
      const int errornum = e.code().value();
      PyErr_SetObject(PyExc_IOError, py::make_tuple(errornum, e.what()).ptr());
    } catch (const tao::pegtl::parse_error &e) {
      PyErr_SetString(PyExc_ValueError, e.what());
    }
  });

  py::module cif = mg.def_submodule("cif", "CIF file format");
  add_cif(cif);
  add_symmetry(mg);
  add_unitcell(mg);
  add_elem(mg);
  add_meta(mg);
  add_mol(mg);
  add_small(mg);
  add_misc(mg);
  add_grid(mg);
  add_recgrid(mg);
  add_ccp4(mg);
  add_sf(mg);
  add_cif_read(cif);
  add_mtz(mg);
  add_hkl(mg);
  add_chemcomp(mg);
  add_monlib(mg);
  add_topo(mg);
  add_alignment(mg);
  add_select(mg);
  add_search(mg);
  add_read_structure(mg);
  add_scaling(mg);
  add_custom(mg);
}
