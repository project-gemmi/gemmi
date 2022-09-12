// Copyright 2017 Global Phasing Ltd.

#include "gemmi/version.hpp"   // for GEMMI_VERSION
#include "gemmi/math.hpp"      // for hc
#include "gemmi/dirwalk.hpp"   // for CifWalk, CoorFileWalk
#include "gemmi/fileutil.hpp"  // for expand_if_pdb_code
#include "gemmi/bessel.hpp"    // for bessel_i1_over_i0
#include "gemmi/third_party/tao/pegtl/parse_error.hpp" // for parse_error

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/numpy.h>    // for vectorize

namespace py = pybind11;

namespace {
template<typename T> int get_max_bin(const T& bins) {
  int max_bin = 0;
  for (int i = 0; i < bins.shape(0); ++i) {
    if (bins(i) < 0)
      throw py::value_error("bins argument must have no negative elements");
    max_bin = std::max(max_bin, bins(i));
  }
  if (max_bin > 1000000)
    throw py::value_error("bin numbers must be smaller than million");
  return max_bin;
}
} // anonymous namespace

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
  m.def("log_bessel_i0", py::vectorize(gemmi::log_bessel_i0));
  m.def("log_cosh", py::vectorize([](double x) {
        // ln(cosh(x)) = ln(e^x + e^-x) - ln(2) = ln(e^x * (1 + e^-2x)) - ln(2)
        x = std::abs(x);
        return x + std::log1p(std::exp(-2 * x)) - std::log(2);
  }));

  // utilities inspired by numpy.bincount()
  m.def("binmean", [](py::array_t<int> bins, py::array_t<double> values) {
      auto bins_ = bins.unchecked<1>();
      auto values_ = values.unchecked<1>();
      auto len = bins_.shape(0);
      if (len != values_.shape(0))
        throw std::domain_error("arrays have different lengths");
      int ret_size = get_max_bin(bins_) + 1;
      py::array_t<double> ret(ret_size);
      double* retp = (double*) ret.request().ptr;
      for (int i = 0; i != ret_size; ++i)
        retp[i] = 0.;
      std::vector<int> counts(ret_size);
      for (int i = 0; i != len; ++i)
        if (!std::isnan(values_(i))) {
          int n = bins_(i);
          counts[n]++;
          retp[n] += values_(i);
        }
      for (int i = 0; i != ret_size; ++i)
        retp[i] /= counts[i];
      return ret;
  }, py::arg("nbins"), py::arg("values"));

  m.def("binrfactor", [](py::array_t<int> bins, py::array_t<double> obs,
                         py::array_t<double> calc, bool riso) {
      auto bins_ = bins.unchecked<1>();
      auto obs_ = obs.unchecked<1>();
      auto calc_ = calc.unchecked<1>();
      auto len = bins_.shape(0);
      if (len != obs_.shape(0) || len != calc_.shape(0))
        throw std::domain_error("arrays have different lengths");
      int ret_size = get_max_bin(bins_) + 1;
      py::array_t<double> ret(ret_size);
      double* retp = (double*) ret.request().ptr;
      for (int i = 0; i != ret_size; ++i)
        retp[i] = 0.;
      std::vector<double> sum_fobs(ret_size);
      for (int i = 0; i != len; ++i)
        if (!std::isnan(obs_(i)) && !std::isnan(calc_(i))) {
          int n = bins_(i);
          retp[n] += std::fabs(obs_(i) - calc_(i));
          sum_fobs[n] += riso ? (obs_(i) + calc_(i)) : obs_(i);
        }
      for (int i = 0; i != ret_size; ++i)
        retp[i] /= (riso ? 0.5 * sum_fobs[i] : sum_fobs[i]);
      return ret;
  }, py::arg("nbins"), py::arg("obs"), py::arg("calc"), py::arg("riso")=false);

  m.def("bincorr", [](py::array_t<int> bins, py::array_t<double> obs,
                                             py::array_t<double> calc) {
      auto bins_ = bins.unchecked<1>();
      auto obs_ = obs.unchecked<1>();
      auto calc_ = calc.unchecked<1>();
      auto len = bins_.shape(0);
      if (len != obs_.shape(0) || len != calc_.shape(0))
        throw std::domain_error("arrays have different lengths");
      int ret_size = get_max_bin(bins_) + 1;
      std::vector<gemmi::Correlation> cor(ret_size);
      for (int i = 0; i != len; ++i)
        if (!std::isnan(obs_(i)) && !std::isnan(calc_(i)))
          cor[bins_(i)].add_point(obs_(i), calc_(i));
      return cor;
  }, py::arg("nbins"), py::arg("obs"), py::arg("calc"));
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
  add_search(mg);
  add_read_structure(mg);
  add_scaling(mg);
  add_custom(mg);
}
