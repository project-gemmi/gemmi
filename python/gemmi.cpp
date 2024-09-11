// Copyright 2017 Global Phasing Ltd.

#include "gemmi/version.hpp"   // for GEMMI_VERSION
#include "gemmi/math.hpp"      // for hc
#include "gemmi/dirwalk.hpp"   // for CifWalk, CoorFileWalk
#include "gemmi/pdb_id.hpp"    // for expand_if_pdb_code
#include "gemmi/bessel.hpp"    // for bessel_i1_over_i0
#include "gemmi/pirfasta.hpp"  // for read_pir_or_fasta
#include "gemmi/seqtools.hpp"  // for calculate_sequence_weight
#include "gemmi/stats.hpp"     // for Correlation
#include "gemmi/third_party/tao/pegtl/parse_error.hpp" // for parse_error

#include "common.h"
#include "array.h"
#include <nanobind/make_iterator.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>  // for calculate_sequence_weight

namespace {
// returns max value in the array bin indices (which are of type int)
template<typename T> size_t get_max_bin(const T& bins) {
  int max_bin = 0;
  for (size_t i = 0; i < bins.shape(0); ++i) {
    if (bins(i) < 0)
      throw nb::value_error("bins argument must have no negative elements");
    max_bin = std::max(max_bin, bins(i));
  }
  if (max_bin > 1000000)
    throw nb::value_error("bin numbers must be smaller than million");
  return (size_t) max_bin;
}
} // anonymous namespace

void add_misc(nb::module_& m) {
  nb::class_<gemmi::CifWalk>(m, "CifWalk")
    .def(nb::init<const char*, char>(), nb::arg("path"), nb::arg("try_pdbid")='\0')
    .def("__iter__", [](gemmi::CifWalk& self) {
        return nb::make_iterator(nb::type<gemmi::CifWalk>(), "iterator", self);
    }, nb::keep_alive<0, 1>());
  nb::class_<gemmi::CoorFileWalk>(m, "CoorFileWalk")
    .def(nb::init<const char*, char>(), nb::arg("path"), nb::arg("try_pdbid")='\0')
    .def("__iter__", [](gemmi::CoorFileWalk& self) {
        return nb::make_iterator(nb::type<gemmi::CoorFileWalk>(), "iterator", self);
    }, nb::keep_alive<0, 1>());
  m.def("is_pdb_code", &gemmi::is_pdb_code);
  m.def("expand_pdb_code_to_path", &gemmi::expand_pdb_code_to_path,
        nb::arg("code"), nb::arg("filetype"), nb::arg("throw_if_unset")=false);
  m.def("expand_if_pdb_code", &gemmi::expand_if_pdb_code,
        nb::arg("code"), nb::arg("filetype")='M');
  m.attr("hc") = nb::float_(gemmi::hc());
  //m.def("bessel_i1_over_i0", nb::vectorize(gemmi::bessel_i1_over_i0));
  //m.def("log_bessel_i0", nb::vectorize(gemmi::log_bessel_i0));
  //m.def("log_cosh", nb::vectorize(gemmi::log_cosh));

  // pirfasta.hpp
  nb::class_<gemmi::FastaSeq>(m, "FastaSeq")
    .def_ro("header", &gemmi::FastaSeq::header)
    .def_ro("seq", &gemmi::FastaSeq::seq)
    ;
  m.def("read_pir_or_fasta", &gemmi::read_pir_or_fasta);

  // seqtools.hpp
  m.def("calculate_sequence_weight", &gemmi::calculate_sequence_weight,
        nb::arg("sequence"), nb::arg("unknown")=0.);
  m.def("one_letter_code", &gemmi::one_letter_code);
  m.def("pdbx_one_letter_code", &gemmi::pdbx_one_letter_code);
  m.def("sequence_kind", &gemmi::sequence_kind);

  // stats.hpp
  nb::class_<gemmi::Correlation>(m, "Correlation")
    .def_ro("n", &gemmi::Correlation::n)
    .def("coefficient", &gemmi::Correlation::coefficient)
    .def("mean_ratio", &gemmi::Correlation::mean_ratio)
    ;

  // utilities inspired by numpy.bincount()
  m.def("binmean", [](cpu_array<int> bins, cpu_array<double> values) {
      auto bins_ = bins.view();
      auto values_ = values.view();
      auto len = bins_.shape(0);
      if (len != values_.shape(0))
        throw std::domain_error("arrays have different lengths");
      size_t ret_size = get_max_bin(bins_) + 1;
      auto ret = make_numpy_array<double>({ret_size});
      double* retp = ret.data();
      for (size_t i = 0; i != ret_size; ++i)
        retp[i] = 0.;
      std::vector<int> counts(ret_size);
      for (size_t i = 0; i != len; ++i)
        if (!std::isnan(values_(i))) {
          int n = bins_(i);
          counts[n]++;
          retp[n] += values_(i);
        }
      for (size_t i = 0; i != ret_size; ++i)
        retp[i] /= counts[i];
      return ret;
  }, nb::arg("nbins"), nb::arg("values"));

  m.def("binrfactor", [](cpu_array<int> bins, cpu_array<double> obs,
                         cpu_array<double> calc, bool riso) {
      auto bins_ = bins.view();
      auto obs_ = obs.view();
      auto calc_ = calc.view();
      auto len = bins_.shape(0);
      if (len != obs_.shape(0) || len != calc_.shape(0))
        throw std::domain_error("arrays have different lengths");
      size_t ret_size = get_max_bin(bins_) + 1;
      auto ret = make_numpy_array<double>({ret_size});
      double* retp = ret.data();
      for (size_t i = 0; i != ret_size; ++i)
        retp[i] = 0.;
      std::vector<double> sum_fobs(ret_size);
      for (size_t i = 0; i != len; ++i)
        if (!std::isnan(obs_(i)) && !std::isnan(calc_(i))) {
          int n = bins_(i);
          retp[n] += std::fabs(obs_(i) - calc_(i));
          sum_fobs[n] += riso ? (obs_(i) + calc_(i)) : obs_(i);
        }
      for (size_t i = 0; i != ret_size; ++i)
        retp[i] /= (riso ? 0.5 * sum_fobs[i] : sum_fobs[i]);
      return ret;
  }, nb::arg("nbins"), nb::arg("obs"), nb::arg("calc"), nb::arg("riso")=false);

  m.def("bincorr", [](cpu_array<int> bins, cpu_array<double> obs, cpu_array<double> calc) {
      auto bins_ = bins.view();
      auto obs_ = obs.view();
      auto calc_ = calc.view();
      auto len = bins_.shape(0);
      if (len != obs_.shape(0) || len != calc_.shape(0))
        throw std::domain_error("arrays have different lengths");
      size_t ret_size = get_max_bin(bins_) + 1;
      std::vector<gemmi::Correlation> cor(ret_size);
      for (size_t i = 0; i != len; ++i)
        if (!std::isnan(obs_(i)) && !std::isnan(calc_(i)))
          cor[bins_(i)].add_point(obs_(i), calc_(i));
      return cor;
  }, nb::arg("nbins"), nb::arg("obs"), nb::arg("calc"));
}

NB_MODULE(gemmi, mg) {
  mg.doc() = "Python bindings to GEMMI - a library used in macromolecular\n"
             "crystallography and related fields";
  mg.attr("__version__") = GEMMI_VERSION;

  nb::register_exception_translator([](const std::exception_ptr& p, void*) {
    try {
      if (p)
        std::rethrow_exception(p);
    } catch (const std::system_error &e) {
      const int errornum = e.code().value();
      PyErr_SetObject(PyExc_IOError, nb::make_tuple(errornum, e.what()).ptr());
    } catch (const tao::pegtl::parse_error &e) {
      PyErr_SetString(PyExc_ValueError, e.what());
    }
  });

  nb::module_ cif = mg.def_submodule("cif", "CIF file format");
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
