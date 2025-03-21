// Copyright 2019 Global Phasing Ltd.

#include "common.h"
#include "array.h"  // numpy_array_from_vector
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include "gemmi/unitcell.hpp"
#include "gemmi/refln.hpp"
#include "gemmi/fourier.hpp"  // for get_size_for_hkl, get_f_phi_on_grid, ...
#include "gemmi/fprime.hpp"   // for cromer_liberman
#include "gemmi/reciproc.hpp" // for count_reflections, make_miller_vector
#include "gemmi/cif2mtz.hpp"  // for CifToMtz
#include "gemmi/mtz2cif.hpp"  // for MtzToCif
#include "gemmi/intensit.hpp" // for Intensities
#include "gemmi/binner.hpp"   // for Binner
#include "gemmi/ecalc.hpp"    // for calculate_amplitude_normalizers

using namespace gemmi;

NB_MAKE_OPAQUE(std::vector<ReflnBlock>)

// HklMatch used to live in C++ library (in binner.hpp),
// but AFAICT it hasn't found good use in C++, so it was moved here.
struct HklMatch {
  std::vector<int> pos;
  size_t hkl_size;

  HklMatch(const Miller* hkl, size_t hkl_size_,
           const Miller* ref, size_t ref_size)
      : pos(ref_size, -1), hkl_size(hkl_size_) {
    // Usually, both datasets are sorted. This make things faster.
    if (std::is_sorted(hkl, hkl + hkl_size) &&
        std::is_sorted(ref, ref + ref_size)) {
      // cf. for_matching_reflections()
      const Miller* a = hkl;
      const Miller* b = ref;
      while (a != hkl + hkl_size && b != ref + ref_size) {
        if (*a == *b)
          pos[b++ - ref] = static_cast<int>(a++ - hkl);
        else if (*a < *b)
          ++a;
        else
          ++b;
      }
    } else {
      std::unordered_map<Miller, int, MillerHash> hkl_index;
      for (int i = 0; i != (int)hkl_size; ++i)
        hkl_index.emplace(hkl[i], i);
      for (size_t i = 0; i != ref_size; ++i) {
        auto it = hkl_index.find(ref[i]);
        if (it != hkl_index.end())
          pos[i] = it->second;
      }
    }
  }

  HklMatch(const std::vector<Miller>& hkl, const std::vector<Miller>& ref)
    : HklMatch(hkl.data(), hkl.size(), ref.data(), ref.size()) {}

  template <typename T> std::vector<T> aligned_(const T* v, size_t size, T nan) {
    if (size != hkl_size)
      fail("HklMatch.aligned(): wrong data, size differs");
    std::vector<T> result(pos.size());
    for (size_t i = 0; i != pos.size(); ++i)
      result[i] = pos[i] >= 0 ? v[pos[i]] : nan;
    return result;
  }
  template <typename T> std::vector<T> aligned(const std::vector<T>& v, T nan) {
    return aligned_(v.data(), v.size(), nan);
  }
};


void add_hkl(nb::module_& m) {
  nb::class_<ReflnBlock> pyReflnBlock(m, "ReflnBlock");
  nb::bind_vector<std::vector<ReflnBlock>, rv_ri>(m, "ReflnBlocks");
  pyReflnBlock
    .def_ro("block", &ReflnBlock::block)
    .def_ro("entry_id", &ReflnBlock::entry_id)
    .def_ro("cell", &ReflnBlock::cell)
    .def_ro("spacegroup", &ReflnBlock::spacegroup)
    .def_ro("wavelength", &ReflnBlock::wavelength)
    .def_ro("default_loop", &ReflnBlock::default_loop)
    .def("column_labels", &ReflnBlock::column_labels)
    .def("make_int_array",
         [](ReflnBlock& self, const std::string& tag, int null) {
           return numpy_array_from_vector(self.make_vector(tag, null));
    }, nb::arg("tag"), nb::arg("null"))
    .def("make_float_array",
         [](ReflnBlock& self, const std::string& tag, double null) {
           return numpy_array_from_vector(self.make_vector(tag, null));
    }, nb::arg("tag"), nb::arg("null")=NAN)
    .def("make_miller_array", [](ReflnBlock& self) {
        return py_array2d_from_vector(self.make_miller_vector());
    })
    .def("make_1_d2_array", [](ReflnBlock& self) {
        return numpy_array_from_vector(self.make_1_d2_vector());
    })
    .def("make_d_array", [](ReflnBlock& self) {
        return numpy_array_from_vector(self.make_d_vector());
    })
    .def("get_size_for_hkl",
         [](const ReflnBlock& self,
            std::array<int,3> min_size, double sample_rate) {
          return get_size_for_hkl(ReflnDataProxy(self), min_size, sample_rate);
    }, nb::arg("min_size")=std::array<int,3>{{0,0,0}},
       nb::arg("sample_rate")=0.)
    .def("data_fits_into", [](const ReflnBlock& self, std::array<int,3> size) {
        return data_fits_into(ReflnDataProxy(self), size);
    }, nb::arg("size"))
    .def("get_f_phi_on_grid", [](const ReflnBlock& self,
                                 const std::string& f_col,
                                 const std::string& phi_col,
                                 std::array<int, 3> size,
                                 bool half_l, AxisOrder order) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        FPhiProxy<ReflnDataProxy> fphi(ReflnDataProxy{self}, f_idx, phi_idx);
        return get_f_phi_on_grid<float>(fphi, size, half_l, order);
    }, nb::arg("f"), nb::arg("phi"), nb::arg("size"),
       nb::arg("half_l")=false, nb::arg("order")=AxisOrder::XYZ)
    .def("get_value_on_grid", [](const ReflnBlock& self,
                                 const std::string& column,
                                 std::array<int, 3> size,
                                 bool half_l, AxisOrder order) {
        size_t col_idx = self.get_column_index(column);
        return get_value_on_grid<float>(ReflnDataProxy(self), col_idx,
                                        size, half_l, order);
    }, nb::arg("column"), nb::arg("size"), nb::arg("half_l")=false,
       nb::arg("order")=AxisOrder::XYZ)
    .def("transform_f_phi_to_map", [](const ReflnBlock& self,
                                      const std::string& f_col,
                                      const std::string& phi_col,
                                      std::array<int, 3> min_size,
                                      std::array<int, 3> exact_size,
                                      double sample_rate,
                                      AxisOrder order) {
        size_t f_idx = self.get_column_index(f_col);
        size_t phi_idx = self.get_column_index(phi_col);
        FPhiProxy<ReflnDataProxy> fphi(ReflnDataProxy{self}, f_idx, phi_idx);
        return transform_f_phi_to_map2<float>(fphi, min_size, sample_rate,
                                              exact_size, order);
    }, nb::arg("f"), nb::arg("phi"),
       nb::arg("min_size")=std::array<int,3>{{0,0,0}},
       nb::arg("exact_size")=std::array<int,3>{{0,0,0}},
       nb::arg("sample_rate")=0.,
       nb::arg("order")=AxisOrder::XYZ)
    .def("get_float", &make_asu_data<float, ReflnBlock>,
         nb::arg("col"), nb::arg("as_is")=false)
    .def("get_int", &make_asu_data<int, ReflnBlock>,
         nb::arg("col"), nb::arg("as_is")=false)
    .def("get_f_phi", [](const ReflnBlock& self, const std::string& f_col,
                                                 const std::string& phi_col,
                                                 bool as_is) {
        return make_asu_data<std::complex<float>, 2>(self, {f_col, phi_col}, as_is);
    }, nb::arg("f"), nb::arg("phi"), nb::arg("as_is")=false)
    .def("get_value_sigma", [](const ReflnBlock& self, const std::string& f_col,
                                                       const std::string& sigma_col,
                                                       bool as_is) {
        return make_asu_data<ValueSigma<float>, 2>(self, {f_col, sigma_col}, as_is);
    }, nb::arg("f"), nb::arg("sigma"), nb::arg("as_is")=false)
    .def("is_merged", &ReflnBlock::is_merged)
    .def("is_unmerged", &ReflnBlock::is_unmerged)  // deprecated
    .def("use_unmerged", &ReflnBlock::use_unmerged)
    .def("__bool__", [](const ReflnBlock& self) { return self.ok(); })
    .def("__repr__", [](const ReflnBlock& self) {
        std::string s = cat("<gemmi.ReflnBlock ", self.block.name, " with ");
        if (self.default_loop)
          cat_to(s, self.default_loop->width(), " x ", self.default_loop->length());
        else
          s += " no ";
        s += " loop>";
        return s;
    });
  m.def("as_refln_blocks",
        [](cif::Document& d) { return as_refln_blocks(std::move(d.blocks)); });
  m.def("hkl_cif_as_refln_block", &hkl_cif_as_refln_block, nb::arg("block"));
  m.def("transform_f_phi_grid_to_map", [](FPhiGrid<float> grid) {
          return transform_f_phi_grid_to_map<float>(std::move(grid));
        }, nb::arg("grid"));
  m.def("transform_map_to_f_phi", &transform_map_to_f_phi<float>,
        nb::arg("map"), nb::arg("half_l")=false, nb::arg("use_scale")=true);
  m.def("cromer_liberman", [](int z, double energy) {
      std::pair<double, double> r;
      r.first = cromer_liberman(z, energy, &r.second);
      return r;
  }, nb::arg("z"), nb::arg("energy"));
  m.def("count_reflections", &count_reflections,
        nb::arg("cell"), nb::arg("spacegroup"), nb::arg("dmin"),
        nb::arg("dmax")=0., nb::arg("unique")=true);
  m.def("make_miller_array", [](const UnitCell& cell, const SpaceGroup* sg,
                                double dmin, double dmax, bool unique) {
      return py_array2d_from_vector(gemmi::make_miller_vector(cell, sg, dmin, dmax, unique));
  }, nb::arg("cell"), nb::arg("spacegroup"), nb::arg("dmin"),
     nb::arg("dmax")=0., nb::arg("unique")=true);

  nb::class_<CifToMtz>(m, "CifToMtz")
    .def(nb::init<>())
    .def_rw("title", &CifToMtz::title)
    .def_rw("history", &CifToMtz::history)
    .def_rw("spec_lines", &CifToMtz::spec_lines)
    .def("convert_block_to_mtz", &CifToMtz::convert_block_to_mtz,
         nb::arg("rblock"), nb::arg("logger")=nb::make_tuple(nb::none(), 6))
    ;

  nb::class_<MtzToCif>(m, "MtzToCif")
    .def(nb::init<>())
    .def_rw("spec_lines", &MtzToCif::spec_lines)
    .def_rw("with_comments", &MtzToCif::with_comments)
    .def_rw("with_history", &MtzToCif::with_history)
    .def_rw("skip_empty", &MtzToCif::skip_empty)
    .def_rw("skip_negative_sigi", &MtzToCif::skip_negative_sigi)
    .def_rw("wavelength", &MtzToCif::wavelength)
    .def_rw("free_flag_value", &MtzToCif::free_flag_value)
    .def("write_cif_to_string", [](MtzToCif& self, const Mtz& mtz, const Mtz* mtz2) {
        std::ostringstream out;
        self.write_cif(mtz, mtz2, nullptr, out);
        return out.str();
    }, nb::arg("mtz"), nb::arg("mtz2")=nb::none())
    ;
  m.def("validate_merged_intensities", &validate_merged_intensities,
        nb::arg(), nb::arg(), nb::arg("relaxed_check")=false, nb::arg("logger"));

  m.def("get_software_from_mtz_history", &get_software_from_mtz_history);

  nb::enum_<DataType>(m, "DataType")
    .value("Unknown", DataType::Unknown)
    .value("Unmerged", DataType::Unmerged)
    .value("Mean", DataType::Mean)
    .value("Anomalous", DataType::Anomalous)
    .value("MergedMA", DataType::MergedMA)
    .value("MergedAM", DataType::MergedAM)
    .value("UAM", DataType::UAM)
    ;
  m.def("check_data_type_under_symmetry", [](const ReflnBlock& data) {
      return check_data_type_under_symmetry(ReflnDataProxy(data));
  });
  m.def("check_data_type_under_symmetry", [](const Mtz& data) {
      return check_data_type_under_symmetry(MtzDataProxy{data});
  });

  nb::class_<MergingStats>(m, "MergingStats")
    .def_ro("all_refl", &MergingStats::all_refl)
    .def_ro("unique_refl", &MergingStats::unique_refl)
    .def_ro("stats_refl", &MergingStats::stats_refl)
    .def("r_merge", &MergingStats::r_merge)
    .def("r_meas", &MergingStats::r_meas)
    .def("r_pim", &MergingStats::r_pim)
    .def("cc_half", &MergingStats::cc_half)
    .def("cc_star", &MergingStats::cc_star)
    ;

  nb::class_<Intensities>(m, "Intensities")
    .def(nb::init<>())
    .def("clone", [](const Intensities& self) {return new Intensities(self); })
    .def("__len__", [](const Intensities& self) { return self.data.size(); })
    .def_rw("spacegroup", &Intensities::spacegroup)
    .def_rw("unit_cell", &Intensities::unit_cell)
    .def_rw("type", &Intensities::type)
    .def("resolution_range", &Intensities::resolution_range)
    .def("remove_systematic_absences", &Intensities::remove_systematic_absences)
    .def("sort", &Intensities::sort)
    .def("merge_in_place", &Intensities::merge_in_place, nb::arg("new_type"))
    .def("calculate_merging_stats", &Intensities::calculate_merging_stats,
         nb::arg("binner").none(), nb::arg("use_weights")='Y')
    .def("prepare_for_merging", &Intensities::prepare_for_merging)
    .def("calculate_correlation", &Intensities::calculate_correlation)
    .def("import_mtz", &Intensities::import_mtz,
         nb::arg(), nb::arg("type")=DataType::Unknown)
    .def("import_xds", &Intensities::import_xds)
    .def("import_refln_block", &Intensities::import_refln_block,
         nb::arg(), nb::arg("type")=DataType::Unknown)
    .def("prepare_merged_mtz", &Intensities::prepare_merged_mtz,
         nb::arg("with_nobs"))
    .def_prop_ro("miller_array", [](Intensities& self) {
      int64_t stride = static_cast<int64_t>(sizeof(Intensities::Refl) / sizeof(int));
      return nb::ndarray<nb::numpy, int, nb::shape<-1,3>>(
          &self.data.data()->hkl[0], {self.data.size(), 3},
          nb::handle(), {stride, 1});
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("value_array", [](Intensities& self) {
      return vector_member_array(self.data, &Intensities::Refl::value);
      //const Intensities::Refl* data = self.data.data();
      //nb::ssize_t stride = (const char*)(data+1) - (const char*)data;
      //return nb::array_t<double>({(nb::ssize_t)self.data.size()}, {stride},
      //                           &data->value, nb::cast(self));
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("sigma_array", [](Intensities& self) {
      return vector_member_array(self.data, &Intensities::Refl::sigma);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("nobs_array", [](Intensities& self) {
      return vector_member_array(self.data, &Intensities::Refl::nobs);
    }, nb::rv_policy::reference_internal)
    .def_prop_ro("isign_array", [](Intensities& self) {
      return vector_member_array(self.data, &Intensities::Refl::isign);
    }, nb::rv_policy::reference_internal)
    .def("set_data", [](Intensities& self,
                        const UnitCell& unit_cell,
                        const SpaceGroup* sg,
                        const cpu_miller_array& hkl,
                        const cpu_array<double>& values,
                        const cpu_array<double>& sigmas) {
      auto h = hkl.view();
      auto v = values.view();
      auto s = sigmas.view();
      if (h.shape(0) != v.shape(0) || h.shape(0) != s.shape(0))
        throw std::domain_error("arrays have different lengths");
      self.unit_cell = unit_cell;
      self.spacegroup = sg;
      self.data.clear();
      self.data.reserve(h.shape(0));
      for (size_t i = 0; i < h.shape(0); ++i)
        self.add_if_valid({h(i, 0), h(i, 1), h(i, 2)}, 0, 0, v(i), s(i));
      self.type = DataType::Unmerged;
      self.switch_to_asu_indices();
    }, nb::arg("cell"), nb::arg("sg").none(false),
       nb::arg("miller_array"), nb::arg("value_array"), nb::arg("sigma_array"))
    ;

  nb::class_<Binner> binner(m, "Binner");
  nb::enum_<Binner::Method>(binner, "Method")
      .value("EqualCount", Binner::Method::EqualCount)
      .value("Dstar", Binner::Method::Dstar)
      .value("Dstar2", Binner::Method::Dstar2)
      .value("Dstar3", Binner::Method::Dstar3)
      ;
  binner
    .def(nb::init<>())
    .def("setup", [](Binner& self, int nbins, Binner::Method method,
                     const Mtz& mtz, const UnitCell* cell) {
        self.setup(nbins, method, MtzDataProxy{mtz}, cell);
    }, nb::arg("nbins"), nb::arg("method"), nb::arg("mtz"), nb::arg("cell")=nb::none())
    .def("setup", [](Binner& self, int nbins, Binner::Method method,
                     const ReflnBlock& r, const UnitCell* cell) {
        self.setup(nbins, method, ReflnDataProxy(r), cell);
    }, nb::arg("nbins"), nb::arg("method"), nb::arg("r"), nb::arg("cell")=nb::none())
    .def("setup", [](Binner& self, int nbins, Binner::Method method,
                     const Intensities& intensities) {
        self.setup(nbins, method, gemmi::IntensitiesDataProxy{intensities});
    }, nb::arg("nbins"), nb::arg("method"), nb::arg("intensities"))
    .def("setup", [](Binner& self, int nbins, Binner::Method method,
                     const cpu_miller_array& hkl, const UnitCell* cell) {
        auto h = hkl.view();
        std::vector<double> inv_d2(h.shape(0));
        if (cell)
          for (size_t i = 0; i < inv_d2.size(); ++i)
            inv_d2[i] = cell->calculate_1_d2_double(h(i, 0), h(i, 1), h(i, 2));
        self.setup_from_1_d2(nbins, method, std::move(inv_d2), cell);
    }, nb::arg("nbins"), nb::arg("method"), nb::arg("hkl"), nb::arg("cell"))
    .def("setup_from_1_d2", [](Binner& self, int nbins, Binner::Method method,
                               const cpu_c_array<double>& inv_d2, const UnitCell* cell) {
        double* ptr = inv_d2.data();
        auto len = inv_d2.shape(0);
        self.setup_from_1_d2(nbins, method, std::vector<double>(ptr, ptr+len), cell);
    }, nb::arg("nbins"), nb::arg("method"), nb::arg("inv_d2"), nb::arg("cell"))
    .def("get_bin", &Binner::get_bin)
    .def("get_bins", [](Binner& self, const Mtz& mtz) {
        return numpy_array_from_vector(self.get_bins(MtzDataProxy{mtz}));
    })
    .def("get_bins", [](Binner& self, const ReflnBlock& r) {
        return numpy_array_from_vector(self.get_bins(ReflnDataProxy(r)));
    })
    .def("get_bins", [](Binner& self, const cpu_miller_array& hkl) {
        if (hkl.stride(1) != 1 || hkl.stride(0) < 3)
          throw std::domain_error("hkl array must be contiguous");
        struct {  // cf. MtzDataProxy
          size_t size_;
          size_t stride_;
          int* data_;
          size_t size() const noexcept { return size_; }
          size_t stride() const noexcept { return stride_; }
          Miller get_hkl(size_t offset) const noexcept {
            return {data_[offset], data_[offset+1], data_[offset+2]};
          }
        } proxy{hkl.size(), (size_t) hkl.stride(0), hkl.data()};
        return numpy_array_from_vector(self.get_bins(proxy));
    })
    .def("get_bins_from_1_d2", [](Binner& self, const cpu_c_array<double>& inv_d2) {
        return numpy_array_from_vector(self.get_bins_from_1_d2(inv_d2.data(), inv_d2.shape(0)));
    })
    .def("dmin_of_bin", &Binner::dmin_of_bin)
    .def("dmax_of_bin", &Binner::dmax_of_bin)
    .def_prop_ro("size", &Binner::size)
    .def_ro("limits", &Binner::limits)
    .def_rw("cell", &Binner::cell)
    .def_ro("min_1_d2", &Binner::min_1_d2)
    .def_ro("max_1_d2", &Binner::max_1_d2)
    ;

  m.def("combine_correlations", &combine_correlations);

  m.def("calculate_amplitude_normalizers",
        [](const Mtz& mtz, const std::string& f_col, const Binner& binner) {
      const Mtz::Column& f = mtz.get_column_with_label(f_col);
      return numpy_array_from_vector(
          calculate_amplitude_normalizers(MtzDataProxy{mtz}, f.idx, binner));
  });

  nb::class_<HklMatch>(m, "HklMatch")
    .def("__init__", [](HklMatch* p, const cpu_c_miller_array& hkl,
                                     const cpu_c_miller_array& ref) {
        static_assert(sizeof(Miller) == 3 * sizeof(int), "sizeof(Miller) problem");
        // NOLINTBEGIN(bugprone-casting-through-void)
        new(p) HklMatch(static_cast<Miller*>((void*)hkl.data()), hkl.shape(0),
                        static_cast<Miller*>((void*)ref.data()), ref.shape(0));
        // NOLINTEND(bugprone-casting-through-void)
    }, nb::arg("hkl"), nb::arg("ref"))
    .def("aligned", [](HklMatch& self, const cpu_c_array<double>& vec) {
        return numpy_array_from_vector(self.aligned_(vec.data(), vec.size(), (double)NAN));
    })
    .def_ro("pos", &HklMatch::pos)
    ;
}
