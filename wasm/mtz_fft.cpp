// Copyright 2019 Global Phasing Ltd.

#include "common.h"
#include "blob_search.h"

#define POCKETFFT_CACHE_SIZE 0
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include <emscripten/val.h>
#include <gemmi/fourier.hpp>
#include <gemmi/grid.hpp>
#include <gemmi/isosurface.hpp>
#include <gemmi/mtz.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using gemmi::Mtz;

namespace {

em::val float_view(const std::vector<float>& v) {
  return em::val(emscripten::typed_memory_view(v.size(), v.data()));
}

em::val uint32_view(const std::vector<uint32_t>& v) {
  return em::val(emscripten::typed_memory_view(v.size(), v.data()));
}

bool calculate_map_grid(const Mtz& mtz,
                        const std::string& buf,
                        const Mtz::Column* f_col,
                        const Mtz::Column* phi_col,
                        gemmi::Grid<float>& out_grid,
                        gemmi::DataStats& out_stats,
                        std::string& error) {
  try {
    const float* raw_data = reinterpret_cast<const float*>(buf.data() + 80);
    gemmi::MtzExternalDataProxy proxy(mtz, raw_data);
    auto size = gemmi::get_size_for_hkl(proxy, {{0, 0, 0}}, 3.);
    gemmi::FPhiProxy<gemmi::MtzExternalDataProxy> fphi(proxy, f_col->idx, phi_col->idx);
    gemmi::FPhiGrid<float> coefs =
      gemmi::get_f_phi_on_grid<float>(fphi, size, /*half_l=*/true, gemmi::AxisOrder::ZYX);

    gemmi::Grid<float> zyx_grid;
    gemmi::transform_f_phi_grid_to_map_(std::move(coefs), zyx_grid);

    out_grid.set_unit_cell(zyx_grid.unit_cell);
    out_grid.set_size(zyx_grid.nw, zyx_grid.nv, zyx_grid.nu);
    for (int x = 0; x < out_grid.nu; ++x)
      for (int y = 0; y < out_grid.nv; ++y)
        for (int z = 0; z < out_grid.nw; ++z)
          out_grid.data[out_grid.index_q(x, y, z)] = zyx_grid.get_value_q(z, y, x);

    out_stats = gemmi::calculate_data_statistics(out_grid.data);
  } catch (std::runtime_error& e) {
    error = e.what();
    return false;
  }
  return true;
}

class MtzMap {
public:
  MtzMap(gemmi::Grid<float>&& grid, const gemmi::DataStats& stats)
    : grid_(std::move(grid)), stats_(stats) {}

  em::val data() const { return float_view(grid_.data); }

  bool extract_isosurface(double radius, double x, double y, double z,
                          double isolevel, const std::string& method) {
    try {
      iso_ = gemmi::extract_isosurface(grid_, {x, y, z}, radius,
                                       isolevel, gemmi::iso_method_from_string(method));
    } catch (std::exception& e) {
      last_error_ = e.what();
      return false;
    }
    return true;
  }

  em::val isosurface_vertices() const { return float_view(iso_.vertices); }
  em::val isosurface_segments() const { return uint32_view(iso_.triangles); }

  int get_nx() const { return grid_.nu; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nw; }
  double get_mean() const { return stats_.dmean; }
  double get_rms() const { return stats_.rms; }
  std::string get_last_error() const { return last_error_; }
  gemmi::UnitCell get_cell() const { return grid_.unit_cell; }
  blob_wasm::BlobSearchResult* find_blobs(double cutoff, double min_volume,
                                          double min_score, double min_peak,
                                          bool negate, gemmi::Structure* st,
                                          int model_index, double mask_radius,
                                          bool mask_waters) const {
    return blob_wasm::find_blobs(grid_, cutoff, min_volume, min_score, min_peak,
                                 negate, st, model_index, mask_radius, mask_waters);
  }

private:
  gemmi::Grid<float> grid_;
  gemmi::DataStats stats_;
  gemmi::IsoSurface iso_;
  std::string last_error_;
};

const Mtz::Column* find_map_columns(const Mtz& mtz, bool diff_map,
                                     const Mtz::Column*& phi_col) {
  static const char* normal_labels[] = {
    "FWT", "PHWT",
    "2FOFCWT", "PH2FOFCWT",
  };
  static const char* diff_labels[] = {
    "DELFWT", "PHDELWT",
    "FOFCWT", "PHFOFCWT",
  };
  const char** labels = diff_map ? diff_labels : normal_labels;
  for (int i = 0; i < 4; i += 2)
    if (const Mtz::Column* f_col = mtz.column_with_label(labels[i]))
      if ((phi_col = mtz.column_with_label(labels[i + 1])))
        return f_col;
  return nullptr;
}

}  // namespace

class MtzFft {
public:
  MtzFft(std::string buf): buf_(std::move(buf)) {}

  bool read() {
    try {
      gemmi::MemoryStream stream(buf_.data(), buf_.size());
      mtz_.read_stream(stream, false);
    } catch (std::runtime_error& e) {
      last_error_ = "Failed to read MTZ file: ";
      last_error_ += e.what();
      return false;
    }
    return true;
  }

  em::val calculate_map_from_columns(const Mtz::Column* f_col,
                                     const Mtz::Column* phi_col) {
    if (!calculate_map_grid(mtz_, buf_, f_col, phi_col, grid_, stats_, last_error_))
      return em::val::null();
    return float_view(grid_.data);
  }

  em::val calculate_map_from_labels(const std::string& f_label,
                                    const std::string& phi_label) {
    if (const Mtz::Column* f_col = mtz_.column_with_label(f_label))
      if (const Mtz::Column* phi_col = mtz_.column_with_label(phi_label))
        return calculate_map_from_columns(f_col, phi_col);
    last_error_ = "Requested labels not found in the MTZ file.";
    return em::val::null();
  }

  em::val calculate_map(bool diff_map) {
    const Mtz::Column* phi_col = nullptr;
    if (const Mtz::Column* f_col = find_map_columns(mtz_, diff_map, phi_col))
      return calculate_map_from_columns(f_col, phi_col);
    last_error_ = "Default map coefficient labels not found";
    return em::val::null();
  }

  MtzMap* calculate_wasm_map_from_columns(const Mtz::Column* f_col,
                                          const Mtz::Column* phi_col) {
    gemmi::Grid<float> grid;
    gemmi::DataStats stats;
    if (!calculate_map_grid(mtz_, buf_, f_col, phi_col, grid, stats, last_error_))
      return nullptr;
    return new MtzMap(std::move(grid), stats);
  }

  MtzMap* calculate_wasm_map_from_labels(const std::string& f_label,
                                         const std::string& phi_label) {
    if (const Mtz::Column* f_col = mtz_.column_with_label(f_label))
      if (const Mtz::Column* phi_col = mtz_.column_with_label(phi_label))
        return calculate_wasm_map_from_columns(f_col, phi_col);
    last_error_ = "Requested labels not found in the MTZ file.";
    return nullptr;
  }

  MtzMap* calculate_wasm_map(bool diff_map) {
    const Mtz::Column* phi_col = nullptr;
    if (const Mtz::Column* f_col = find_map_columns(mtz_, diff_map, phi_col))
      return calculate_wasm_map_from_columns(f_col, phi_col);
    last_error_ = "Default map coefficient labels not found";
    return nullptr;
  }

  int get_nx() const { return grid_.nu; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nw; }
  double get_rmsd() const { return stats_.rms; }
  std::string get_last_error() const { return last_error_; }
  gemmi::UnitCell get_cell() const { return mtz_.cell; }

private:
  gemmi::Mtz mtz_;
  std::string buf_;
  gemmi::Grid<float> grid_;
  gemmi::DataStats stats_;
  std::string last_error_;
};

void add_mtz_fft() {
  em::class_<MtzMap>("MtzMap")
    .function("data", &MtzMap::data)
    .function("extract_isosurface", &MtzMap::extract_isosurface)
    .function("isosurface_vertices", &MtzMap::isosurface_vertices)
    .function("isosurface_segments", &MtzMap::isosurface_segments)
    .function("find_blobs", &MtzMap::find_blobs, em::allow_raw_pointers())
    .property("nx", &MtzMap::get_nx)
    .property("ny", &MtzMap::get_ny)
    .property("nz", &MtzMap::get_nz)
    .property("mean", &MtzMap::get_mean)
    .property("rms", &MtzMap::get_rms)
    .property("last_error", &MtzMap::get_last_error)
    .property("cell", &MtzMap::get_cell)
    ;

  em::class_<MtzFft>("Mtz")
    .constructor<std::string>()
    .function("read", &MtzFft::read)
    .function("calculate_map", &MtzFft::calculate_map)
    .function("calculate_map_from_labels", &MtzFft::calculate_map_from_labels)
    .function("calculate_wasm_map", &MtzFft::calculate_wasm_map, em::allow_raw_pointers())
    .function("calculate_wasm_map_from_labels", &MtzFft::calculate_wasm_map_from_labels,
              em::allow_raw_pointers())
    .property("nx", &MtzFft::get_nx)
    .property("ny", &MtzFft::get_ny)
    .property("nz", &MtzFft::get_nz)
    .property("rmsd", &MtzFft::get_rmsd)
    .property("last_error", &MtzFft::get_last_error)
    .property("cell", &MtzFft::get_cell)
    ;
}

#ifdef STANDALONE_MTZ
EMSCRIPTEN_BINDINGS(GemmiMtz) {
  add_cell();
  add_mtz_fft();
}
#endif
